/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2008-2011 Ames Laboratory
Iowa State University and The Regents of the University of Michigan All rights
reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

You may redistribute, use, and create derivate works of HOOMD-blue, in source
and binary forms, provided you abide by the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer both in the code and
prominently in any materials provided with the distribution.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* All publications and presentations based on HOOMD-blue, including any reports
or published results obtained, in whole or in part, with HOOMD-blue, will
acknowledge its use according to the terms posted at the time of submission on:
http://codeblue.umich.edu/hoomd-blue/citations.html

* Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
http://codeblue.umich.edu/hoomd-blue/

* Apart from the above required attributions, neither the name of the copyright
holder nor the names of HOOMD-blue's contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Maintainer: joaander

#include <boost/python.hpp>
using namespace boost::python;

#include "TwoStepBerendsen.h"
#ifdef ENABLE_CUDA
#include "TwoStepBerendsenGPU.cuh"
#endif

/*! \file TwoStepBerendsen.cc
    \brief Definition of Berendsen thermostat
*/

// ********************************
// here follows the code for Berendsen on the CPU

/*! \param sysdef System to zero the velocities of
    \param group Group of particles on which this method will act
    \param thermo compute for thermodynamic quantities
    \param tau Berendsen time constant
    \param T Temperature set point
*/
TwoStepBerendsen::TwoStepBerendsen(boost::shared_ptr<SystemDefinition> sysdef,
                                   boost::shared_ptr<ParticleGroup> group,
                                   boost::shared_ptr<ComputeThermo> thermo,
                                   Scalar tau,
                                   boost::shared_ptr<Variant> T)
    : IntegrationMethodTwoStep(sysdef, group), m_thermo(thermo), m_tau(tau), m_T(T)
    {
    m_exec_conf->msg->notice(5) << "Constructing TwoStepBerendsen" << endl;

    if (m_tau <= 0.0)
        m_exec_conf->msg->warning() << "integrate.berendsen: tau set less than 0.0" << endl;
    }

TwoStepBerendsen::~TwoStepBerendsen()
    {
    m_exec_conf->msg->notice(5) << "Destroying TwoStepBerendsen" << endl;
    }

/*! Perform the needed calculations to zero the system's velocity
    \param timestep Current time step of the simulation
*/
void TwoStepBerendsen::integrateStepOne(unsigned int timestep)
    {
    unsigned int group_size = m_group->getNumMembers();
    if (group_size == 0)
        return;

    // profile this step
    if (m_prof)
        m_prof->push("Berendsen step 1");

    // compute the current thermodynamic properties and get the temperature
    m_thermo->compute(timestep);
    Scalar curr_T = m_thermo->getTemperature();

    // compute the value of lambda for the current timestep
    Scalar lambda = sqrt(Scalar(1.0) + m_deltaT / m_tau * (m_T->getValue(timestep) / curr_T - Scalar(1.0)));

    // access the particle data for writing on the CPU
    assert(m_pdata);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);


    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // advance velocity forward by half a timestep and position forward by a full timestep
        h_vel.data[j].x = lambda * (h_vel.data[j].x + h_accel.data[j].x * m_deltaT * Scalar(1.0 / 2.0));
        h_pos.data[j].x += h_vel.data[j].x * m_deltaT;

        h_vel.data[j].y = lambda * (h_vel.data[j].y + h_accel.data[j].y * m_deltaT * Scalar(1.0 / 2.0));
        h_pos.data[j].y += h_vel.data[j].y * m_deltaT;

        h_vel.data[j].z = lambda * (h_vel.data[j].z + h_accel.data[j].z * m_deltaT * Scalar(1.0 / 2.0));
        h_pos.data[j].z += h_vel.data[j].z * m_deltaT;
        }

    /* particles may have been moved slightly outside the box by the above steps so we should wrap
        them back into place */
    const BoxDim& box = m_pdata->getBox();

    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);

    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);
        box.wrap(h_pos.data[j], h_image.data[j]);
        }

    if (m_prof)
        m_prof->pop();
    }

/*! \param timestep Current timestep
    \post particle velocities are moved forward to timestep+1
*/
void TwoStepBerendsen::integrateStepTwo(unsigned int timestep)
    {
    unsigned int group_size = m_group->getNumMembers();
    if (group_size == 0)
        return;

    // access the particle data for writing on the CPU
    assert(m_pdata);
    ArrayHandle<Scalar4> h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(m_pdata->getAccelerations(), access_location::host, access_mode::readwrite);

    // access the force data
    const GPUArray< Scalar4 >& net_force = m_pdata->getNetForce();
    ArrayHandle< Scalar4 > h_net_force(net_force, access_location::host, access_mode::read);

    // profile this step
    if (m_prof)
        m_prof->push("Berendsen step 2");

    // integrate the particle velocities to timestep+1
    for (unsigned int group_idx = 0; group_idx < group_size; group_idx++)
        {
        unsigned int j = m_group->getMemberIndex(group_idx);

        // calculate the acceleration from the net force
        Scalar minv = Scalar(1.0) / h_vel.data[j].w;
        h_accel.data[j].x = h_net_force.data[j].x * minv;
        h_accel.data[j].y = h_net_force.data[j].y * minv;
        h_accel.data[j].z = h_net_force.data[j].z * minv;

        // update the velocity
        h_vel.data[j].x += h_accel.data[j].x * m_deltaT / Scalar(2.0);
        h_vel.data[j].y += h_accel.data[j].y * m_deltaT / Scalar(2.0);
        h_vel.data[j].z += h_accel.data[j].z * m_deltaT / Scalar(2.0);
        }

    }

void export_Berendsen()
    {
    class_<TwoStepBerendsen, boost::shared_ptr<TwoStepBerendsen>, bases<IntegrationMethodTwoStep>, boost::noncopyable>
    ("TwoStepBerendsen", init< boost::shared_ptr<SystemDefinition>,
                         boost::shared_ptr<ParticleGroup>,
                         boost::shared_ptr<ComputeThermo>,
                         Scalar,
                         boost::shared_ptr<Variant>
                         >())
        .def("setT", &TwoStepBerendsen::setT)
        .def("setTau", &TwoStepBerendsen::setTau)
        ;
    }

