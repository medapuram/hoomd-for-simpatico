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

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif

#include <boost/python.hpp>
using namespace boost::python;

#include "IntegrationMethodTwoStep.h"

#ifdef ENABLE_MPI
#include "Communicator.h"
#endif

/*! \file IntegrationMethodTwoStep.h
    \brief Contains code for the IntegrationMethodTwoStep class
*/

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
    \post The method is constructed with the given particle data and a NULL profiler.
*/
IntegrationMethodTwoStep::IntegrationMethodTwoStep(boost::shared_ptr<SystemDefinition> sysdef,
                                                   boost::shared_ptr<ParticleGroup> group)
    : m_sysdef(sysdef), m_group(group), m_pdata(m_sysdef->getParticleData()), exec_conf(m_pdata->getExecConf()), 
      m_deltaT(Scalar(0.0)), m_valid_restart(false)
    {
    // sanity check
    assert(m_sysdef);
    assert(m_pdata);
    assert(m_group);
    m_exec_conf = exec_conf;

    m_integrator_id = m_sysdef->getIntegratorData()->registerIntegrator();
    }

/*! It is useful for the user to know where computation time is spent, so all integration methods
    should profile themselves. This method sets the profiler for them to use.
    This method does not need to be called, as Computes will not profile themselves
    on a NULL profiler
    \param prof Pointer to a profiler for the compute to use. Set to NULL
        (boost::shared_ptr<Profiler>()) to stop the
        analyzer from profiling itself.
    \note Derived classes MUST check if m_prof is set before calling any profiler methods.
*/
void IntegrationMethodTwoStep::setProfiler(boost::shared_ptr<Profiler> prof)
    {
    m_prof = prof;
    }

/*! \param deltaT New time step to set
*/
void IntegrationMethodTwoStep::setDeltaT(Scalar deltaT)
    {
    m_deltaT = deltaT;
    }


/*! \param v is the restart variables for the current integrator
    \param type is the type of expected integrator type
    \param nvariables is the expected number of variables
    
    If the either the integrator type or number of variables does not match the
    expected values, this function throws the appropriate warning and returns
    "false."  Otherwise, the function returns true.
*/
bool IntegrationMethodTwoStep::restartInfoTestValid(IntegratorVariables& v, std::string type, unsigned int nvariables)
    {
    bool good = true;
    if (v.type == "")
        good = false;
    else if (v.type != type && v.type != "")
        {
        m_exec_conf->msg->warning() << "Integrator #"<<  m_integrator_id <<" type "<< type <<" does not match type ";
        m_exec_conf->msg->warning() << v.type << " found in restart file. " << endl;
        m_exec_conf->msg->warning() << "Ensure that the integrator order is consistent for restarted simulations. " << endl;
        m_exec_conf->msg->warning() << "Continuing while ignoring restart information..." << endl;
        good = false;
        }
    else if (v.type == type)
        {
        if (v.variable.size() != nvariables)
            {
            m_exec_conf->msg->warning() << "Integrator #"<< m_integrator_id <<" type "<< type << endl;
            m_exec_conf->msg->warning() << "appears to contain bad or incomplete restart information. " << endl;
            m_exec_conf->msg->warning() << "Continuing while ignoring restart information..." << endl;
            good = false;
            }
        }
    return good;
    }

/*! \param query_group Group over which to count degrees of freedom.
    A majority of the integration methods add D degrees of freedom per particle in \a query_group that is also in the
    group assigned to the method. Hence, the base class IntegrationMethodTwoStep will implement that counting.
    Derived classes can ovveride if needed.
*/
unsigned int IntegrationMethodTwoStep::getNDOF(boost::shared_ptr<ParticleGroup> query_group)
    {
    // get the size of the intersecion between query_group and m_group
    unsigned int intersect_size = ParticleGroup::groupIntersection(query_group, m_group)->getNumMembersGlobal();
    
    return m_sysdef->getNDimensions() * intersect_size;
    }

/*! Checks that every particle in the group is valid. This method may be called by anyone wishing to make this
    error check.

    The base class defines a valid particle as one that does not belong to a rigid body (as this is the common case).
    Derived classes may override this method to perform custom checks.
*/
void IntegrationMethodTwoStep::validateGroup()
    {
    for (unsigned int gidx = 0; gidx < m_group->getNumMembersGlobal(); gidx++)
        {
        unsigned int tag = m_group->getMemberTag(gidx);
        if (m_pdata->isParticleLocal(tag))
            {
            ArrayHandle<unsigned int> h_body(m_pdata->getBodies(), access_location::host, access_mode::read);
            ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);

            unsigned int body = h_body.data[h_rtag.data[tag]];

            if (body != NO_BODY)
                {
                m_exec_conf->msg->error() << "Particle " << tag << " belongs to a rigid body. "
                     << "This integration method does not operate on rigid bodies" << endl;
                    
                throw std::runtime_error("Error initializing integration method");
                }
            }
        }
    }


void export_IntegrationMethodTwoStep()
    {
    class_<IntegrationMethodTwoStep, boost::shared_ptr<IntegrationMethodTwoStep>, boost::noncopyable>
        ("IntegrationMethodTwoStep", init< boost::shared_ptr<SystemDefinition>, boost::shared_ptr<ParticleGroup> >())
        .def("validateGroup", &IntegrationMethodTwoStep::validateGroup)
#ifdef ENABLE_MPI
        .def("setCommunicator", &IntegrationMethodTwoStep::setCommunicator)
#endif
        ;
    }

#ifdef WIN32
#pragma warning( pop )
#endif

