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
#include <boost/bind.hpp>
using namespace boost;

#include "TwoStepBDNVTGPU.h"
#include "TwoStepNVEGPU.cuh"
#include "TwoStepBDNVTGPU.cuh"

#ifdef ENABLE_MPI
#include "HOOMDMPI.h"
#endif

/*! \file TwoStepBDNVTGPU.h
    \brief Contains code for the TwoStepBDNVTGPU class
*/

/*! \param sysdef SystemDefinition this method will act on. Must not be NULL.
    \param group The group of particles this integration method is to work on
    \param T Temperature set point as a function of time
    \param seed Random seed to use in generating random numbers
    \param gamma_diam Set gamma to the particle diameter of each particle if true, otherwise use a per-type
                      gamma via setGamma()
    \param suffix Suffix to attach to the end of log quantity names
*/
TwoStepBDNVTGPU::TwoStepBDNVTGPU(boost::shared_ptr<SystemDefinition> sysdef,
                                 boost::shared_ptr<ParticleGroup> group,
                                 boost::shared_ptr<Variant> T,
                                 unsigned int seed,
                                 bool gamma_diam,
                                 const std::string& suffix)
    : TwoStepBDNVT(sysdef, group, T, seed, gamma_diam, suffix)
    {
    // only one GPU is supported
    if (!exec_conf->isCUDAEnabled())
        {
        m_exec_conf->msg->error() << "Creating a TwoStepNVEGPU what CUDA is disabled" << endl;
        throw std::runtime_error("Error initializing TwoStepNVEGPU");
        }
        
    // allocate the sum arrays
    GPUArray<float> sum(1, exec_conf);
    m_sum.swap(sum);
    
    // initialize the partial sum array
    m_block_size = 256; 
    unsigned int group_size = m_group->getNumMembers();
    m_num_blocks = group_size / m_block_size + 1;
    GPUArray<float> partial_sum1(m_num_blocks, exec_conf);
    m_partial_sum1.swap(partial_sum1);          
    }

/*! \param timestep Current time step
    \post Particle positions are moved forward to timestep+1 and velocities to timestep+1/2 per the velocity verlet
          method.
    
    This method is copied directoy from TwoStepNVEGPU::integrateStepOne() and reimplemented here to avoid multiple
    inheritance.
*/
void TwoStepBDNVTGPU::integrateStepOne(unsigned int timestep)
    {
    // profile this step
    if (m_prof)
        m_prof->push(exec_conf, "NVE step 1");
    
    // access all the needed data
    BoxDim box = m_pdata->getBox();
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
    unsigned int group_size = m_group->getNumMembers();

    ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
    ArrayHandle<Scalar3> d_accel(m_pdata->getAccelerations(), access_location::device, access_mode::readwrite);
    ArrayHandle<int3> d_image(m_pdata->getImages(), access_location::device, access_mode::readwrite);

    // perform the update on the GPU
    gpu_nve_step_one(d_pos.data,
                     d_vel.data,
                     d_accel.data,
                     d_image.data,
                     d_index_array.data,
                     group_size,
                     box,
                     m_deltaT,
                     m_limit,
                     m_limit_val,
                     m_zero_force);

    if (exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    
    // done profiling
    if (m_prof)
        m_prof->pop(exec_conf);
    }
        
/*! \param timestep Current time step
    \post particle velocities are moved forward to timestep+1 on the GPU
*/
void TwoStepBDNVTGPU::integrateStepTwo(unsigned int timestep)
    {
    const GPUArray< Scalar4 >& net_force = m_pdata->getNetForce();
    
    // profile this step
    if (m_prof)
        m_prof->push(exec_conf, "NVE step 2");
    
    // get the dimensionality of the system
    const Scalar D = Scalar(m_sysdef->getNDimensions());
    
    ArrayHandle<Scalar4> d_net_force(net_force, access_location::device, access_mode::read);
    ArrayHandle<Scalar> d_gamma(m_gamma, access_location::device, access_mode::read);
    ArrayHandle< unsigned int > d_index_array(m_group->getIndexArray(), access_location::device, access_mode::read);
 
        {
        ArrayHandle<float> d_partial_sumBD(m_partial_sum1, access_location::device, access_mode::overwrite);
        ArrayHandle<float> d_sumBD(m_sum, access_location::device, access_mode::overwrite);
        ArrayHandle<Scalar4> d_pos(m_pdata->getPositions(), access_location::device, access_mode::read);
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);
        ArrayHandle<Scalar3> d_accel(m_pdata->getAccelerations(), access_location::device, access_mode::readwrite);
        ArrayHandle<Scalar> d_diameter(m_pdata->getDiameters(), access_location::device, access_mode::read);
        ArrayHandle<unsigned int> d_tag(m_pdata->getTags(), access_location::device, access_mode::read);
        
        unsigned int group_size = m_group->getNumMembers();
        m_num_blocks = group_size / m_block_size + 1;

        // perform the update on the GPU
        bdnvt_step_two_args args;
        args.d_gamma = d_gamma.data;
        args.n_types = m_gamma.getNumElements();
        args.gamma_diam = m_gamma_diam;
        args.T = m_T->getValue(timestep);
        args.timestep = timestep;
        args.seed = m_seed;
        args.d_sum_bdenergy = d_sumBD.data;
        args.d_partial_sum_bdenergy = d_partial_sumBD.data;
        args.block_size = m_block_size;
        args.num_blocks = m_num_blocks;
        args.tally = m_tally;
        
        gpu_bdnvt_step_two(d_pos.data,
                           d_vel.data,
                           d_accel.data,
                           d_diameter.data,
                           d_tag.data,
                           d_index_array.data,
                           group_size,
                           d_net_force.data,
                           args,
                           m_deltaT,
                           D,
                           m_limit,
                           m_limit_val);

        if (exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        
        }
 
    if (m_tally)
        {
        ArrayHandle<float> h_sumBD(m_sum, access_location::host, access_mode::read);   
        #ifdef ENABLE_MPI
        if (m_comm)
            {
            MPI_Allreduce(MPI_IN_PLACE, &h_sumBD.data[0], 1, MPI_HOOMD_SCALAR, MPI_SUM, m_exec_conf->getMPICommunicator()); 
            } 
        #endif
        m_reservoir_energy -= h_sumBD.data[0]*m_deltaT;
        m_extra_energy_overdeltaT= 0.5*h_sumBD.data[0];
        }
    // done profiling
    if (m_prof)
        m_prof->pop(exec_conf);
    }

void export_TwoStepBDNVTGPU()
    {
    class_<TwoStepBDNVTGPU, boost::shared_ptr<TwoStepBDNVTGPU>, bases<TwoStepBDNVT>, boost::noncopyable>
        ("TwoStepBDNVTGPU", init< boost::shared_ptr<SystemDefinition>,
                         boost::shared_ptr<ParticleGroup>,
                         boost::shared_ptr<Variant>,
                         unsigned int,
                         bool,
                         const std::string&
                         >())
        ;
    }

#ifdef WIN32
#pragma warning( pop )
#endif

