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

// Maintainer: joaander, grva, baschult

/*! \file ForceCompute.cc
    \brief Defines the ForceCompute class
*/

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4244 )
#endif

#include "ForceCompute.h"
#include <iostream>
using namespace std;

#include <boost/python.hpp>
using namespace boost::python;

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
using namespace boost;

#ifdef ENABLE_MPI
#include "Communicator.h"
#endif

/*! \param sysdef System to compute forces on
    \post The Compute is initialized and all memory needed for the forces is allocated
    \post \c force and \c virial GPUarrays are initialized
    \post All forces are initialized to 0
*/
ForceCompute::ForceCompute(boost::shared_ptr<SystemDefinition> sysdef) : Compute(sysdef), m_particles_sorted(false),
    m_index_thread_partial(0)
    {
    assert(m_pdata);
    assert(m_pdata->getMaxN() > 0);
    
    // allocate data on the host
    unsigned int max_num_particles = m_pdata->getMaxN();
    GPUArray<Scalar4>  force(max_num_particles,exec_conf);
    GPUArray<Scalar>   virial(max_num_particles,6,exec_conf);
    GPUArray<Scalar4>  torque(max_num_particles,exec_conf);
    m_force.swap(force);
    m_virial.swap(virial);
    m_torque.swap(torque);

    m_virial_pitch = m_virial.getPitch();
    m_fdata_partial = NULL;
    m_virial_partial = NULL;
    m_torque_partial = NULL;
  
    // connect to the ParticleData to recieve notifications when particles change order in memory
    m_sort_connection = m_pdata->connectParticleSort(bind(&ForceCompute::setParticlesSorted, this));

    // connect to the ParticleData to receive notifications when the maximum number of particles changes
    m_max_particle_num_change_connection = m_pdata->connectMaxParticleNumberChange(bind(&ForceCompute::reallocate, this));

    // reset external virial
    for (unsigned int i = 0; i < 6; ++i)
        m_external_virial[i] = Scalar(0.0);

    }

/*! \post m_fdata and virial _partial are both allocated, and m_index_thread_partial is intiialized for indexing them
*/
void ForceCompute::allocateThreadPartial()
    {
    assert(exec_conf->n_cpu >= 1);
    m_index_thread_partial = Index2D(m_pdata->getMaxN(), exec_conf->n_cpu);
    //Don't use GPU arrays here, *_partial's only used on CPU
    m_fdata_partial = new Scalar4[m_index_thread_partial.getNumElements()];
    m_virial_partial = new Scalar[6*m_index_thread_partial.getNumElements()];
    m_torque_partial = new Scalar4[m_index_thread_partial.getNumElements()];
    }

/*! \post m_force, m_virial and m_torque are resized to the current maximum particle number
 */
void ForceCompute::reallocate()
    {
    m_force.resize(m_pdata->getMaxN());
    m_virial.resize(m_pdata->getMaxN(),6);
    m_torque.resize(m_pdata->getMaxN());

    // the pitch of the virial array may have changed
    m_virial_pitch = m_virial.getPitch();

    reallocateThreadPartial();
    }

/*! \post m_fdata and virial _partial are both reallocated, and m_index_thread_partial is intiialized for indexing them
*/
void ForceCompute::reallocateThreadPartial()
    {
    assert(exec_conf->n_cpu >= 1);

    // never allocated ? do nothing.
    if (! m_fdata_partial || ! m_virial_partial || ! m_torque_partial) return;

    delete[] m_fdata_partial;
    delete[] m_virial_partial;
    delete[] m_torque_partial;
    allocateThreadPartial();
    }

/*! Frees allocated memory
*/
ForceCompute::~ForceCompute()
    {
    if (m_fdata_partial)
        {
        delete[] m_fdata_partial;
        m_fdata_partial = NULL;
        delete[] m_virial_partial;
        m_virial_partial = NULL;
        delete[] m_torque_partial;
        m_torque_partial=NULL;
        }
    m_sort_connection.disconnect();
    m_max_particle_num_change_connection.disconnect();
    }

/*! Sums the total potential energy calculated by the last call to compute() and returns it.
*/
Scalar ForceCompute::calcEnergySum()
    {
    ArrayHandle<Scalar4> h_force(m_force,access_location::host,access_mode::read);   
    // always perform the sum in double precision for better accuracy
    // this is cheating and is really just a temporary hack to get logging up and running
    // the potential accuracy loss in simulations needs to be evaluated here and a proper
    // summation algorithm put in place
    double pe_total = 0.0;
    for (unsigned int i=0; i < m_pdata->getN(); i++)
        {
        pe_total += (double)h_force.data[i].w;
        }
#ifdef ENABLE_MPI
    if (m_comm)
        {
        // reduce potential energy on all processors
        MPI_Allreduce(MPI_IN_PLACE, &pe_total, 1, MPI_DOUBLE, MPI_SUM, m_exec_conf->getMPICommunicator());
        }
#endif
    return Scalar(pe_total);
    }

/*! Performs the force computation.
    \param timestep Current Timestep
    \note If compute() has previously been called with a value of timestep equal to
        the current value, the forces are assumed to already have been computed and nothing will
        be done
*/

void ForceCompute::compute(unsigned int timestep)
    {
    // skip if we shouldn't compute this step
    if (!m_particles_sorted && !shouldCompute(timestep))
        return;
        
    computeForces(timestep);
    m_particles_sorted = false;
    }

/*! \param num_iters Number of iterations to average for the benchmark
    \returns Milliseconds of execution time per calculation

    Calls computeForces repeatedly to benchmark the force compute.
*/

double ForceCompute::benchmark(unsigned int num_iters)
    {
    ClockSource t;
    // warm up run
    computeForces(0);
    
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        {
        cudaThreadSynchronize();
        CHECK_CUDA_ERROR();
        }
#endif
    
    // benchmark
    uint64_t start_time = t.getTime();
    for (unsigned int i = 0; i < num_iters; i++)
        computeForces(0);
        
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        cudaThreadSynchronize();
#endif
    uint64_t total_time_ns = t.getTime() - start_time;
    
    // convert the run time to milliseconds
    return double(total_time_ns) / 1e6 / double(num_iters);
    }

/*! \param tag Global particle tag
    \returns Torque of particle referenced by tag
 */
Scalar4 ForceCompute::getTorque(unsigned int tag)
    {
    unsigned int i = m_pdata->getRTag(tag);
    bool found = (i < m_pdata->getN());
    Scalar4 result = make_scalar4(0.0,0.0,0.0,0.0);
    if (found)
        {
        ArrayHandle<Scalar4> h_torque(m_torque, access_location::host, access_mode::read);
        result = h_torque.data[i];
        }
    #ifdef ENABLE_MPI
    if (m_pdata->getDomainDecomposition())
        {
        unsigned int owner_rank = m_pdata->getOwnerRank(tag);
        MPI_Bcast(&result, sizeof(Scalar4), MPI_BYTE, owner_rank, m_exec_conf->getMPICommunicator());
        }
    #endif
    return result;
    }

/*! \param tag Global particle tag
    \returns Force of particle referenced by tag
 */
Scalar3 ForceCompute::getForce(unsigned int tag)
    {
    unsigned int i = m_pdata->getRTag(tag);
    bool found = (i < m_pdata->getN());
    Scalar3 result = make_scalar3(0.0,0.0,0.0);
    if (found)
        {
        ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::read);
        result = make_scalar3(h_force.data[i].x, h_force.data[i].y, h_force.data[i].z);
        }
    #ifdef ENABLE_MPI
    if (m_pdata->getDomainDecomposition())
        {
        unsigned int owner_rank = m_pdata->getOwnerRank(tag);
        MPI_Bcast(&result, sizeof(Scalar3), MPI_BYTE, owner_rank, m_exec_conf->getMPICommunicator());
        }
    #endif
    return result;
    }

/*! \param tag Global particle tag
    \param component Virial component (0=xx, 1=xy, 2=xz, 3=yy, 4=yz, 5=zz)
    \returns Force of particle referenced by tag
 */
Scalar ForceCompute::getVirial(unsigned int tag, unsigned int component)
    {
    unsigned int i = m_pdata->getRTag(tag);
    bool found = (i < m_pdata->getN());
    Scalar result = Scalar(0.0);
    if (found)
        {
        ArrayHandle<Scalar> h_virial(m_virial, access_location::host, access_mode::read);
        result = h_virial.data[m_virial_pitch*component+i];
        }
    #ifdef ENABLE_MPI
    if (m_pdata->getDomainDecomposition())
        {
        unsigned int owner_rank = m_pdata->getOwnerRank(tag);
        MPI_Bcast(&result, sizeof(Scalar), MPI_BYTE, owner_rank, m_exec_conf->getMPICommunicator());
        }
    #endif
    return result;
    }

/*! \param tag Global particle tag
    \returns Energy of particle referenced by tag
 */
Scalar ForceCompute::getEnergy(unsigned int tag)
    {
    unsigned int i = m_pdata->getRTag(tag);
    bool found = (i < m_pdata->getN());
    Scalar result = Scalar(0.0);
    if (found)
        {
        ArrayHandle<Scalar4> h_force(m_force, access_location::host, access_mode::read);
        result = h_force.data[i].w;
        }
    #ifdef ENABLE_MPI
    if (m_pdata->getDomainDecomposition())
        {
        unsigned int owner_rank = m_pdata->getOwnerRank(tag);
        MPI_Bcast(&result, sizeof(Scalar), MPI_BYTE, owner_rank, m_exec_conf->getMPICommunicator());
        }
    #endif
    return result;
    }


//! Wrapper class for wrapping pure virtual methodos of ForceCompute in python
class ForceComputeWrap : public ForceCompute, public wrapper<ForceCompute>
    {
    public:
        //! Constructor
        /*! \param sysdef Particle data passed to the base class */
        ForceComputeWrap(shared_ptr<SystemDefinition> sysdef) : ForceCompute(sysdef) { }
    protected:
        //! Calls the overidden ForceCompute::computeForces()
        /*! \param timestep parameter to pass on to the overidden method 
         */
        void computeForces(unsigned int timestep)
            {
            this->get_override("computeForces")(timestep);
            }
    };

void export_ForceCompute()
    {
    class_< ForceComputeWrap, boost::shared_ptr<ForceComputeWrap>, bases<Compute>, boost::noncopyable >
    ("ForceCompute", init< boost::shared_ptr<SystemDefinition> >())
    .def("getForce", &ForceCompute::getForce)
    .def("getTorque", &ForceCompute::getTorque)
    .def("getVirial", &ForceCompute::getVirial)
    .def("getEnergy", &ForceCompute::getEnergy)
    ;
    }

#ifdef WIN32
#pragma warning( pop )
#endif

