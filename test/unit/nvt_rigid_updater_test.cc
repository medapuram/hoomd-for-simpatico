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


#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 4244 )
#endif

#include <iostream>

//! name the boost unit test module
#define BOOST_TEST_MODULE TwoStepNVTRigidTests
#include "boost_utf_configure.h"

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "TwoStepNVTRigid.h"
#ifdef ENABLE_CUDA
#include "TwoStepNVTRigidGPU.h"
#endif

#include "IntegratorTwoStep.h"

#include "AllPairPotentials.h"
#include "BinnedNeighborList.h"
#include "Initializers.h"

#ifdef ENABLE_CUDA
#include "BinnedNeighborListGPU.h"
#endif

#include "saruprng.h"
#include <math.h>
#include <time.h>

using namespace std;
using namespace boost;

/*! \file nvt_rigid_updater_test.cc
    \brief Implements unit tests for TwoStepNVTRigid
    \ingroup unit_tests
*/


//! Tolerance for floating point comparisons
#ifdef SINGLE_PRECISION
const Scalar tol = Scalar(1e-2);
#else
const Scalar tol = 1e-3;
#endif

//! Typedef'd TwoStepNVTRigid class factory
typedef boost::function<shared_ptr<TwoStepNVTRigid> (shared_ptr<SystemDefinition> sysdef, 
                                                shared_ptr<ParticleGroup> group, Scalar T)> nvtup_creator;

void nvt_updater_energy_tests(nvtup_creator nvt_creator, const ExecutionConfiguration& exec_conf)
    {
#ifdef ENABLE_CUDA
    g_gpu_error_checking = true;
#endif
    
    // check that the nve updater can actually integrate particle positions and velocities correctly
    // start with a 2 particle system to keep things simple: also put everything in a huge box so boundary conditions
    // don't come into play
    unsigned int nbodies = 3000;
    unsigned int nparticlesperbody = 5;
    unsigned int N = nbodies * nparticlesperbody;
    Scalar box_length = 60.0;
    
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(N, BoxDim(box_length), 1, 0, 0, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata = sysdef->getParticleData();
    shared_ptr<ParticleSelector> selector_all(new ParticleSelectorTag(sysdef, 0, pdata->getN()-1));
    shared_ptr<ParticleGroup> group_all(new ParticleGroup(sysdef, selector_all));
    BoxDim box = pdata->getBox();
    
    Scalar temperature = 2.5;
    unsigned int steps = 10000;
    unsigned int sampling = 1000;
    
    // setup a simple initial state
    unsigned int ibody = 0;
    unsigned int iparticle = 0;
    Scalar x0 = box.xlo + 0.01;
    Scalar y0 = box.ylo + 0.01;
    Scalar z0 = box.zlo + 0.01;
    Scalar xspacing = 6.0f;
    Scalar yspacing = 2.0f;
    Scalar zspacing = 2.0f;
    
    unsigned int seed = 10483;
    boost::shared_ptr<Saru> random = boost::shared_ptr<Saru>(new Saru(seed));
    Scalar KE = Scalar(0.0);

    {
    ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_body(pdata->getVelocities(), access_location::host, access_mode::readwrite);
    
    // initialize bodies in a cubic lattice with some velocity profile
    for (unsigned int i = 0; i < nbodies; i++)
        {
        for (unsigned int j = 0; j < nparticlesperbody; j++)
            {
            h_pos.data[iparticle].x = x0 + 1.0 * j;
            h_pos.data[iparticle].y = y0 + 0.0;
            h_pos.data[iparticle].z = z0 + 0.0;
            
            h_vel.data[iparticle].x = random->d();
            h_vel.data[iparticle].y = random->d();
            h_vel.data[iparticle].z = random->d();
            
            KE += Scalar(0.5) * (h_vel.data[iparticle].x*h_vel.data[iparticle].x + h_vel.data[iparticle].x*h_vel.data[iparticle].y + h_vel.data[iparticle].z*h_vel.data[iparticle].z);
            
            h_body.data[iparticle] = ibody;
            
            iparticle++;
            }
            
        x0 += xspacing;
        if (x0 + xspacing >= box.xhi)
            {
            x0 = box.xlo + 0.01;
            
            y0 += yspacing;
            if (y0 + yspacing >= box.yhi)
                {
                y0 = box.ylo + 0.01;
                
                z0 += zspacing;
                if (z0 + zspacing >= box.zhi)
                    z0 = box.zlo + 0.01;
                }
            }
            
        ibody++;
        }
        
    assert(iparticle == N);
    
    }
    
    Scalar deltaT = Scalar(0.001);
    Scalar Q = Scalar(2.0);
    Scalar tau = sqrt(Q / (Scalar(3.0) * temperature));
    
    shared_ptr<TwoStepNVTRigid> two_step_nvt = nvt_creator(sysdef, group_all, temperature);
    shared_ptr<IntegratorTwoStep> nvt_up(new IntegratorTwoStep(sysdef, deltaT));
    nvt_up->addIntegrationMethod(two_step_nvt);

    shared_ptr<NeighborList> nlist(new NeighborList(sysdef, Scalar(2.5), Scalar(0.8)));
    shared_ptr<PotentialPairLJ> fc(new PotentialPairLJ(sysdef, nlist));
    fc->setRcut(0, 0, Scalar(1.122));
    
    // setup some values for alpha and sigma
    Scalar epsilon = Scalar(1.0);
    Scalar sigma = Scalar(1.0);
    Scalar alpha = Scalar(1.0);
    Scalar lj1 = Scalar(4.0) * epsilon * pow(sigma, Scalar(12.0));
    Scalar lj2 = alpha * Scalar(4.0) * epsilon * pow(sigma, Scalar(6.0));
    
    // specify the force parameters
    fc->setParams(0,0,make_scalar2(lj1,lj2));
    
    nvt_up->addForceCompute(fc);
    
    // initialize the rigid bodies
    sysdef->init();
    
    Scalar PE;
    
    shared_ptr<RigidData> rdata = sysdef->getRigidData();
    unsigned int nrigid_dof = rdata->getNumDOF();
    
    // Rescale particle velocities to desired temperature:
    Scalar current_temp = 2.0 * KE / nrigid_dof;
    Scalar factor = sqrt(temperature / current_temp);
    
    {
    ArrayHandle< Scalar4 > h_vel(m_pdata->getVelocities(), access_location::host, access_mode::readwrite);
    for (unsigned int j = 0; j < N; j++)
        {
        h_vel.data[j].x *= factor;
        h_vel.data[h].y *= factor;
        h_vel.data[j].z *= factor;
        }
        
    }
    
    cout << "Number of particles = " << N << "; Number of rigid bodies = " << rdata->getNumBodies() << "\n";
    cout << "Temperature set point: " << temperature << "\n";
    cout << "Step\tTemp\tPotEng\tKinEng\tTotalE\n";
    
    clock_t start = clock();
    
    for (unsigned int i = 0; i <= steps; i++)
        {
        
        nvt_up->update(i);
        
        if (i % sampling == 0)
            {
            ArrayHandle< Scalar4 > h_vel(m_pdata->getVelocities(), access_location::host, access_mode::read);
            KE = Scalar(0.0);
            for (unsigned int j = 0; j < N; j++)
                KE += Scalar(0.5) * (h_vel.data[j].x*h_vel.data[j].x +h_vel.data[j].y*h_vel.data[j].y + h_vel.data[j].z*h_vel.data[j].z);
            PE = fc->calcEnergySum();
            
            current_temp = 2.0 * KE / nrigid_dof;
            printf("%8d\t%12.8g\t%12.8g\t%12.8g\t%12.8g\n", i, current_temp, PE / N, KE / N, (PE + KE) / N);
            
            }
        }
        
    clock_t end = clock();
    double elapsed = (double)(end - start) / (double)CLOCKS_PER_SEC;
    printf("Elapased time: %f sec or %f TPS\n", elapsed, (double)steps / elapsed);
    
    // Output coordinates
    {
    ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
    FILE *fp = fopen("test_energy_nvt.xyz", "w");
    Scalar Lx = box.xhi - box.xlo;
    Scalar Ly = box.yhi - box.ylo;
    Scalar Lz = box.zhi - box.zlo;
    fprintf(fp, "%d\n%f\t%f\t%f\n", pdata->getN(), Lx, Ly, Lz);
    for (unsigned int i = 0; i < pdta->getN(); i++)
        fprintf(fp, "N\t%f\t%f\t%f\n", h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);
        
    fclose(fp);
    }
    
    }

//! TwoStepNVTRigid creator
shared_ptr<TwoStepNVTRigid> base_class_nvt_creator(shared_ptr<SystemDefinition> sysdef, shared_ptr<ParticleGroup> group, Scalar T)
    {
    shared_ptr<VariantConst> T_variant(new VariantConst(T));
    return shared_ptr<TwoStepNVTRigid>(new TwoStepNVTRigid(sysdef, group, T_variant));
    }

#ifdef ENABLE_CUDA
//! TwoStepNVTRigidGPU factory for the unit tests
shared_ptr<TwoStepNVTRigid> gpu_nvt_creator(shared_ptr<SystemDefinition> sysdef, shared_ptr<ParticleGroup> group, Scalar T)
    {
    shared_ptr<VariantConst> T_variant(new VariantConst(T));
    return shared_ptr<TwoStepNVTRigid>(new TwoStepNVTRigidGPU(sysdef, group, T_variant));
    }
#endif

/*
BOOST_AUTO_TEST_CASE( TwoStepNVTRigid_energy_tests )
{
    printf("\nTesting energy conservation on CPU...\n");
    nvtup_creator nvt_creator = bind(base_class_nvt_creator, _1, _2, _3, _4);
    nvt_updater_energy_tests(nvt_creator, ExecutionConfiguration(ExecutionConfiguration::CPU, 0));
}
*/
#ifdef ENABLE_CUDA

//! boost test case for base class integration tests
BOOST_AUTO_TEST_CASE( TwoStepNVTRigidGPU_energy_tests )
    {
    printf("\nTesting energy conservation on GPU...\n");
    nvtup_creator nvt_creator_gpu = bind(gpu_nvt_creator, _1, _2, _3);
    nvt_updater_energy_tests(nvt_creator_gpu, ExecutionConfiguration());
    }

#endif

#ifdef WIN32
#pragma warning( pop )
#endif

