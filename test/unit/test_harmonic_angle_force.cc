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

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "HarmonicAngleForceCompute.h"
#include "ConstForceCompute.h"
#ifdef ENABLE_CUDA
#include "HarmonicAngleForceComputeGPU.h"
#endif

#include <stdio.h>

#include "Initializers.h"
#include "SnapshotSystemData.h"

using namespace std;
using namespace boost;

//! Name the boost unit test module
#define BOOST_TEST_MODULE AngleForceTests
#include "boost_utf_configure.h"

//! Typedef to make using the boost::function factory easier
typedef boost::function<shared_ptr<HarmonicAngleForceCompute>  (shared_ptr<SystemDefinition> sysdef)> angleforce_creator;

//! Perform some simple functionality tests of any BondForceCompute
void angle_force_basic_tests(angleforce_creator af_creator, boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    /////////////////////////////////////////////////////////
    // start with the simplest possible test: 3 particles in a huge box with only one bond type !!!! NO ANGLES
    shared_ptr<SystemDefinition> sysdef_3(new SystemDefinition(3, BoxDim(1000.0), 1, 1, 1, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_3 = sysdef_3->getParticleData();

    pdata_3->setPosition(0,make_scalar3(-1.23,2.0,0.1));
    pdata_3->setPosition(1,make_scalar3(1.0,1.0,1.0));
    pdata_3->setPosition(2,make_scalar3(1.0,0.0,0.5));
    
    //printf(" Particle 1: x = %f  y = %f  z = %f \n", arrays.x[0], arrays.y[0], arrays.z[0]);
    //printf(" Particle 2: x = %f  y = %f  z = %f \n", arrays.x[1], arrays.y[1], arrays.z[1]);
    //printf(" Particle 3: x = %f  y = %f  z = %f \n", arrays.x[2], arrays.y[2], arrays.z[2]);
    //printf("\n");
    
    // create the angle force compute to check
    shared_ptr<HarmonicAngleForceCompute> fc_3 = af_creator(sysdef_3);
    fc_3->setParams(0, Scalar(1.0), Scalar(0.785398)); // type=0, K=1.0,theta_0=pi/4=0.785398
    
    // compute the force and check the results
    fc_3->compute(0);
    {
    GPUArray<Scalar4>& force_array_1 =  fc_3->getForceArray();
    GPUArray<Scalar>& virial_array_1 =  fc_3->getVirialArray();
    unsigned int pitch = virial_array_1.getPitch();
    ArrayHandle<Scalar4> h_force_1(force_array_1,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_1(virial_array_1,access_location::host,access_mode::read);
    
    // check that the force is correct, it should be 0 since we haven't created any angles yet
    MY_BOOST_CHECK_SMALL(h_force_1.data[0].x, tol);
    MY_BOOST_CHECK_SMALL(h_force_1.data[0].y, tol);
    MY_BOOST_CHECK_SMALL(h_force_1.data[0].z, tol);
    MY_BOOST_CHECK_SMALL(h_force_1.data[0].w, tol);
    MY_BOOST_CHECK_SMALL(h_virial_1.data[0*pitch+0], tol);
    MY_BOOST_CHECK_SMALL(h_virial_1.data[1*pitch+0], tol);
    MY_BOOST_CHECK_SMALL(h_virial_1.data[2*pitch+0], tol);
    MY_BOOST_CHECK_SMALL(h_virial_1.data[3*pitch+0], tol);
    MY_BOOST_CHECK_SMALL(h_virial_1.data[4*pitch+0], tol);
    MY_BOOST_CHECK_SMALL(h_virial_1.data[5*pitch+0], tol);
    }

    // add an angle and check again
    sysdef_3->getAngleData()->addAngle(Angle(0,0,1,2)); // add type 0 bewtween angle formed by atom 0-1-2
    fc_3->compute(1);
    
    // this time there should be a force
    {
    GPUArray<Scalar4>& force_array_2 =  fc_3->getForceArray();
    GPUArray<Scalar>& virial_array_2 =  fc_3->getVirialArray();
    unsigned int pitch = virial_array_2.getPitch();
    ArrayHandle<Scalar4> h_force_2(force_array_2,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_2(virial_array_2,access_location::host,access_mode::read);
    MY_BOOST_CHECK_CLOSE(h_force_2.data[0].x, -0.061684, tol);
    MY_BOOST_CHECK_CLOSE(h_force_2.data[0].y, -0.313469, tol);
    MY_BOOST_CHECK_CLOSE(h_force_2.data[0].z, -0.195460, tol);
    MY_BOOST_CHECK_CLOSE(h_force_2.data[0].w, 0.158576, tol);
    MY_BOOST_CHECK_SMALL(h_virial_2.data[0*pitch+0]
                        +h_virial_2.data[3*pitch+0]
                        +h_virial_2.data[5*pitch+0], tol);
    
    //MY_BOOST_CHECK_SMALL(h_force_2.data[0].y, tol);
    //MY_BOOST_CHECK_CLOSE(h_force_2.data[0].z, 0.564651,tol);
    //MY_BOOST_CHECK_CLOSE(h_force_2.data[0].w, 0.298813, tol);
    //MY_BOOST_CHECK_CLOSE(h_virial_2.data[0], 0.0000001, tol);
    }
    /*
        printf("\n");
    */
    
    // rearrange the two particles in memory and see if they are properly updated
    {
    ArrayHandle<Scalar4> h_pos(pdata_3->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_tag(pdata_3->getTags(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_rtag(pdata_3->getRTags(), access_location::host, access_mode::readwrite);

    h_pos.data[1].x = Scalar(-1.23); // put atom a at (-1,0,0.1)
    h_pos.data[1].y = Scalar(2.0);
    h_pos.data[1].z = Scalar(0.1);
    
    h_pos.data[0].x = h_pos.data[0].y = h_pos.data[0].z = Scalar(1.0); // put atom b at (0,0,0)
    
    h_tag.data[0] = 1;
    h_tag.data[1] = 0;
    h_rtag.data[0] = 1;
    h_rtag.data[1] = 0;
    }
    
    // notify that we made the sort
    pdata_3->notifyParticleSort();
    // recompute at the same timestep, the forces should still be updated
    fc_3->compute(1);
    
    {
    GPUArray<Scalar4>& force_array_3 =  fc_3->getForceArray();
    GPUArray<Scalar>& virial_array_3 =  fc_3->getVirialArray();
    unsigned int pitch = virial_array_3.getPitch();
    ArrayHandle<Scalar4> h_force_3(force_array_3,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_3(virial_array_3,access_location::host,access_mode::read);
    
    MY_BOOST_CHECK_CLOSE(h_force_3.data[1].x, -0.061684, tol);
    MY_BOOST_CHECK_CLOSE(h_force_3.data[1].y, -0.3134695, tol);
    MY_BOOST_CHECK_CLOSE(h_force_3.data[1].z, -0.195460, tol);
    MY_BOOST_CHECK_CLOSE(h_force_3.data[1].w, 0.158576, tol);
    MY_BOOST_CHECK_SMALL(h_virial_3.data[0*pitch+1]
                        +h_virial_3.data[3*pitch+1]
                        +h_virial_3.data[5*pitch+1], tol);
    }

    ////////////////////////////////////////////////////////////////////
    // now, lets do a more thorough test and include boundary conditions
    // there are way too many permutations to test here, so I will simply
    // test +x, -x, +y, -y, +z, and -z independantly
    // build a 6 particle system with particles across each boundary
    // also test more than one type of bond
    unsigned int num_angles_to_test = 3;
    shared_ptr<SystemDefinition> sysdef_6(new SystemDefinition(6, BoxDim(20.0, 40.0, 60.0), 1, 1, num_angles_to_test, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_6 = sysdef_6->getParticleData();
    
    pdata_6->setPosition(0, make_scalar3(-9.6,0.0,0.0));
    pdata_6->setPosition(1, make_scalar3(9.6, 0.0,0.0));
    pdata_6->setPosition(2, make_scalar3(0.0,-19.6,0.0));
    pdata_6->setPosition(3, make_scalar3(0.0,19.6,0.0));
    pdata_6->setPosition(4, make_scalar3(0.0,0.0,-29.6));
    pdata_6->setPosition(5, make_scalar3(0.0,0.0,29.6));

    shared_ptr<HarmonicAngleForceCompute> fc_6 = af_creator(sysdef_6);
    fc_6->setParams(0, Scalar(1.0), Scalar(0.785398));
    fc_6->setParams(1, Scalar(2.0), Scalar(1.46));
    //fc_6->setParams(2, 1.5, 1.68);
    
    sysdef_6->getAngleData()->addAngle(Angle(0, 0,1,2));
    sysdef_6->getAngleData()->addAngle(Angle(1, 3,4,5));
    //pdata_6->getAngleData()->addAngle(Angle(2, 3,4,5));
    
    fc_6->compute(0);
    
    {
    // check that the forces are correctly computed
    GPUArray<Scalar4>& force_array_4 =  fc_6->getForceArray();
    GPUArray<Scalar>& virial_array_4 =  fc_6->getVirialArray();
    unsigned int pitch = virial_array_4.getPitch();
    ArrayHandle<Scalar4> h_force_4(force_array_4,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_4(virial_array_4,access_location::host,access_mode::read);
    
    
    //printf("\n");
    
    
    MY_BOOST_CHECK_SMALL(h_force_4.data[0].x, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[0].y, -1.55106342,tol);
    MY_BOOST_CHECK_SMALL(h_force_4.data[0].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[0].w, 0.256618, tol);
    MY_BOOST_CHECK_SMALL(h_virial_4.data[0*pitch+0]
                        +h_virial_4.data[3*pitch+0]
                        +h_virial_4.data[5*pitch+0], tol);
    
    MY_BOOST_CHECK_CLOSE(h_force_4.data[1].x, -0.0510595, loose_tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[1].y, 1.5760721,tol);
    MY_BOOST_CHECK_SMALL(h_force_4.data[1].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[1].w, 0.256618, tol);
    MY_BOOST_CHECK_SMALL(h_virial_4.data[0*pitch+1]
                        +h_virial_4.data[3*pitch+1]
                        +h_virial_4.data[5*pitch+1], tol);
    
    MY_BOOST_CHECK_CLOSE(h_force_4.data[2].x, 0.0510595,tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[2].y, -0.0250087, tol);
    MY_BOOST_CHECK_SMALL(h_force_4.data[2].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[2].w, 0.256618, tol);
    MY_BOOST_CHECK_SMALL(h_virial_4.data[0*pitch+2]
                        +h_virial_4.data[3*pitch+2]
                        +h_virial_4.data[5*pitch+2], tol);
    
    MY_BOOST_CHECK_SMALL(h_force_4.data[3].x, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[3].y, 0.05151510, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[3].z, -0.03411135,tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[3].w, 0.400928, tol);
    MY_BOOST_CHECK_SMALL(h_virial_4.data[0*pitch+3]
                        +h_virial_4.data[3*pitch+3]
                        +h_virial_4.data[5*pitch+3], tol);
    
    MY_BOOST_CHECK_SMALL(h_force_4.data[4].x, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[4].y, -2.79330492,tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[4].z, 0.034110874, loose_tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[4].w, 0.400928, tol);
    MY_BOOST_CHECK_SMALL(h_virial_4.data[0*pitch+4]
                        +h_virial_4.data[3*pitch+4]
                        +h_virial_4.data[5*pitch+4], tol);
    
    MY_BOOST_CHECK_SMALL(h_force_4.data[5].x, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[5].y, 2.74179,tol);
    MY_BOOST_CHECK_SMALL(h_force_4.data[5].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_4.data[5].w, 0.400928, tol);
    MY_BOOST_CHECK_SMALL(h_virial_4.data[0*pitch+5]
                        +h_virial_4.data[3*pitch+5]
                        +h_virial_4.data[5*pitch+5], tol);
    }
    
    
    // one more test: this one will test two things:
    // 1) That the forces are computed correctly even if the particles are rearranged in memory
    // and 2) That two forces can add to the same particle
    shared_ptr<SystemDefinition> sysdef_4(new SystemDefinition(4, BoxDim(100.0, 100.0, 100.0), 1, 1, 1, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_4 = sysdef_4->getParticleData();

    // make a square of particles
    {
    ArrayHandle<Scalar4> h_pos(pdata_4->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_tag(pdata_4->getTags(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_rtag(pdata_4->getRTags(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = 0.0; h_pos.data[0].y = 0.0; h_pos.data[0].z = 0.0;
    h_pos.data[1].x = 1.0; h_pos.data[1].y = 0; h_pos.data[1].z = 0.0;
    h_pos.data[2].x = 0; h_pos.data[2].y = 1.0; h_pos.data[2].z = 0.0;
    h_pos.data[3].x = 1.0; h_pos.data[3].y = 1.0; h_pos.data[3].z = 0.0;

    h_tag.data[0] = 2;
    h_tag.data[1] = 3;
    h_tag.data[2] = 0;
    h_tag.data[3] = 1;
    h_rtag.data[h_tag.data[0]] = 0;
    h_rtag.data[h_tag.data[1]] = 1;
    h_rtag.data[h_tag.data[2]] = 2;
    h_rtag.data[h_tag.data[3]] = 3;
    }
    
    // build the bond force compute and try it out
    shared_ptr<HarmonicAngleForceCompute> fc_4 = af_creator(sysdef_4);
    fc_4->setParams(0, 1.5, 1.75);
    // only add bonds on the left, top, and bottom of the square
    sysdef_4->getAngleData()->addAngle(Angle(0, 0,1,2));
    sysdef_4->getAngleData()->addAngle(Angle(0, 1,2,3));
    sysdef_4->getAngleData()->addAngle(Angle(0, 0,1,3));
    
    fc_4->compute(0);
    
    {
    GPUArray<Scalar4>& force_array_5 =  fc_4->getForceArray();
    GPUArray<Scalar>& virial_array_5 =  fc_4->getVirialArray();
    unsigned int pitch = virial_array_5.getPitch();
    ArrayHandle<Scalar4> h_force_5(force_array_5,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_5(virial_array_5,access_location::host,access_mode::read);
    
    
    // the first particles shoul only have a force pulling them right
    MY_BOOST_CHECK_CLOSE(h_force_5.data[0].x, 1.446903, tol);
    MY_BOOST_CHECK_SMALL(h_force_5.data[0].y, tol);
    MY_BOOST_CHECK_SMALL(h_force_5.data[0].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[0].w, 0.465228, tol);
    MY_BOOST_CHECK_SMALL(h_virial_5.data[0*pitch+0]
                        +h_virial_5.data[3*pitch+0]
                        +h_virial_5.data[5*pitch+0], tol);
    
    // and the bottom left particle should have a force pulling up and to the right
    MY_BOOST_CHECK_CLOSE(h_force_5.data[1].x, 0.2688054, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[1].y, -1.446902,tol);
    MY_BOOST_CHECK_SMALL(h_force_5.data[1].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[1].w, 0.240643, tol);
    MY_BOOST_CHECK_SMALL(h_virial_5.data[0*pitch+1]
                        +h_virial_5.data[3*pitch+1]
                        +h_virial_5.data[5*pitch+1], tol);
    
    // the bottom left particle should have a force pulling down
    MY_BOOST_CHECK_SMALL(h_force_5.data[2].x, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[2].y, 1.715708, tol);
    MY_BOOST_CHECK_SMALL(h_force_5.data[2].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[2].w, 0.240643,tol);
    MY_BOOST_CHECK_SMALL(h_virial_5.data[0*pitch+2]
                        +h_virial_5.data[3*pitch+2]
                        +h_virial_5.data[5*pitch+2], tol);
    
    // and the top left particle should have a force pulling up and to the left
    MY_BOOST_CHECK_CLOSE(h_force_5.data[3].x, -1.715708, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[3].y, -0.268805, tol);
    MY_BOOST_CHECK_SMALL(h_force_5.data[3].z, tol);
    MY_BOOST_CHECK_CLOSE(h_force_5.data[3].w, 0.473257, tol);
    MY_BOOST_CHECK_SMALL(h_virial_5.data[0*pitch+3]
                        +h_virial_5.data[3*pitch+3]
                        +h_virial_5.data[5*pitch+3], tol);
    }
    }





//! Compares the output of two HarmonicAngleForceComputes
void angle_force_comparison_tests(angleforce_creator af_creator1, angleforce_creator af_creator2, boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    const unsigned int N = 1000;
    
    // create a particle system to sum forces on
    // just randomly place particles. We don't really care how huge the bond forces get: this is just a unit test
    RandomInitializer rand_init(N, Scalar(0.2), Scalar(0.9), "A");
    boost::shared_ptr<SnapshotSystemData> snap =  rand_init.getSnapshot();
    snap->angle_data.type_mapping.push_back("A");
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snap, exec_conf));
    
    shared_ptr<HarmonicAngleForceCompute> fc1 = af_creator1(sysdef);
    shared_ptr<HarmonicAngleForceCompute> fc2 = af_creator2(sysdef);
    fc1->setParams(0, Scalar(1.0), Scalar(1.348));
    fc2->setParams(0, Scalar(1.0), Scalar(1.348));
    
    // add angles
    for (unsigned int i = 0; i < N-2; i++)
        {
        sysdef->getAngleData()->addAngle(Angle(0, i, i+1,i+2));
        }
        
    // compute the forces
    fc1->compute(0);
    fc2->compute(0);
    
    {
    GPUArray<Scalar4>& force_array_7 =  fc1->getForceArray();
    GPUArray<Scalar>& virial_array_7 =  fc1->getVirialArray();
    unsigned int pitch = virial_array_7.getPitch();
    ArrayHandle<Scalar4> h_force_7(force_array_7,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_7(virial_array_7,access_location::host,access_mode::read);
    GPUArray<Scalar4>& force_array_8 =  fc2->getForceArray();
    GPUArray<Scalar>& virial_array_8 =  fc2->getVirialArray();
    ArrayHandle<Scalar4> h_force_8(force_array_8,access_location::host,access_mode::read);
    ArrayHandle<Scalar> h_virial_8(virial_array_8,access_location::host,access_mode::read);
    
    // compare average deviation between the two computes
    double deltaf2 = 0.0;
    double deltape2 = 0.0;
    double deltav2[6];
    for (unsigned int i = 0; i < 6; i++)
        deltav2[i] = 0;

        
    for (unsigned int i = 0; i < N; i++)
        {
        deltaf2 += double(h_force_8.data[i].x - h_force_7.data[i].x) * double(h_force_8.data[i].x - h_force_7.data[i].x);
        deltaf2 += double(h_force_8.data[i].y - h_force_7.data[i].y) * double(h_force_8.data[i].y - h_force_7.data[i].y);
        deltaf2 += double(h_force_8.data[i].z - h_force_7.data[i].z) * double(h_force_8.data[i].z - h_force_7.data[i].z);
        deltape2 += double(h_force_8.data[i].w - h_force_7.data[i].w) * double(h_force_8.data[i].w - h_force_7.data[i].w);
        for (unsigned int j =0; j < 6; j++)
            deltav2[j] += double(h_virial_8.data[j*pitch+i] - h_virial_7.data[j*pitch+i]) * double(h_virial_8.data[j*pitch+i] - h_virial_7.data[j*pitch+i]);

        // also check that each individual calculation is somewhat close
        }
    deltaf2 /= double(N);
    deltape2 /= double(N);
    for (unsigned int i = 0; i < 6; i++)
        deltav2[i] /= double(N);
    BOOST_CHECK_SMALL(deltaf2, double(tol_small));
    BOOST_CHECK_SMALL(deltape2, double(tol_small));
    BOOST_CHECK_SMALL(deltav2[0], double(tol_small));
    BOOST_CHECK_SMALL(deltav2[1], double(tol_small));
    BOOST_CHECK_SMALL(deltav2[2], double(tol_small));
    BOOST_CHECK_SMALL(deltav2[3], double(tol_small));
    BOOST_CHECK_SMALL(deltav2[4], double(tol_small));
    BOOST_CHECK_SMALL(deltav2[5], double(tol_small));
    }
    }

//! HarmonicAngleForceCompute creator for angle_force_basic_tests()
shared_ptr<HarmonicAngleForceCompute> base_class_af_creator(shared_ptr<SystemDefinition> sysdef)
    {
    return shared_ptr<HarmonicAngleForceCompute>(new HarmonicAngleForceCompute(sysdef));
    }

#ifdef ENABLE_CUDA
//! AngleForceCompute creator for bond_force_basic_tests()
shared_ptr<HarmonicAngleForceCompute> gpu_af_creator(shared_ptr<SystemDefinition> sysdef)
    {
    return shared_ptr<HarmonicAngleForceCompute>(new HarmonicAngleForceComputeGPU(sysdef));
    }
#endif

//! boost test case for angle forces on the CPU
BOOST_AUTO_TEST_CASE( HarmonicAngleForceCompute_basic )
    {
    printf(" IN BOOST_AUTO_TEST_CASE: CPU \n");
    //cout << " IN BOOST_AUTO_TEST_CASE: CPU \n";
    angleforce_creator af_creator = bind(base_class_af_creator, _1);
    angle_force_basic_tests(af_creator, boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }

#ifdef ENABLE_CUDA
//! boost test case for angle forces on the GPU
BOOST_AUTO_TEST_CASE( HarmonicAngleForceComputeGPU_basic )
    {
    printf(" IN BOOST_AUTO_TEST_CASE: GPU \n");
    cout << " IN BOOST_AUTO_TEST_CASE: GPU \n";
    angleforce_creator af_creator = bind(gpu_af_creator, _1);
    angle_force_basic_tests(af_creator, boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }


//! boost test case for comparing bond GPU and CPU BondForceComputes
BOOST_AUTO_TEST_CASE( HarmonicAngleForceComputeGPU_compare )
    {
    angleforce_creator af_creator_gpu = bind(gpu_af_creator, _1);
    angleforce_creator af_creator = bind(base_class_af_creator, _1);
    angle_force_comparison_tests(af_creator, af_creator_gpu, boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }


#endif

