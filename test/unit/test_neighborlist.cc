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
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "NeighborList.h"
#include "NeighborListBinned.h"
#include "Initializers.h"

#ifdef ENABLE_CUDA
#include "NeighborListGPU.h"
#include "NeighborListGPUBinned.h"
#endif

using namespace std;
using namespace boost;

//! Define the name of the boost test module
#define BOOST_TEST_MODULE NeighborListTest
#include "boost_utf_configure.h"

//! Performs basic functionality tests on a neighbor list
template <class NL>
void neighborlist_basic_tests(boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    /////////////////////////////////////////////////////////
    // start with the simplest possible test: 2 particles in a huge box
    shared_ptr<SystemDefinition> sysdef_2(new SystemDefinition(2, BoxDim(25.0), 1, 0, 0, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_2 = sysdef_2->getParticleData();

    {
    ArrayHandle<Scalar4> h_pos(pdata_2->getPositions(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = h_pos.data[0].y = h_pos.data[0].z = 0.0;
    h_pos.data[1].x = h_pos.data[1].y = h_pos.data[1].z = 3.25;
    }
    
    // test construction of the neighborlist
    shared_ptr<NeighborList> nlist_2(new NL(sysdef_2, 3.0, 0.25));
    nlist_2->compute(1);
    
    // with the given radius, there should be no neighbors: check that
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[0], 0);
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[1], 0);
        }
    
    // adjust the radius to include the particles and see if we get some now
    nlist_2->setRCut(5.5, 0.5);
    nlist_2->compute(2);
    // some neighbor lists default to full because they don't support half: ignore them
    if (nlist_2->getStorageMode() == NeighborList::half)
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_2->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_2->getNListIndexer();
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        // since this is a half list, only 0 stores 1 as a neighbor
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[1], 0);
        }
        
    // change to full mode to check that
    nlist_2->setStorageMode(NeighborList::full);
    nlist_2->compute(3);
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_2->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_2->getNListIndexer();
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);
        }
    
    
    ////////////////////////////////////////////////////////////////////
    // now, lets do a more thorough test and include boundary conditions
    // there are way too many permutations to test here, so I will simply
    // test +x, -x, +y, -y, +z, and -z independantly
    // build a 6 particle system with particles across each boundary
    
    shared_ptr<SystemDefinition> sysdef_6(new SystemDefinition(6, BoxDim(20.0, 40.0, 60.0), 1, 0, 0, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_6 = sysdef_6->getParticleData();
    
    {
    ArrayHandle<Scalar4> h_pos(pdata_6->getPositions(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = Scalar(-9.6); h_pos.data[0].y = 0; h_pos.data[0].z = 0.0;
    h_pos.data[1].x =  Scalar(9.6); h_pos.data[1].y = 0; h_pos.data[1].z = 0.0;
    h_pos.data[2].x = 0; h_pos.data[2].y = Scalar(-19.6); h_pos.data[2].z = 0.0;
    h_pos.data[3].x = 0; h_pos.data[3].y = Scalar(19.6); h_pos.data[3].z = 0.0;
    h_pos.data[4].x = 0; h_pos.data[4].y = 0; h_pos.data[4].z = Scalar(-29.6);
    h_pos.data[5].x = 0; h_pos.data[5].y = 0; h_pos.data[5].z =  Scalar(29.6);
    }
    
    shared_ptr<NeighborList> nlist_6(new NL(sysdef_6, 3.0, 0.25));
    nlist_6->setStorageMode(NeighborList::full);
    nlist_6->compute(0);
    // verify the neighbor list
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_6->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_6->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_6->getNListIndexer();
        
        BOOST_REQUIRE(nli.getW() >= 6);
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 3);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[3], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,0)], 2);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[4], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,0)], 5);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[5], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,0)], 4);
        }
    
    // swap the order of the particles around to look for subtle directional bugs
    {
    ArrayHandle<Scalar4> h_pos(pdata_6->getPositions(), access_location::host, access_mode::readwrite);

    h_pos.data[1].x = Scalar(-9.6); h_pos.data[1].y = 0; h_pos.data[1].z = 0.0;
    h_pos.data[0].x =  Scalar(9.6); h_pos.data[0].y = 0; h_pos.data[0].z = 0.0;
    h_pos.data[3].x = 0; h_pos.data[3].y = Scalar(-19.6); h_pos.data[3].z = 0.0;
    h_pos.data[2].x = 0; h_pos.data[2].y = Scalar(19.6); h_pos.data[2].z = 0.0;
    h_pos.data[5].x = 0; h_pos.data[5].y = 0; h_pos.data[5].z = Scalar(-29.6);
    h_pos.data[4].x = 0; h_pos.data[4].y = 0; h_pos.data[4].z =  Scalar(29.6);
    }
    
    // verify the neighbor list
    nlist_6->compute(1);
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_6->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_6->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_6->getNListIndexer();
        
        BOOST_REQUIRE(nli.getW() >= 6);
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 3);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[3], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,0)], 2);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[4], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,0)], 5);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[5], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,0)], 4);
        }
    
    // one last test, we should check that more than one neighbor can be generated
    {
    ArrayHandle<Scalar4> h_pos(pdata_6->getPositions(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = 0; h_pos.data[0].y = 0; h_pos.data[0].z = 0.0;
    h_pos.data[1].x = 0; h_pos.data[1].y = 0; h_pos.data[1].z = 0.0;
    h_pos.data[2].x = 0; h_pos.data[2].y = Scalar(-19.6); h_pos.data[2].z = 0.0;
    h_pos.data[3].x = 0; h_pos.data[3].y = Scalar(19.6); h_pos.data[3].z = 0.0;
    h_pos.data[4].x = 0; h_pos.data[4].y = 0; h_pos.data[4].z = 0;
    h_pos.data[5].x = 0; h_pos.data[5].y = 0; h_pos.data[5].z =  0;
    }
    
    nlist_6->compute(20);
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_6->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_6->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_6->getNListIndexer();
        
        BOOST_REQUIRE(nli.getW() >= 6);
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,1)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,2)], 5);
        }
    }

//! Tests the ability of the neighbor list to exclude particle pairs
template <class NL>
void neighborlist_exclusion_tests(boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    shared_ptr<SystemDefinition> sysdef_6(new SystemDefinition(6, BoxDim(20.0, 40.0, 60.0), 1, 0, 0, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_6 = sysdef_6->getParticleData();
    
    // lets make this test simple: put all 6 particles on top of each other and
    // see if the exclusion code can ignore 4 of the particles
    {
    ArrayHandle<Scalar4> h_pos(pdata_6->getPositions(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = 0; h_pos.data[0].y = 0; h_pos.data[0].z = 0.0;
    h_pos.data[1].x = 0; h_pos.data[1].y = 0; h_pos.data[1].z = 0.0;
    h_pos.data[2].x = 0; h_pos.data[2].y = 0; h_pos.data[2].z = 0.0;
    h_pos.data[3].x = 0; h_pos.data[3].y = 0; h_pos.data[3].z = 0.0;
    h_pos.data[4].x = 0; h_pos.data[4].y = 0; h_pos.data[4].z = 0;
    h_pos.data[5].x = 0; h_pos.data[5].y = 0; h_pos.data[5].z =  0;
    }

    shared_ptr<NeighborList> nlist_6(new NL(sysdef_6, 3.0, 0.25));
    nlist_6->setStorageMode(NeighborList::full);
    nlist_6->addExclusion(0,1);
    nlist_6->addExclusion(0,2);
    nlist_6->addExclusion(0,3);
    nlist_6->addExclusion(0,4);
    
    nlist_6->compute(0);
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_6->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_6->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_6->getNListIndexer();
        
        BOOST_REQUIRE(nli.getW() >= 6);
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 5);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,1)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,2)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,3)], 5);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,1)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,2)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,3)], 5);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[3], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,2)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,3)], 5);
    
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[4], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,2)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,3)], 5);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[5], 5);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,1)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,2)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,3)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,4)], 4);
        }
    }

//! Tests the ability of the neighbor list to exclude particles from the same body
template <class NL>
void neighborlist_body_filter_tests(boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    shared_ptr<SystemDefinition> sysdef_6(new SystemDefinition(6, BoxDim(20.0, 40.0, 60.0), 1, 0, 0, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_6 = sysdef_6->getParticleData();
    
    // lets make this test simple: put all 6 particles on top of each other and
    // see if the exclusion code can ignore 4 of the particles
    {
    ArrayHandle<Scalar4> h_pos(pdata_6->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_body(pdata_6->getBodies(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = 0; h_pos.data[0].y = 0; h_pos.data[0].z = 0; h_body.data[0] = NO_BODY;
    h_pos.data[1].x = 0; h_pos.data[1].y = 0; h_pos.data[1].z = 0; h_body.data[1] = 0;
    h_pos.data[2].x = 0; h_pos.data[2].y = 0; h_pos.data[2].z = 0; h_body.data[2] = 1;
    h_pos.data[3].x = 0; h_pos.data[3].y = 0; h_pos.data[3].z = 0; h_body.data[3] = 0;
    h_pos.data[4].x = 0; h_pos.data[4].y = 0; h_pos.data[4].z = 0; h_body.data[4] = 1;
    h_pos.data[5].x = 0; h_pos.data[5].y = 0; h_pos.data[5].z = 0; h_body.data[5] = NO_BODY;
    }
    
    // this test uses rigid bodies, initialize them
    sysdef_6->getRigidData()->initializeData();

    shared_ptr<NeighborList> nlist_6(new NL(sysdef_6, 3.0, 0.25));
    nlist_6->setFilterBody(true);
    nlist_6->setStorageMode(NeighborList::full);
    
    nlist_6->compute(0);
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_6->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_6->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_6->getNListIndexer();
        
        BOOST_REQUIRE(nli.getW() >= 6);
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 5);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,2)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,3)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,4)], 5);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,2)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,3)], 5);
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,1)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,2)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,3)], 5);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[3], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,2)], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(3,3)], 5);
    
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[4], 4);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,1)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,2)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(4,3)], 5);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[5], 5);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,1)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,2)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,3)], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(5,4)], 4);
        }
    }

//! Tests the ability of the neighbor list to filter by diameter
template <class NL>
void neighborlist_diameter_filter_tests(boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    /////////////////////////////////////////////////////////
    // start with the simplest possible test: 3 particles in a huge box
    shared_ptr<SystemDefinition> sysdef_3(new SystemDefinition(4, BoxDim(25.0), 1, 0, 0, 0, 0, exec_conf));
    shared_ptr<ParticleData> pdata_3 = sysdef_3->getParticleData();

    {
    ArrayHandle<Scalar4> h_pos(pdata_3->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_diameter(pdata_3->getDiameters(), access_location::host, access_mode::readwrite);

    h_pos.data[0].x = 0; h_pos.data[0].y = 0; h_pos.data[0].z = 0.0; h_diameter.data[0] = 3.0;
    h_pos.data[2].x = 0; h_pos.data[2].y = 0; h_pos.data[2].z = 2.5; h_diameter.data[2] = 2.0;
    h_pos.data[1].x = 0; h_pos.data[1].y = 0; h_pos.data[1].z = -3.0; h_diameter.data[1] = 1.0;
    h_pos.data[3].x = 0; h_pos.data[3].y = 2.51; h_pos.data[3].z = 0; h_diameter.data[3] = 0;
    }
    
    // test construction of the neighborlist
    shared_ptr<NeighborList> nlist_2(new NL(sysdef_3, 1.5, 0.5));
    nlist_2->compute(1);
    nlist_2->setStorageMode(NeighborList::full);

    // with the given settings, there should be no neighbors: check that
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[0], 0);
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[1], 0);
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[2], 0);
        }
    
    // set a test maximum diameter of 2.0
    nlist_2->setMaximumDiameter(2.0);
    nlist_2->compute(2);
    
    // 0 and 1 should be neighbors now, as well as 0 and 2
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_2->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_2->getNListIndexer();
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,2)], 3);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 0);
        }
    
    // bump it up to 3.0
    nlist_2->setMaximumDiameter(3.0);
    nlist_2->compute(3);
    
    // should be the same as above
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_2->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_2->getNListIndexer();
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 3);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,1)], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,2)], 3);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,1)], 3);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 0);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,1)], 3);
        }
    
    // enable diameter filtering and verify the result is still correct
    nlist_2->setFilterDiameter(true);
    nlist_2->compute(4);

    // the particle 0 should now be neighbors with 1 and 2
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist_2->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist_2->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist_2->getNListIndexer();
        
        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[0], 2);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,0)], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(0,1)], 2);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[1], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(1,0)], 0);

        BOOST_REQUIRE_EQUAL_UINT(h_n_neigh.data[2], 1);
        BOOST_CHECK_EQUAL_UINT(h_nlist.data[nli(2,0)], 0);
        }
    }

//! Test two implementations of NeighborList and verify that the output is identical
template <class NLA, class NLB>
void neighborlist_comparison_test(boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    // construct the particle system
    RandomInitializer init(1000, Scalar(0.016778), Scalar(0.9), "A");
    boost::shared_ptr<SnapshotSystemData> snap = init.getSnapshot();
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snap, exec_conf));
    shared_ptr<ParticleData> pdata = sysdef->getParticleData();
    
    shared_ptr<NeighborList> nlist1(new NLA(sysdef, Scalar(3.0), Scalar(0.4)));
    nlist1->setStorageMode(NeighborList::full);
    
    shared_ptr<NeighborList> nlist2(new NLB(sysdef, Scalar(3.0), Scalar(0.4)));
    nlist2->setStorageMode(NeighborList::full);
    
    // setup some exclusions: try to fill out all four exclusions for each particle
    for (unsigned int i=0; i < pdata->getN()-2; i++)
        {
        nlist1->addExclusion(i,i+1);
        nlist1->addExclusion(i,i+2);
        
        nlist2->addExclusion(i,i+1);
        nlist2->addExclusion(i,i+2);
        }
        
    // compute each of the lists
    nlist1->compute(0);
    nlist2->compute(0);
    
    // verify that both new ones match the basic
    ArrayHandle<unsigned int> h_n_neigh1(nlist1->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist1(nlist1->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_n_neigh2(nlist2->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist2(nlist2->getNListArray(), access_location::host, access_mode::read);
    Index2D nli = nlist1->getNListIndexer();
    
    // temporary vectors for holding the lists: they will be sorted for comparison
    std::vector<unsigned int> tmp_list1;
    std::vector<unsigned int> tmp_list2;
    
    // check to make sure that every neighbor matches
    for (unsigned int i = 0; i < pdata->getN(); i++)
        {
        BOOST_REQUIRE_EQUAL(h_n_neigh1.data[i], h_n_neigh2.data[i]);
        
        tmp_list1.resize(h_n_neigh1.data[i]);
        tmp_list2.resize(h_n_neigh1.data[i]);
        
        for (unsigned int j = 0; j < h_n_neigh1.data[i]; j++)
            {
            tmp_list1[j] = h_nlist1.data[nli(i,j)];
            tmp_list2[j] = h_nlist2.data[nli(i,j)];
            }
        
        sort(tmp_list1.begin(), tmp_list1.end());
        sort(tmp_list2.begin(), tmp_list2.end());
        
        for (unsigned int j = 0; j < tmp_list1.size(); j++)
            {
            BOOST_CHECK_EQUAL(tmp_list1[j], tmp_list2[j]);
            }
        }
    }

//! Test that a NeighborList can successfully exclude a ridiculously large number of particles
template <class NL>
void neighborlist_large_ex_tests(boost::shared_ptr<ExecutionConfiguration> exec_conf)
    {
    // construct the particle system
    RandomInitializer init(1000, Scalar(0.016778), Scalar(0.9), "A");
    boost::shared_ptr<SnapshotSystemData> snap = init.getSnapshot();
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snap, exec_conf));
    shared_ptr<ParticleData> pdata = sysdef->getParticleData();
    
    shared_ptr<NeighborList> nlist(new NL(sysdef, Scalar(8.0), Scalar(0.4)));
    nlist->setStorageMode(NeighborList::full);
    
    // add every single neighbor as an exclusion
    nlist->compute(0);
        {
        ArrayHandle<unsigned int> h_n_neigh(nlist->getNNeighArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_nlist(nlist->getNListArray(), access_location::host, access_mode::read);
        Index2D nli = nlist->getNListIndexer();
        
        for (unsigned int i = 0; i < pdata->getN(); i++)
            {
            for (unsigned int neigh = 0; neigh < h_n_neigh.data[i]; neigh++)
                {
                unsigned int j = h_nlist.data[nli(i, neigh)];
                nlist->addExclusion(i,j);
                }
            }
        }
    
    // compute the nlist again
    nlist->compute(0);
    
    // verify that there are now 0 neighbors for each particle
    ArrayHandle<unsigned int> h_n_neigh(nlist->getNNeighArray(), access_location::host, access_mode::read);
    
    // check to make sure that every neighbor matches
    for (unsigned int i = 0; i < pdata->getN(); i++)
        {
        BOOST_CHECK_EQUAL_UINT(h_n_neigh.data[i], 0);
        }
    }

//! basic test case for base class
BOOST_AUTO_TEST_CASE( NeighborList_basic )
    {
    neighborlist_basic_tests<NeighborList>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! exclusion test case for base class
BOOST_AUTO_TEST_CASE( NeighborList_exclusion )
    {
    neighborlist_exclusion_tests<NeighborList>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! large exclusion test case for base class
BOOST_AUTO_TEST_CASE( NeighborList_large_ex )
    {
    neighborlist_large_ex_tests<NeighborList>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! body filter test case for base class
BOOST_AUTO_TEST_CASE( NeighborList_body_filter)
    {
    neighborlist_body_filter_tests<NeighborList>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! diameter filter test case for base class
BOOST_AUTO_TEST_CASE( NeighborList_diameter_filter )
    {
    neighborlist_diameter_filter_tests<NeighborList>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }

//! basic test case for binned class
BOOST_AUTO_TEST_CASE( NeighborListBinned_basic )
    {
    neighborlist_basic_tests<NeighborListBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! exclusion test case for binned class
BOOST_AUTO_TEST_CASE( NeighborListBinned_exclusion )
    {
    neighborlist_exclusion_tests<NeighborListBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! large exclusion test case for binned class
BOOST_AUTO_TEST_CASE( NeighborListBinned_large_ex )
    {
    neighborlist_large_ex_tests<NeighborListBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! body filter test case for binned class
BOOST_AUTO_TEST_CASE( NeighborListBinned_body_filter)
    {
    neighborlist_body_filter_tests<NeighborListBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! diameter filter test case for binned class
BOOST_AUTO_TEST_CASE( NeighborListBinned_diameter_filter )
    {
    neighborlist_diameter_filter_tests<NeighborListBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }
//! comparison test case for binned class
BOOST_AUTO_TEST_CASE( NeighborListBinned_comparison )
    {
    neighborlist_comparison_test<NeighborList, NeighborListBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::CPU)));
    }

#ifdef ENABLE_CUDA

//! basic test case for GPU class
BOOST_AUTO_TEST_CASE( NeighborListGPU_basic )
    {
    neighborlist_basic_tests<NeighborListGPU>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! exclusion test case for GPU class
BOOST_AUTO_TEST_CASE( NeighborListGPU_exclusion )
    {
    neighborlist_exclusion_tests<NeighborListGPU>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! large exclusion test case for GPU class
BOOST_AUTO_TEST_CASE( NeighborListGPU_large_ex )
    {
    neighborlist_large_ex_tests<NeighborListGPU>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }

// disabled as NeighborListGPU doesn't support these filters yet
/*
//! body filter test case for GPU class
BOOST_AUTO_TEST_CASE( NeighborListGPU_body_filter)
    {
    neighborlist_body_filter_tests<NeighborListGPU>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! diameter filter test case for GPU class
BOOST_AUTO_TEST_CASE( NeighborListGPU_diameter_filter )
    {
    neighborlist_diameter_filter_tests<NeighborListGPU>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
*/
//! comparison test case for GPU class
BOOST_AUTO_TEST_CASE( NeighborListGPU_comparison )
    {
    neighborlist_comparison_test<NeighborList, NeighborListGPU>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }

//! basic test case for GPUBinned class
BOOST_AUTO_TEST_CASE( NeighborListGPUBinned_basic )
    {
    neighborlist_basic_tests<NeighborListGPUBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! exclusion test case for GPUBinned class
BOOST_AUTO_TEST_CASE( NeighborListGPUBinned_exclusion )
    {
    neighborlist_exclusion_tests<NeighborListGPUBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! large exclusion test case for GPUBinned class
BOOST_AUTO_TEST_CASE( NeighborListGPUBinned_large_ex )
    {
    neighborlist_large_ex_tests<NeighborListGPUBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! body filter test case for GPUBinned class
BOOST_AUTO_TEST_CASE( NeighborListGPUBinned_body_filter)
    {
    neighborlist_body_filter_tests<NeighborListGPUBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! diameter filter test case for GPUBinned class
BOOST_AUTO_TEST_CASE( NeighborListGPUBinned_diameter_filter )
    {
    neighborlist_diameter_filter_tests<NeighborListGPUBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }
//! comparison test case for GPUBinned class
BOOST_AUTO_TEST_CASE( NeighborListGPUBinned_comparison )
    {
    neighborlist_comparison_test<NeighborList, NeighborListGPUBinned>(boost::shared_ptr<ExecutionConfiguration>(new ExecutionConfiguration(ExecutionConfiguration::GPU)));
    }

#endif

#ifdef WIN32
#pragma warning( pop )
#endif




