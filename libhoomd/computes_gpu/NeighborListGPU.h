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

#include "NeighborList.h"
#include "GPUFlags.h"

/*! \file NeighborListGPU.h
    \brief Declares the NeighborListGPU class
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#ifndef __NEIGHBORLISTGPU_H__
#define __NEIGHBORLISTGPU_H__

//! Neighbor list build on the GPU
/*! Implements the O(N^2) neighbor list build on the GPU. Also implements common functions (like distance check)
    on the GPU for use by other GPU nlist classes derived from NeighborListGPU.
    
    GPU kernel methods are defined in NeighborListGPU.cuh and defined in NeighborListGPU.cu.
    
    \ingroup computes
*/
class NeighborListGPU : public NeighborList
    {
    public:
        //! Constructs the compute
        NeighborListGPU(boost::shared_ptr<SystemDefinition> sysdef, Scalar r_cut, Scalar r_buff)
            : NeighborList(sysdef, r_cut, r_buff)
            {
            GPUFlags<uint2> flags(exec_conf);
            m_flags.swap(flags);
            m_flags.resetFlags(make_uint2(0,0));
            // default to full mode
            m_storage_mode = full;
            m_block_size_filter = 192;
            m_checkn = 1;
            }
        
        //! Destructor
        virtual ~NeighborListGPU()
            {
            }

        //! Set block size for filter kernel
        void setBlockSizeFilter(unsigned int block_size)
            {
            m_block_size_filter = block_size;
            }
        
        //! Benchmark the filter kernel
        double benchmarkFilter(unsigned int num_iters);

        //! Update the exclusion list on the GPU
        virtual void updateExListIdx();

    protected:
        GPUFlags<uint2> m_flags;     //!< Storage for device flags on the GPU

        //! Builds the neighbor list
        virtual void buildNlist(unsigned int timestep);

        //! Perform the nlist distance check on the GPU
        virtual bool distanceCheck();
        
        //! GPU nlists set their last updated pos in the compute kernel, this call only resets the last box length
        virtual void setLastUpdatedPos()
            {
            m_last_L = m_pdata->getGlobalBox().getL(); 
            }
        
        //! Filter the neighbor list of excluded particles
        virtual void filterNlist();

    private:
        unsigned int m_block_size_filter;   //!< Block size for the filter kernel
        unsigned int m_checkn;              //!< Internal counter to assign when checking if the nlist needs an update
    };

//! Exports NeighborListGPU to python
void export_NeighborListGPU();

#endif

