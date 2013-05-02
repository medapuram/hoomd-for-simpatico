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

#include "Enforce2DUpdaterGPU.cuh"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

#include <stdio.h>

/*! \file Enforce2DUpdaterGPU.cu
    \brief Defines GPU kernel code for constraining systems to a 2D plane on 
    the GPU. Used by Enforce2DUpdaterGPU.
*/

//! Constrains partcles to the xy plane on the GPU
/*! \param N number of particles in system
    \param d_vel Particle velocities to constrain to xy plane
    \param d_accel Particle accelerations to constrain to xy plane
*/
extern "C" __global__ 
void gpu_enforce2d_kernel(const unsigned int N,
                          Scalar4 *d_vel,
                          Scalar3 *d_accel)
    {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < N)
        {        
        // read the particle's velocity and acceleration (MEM TRANSFER: 32 bytes)
        Scalar4 vel = d_vel[idx];
        Scalar3 accel = d_accel[idx];
                
        // zero the z-velocity and z-acceleration(FLOPS: ?)
        vel.z = 0.0f;
        accel.z = 0.0f;
                
        // write out the results (MEM_TRANSFER: 32 bytes)
        d_vel[idx] = vel;
        d_accel[idx] = accel;
        }
    }

/*! \param N number of particles in system
    \param d_vel Particle velocities to constrain to xy plane
    \param d_accel Particle accelerations to constrain to xy plane
*/
cudaError_t gpu_enforce2d(const unsigned int N,
                          Scalar4 *d_vel,
                          Scalar3 *d_accel)
    {
    // setup the grid to run the kernel
    int block_size = 256;
    dim3 grid( (N/block_size) + 1, 1, 1);
    dim3 threads(block_size, 1, 1);
            
    // run the kernel
    gpu_enforce2d_kernel<<< grid, threads >>>(N, d_vel, d_accel);
    
    return cudaSuccess;
    }

