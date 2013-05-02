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

#include "ExampleUpdater.cuh"

/*! \file ExampleUpdater.cu
    \brief CUDA kernels for ExampleUpdater
*/

// First, the kernel code for zeroing the velocities on the GPU
//! Kernel that zeroes velocities on the GPU
/*! \param d_vel Velocity-mass array from the ParticleData
    \param N Number of particles
    
    This kernel executes one thread per particle and zeros the velocity of each. It can be run with any 1D block size
    as long as block_size * num_blocks is >= the number of particles.
*/
extern "C" __global__ 
void gpu_zero_velocities_kernel(float4 *d_vel, unsigned int N)
    {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < N)
        {
        // vel.w is the mass, don't want to modify that
        float4 vel = d_vel[idx];
        vel.x = vel.y = vel.z = 0.0f;
        d_vel[idx] = vel;
        }
    }

/*! \param d_vel Velocity-mass array from the ParticleData
    \param N Number of particles
    This is just a driver for gpu_zero_velocities_kernel(), see it for the details
*/
cudaError_t gpu_zero_velocities(float4 *d_vel, unsigned int N)
    {
    // setup the grid to run the kernel
    int block_size = 256;
    dim3 grid( (int)ceil((double)N / (double)block_size), 1, 1);
    dim3 threads(block_size, 1, 1);
    
    // run the kernel
    gpu_zero_velocities_kernel<<< grid, threads >>>(d_vel, N);
    
    // this method always succeds. If you had a cuda* call in this driver, you could return its error code if not
    // cudaSuccess
    return cudaSuccess;
    }

