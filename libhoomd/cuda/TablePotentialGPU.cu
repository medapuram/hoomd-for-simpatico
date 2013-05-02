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

#include "TablePotentialGPU.cuh"

#include "Index1D.h"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

/*! \file TablePotentialGPU.cu
    \brief Defines GPU kernel code for calculating the table pair forces. Used by TablePotentialGPU.
*/

//! Texture for reading particle positions
texture<float4, 1, cudaReadModeElementType> pdata_pos_tex;

//! Texture for reading table values
texture<float2, 1, cudaReadModeElementType> tables_tex;

/*!  This kernel is called to calculate the table pair forces on all N particles

    \param d_force Device memory to write computed forces
    \param d_virial Device memory to write computed virials
    \param virial_pitch Pitch of 2D virial array
    \param N number of particles in system
    \param d_pos device array of particle positions
    \param box Box dimensions used to implement periodic boundary conditions
    \param d_n_neigh Device memory array listing the number of neighbors for each particle
    \param d_nlist Device memory array containing the neighbor list contents
    \param nli Indexer for indexing \a d_nlist
    \param d_params Parameters for each table associated with a type pair
    \param ntypes Number of particle types in the system
    \param table_width Number of points in each table

    See TablePotential for information on the memory layout.

    \b Details:
    * Table entries are read from tables_tex. Note that currently this is bound to a 1D memory region. Performance tests
      at a later date may result in this changing.
*/
__global__ void gpu_compute_table_forces_kernel(float4* d_force,
                                                float* d_virial,
                                                const unsigned virial_pitch,
                                                const unsigned int N,
                                                const Scalar4 *d_pos,
                                                const BoxDim box,
                                                const unsigned int *d_n_neigh,
                                                const unsigned int *d_nlist,
                                                const Index2D nli,
                                                const float4 *d_params,
                                                const unsigned int ntypes,
                                                const unsigned int table_width)
    {
    // index calculation helpers
    Index2DUpperTriangular table_index(ntypes);
    Index2D table_value(table_width);

    // read in params for easy and fast access in the kernel
    extern __shared__ float4 s_params[];
    for (unsigned int cur_offset = 0; cur_offset < table_index.getNumElements(); cur_offset += blockDim.x)
        {
        if (cur_offset + threadIdx.x < table_index.getNumElements())
            s_params[cur_offset + threadIdx.x] = d_params[cur_offset + threadIdx.x];
        }
    __syncthreads();

    // start by identifying which particle we are to handle
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= N)
        return;

    // load in the length of the list
    unsigned int n_neigh = d_n_neigh[idx];

    // read in the position of our particle. Texture reads of float4's are faster than global reads on compute 1.0 hardware
    Scalar4 postype = tex1Dfetch(pdata_pos_tex, idx);
    Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);
    unsigned int typei = __float_as_int(postype.w);

    // initialize the force to 0
    float4 force = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    float virialxx = 0.0f;
    float virialxy = 0.0f;
    float virialxz = 0.0f;
    float virialyy = 0.0f;
    float virialyz = 0.0f;
    float virialzz = 0.0f;

    // prefetch neighbor index
    unsigned int cur_neigh = 0;
    unsigned int next_neigh = d_nlist[nli(idx, 0)];

    // loop over neighbors
    // on pre Fermi hardware, there is a bug that causes rare and random ULFs when simply looping over n_neigh
    // the workaround (activated via the template paramter) is to loop over nlist.height and put an if (i < n_neigh)
    // inside the loop
    #if (__CUDA_ARCH__ < 200)
    for (int neigh_idx = 0; neigh_idx < nli.getH(); neigh_idx++)
    #else
    for (int neigh_idx = 0; neigh_idx < n_neigh; neigh_idx++)
    #endif
        {
        #if (__CUDA_ARCH__ < 200)
        if (neigh_idx < n_neigh)
        #endif
            {
            // read the current neighbor index
            // prefetch the next value and set the current one
            cur_neigh = next_neigh;
            next_neigh = d_nlist[nli(idx, (neigh_idx+1))];

            // get the neighbor's position
            float4 neigh_postype = tex1Dfetch(pdata_pos_tex, cur_neigh);
            float3 neigh_pos = make_float3(neigh_postype.x, neigh_postype.y, neigh_postype.z);

            // calculate dr (with periodic boundary conditions)
            float3 dx = pos - neigh_pos;

            // apply periodic boundary conditions
            dx = box.minImage(dx);

            // access needed parameters
            unsigned int typej = __float_as_int(neigh_postype.w);
            unsigned int cur_table_index = table_index(typei, typej);
            float4 params = s_params[cur_table_index];
            float rmin = params.x;
            float rmax = params.y;
            float delta_r = params.z;

            // calculate r
            float rsq = dot(dx, dx);
            float r = sqrtf(rsq);

            if (r < rmax && r >= rmin)
                {
                // precomputed term
                float value_f = (r - rmin) / delta_r;

                // compute index into the table and read in values
                unsigned int value_i = floor(value_f);
                float2 VF0 = tex1Dfetch(tables_tex, table_value(value_i, cur_table_index));
                float2 VF1 = tex1Dfetch(tables_tex, table_value(value_i+1, cur_table_index));
                // unpack the data
                float V0 = VF0.x;
                float V1 = VF1.x;
                float F0 = VF0.y;
                float F1 = VF1.y;

                // compute the linear interpolation coefficient
                float f = value_f - float(value_i);

                // interpolate to get V and F;
                float V = V0 + f * (V1 - V0);
                float F = F0 + f * (F1 - F0);

                // convert to standard variables used by the other pair computes in HOOMD-blue
                float forcemag_divr = 0.0f;
                if (r > 0.0f)
                    forcemag_divr = F / r;
                float pair_eng = V;
                // calculate the virial
                float force_div2r = float(0.5) * forcemag_divr;
                virialxx +=  dx.x * dx.x * force_div2r;
                virialxy +=  dx.x * dx.y * force_div2r;
                virialxz +=  dx.x * dx.z * force_div2r;
                virialyy +=  dx.y * dx.y * force_div2r;
                virialyz +=  dx.y * dx.z * force_div2r;
                virialzz +=  dx.z * dx.z * force_div2r;

                // add up the force vector components (FLOPS: 7)
                force.x += dx.x * forcemag_divr;
                force.y += dx.y * forcemag_divr;
                force.z += dx.z * forcemag_divr;
                force.w += pair_eng;
                }
            }
        }

    // potential energy per particle must be halved
    force.w *= 0.5f;
    // now that the force calculation is complete, write out the result
    d_force[idx] = force;
    d_virial[0*virial_pitch+idx] = virialxx;
    d_virial[1*virial_pitch+idx] = virialxy;
    d_virial[2*virial_pitch+idx] = virialxz;
    d_virial[3*virial_pitch+idx] = virialyy;
    d_virial[4*virial_pitch+idx] = virialyz;
    d_virial[5*virial_pitch+idx] = virialzz;
    }

/*! \param d_force Device memory to write computed forces
    \param d_virial Device memory to write computed virials
    \param virial_pitch pitch of 2D virial array
    \param N number of particles
    \param n_ghost number of ghost particles
    \param d_pos particle positions on the device
    \param box Box dimensions used to implement periodic boundary conditions
    \param d_n_neigh Device memory array listing the number of neighbors for each particle
    \param d_nlist Device memory array containing the neighbor list contents
    \param nli Indexer for indexing \a d_nlist
    \param d_tables Tables of the potential and force
    \param d_params Parameters for each table associated with a type pair
    \param ntypes Number of particle types in the system
    \param table_width Number of points in each table
    \param block_size Block size at which to run the kernel

    \note This is just a kernel driver. See gpu_compute_table_forces_kernel for full documentation.
*/
cudaError_t gpu_compute_table_forces(float4* d_force,
                                     float* d_virial,
                                     const unsigned int virial_pitch,
                                     const unsigned int N,
                                     const unsigned int n_ghost,
                                     const Scalar4 *d_pos,
                                     const BoxDim& box,
                                     const unsigned int *d_n_neigh,
                                     const unsigned int *d_nlist,
                                     const Index2D& nli,
                                     const float2 *d_tables,
                                     const float4 *d_params,
                                     const unsigned int ntypes,
                                     const unsigned int table_width,
                                     const unsigned int block_size)
    {
    assert(d_params);
    assert(d_tables);
    assert(ntypes > 0);
    assert(table_width > 1);

    // index calculation helper
    Index2DUpperTriangular table_index(ntypes);

    // setup the grid to run the kernel
    dim3 grid( (int)ceil((double)N / (double)block_size), 1, 1);
    dim3 threads(block_size, 1, 1);

    // bind the pdata position texture
    pdata_pos_tex.normalized = false;
    pdata_pos_tex.filterMode = cudaFilterModePoint;
    cudaError_t error = cudaBindTexture(0, pdata_pos_tex, d_pos, sizeof(float4) * (N+n_ghost));
    if (error != cudaSuccess)
        return error;

    // bind the tables texture
    tables_tex.normalized = false;
    tables_tex.filterMode = cudaFilterModePoint;
    error = cudaBindTexture(0, tables_tex, d_tables, sizeof(float2) * table_width * table_index.getNumElements());
    if (error != cudaSuccess)
        return error;

    gpu_compute_table_forces_kernel<<< grid, threads, sizeof(float4)*table_index.getNumElements() >>>
            (d_force, d_virial, virial_pitch, N, d_pos, box, d_n_neigh, d_nlist, nli, d_params, ntypes, table_width);

    return cudaSuccess;
    }

// vim:syntax=cpp

