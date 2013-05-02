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

#include "BondTablePotentialGPU.cuh"


#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

/*! \file BondTablePotentialGPU.cu
    \brief Defines GPU kernel code for calculating the table bond forces. Used by BondTablePotentialGPU.
*/


//! Texture for reading table values
texture<float2, 1, cudaReadModeElementType> tables_tex;

/*!  This kernel is called to calculate the table pair forces on all N particles

    \param d_force Device memory to write computed forces
    \param d_virial Device memory to write computed virials
    \param virial_pitch Pitch of 2D virial array
    \param N number of particles in system
    \param d_pos device array of particle positions
    \param box Box dimensions used to implement periodic boundary conditions
    \param blist List of bonds stored on the GPU
    \param pitch Pitch of 2D bond list
    \param n_bonds_list List of numbers of bonds stored on the GPU
    \param n_bond_type number of bond types
    \param d_params Parameters for each table associated with a type pair
    \param table_value index helper function
    \param d_flags Flag allocated on the device for use in checking for bonds that cannot be evaluated

    See BondTablePotential for information on the memory layout.

    \b Details:
    * Table entries are read from tables_tex. Note that currently this is bound to a 1D memory region. Performance tests
      at a later date may result in this changing.
*/
__global__ void gpu_compute_bondtable_forces_kernel(float4* d_force,
                                     float* d_virial,
                                     const unsigned int virial_pitch,
                                     const unsigned int N,
                                     const Scalar4 *d_pos,
                                     const BoxDim box,
                                     const uint2 *blist,
                                     const unsigned int pitch,
                                     const unsigned int *n_bonds_list,
                                     const unsigned int n_bond_type,
                                     const float4 *d_params,
                                     const Index2D table_value,
                                     unsigned int *d_flags)
    {

    
    // read in params for easy and fast access in the kernel
    extern __shared__ float4 s_params[];
    for (unsigned int cur_offset = 0; cur_offset < n_bond_type; cur_offset += blockDim.x)
        {
        if (cur_offset + threadIdx.x < n_bond_type)
            s_params[cur_offset + threadIdx.x] = d_params[cur_offset + threadIdx.x];
        }
    __syncthreads();


    // start by identifying which particle we are to handle
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= N)
        return;

    // load in the length of the list for this thread (MEM TRANSFER: 4 bytes)
    int n_bonds =n_bonds_list[idx];

    // read in the position of our particle.
    Scalar4 postype = d_pos[idx];
    Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);

    // initialize the force to 0
    float4 force = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    // initialize the virial tensor to 0
    float virial[6];
    for (unsigned int i = 0; i < 6; i++)
        virial[i] = 0;

    // loop over neighbors
    for (int bond_idx = 0; bond_idx < n_bonds; bond_idx++)
        {
        // MEM TRANSFER: 8 bytes
        uint2 cur_bond = blist[pitch*bond_idx + idx];

        int cur_bond_idx = cur_bond.x;
        int cur_bond_type = cur_bond.y;

        // get the bonded particle's position (MEM_TRANSFER: 16 bytes)
        Scalar4 neigh_postype = d_pos[cur_bond_idx];
        Scalar3 neigh_pos = make_scalar3(neigh_postype.x, neigh_postype.y, neigh_postype.z);

        // calculate dr (FLOPS: 3)
        float3 dx = pos - neigh_pos;

        // apply periodic boundary conditions (FLOPS: 12)
        dx = box.minImage(dx);

        // access needed parameters
        float4 params = s_params[cur_bond_type];
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
            float2 VF0 = tex1Dfetch(tables_tex, table_value(value_i, cur_bond_type));
            float2 VF1 = tex1Dfetch(tables_tex, table_value(value_i+1, cur_bond_type));
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
            float bond_eng = V;
            // calculate the virial
            float force_div2r = float(0.5) * forcemag_divr;
            virial[0] += dx.x * dx.x * force_div2r; // xx
            virial[1] += dx.x * dx.y * force_div2r; // xy
            virial[2] += dx.x * dx.z * force_div2r; // xz
            virial[3] += dx.y * dx.y * force_div2r; // yy
            virial[4] += dx.y * dx.z * force_div2r; // yz
            virial[5] += dx.z * dx.z * force_div2r; // zz

            // add up the force vector components (FLOPS: 7)
            force.x += dx.x * forcemag_divr;
            force.y += dx.y * forcemag_divr;
            force.z += dx.z * forcemag_divr;
            force.w += bond_eng * 0.5f;
            }
        else
            {
            *d_flags = 1;
            }
        }


    // now that the force calculation is complete, write out the result (MEM TRANSFER: 20 bytes);
    d_force[idx] = force;
    for (unsigned int i = 0; i < 6 ; i++)
        d_virial[i*virial_pitch + idx] = virial[i];
    }


/*! \param d_force Device memory to write computed forces
    \param d_virial Device memory to write computed virials
    \param virial_pitch pitch of 2D virial array
    \param N number of particles
    \param d_pos particle positions on the device
    \param box Box dimensions used to implement periodic boundary conditions
    \param blist List of bonds stored on the GPU
    \param pitch Pitch of 2D bond list
    \param n_bonds_list List of numbers of bonds stored on the GPU
    \param n_bond_type number of bond types
    \param d_tables Tables of the potential and force
    \param d_params Parameters for each table associated with a type pair
    \param table_width Number of entries in the table
    \param table_value indexer helper
    \param d_flags flags on the device - a 1 will be written if evaluation
                   of forces failed for any bond
    \param block_size Block size at which to run the kernel

    \note This is just a kernel driver. See gpu_compute_bondtable_forces_kernel for full documentation.
*/
cudaError_t gpu_compute_bondtable_forces(float4* d_force,
                                     float* d_virial,
                                     const unsigned int virial_pitch,
                                     const unsigned int N,
                                     const Scalar4 *d_pos,
                                     const BoxDim &box,
                                     const uint2 *blist,
                                     const unsigned int pitch,
                                     const unsigned int *n_bonds_list,
                                     const unsigned int n_bond_type,
                                     const float2 *d_tables,
                                     const float4 *d_params,
                                     const unsigned int table_width,
                                     const Index2D &table_value,
                                     unsigned int *d_flags,
                                     const unsigned int block_size)
    {
    assert(d_params);
    assert(d_tables);
    assert(n_bond_type > 0);
    assert(table_width > 1);


    // setup the grid to run the kernel
    dim3 grid( (int)ceil((double)N / (double)block_size), 1, 1);
    dim3 threads(block_size, 1, 1);


    // bind the tables texture
    tables_tex.normalized = false;
    tables_tex.filterMode = cudaFilterModePoint;
    cudaError_t error = cudaBindTexture(0, tables_tex, d_tables, sizeof(float2) * table_value.getNumElements());
    if (error != cudaSuccess)
        return error;

    gpu_compute_bondtable_forces_kernel<<< grid, threads, sizeof(float4)*n_bond_type >>>
            (d_force,
             d_virial,
             virial_pitch,
             N,
             d_pos,
             box,
             blist,
             pitch,
             n_bonds_list,
             n_bond_type,
             d_params,
             table_value,
             d_flags);

    return cudaSuccess;
    }

// vim:syntax=cpp
