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

#include "ComputeThermoGPU.cuh"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

//! Shared memory used in reducing the sums
extern __shared__ float3 compute_thermo_sdata[];
//! Shared memory used in reducing the sums of the pressure tensor
extern __shared__ float compute_pressure_tensor_sdata[];

/*! \file ComputeThermoGPU.cu
    \brief Defines GPU kernel code for computing thermodynamic properties on the GPU. Used by ComputeThermoGPU.
*/

//! Perform partial sums of the thermo properties on the GPU
/*! \param d_scratch Scratch space to hold partial sums. One element is written per block
    \param d_net_force Net force / pe array from ParticleData
    \param d_net_virial Net virial array from ParticleData
    \param virial_pitch pitch of 2D virial array
    \param d_velocity Particle velocity and mass array from ParticleData
    \param d_group_members List of group members for which to sum properties
    \param group_size Number of particles in the group
    
    All partial sums are packaged up in a float4 to keep pointer management down.
     - 2*Kinetic energy is summed in .x
     - Potential energy is summed in .y
     - W is summed in .z
    
    One thread is executed per group member. That thread reads in the values for its member into shared memory
    and then the block performs a reduction in parallel to produce a partial sum output for the block. These
    partial sums are written to d_scratch[blockIdx.x]. sizeof(float3)*block_size of dynamic shared memory are needed
    for this kernel to run.
*/

__global__ void gpu_compute_thermo_partial_sums(float4 *d_scratch,
                                                float4 *d_net_force,
                                                float *d_net_virial,
                                                const unsigned int virial_pitch,
                                                Scalar4 *d_velocity,
                                                unsigned int *d_group_members,
                                                unsigned int group_size)
    {
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    float3 my_element; // element of scratch space read in
    if (group_idx < group_size)
        {
        unsigned int idx = d_group_members[group_idx];
   
        // update positions to the next timestep and update velocities to the next half step
        float4 net_force = d_net_force[idx];
        float net_isotropic_virial;
        // (1/3)*trace of virial tensor
        net_isotropic_virial = Scalar(1.0/3.0)*
                               (d_net_virial[0*virial_pitch+idx]   // xx
                               +d_net_virial[3*virial_pitch+idx]   // yy
                               +d_net_virial[5*virial_pitch+idx]); // zz
        float4 vel = d_velocity[idx];
        float mass = vel.w;
        
        // compute our contribution to the sum
        my_element.x = mass * (vel.x*vel.x + vel.y*vel.y + vel.z*vel.z);
        my_element.y = net_force.w;
        my_element.z = net_isotropic_virial;
        }
    else
        {
        // non-participating thread: contribute 0 to the sum
        my_element = make_float3(0.0f, 0.0f, 0.0f);
        }
        
    compute_thermo_sdata[threadIdx.x] = my_element;
    __syncthreads();
    
    // reduce the sum in parallel
    int offs = blockDim.x >> 1;
    while (offs > 0)
        {
        if (threadIdx.x < offs)
            {
            compute_thermo_sdata[threadIdx.x].x += compute_thermo_sdata[threadIdx.x + offs].x;
            compute_thermo_sdata[threadIdx.x].y += compute_thermo_sdata[threadIdx.x + offs].y;
            compute_thermo_sdata[threadIdx.x].z += compute_thermo_sdata[threadIdx.x + offs].z;
            }
        offs >>= 1;
        __syncthreads();
        }
        
    // write out our partial sum
    if (threadIdx.x == 0)
        {
        float3 res = compute_thermo_sdata[0];
        d_scratch[blockIdx.x] = make_float4(res.x, res.y, res.z, 0.0f);
        }
    }

//! Perform partial sums of the pressure tensor on the GPU
/*! \param d_scratch Scratch space to hold partial sums. One element is written per block
    \param d_net_force Net force / pe array from ParticleData
    \param d_net_virial Net virial array from ParticleData
    \param virial_pitch pitch of 2D virial array
    \param d_velocity Particle velocity and mass array from ParticleData
    \param d_group_members List of group members for which to sum properties
    \param group_size Number of particles in the group

    One thread is executed per group member. That thread reads in the six values (components of the presure tensor)
    for its member into shared memory and then the block performs a reduction in parallel to produce a partial sum output for the block.
    These partial sums are written to d_scratch[i*gridDim.x + blockIdx.x], where i=0..5 is the index of the component.
    For this kernel to run, 6*sizeof(float)*block_size of dynamic shared memory are needed.
*/

__global__ void gpu_compute_pressure_tensor_partial_sums(float *d_scratch,
                                                float4 *d_net_force,
                                                float *d_net_virial,
                                                const unsigned int virial_pitch,
                                                float4 *d_velocity,
                                                unsigned int *d_group_members,
                                                unsigned int group_size)
    {
    // determine which particle this thread works on
    int group_idx = blockIdx.x * blockDim.x + threadIdx.x;

    float my_element[6]; // element of scratch space read in
    if (group_idx < group_size)
        {
        unsigned int idx = d_group_members[group_idx];

        // compute contribution to pressure tensor and store it in my_element
        float4 vel = d_velocity[idx];
        float mass = vel.w;
        my_element[0] = mass*vel.x*vel.x + d_net_virial[0*virial_pitch+idx];   // xx
        my_element[1] = mass*vel.x*vel.y + d_net_virial[1*virial_pitch+idx];   // xy
        my_element[2] = mass*vel.x*vel.z + d_net_virial[2*virial_pitch+idx];   // xz
        my_element[3] = mass*vel.y*vel.y + d_net_virial[3*virial_pitch+idx];   // yy
        my_element[4] = mass*vel.y*vel.z + d_net_virial[4*virial_pitch+idx];   // yz
        my_element[5] = mass*vel.z*vel.z + d_net_virial[5*virial_pitch+idx];   // zz
        }
    else
        {
        // non-participating thread: contribute 0 to the sum
        my_element[0] = 0.0f;
        my_element[1] = 0.0f;
        my_element[2] = 0.0f;
        my_element[3] = 0.0f;
        my_element[4] = 0.0f;
        my_element[5] = 0.0f;
        }

    for (unsigned int i = 0; i < 6; i++)
        compute_pressure_tensor_sdata[i*blockDim.x+threadIdx.x] = my_element[i];

    __syncthreads();

    // reduce the sum in parallel
    int offs = blockDim.x >> 1;
    while (offs > 0)
        {
        if (threadIdx.x < offs)
            {
            for (unsigned int i = 0; i < 6; i++)
                compute_pressure_tensor_sdata[i*blockDim.x+threadIdx.x] += compute_pressure_tensor_sdata[i*blockDim.x + threadIdx.x + offs];
            }
        offs >>= 1;
        __syncthreads();
        }

    // write out our partial sum
    if (threadIdx.x == 0)
        {
        for (unsigned int i = 0; i < 6; i++)
            d_scratch[gridDim.x * i + blockIdx.x] = compute_pressure_tensor_sdata[i*blockDim.x];
        }
    }

//! Complete partial sums and compute final thermodynamic quantities (for pressure, only isotropic contribution)
/*! \param d_properties Property array to write final values
    \param d_scratch Partial sums
    \param ndof Number of degrees of freedom this group posesses
    \param box Box the particles are in
    \param D Dimensionality of the system
    \param group_size Number of particles in the group
    \param num_partial_sums Number of partial sums in \a d_scratch
    \param external_virial External contribution to virial (1/3 trace)

   
    Only one block is executed. In that block, the partial sums are read in and reduced to final values. From the final
    sums, the thermodynamic properties are computed and written to d_properties.
    
    sizeof(float3)*block_size bytes of shared memory are needed for this kernel to run.
*/
__global__ void gpu_compute_thermo_final_sums(float *d_properties,
                                              float4 *d_scratch,
                                              unsigned int ndof,
                                              BoxDim box,
                                              unsigned int D,
                                              unsigned int group_size,
                                              unsigned int num_partial_sums,
                                              Scalar external_virial
                                              )
    {
    float3 final_sum = make_float3(0.0f, 0.0f, 0.0f);
    
    // sum up the values in the partial sum via a sliding window
    for (int start = 0; start < num_partial_sums; start += blockDim.x)
        {
        __syncthreads();
        if (start + threadIdx.x < num_partial_sums)
            {
            float4 scratch = d_scratch[start + threadIdx.x];
            compute_thermo_sdata[threadIdx.x] = make_float3(scratch.x, scratch.y, scratch.z);
            }
        else
            compute_thermo_sdata[threadIdx.x] = make_float3(0.0f, 0.0f, 0.0f);
        __syncthreads();
        
        // reduce the sum in parallel
        int offs = blockDim.x >> 1;
        while (offs > 0)
            {
            if (threadIdx.x < offs)
                {
                compute_thermo_sdata[threadIdx.x].x += compute_thermo_sdata[threadIdx.x + offs].x;
                compute_thermo_sdata[threadIdx.x].y += compute_thermo_sdata[threadIdx.x + offs].y;
                compute_thermo_sdata[threadIdx.x].z += compute_thermo_sdata[threadIdx.x + offs].z;
                }
            offs >>= 1;
            __syncthreads();
            }
        
        if (threadIdx.x == 0)
            {
            final_sum.x += compute_thermo_sdata[0].x;
            final_sum.y += compute_thermo_sdata[0].y;
            final_sum.z += compute_thermo_sdata[0].z;
            }
        }
        
    if (threadIdx.x == 0)
        {
        // compute final quantities
        float ke_total = final_sum.x * 0.5f;
        float pe_total = final_sum.y;
        float W = final_sum.z + external_virial;
        
        // compute the temperature
        Scalar temperature = Scalar(2.0) * Scalar(ke_total) / Scalar(ndof);

        // compute the pressure
        // volume/area & other 2D stuff needed

        Scalar volume;
        Scalar3 L = box.getL();

        if (D == 2)
            {
            // "volume" is area in 2D
            volume = L.x * L.y;
            // W needs to be corrected since the 1/3 factor is built in
            W *= Scalar(3.0)/Scalar(2.0);
            }
        else
            {
            volume = L.x * L.y * L.z;
            }

        // pressure: P = (N * K_B * T + W)/V
        Scalar pressure =  (Scalar(2.0) * ke_total / Scalar(D) + W) / volume;

        // fill out the GPUArray
        d_properties[thermo_index::temperature] = temperature;
        d_properties[thermo_index::pressure] = pressure;
        d_properties[thermo_index::kinetic_energy] = Scalar(ke_total);
        d_properties[thermo_index::potential_energy] = Scalar(pe_total);
        }
    }

//! Complete partial sums and compute final pressure tensor
/*! \param d_properties Property array to write final values
    \param d_scratch Partial sums
    \param box Box the particles are in
    \param group_size Number of particles in the group
    \param num_partial_sums Number of partial sums in \a d_scratch

    \param external_virial_xx External contribution to virial (xx component)
    \param external_virial_xy External contribution to virial (xy component)
    \param external_virial_xz External contribution to virial (xz component)
    \param external_virial_yy External contribution to virial (yy component)
    \param external_virial_yz External contribution to virial (yz component)
    \param external_virial_zz External contribution to virial (zz component)

    Only one block is executed. In that block, the partial sums are read in and reduced to final values. From the final
    sums, the thermodynamic properties are computed and written to d_properties.

    6*sizeof(float)*block_size bytes of shared memory are needed for this kernel to run.
*/
__global__ void gpu_compute_pressure_tensor_final_sums(float *d_properties,
                                              float *d_scratch,
                                              BoxDim box,
                                              unsigned int group_size,
                                              unsigned int num_partial_sums,
                                              Scalar external_virial_xx,
                                              Scalar external_virial_xy,
                                              Scalar external_virial_xz,
                                              Scalar external_virial_yy,
                                              Scalar external_virial_yz,
                                              Scalar external_virial_zz)
    {
    float final_sum[6];

    final_sum[0] = external_virial_xx;
    final_sum[1] = external_virial_xy;
    final_sum[2] = external_virial_xz;
    final_sum[3] = external_virial_yy;
    final_sum[4] = external_virial_yz;
    final_sum[5] = external_virial_zz;

    // sum up the values in the partial sum via a sliding window
    for (int start = 0; start < num_partial_sums; start += blockDim.x)
        {
        __syncthreads();
        if (start + threadIdx.x < num_partial_sums)
            {
            for (unsigned int i = 0; i < 6; i++)
                compute_pressure_tensor_sdata[i * blockDim.x + threadIdx.x] = d_scratch[i*num_partial_sums + start + threadIdx.x];
            }
        else
            for (unsigned int i = 0; i < 6; i++)
                compute_pressure_tensor_sdata[i * blockDim.x + threadIdx.x] = 0.0f;
        __syncthreads();

        // reduce the sum in parallel
        int offs = blockDim.x >> 1;
        while (offs > 0)
            {
            if (threadIdx.x < offs)
                {
                for (unsigned int i = 0; i < 6; i++)
                    compute_pressure_tensor_sdata[i*blockDim.x + threadIdx.x] += compute_pressure_tensor_sdata[i*blockDim.x + threadIdx.x + offs];
                }
            offs >>= 1;
            __syncthreads();
            }

        if (threadIdx.x == 0)
            {
            for (unsigned int i = 0; i < 6; i++)
                final_sum[i] += compute_pressure_tensor_sdata[i*blockDim.x];
            }
        }

    if (threadIdx.x == 0)
        {
        // fill out the GPUArray
        // we have thus far calculated the sum of the kinetic part of the pressure tensor
        // and the virial part, the definition includes an inverse factor of the box volume
        Scalar3 L = box.getL();
        float V = L.x * L.y * L.z;

        d_properties[thermo_index::pressure_xx] = final_sum[0]/V;
        d_properties[thermo_index::pressure_xy] = final_sum[1]/V;
        d_properties[thermo_index::pressure_xz] = final_sum[2]/V;
        d_properties[thermo_index::pressure_yy] = final_sum[3]/V;
        d_properties[thermo_index::pressure_yz] = final_sum[4]/V;
        d_properties[thermo_index::pressure_zz] = final_sum[5]/V;
        }
    }

//! Compute thermodynamic properties of a group on the GPU
/*! \param d_properties Array to write computed properties
    \param d_vel particle velocities and masses on the GPU
    \param d_group_members List of group members
    \param group_size Number of group members
    \param box Box the particles are in
    \param args Additional arguments
    \param compute_pressure_tensor whether to compute the full pressure tensor
    
    This function drives gpu_compute_thermo_partial_sums and gpu_compute_thermo_final_sums, see them for details.
*/

cudaError_t gpu_compute_thermo(float *d_properties,
                               Scalar4 *d_vel,
                               unsigned int *d_group_members,
                               unsigned int group_size,
                               const BoxDim& box,
                               const compute_thermo_args& args,
                               const bool compute_pressure_tensor
                               )
    {
    assert(d_properties);
    assert(d_vel);
    assert(d_group_members);
    assert(args.d_net_force);
    assert(args.d_net_virial);
    assert(args.d_scratch);
    
    dim3 grid(args.n_blocks, 1, 1);
    dim3 threads(args.block_size, 1, 1);
    unsigned int shared_bytes = sizeof(float3)*args.block_size;
   
    Scalar external_virial = Scalar(1.0/3.0)*(args.external_virial_xx
                             + args.external_virial_yy
                             + args.external_virial_zz);

    gpu_compute_thermo_partial_sums<<<grid,threads, shared_bytes>>>(args.d_scratch,
                                                                    args.d_net_force,
                                                                    args.d_net_virial,
                                                                    args.virial_pitch,
                                                                    d_vel,
                                                                    d_group_members,
                                                                    group_size);


    if (compute_pressure_tensor)
        {
        assert(args.d_scratch_pressure_tensor);

        shared_bytes = 6 * sizeof(float) * args.block_size;
        // run the kernel
        gpu_compute_pressure_tensor_partial_sums<<<grid, threads, shared_bytes>>>(args.d_scratch_pressure_tensor,
                                                                                  args.d_net_force,
                                                                                  args.d_net_virial,
                                                                                  args.virial_pitch,
                                                                                  d_vel,
                                                                                  d_group_members,
                                                                                  group_size);
        }

    // setup the grid to run the final kernel
    int final_block_size = 512;
    grid = dim3(1, 1, 1);
    threads = dim3(final_block_size, 1, 1);
    shared_bytes = sizeof(float3)*final_block_size;
    
    // run the kernel
    gpu_compute_thermo_final_sums<<<grid, threads, shared_bytes>>>(d_properties,
                                                                   args.d_scratch,
                                                                   args.ndof,
                                                                   box,
                                                                   args.D,
                                                                   group_size,
                                                                   args.n_blocks,
                                                                   external_virial);
    
    if (compute_pressure_tensor)
        {
        shared_bytes = 6 * sizeof(float) * final_block_size;
        // run the kernel
        gpu_compute_pressure_tensor_final_sums<<<grid, threads, shared_bytes>>>(d_properties,
                                                                               args.d_scratch_pressure_tensor,
                                                                               box,
                                                                               group_size,
                                                                               args.n_blocks,
                                                                               args.external_virial_xx,
                                                                               args.external_virial_xy,
                                                                               args.external_virial_xz,
                                                                               args.external_virial_yy,
                                                                               args.external_virial_yz,
                                                                               args.external_virial_zz
                                                                               );
        }

    return cudaSuccess;
    }

