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

#include "HOOMDMath.h"
#include "ParticleData.cuh"
#include "Index1D.h"

#ifdef WIN32
#include <cassert>
#else
#include <assert.h>
#endif

/*! \file PotentialBondGPU.cuh
    \brief Defines templated GPU kernel code for calculating the bond forces.
*/

#ifndef __POTENTIAL_BOND_GPU_CUH__
#define __POTENTIAL_BOND_GPU_CUH__

//! Wrapps arguments to gpu_cgbf
struct bond_args_t
    {
    //! Construct a bond_args_t
    bond_args_t(float4 *_d_force,
              float *_d_virial,
              const unsigned int _virial_pitch,
              const unsigned int _N,
              const unsigned int _n_ghost,
              const Scalar4 *_d_pos,
              const Scalar *_d_charge,
              const Scalar *_d_diameter,
              const BoxDim& _box,
              const uint2 *_d_gpu_bondlist,
              const unsigned int _pitch,
              const unsigned int *_d_gpu_n_bonds,
              const unsigned int _n_bond_types,
              const unsigned int _block_size)
                : d_force(_d_force),
                  d_virial(_d_virial),
                  virial_pitch(_virial_pitch),
                  N(_N),
                  n_ghost(_n_ghost),
                  d_pos(_d_pos),
                  d_charge(_d_charge),
                  d_diameter(_d_diameter),
                  box(_box),
                  d_gpu_bondlist(_d_gpu_bondlist),
                  pitch(_pitch),
                  d_gpu_n_bonds(_d_gpu_n_bonds),
                  n_bond_types(_n_bond_types),
                  block_size(_block_size)
        {
        };

    float4 *d_force;                   //!< Force to write out
    float *d_virial;                   //!< Virial to write out
    const unsigned int virial_pitch;   //!< pitch of 2D array of virial matrix elements
    unsigned int N;                    //!< number of particles
    unsigned int n_ghost;              //!< number of ghost particles
    const Scalar4 *d_pos;              //!< particle positions
    const Scalar *d_charge;            //!< particle charges
    const Scalar *d_diameter;          //!< particle diameters
    const BoxDim& box;            //!< Simulation box in GPU format
    const uint2 *d_gpu_bondlist;       //!< List of bonds stored on the GPU
    const unsigned int pitch;          //!< Pitch of 2D bond list
    const unsigned int *d_gpu_n_bonds; //!< List of number of bonds stored on the GPU
    const unsigned int n_bond_types;   //!< Number of bond types in the simulation
    const unsigned int block_size;     //!< Block size to execute
    };

#ifdef NVCC
//! Texture for reading particle positions
texture<float4, 1, cudaReadModeElementType> pdata_pos_tex;

//! Texture for reading particle diameters
texture<float, 1, cudaReadModeElementType> pdata_diam_tex;

//! Texture for reading particle charges
texture<float, 1, cudaReadModeElementType> pdata_charge_tex;

//! Kernel for calculating bond forces
/*! This kernel is called to calculate the bond forces on all N particles. Actual evaluation of the potentials and
    forces for each bond is handled via the template class \a evaluator.

    \param d_force Device memory to write computed forces
    \param d_virial Device memory to write computed virials
    \param virial_pitch pitch of 2D virial array
    \param N Number of particles in the system
    \param d_pos particle positions on the GPU
    \param d_charge particle charges
    \param d_diameter particle diameters
    \param box Box dimensions used to implement periodic boundary conditions
    \param blist List of bonds stored on the GPU
    \param pitch Pitch of 2D bond list
    \param n_bonds_list List of numbers of bonds stored on the GPU
    \param n_bond_type number of bond types
    \param d_params Parameters for the potential, stored per bond type
    \param d_flags Flag allocated on the device for use in checking for bonds that cannot be evaluated


    Certain options are controlled via template parameters to avoid the performance hit when they are not enabled.
    \tparam evaluator EvaluatorBond class to evualuate V(r) and -delta V(r)/r

*/
template< class evaluator >
__global__ void gpu_compute_bond_forces_kernel(float4 *d_force,
                                               float *d_virial,
                                               const unsigned int virial_pitch,
                                               const unsigned int N,
                                               const Scalar4 *d_pos,
                                               const Scalar *d_charge,
                                               const Scalar *d_diameter,
                                               const BoxDim box,
                                               const uint2 *blist,
                                               const unsigned int pitch,
                                               const unsigned int *n_bonds_list,
                                               const unsigned int n_bond_type,
                                               const typename evaluator::param_type *d_params,
                                               unsigned int *d_flags)
    {
    // start by identifying which particle we are to handle
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= N)
        return;

    // load in the length of the list for this thread (MEM TRANSFER: 4 bytes)
    int n_bonds =n_bonds_list[idx];

    // shared array for per bond type parameters
    extern __shared__ char s_data[];
    typename evaluator::param_type *s_params =
        (typename evaluator::param_type *)(&s_data[0]);

    // load in per bond type parameters
    if (threadIdx.x < n_bond_type)
       s_params[threadIdx.x] = d_params[threadIdx.x];
    __syncthreads();

    // read in the position of our particle. (MEM TRANSFER: 16 bytes)
    float4 postype = tex1Dfetch(pdata_pos_tex, idx);
    Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);

    // read in the diameter of our particle if needed
    Scalar diam(0.0);
    if (evaluator::needsDiameter())
        diam = tex1Dfetch(pdata_diam_tex, idx);
    else
        diam += 1.0f; // shutup compiler warning

    Scalar q(0.0);
    if (evaluator::needsCharge())
        q = tex1Dfetch(pdata_charge_tex, idx);
    else
        q += 0.0f; // shutup compiler warning

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
        // the volatile is needed to force the compiler to load the uint2 coalesced
        volatile uint2 cur_bond = blist[pitch*bond_idx + idx];

        int cur_bond_idx = cur_bond.x;
        int cur_bond_type = cur_bond.y;

        // get the bonded particle's position (MEM_TRANSFER: 16 bytes)
        float4 neigh_postypej = tex1Dfetch(pdata_pos_tex, cur_bond_idx);
        float3 neigh_pos= make_float3(neigh_postypej.x, neigh_postypej.y, neigh_postypej.z);

        // calculate dr (FLOPS: 3)
        float3 dx = pos - neigh_pos;

        // apply periodic boundary conditions (FLOPS: 12)
        dx = box.minImage(dx);

        // get the bond parameters (MEM TRANSFER: 8 bytes)
        typename evaluator::param_type param = s_params[cur_bond_type];

        float rsq = dot(dx, dx);

        // evaluate the potential
        float force_divr = 0.0f;
        float bond_eng = 0.0f;

        evaluator eval(rsq, param);

        // get the bonded particle's diameter if needed
        if (evaluator::needsDiameter())
            {
            float neigh_diam = tex1Dfetch(pdata_diam_tex, cur_bond_idx);
            eval.setDiameter(diam, neigh_diam);
            }
        if (evaluator::needsCharge())
            {
            float neigh_q = tex1Dfetch(pdata_charge_tex, cur_bond_idx);
            eval.setCharge(q, neigh_q);
            }

        bool evaluated = eval.evalForceAndEnergy(force_divr, bond_eng);

        if (evaluated)
            {
            // add up the virial (double counting, multiply by 0.5)
            float force_div2r = force_divr/2.0f;
            virial[0] += dx.x * dx.x * force_div2r; // xx
            virial[1] += dx.x * dx.y * force_div2r; // xy
            virial[2] += dx.x * dx.z * force_div2r; // xz
            virial[3] += dx.y * dx.y * force_div2r; // yy
            virial[4] += dx.y * dx.z * force_div2r; // yz
            virial[5] += dx.z * dx.z * force_div2r; // zz

            // add up the forces
            force.x += dx.x * force_divr;
            force.y += dx.y * force_divr;
            force.z += dx.z * force_divr;
            // energy is double counted: multiply by 0.5
            force.w += bond_eng * 0.5f;
            }
        else
            {
            *d_flags = 1;
            return;
            }
        }

    // now that the force calculation is complete, write out the result (MEM TRANSFER: 20 bytes);
    d_force[idx] = force;
    for (unsigned int i = 0; i < 6 ; i++)
        d_virial[i*virial_pitch + idx] = virial[i];
    }

//! Kernel driver that computes lj forces on the GPU for LJForceComputeGPU
/*! \param bond_args Other arugments to pass onto the kernel
    \param d_params Parameters for the potential, stored per bond type
    \param d_flags flags on the device - a 1 will be written if evaluation
                   of forces failed for any bond

    This is just a driver function for gpu_compute_bond_forces_kernel(), see it for details.
*/
template< class evaluator >
cudaError_t gpu_compute_bond_forces(const bond_args_t& bond_args,
                                    const typename evaluator::param_type *d_params,
                                    unsigned int *d_flags)
    {
    assert(d_params);
    assert(bond_args.n_bond_types > 0);

    // chck that block_size is valid
    assert(bond_args.block_size != 0);

    // setup the grid to run the kernel
    dim3 grid( bond_args.N / bond_args.block_size + 1, 1, 1);
    dim3 threads(bond_args.block_size, 1, 1);

    // bind the position texture
    pdata_pos_tex.normalized = false;
    pdata_pos_tex.filterMode = cudaFilterModePoint;
    cudaError_t error = cudaBindTexture(0, pdata_pos_tex, bond_args.d_pos, sizeof(Scalar4)*(bond_args.N+bond_args.n_ghost));
    if (error != cudaSuccess)
        return error;

    // bind the diamter texture
    pdata_diam_tex.normalized = false;
    pdata_diam_tex.filterMode = cudaFilterModePoint;
    error = cudaBindTexture(0, pdata_diam_tex, bond_args.d_diameter, sizeof(Scalar) *(bond_args.N+bond_args.n_ghost));
    if (error != cudaSuccess)
        return error;
    
    pdata_charge_tex.normalized = false;
    pdata_charge_tex.filterMode = cudaFilterModePoint;
    error = cudaBindTexture(0, pdata_charge_tex, bond_args.d_charge, sizeof(Scalar) * (bond_args.N+bond_args.n_ghost));
    if (error != cudaSuccess)
        return error;

    unsigned int shared_bytes = sizeof(typename evaluator::param_type) *
                                bond_args.n_bond_types;

    // run the kernel
    gpu_compute_bond_forces_kernel<evaluator><<<grid, threads, shared_bytes>>>(
        bond_args.d_force, bond_args.d_virial, bond_args.virial_pitch, bond_args.N,
        bond_args.d_pos, bond_args.d_charge, bond_args.d_diameter, bond_args.box, bond_args.d_gpu_bondlist,
        bond_args.pitch, bond_args.d_gpu_n_bonds, bond_args.n_bond_types, d_params, d_flags);

    return cudaSuccess;
    }
#endif

#endif // __POTENTIAL_BOND_GPU_CUH__

