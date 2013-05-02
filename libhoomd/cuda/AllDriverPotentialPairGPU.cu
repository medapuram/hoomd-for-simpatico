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

// Maintainer: joaander / Everyone is free to add additional potentials

/*! \file AllDriverPotentialPairGPU.cu
    \brief Defines the driver functions for computing all types of pair forces on the GPU
*/

#include "EvaluatorPairLJ.h"
#include "EvaluatorPairGauss.h"
#include "EvaluatorPairSLJ.h"
#include "EvaluatorPairYukawa.h"
#include "EvaluatorPairMorse.h"
#include "PotentialPairDPDThermoGPU.cuh"
#include "EvaluatorPairDPDThermo.h"
#include "AllDriverPotentialPairGPU.cuh"
#include "EvaluatorPairEwald.h"
#include "EvaluatorPairDPDLJThermo.h"
#include "EvaluatorPairForceShiftedLJ.h"

cudaError_t gpu_compute_ljtemp_forces(const pair_args_t& pair_args,
                                      const float2 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairLJ>(pair_args,
                                                    d_params);
    }

cudaError_t gpu_compute_gauss_forces(const pair_args_t& pair_args,
                                     const float2 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairGauss>(pair_args,
                                                       d_params);
    }

cudaError_t gpu_compute_slj_forces(const pair_args_t& pair_args,
                                   const float2 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairSLJ>(pair_args,
                                                     d_params);
    }

cudaError_t gpu_compute_yukawa_forces(const pair_args_t& pair_args,
                                      const float2 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairYukawa>(pair_args,
                                                        d_params);
    }


cudaError_t gpu_compute_morse_forces(const pair_args_t& pair_args,
                                      const float4 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairMorse>(pair_args,
                                                       d_params);
    }

cudaError_t gpu_compute_dpdthermodpd_forces(const dpd_pair_args_t& args,
                                            const float2 *d_params)
    {
    return gpu_compute_dpd_forces<EvaluatorPairDPDThermo>(args,
                                                          d_params);
    }


cudaError_t gpu_compute_dpdthermo_forces(const pair_args_t& pair_args,
                                         const float2 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairDPDThermo>(pair_args,
                                                           d_params);
    }


cudaError_t gpu_compute_ewald_forces(const pair_args_t& pair_args,
                                     const float *d_params)
    {
    return  gpu_compute_pair_forces<EvaluatorPairEwald>(pair_args,
                                                        d_params);
    }


cudaError_t gpu_compute_dpdljthermodpd_forces(const dpd_pair_args_t& args,
                                              const float4 *d_params)
    {
    return gpu_compute_dpd_forces<EvaluatorPairDPDLJThermo>(args,
                                                            d_params);
    }


cudaError_t gpu_compute_dpdljthermo_forces(const pair_args_t& args,
                                           const float4 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairDPDLJThermo>(args,
                                                             d_params);
    }

cudaError_t gpu_compute_force_shifted_lj_forces(const pair_args_t & args,
                                                const float2 *d_params)
    {
    return gpu_compute_pair_forces<EvaluatorPairForceShiftedLJ>(args,
                                                                d_params);
    }
