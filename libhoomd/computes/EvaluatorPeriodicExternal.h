/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2008, 2009 Ames Laboratory
Iowa State University and The Regents of the University of Michigan All rights
reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

Redistribution and use of HOOMD-blue, in source and binary forms, with or
without modification, are permitted, provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of HOOMD-blue's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR
ANY WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// Maintainer: jglaser

#ifndef __EVALUATOR_PERIODIC_EXTERNAL_H__
#define __EVALUATOR_PERIODIC_EXTERNAL_H__

#ifndef NVCC
#include <string>
#endif

#include <math.h>
#include "HOOMDMath.h"
#include "BoxDim.h"
#include "PeriodicExternalParams.h"

/*! \file EvaluatorPeriodicExternal.h
    \brief Defines the external potential evaluator to induce a periodic ordered phase
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

// SCALARASINT resolves to __float_as_int on the device and to __scalar_as_int on the host
#ifdef NVCC
#define SCALARASINT(x) __float_as_int(x)
#else
#define SCALARASINT(x) __scalar_as_int(x)
#endif

//! Class for evaluating sphere constraints
/*! <b>General Overview</b>
    EvaluatorPeriodicExternal is an evaluator to induce a periodic modulation on the concentration profile
    in the system, e.g. to generate a periodic phase in a system of diblock copolymers.

    The external potential \f$V(\vec{r}) \f$ is implemented using the following formula:

    \f[
    V(\vec{r}) = A * \tanh\left[\frac{1}{2 \pi p w} \cos\left(p \vec{b}_i\cdot\vec{r}\right)\right]
    \f]

    where \f$A\f$ is the ordering parameter, \f$\vec{b}_i\f$ is the reciprocal lattice vector direction
    \f$i=0..2\f$, \f$p\f$ the periodicity and \f$w\f$ the interface width
    (relative to the distance \f$2\pi/|\mathbf{b_i}|\f$ between planes in the \f$i\f$-direction).
    The modulation is one-dimensional. It extends along the lattice vector \f$\mathbf{a}_i\f$ of the
    simulation cell.
*/
class EvaluatorPeriodicExternal
    {
    public:

        //! type of parameters this external potential accepts
        typedef PeriodicExternalParams param_type;

        //! Constructs the constraint evaluator
        /*! \param X position of particle
            \param box box dimensions
            \param params per-type parameters of external potential
        */
        DEVICE EvaluatorPeriodicExternal(Scalar3 X, const BoxDim& box, const param_type& params)
            : m_pos(X),
              m_box(box)
            {
            m_order_parameter = params.order_parameter;
            m_lattice_vector_1 = params.lattice_vector_1;
            m_lattice_vector_2 = params.lattice_vector_2;
            m_lattice_vector_3 = params.lattice_vector_3;
            m_interface_width = params.interface_width;
            m_periodicity = params.periodicity;
            }

        //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial)
            {
            F.x = Scalar(0.0);
            F.y = Scalar(0.0);
            F.z = Scalar(0.0);
            energy = Scalar(0.0);

            // For this potential, since it uses scaled positions, the virial is always zero.
            for (unsigned int i = 0; i < 6; i++)
                virial[i] = Scalar(0.0);

            Scalar3 L = m_box.getL();

            Scalar cosine = Scalar(0.0);
            Scalar3 deriv = make_scalar3(0.0,0.0,0.0);
            Scalar clip_parameter = Scalar(1.0)/(Scalar(2.0*M_PI)*m_interface_width*m_periodicity);

            Scalar3 q1 = make_scalar3(2.0*M_PI*m_lattice_vector_1.x/L.x,
                                      2.0*M_PI*m_lattice_vector_1.y/L.y,
                                      2.0*M_PI*m_lattice_vector_1.z/L.z);
            Scalar3 q2 = make_scalar3(2.0*M_PI*m_lattice_vector_2.x/L.x,
                                      2.0*M_PI*m_lattice_vector_2.y/L.y,
                                      2.0*M_PI*m_lattice_vector_2.z/L.z);
            Scalar3 q3 = make_scalar3(2.0*M_PI*m_lattice_vector_3.x/L.x,
                                      2.0*M_PI*m_lattice_vector_3.y/L.y,
                                      2.0*M_PI*m_lattice_vector_3.z/L.z);
            Scalar arg1, sine1, arg2, sine2, arg3, sine3;
            arg1 = dot(m_pos, q1);
            arg2 = dot(m_pos, q2);
            arg3 = dot(m_pos, q3);
            sine1 = -1.0*sinf(arg1);
            sine2 = -1.0*sinf(arg2);
            sine3 = -1.0*sinf(arg3);

            deriv = deriv + sine1*q1 + sine2*q2 + sine3*q3;
            cosine = cosf(arg1) + cosf(arg2) + cosf(arg3);
            cosine *= clip_parameter;
            deriv *= clip_parameter;
            Scalar tanH = tanhf(cosine);
            energy = m_order_parameter*tanH;

            Scalar sechSq = (1.0 - tanH*tanH);
            Scalar f = m_order_parameter*sechSq;
            F = f*deriv;
            }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Must be short and all lowercase, as this is the name energies will be logged as
            via analyze.log.
        */
        static std::string getName()
            {
            return std::string("periodic");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        Scalar m_order_parameter;
        int3  m_lattice_vector_1;
        int3  m_lattice_vector_2;
        int3  m_lattice_vector_3;
        Scalar m_interface_width;      //!< width of interface between lamellae (relative to box length)
        unsigned int m_periodicity;   //!< number of lamellae of each type

   };


#endif // __EVALUATOR_EXTERNAL_LAMELLAR_H__
