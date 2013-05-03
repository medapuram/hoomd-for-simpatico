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

#ifndef __EVALUATOR_LOCAL_EXTERNAL_H__
#define __EVALUATOR_LOCAL_EXTERNAL_H__

#ifndef NVCC
#include <string>
#endif

#include <math.h>
#include "HOOMDMath.h"
#include "BoxDim.h"
#include "LocalExternalParams.h"

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
class EvaluatorLocalExternal
    {
    public:

        typedef LocalExternalParams param_type;

        //! Constructs the constraint evaluator
        /*! \param X position of particle
            \param box box dimensions
            \param params per-type parameters of external potential
        */
        DEVICE EvaluatorLocalExternal(Scalar3 X, const BoxDim& box, const param_type& params)
            : m_pos(X),
              m_box(box)
            {
            m_perp_index =  params.perp_index;
            m_parallel_index = params.parallel_index;
            m_fraction = params.fraction_length;
            m_orderParameter = params.order_parameter;
            m_interfaceWidth = params.interface_width;
            m_periodicity = params.periodicity;
            }

        //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial)
            {
            Scalar3 a2, a3;
            Scalar3 b2, b3;

            F.x = Scalar(0.0);
            F.y = Scalar(0.0);
            F.z = Scalar(0.0);
            energy = Scalar(0.0);

            // For this potential, since it uses scaled positions, the virial is always zero.
            for (unsigned int i = 0; i < 6; i++)
                virial[i] = Scalar(0.0);

            Scalar V_box = m_box.getVolume();
            // compute the vector pointing from P to V
            if (m_perp_index == 0)
                {
                a2 = m_box.getLatticeVector(1);
                a3 = m_box.getLatticeVector(2);
                }
            else if (m_perp_index == 1)
                {
                a2 = m_box.getLatticeVector(2);
                a3 = m_box.getLatticeVector(0);
                }
            else if (m_perp_index == 2)
                {
                a2 = m_box.getLatticeVector(0);
                a3 = m_box.getLatticeVector(1);
                }
            if (m_parallel_index == 0)
                {
                b2 = m_box.getLatticeVector(1);
                b3 = m_box.getLatticeVector(2);
                }
           else if (m_parallel_index == 1)
                {
                b2 = m_box.getLatticeVector(2);
                b3 = m_box.getLatticeVector(0);
                }
            else if (m_parallel_index == 2)
                {
                b2 = m_box.getLatticeVector(0);
                b3 = m_box.getLatticeVector(1);
                }

            Scalar3 a = Scalar(2.0*M_PI)*make_scalar3(a2.y*a3.z-a2.z*a3.y,
                                                      a2.z*a3.x-a2.x*a3.z,
                                                      a2.x*a3.y-a2.y*a3.x)/V_box;

            Scalar3 b = make_scalar3(b2.y*b3.z-b2.z*b3.y,
                                     b2.z*b3.x-b2.x*b3.z,
                                     b2.x*b3.y-b2.y*b3.x)/V_box;

            Scalar clipParameter, arg, clipcos, tanH, sechSq, parallelFactor1, parallelFactor2, parallelFactor;
            Scalar parallelSechSq1, parallelSechSq2;

            Scalar3 q = a*m_periodicity;
            clipParameter   = Scalar(1.0)/Scalar(2.0*M_PI)/(m_periodicity*m_interfaceWidth);
            arg = dot(m_pos,q);
            clipcos = clipParameter*cosf(arg);
            tanH = tanhf(clipcos);
            parallelFactor1 = tanhf( Scalar(2.0*M_PI)*((dot(m_pos,b) - Scalar(0.5) + (Scalar(0.5)*m_fraction))/m_fraction) );
            parallelFactor2 = tanhf( Scalar(2.0*M_PI)*((-dot(m_pos,b) + Scalar(0.5) + (Scalar(0.5)*m_fraction))/m_fraction) );
            parallelFactor = (Scalar(1.0) + parallelFactor1*parallelFactor2)/2;
            sechSq = (Scalar(1.0) - tanH*tanH);
            parallelSechSq1 = ((Scalar(1.0) - parallelFactor1*parallelFactor1)*parallelFactor2*Scalar(-0.5))/m_fraction;
            parallelSechSq2 = ((Scalar(1.0) - parallelFactor2*parallelFactor2)*parallelFactor1*Scalar(0.5))/m_fraction;

            F = m_orderParameter*sechSq*clipParameter*sinf(arg)*parallelFactor*q + m_orderParameter*tanH*(parallelSechSq1 + parallelSechSq2)*b;
            energy = m_orderParameter*tanH*parallelFactor;
            }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Must be short and all lowercase, as this is the name energies will be logged as
            via analyze.log.
        */
        static std::string getName()
            {
            return std::string("local");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        unsigned int m_perp_index;         //!< cartesian index of direction along which the lammellae should be orientied
        unsigned int m_parallel_index;         //!< cartesian index of direction along which the lammellae should be orientied
        Scalar m_fraction;      //!< ordering parameter
        Scalar m_orderParameter;      //!< ordering parameter
        Scalar m_interfaceWidth;      //!< width of interface between lamellae (relative to box length)
        unsigned int m_periodicity;   //!< number of lamellae of each type
   };


#endif // __EVALUATOR_LOCAL_EXTERNAL_H__
