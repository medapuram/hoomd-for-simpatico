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

#ifndef __EVALUATOR_CONSTRAINT_H__
#define __EVALUATOR_CONSTRAINT_H__

#include "HOOMDMath.h"

/*! \file EvaluatorConstraint.h
    \brief Defines basic evaluation methods common to all constraint force implementations
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for evaluating constraint forces
/*! <b>General Overview</b>
    EvaluatorConstraint is a low level computation class for use on the CPU and GPU. It provides basic functionality
    needed by all constraint forces.
*/
class EvaluatorConstraint
    {
    public:
        //! Constructs the constraint evaluator
        /*! \param _X Current position at the time of the constraint force calculation
            \param V Current velocity at the time of the constraint force calculation
            \param F Current net force at the time of the constraint force calculation
            \param _m Mass of the particle
            \param _deltaT Step size delta t
        */
        DEVICE EvaluatorConstraint(Scalar3 _X, Scalar3 V, Scalar3 F, Scalar _m, Scalar _deltaT)
            : X(_X), m(_m), deltaT(_deltaT)
            {
            // perform step 2 of this velocity verlet update and step 1 of the next to get
            // U = X(t+2deltaT) given X = X(t+deltaT)
            Scalar minv = Scalar(1.0)/m;
            Scalar dtsqdivm = deltaT*deltaT * minv;
            U.x = X.x + V.x * deltaT + F.x * dtsqdivm;
            U.y = X.y + V.y * deltaT + F.y * dtsqdivm;
            U.z = X.z + V.z * deltaT + F.z * dtsqdivm;
            }
        
        //! Evaluate the unconstrained position update U
        /*! \returns The unconstrained position update U
        */
        DEVICE Scalar3 evalU()
            {
            return U;
            }
        
        //! Evaluate the additional constraint force
        /*! \param FC output parameter where the computed force is written
            \param virial array of six scalars the computed virial tensor is written
            \param C constrained position particle will be moved to at the next step
            \return Additional force \a F needed to satisfy the constraint
        */
        DEVICE void evalConstraintForce(Scalar3& FC, Scalar *virial, const Scalar3& C)
            {
            // subtract a constrained update from U and get that F = (C-U)*m/dt^2
            Scalar moverdtsq = m / (deltaT * deltaT);
            FC.x = (C.x - U.x) * moverdtsq;
            FC.y = (C.y - U.y) * moverdtsq;
            FC.z = (C.z - U.z) * moverdtsq;
            
            // compute virial
            virial[0] = FC.x * X.x;
            virial[1] = Scalar(1./2.)*(FC.y * X.x + FC.x * X.y);
            virial[2] = Scalar(1./2.)*(FC.z * X.x + FC.x * X.z);
            virial[3] = FC.y * X.y;
            virial[4] = Scalar(1./2.)*(FC.z * X.y + FC.y * X.z);
            virial[5] = FC.z * X.z;
            }
        
    protected:
        Scalar3 U;      //!< Unconstrained position update
        Scalar3 X;      //!< Current particle position
        Scalar m;       //!< Saved mass value
        Scalar deltaT;  //!< Saved delta T value

    };


#endif // __PAIR_EVALUATOR_LJ_H__

