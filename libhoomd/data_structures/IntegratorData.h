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

/*! \file IntegratorData.h
    \brief Contains declarations for IntegratorData.
 */

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#ifndef __INTEGRATORDATA_H__
#define __INTEGRATORDATA_H__

#include "ParticleData.h"
#include <string>


//! Stores integrator variables
/*! The integration state is necessary for exact restarts.  Extended systems 
    integrators in the spirit of Nose-Hoover store the positions, velocities, 
    etc. of the fictitious variables.  Other integrators store a random number 
    seed.
    \ingroup data_structs
*/
struct IntegratorVariables
    {
    std::string type;                   //!<The type of integrator (NVT, NPT, etc.)
    std::vector<Scalar> variable;       //!<Variables that define the integration state
    };

//! Stores all integrator variables in the simulation
/*! IntegratorData keeps track of the parameters for all of the integrators 
    defined in the simulation, so that they can be saved and reloaded from data files.
    
    Each integrator must register with IntegratorData by calling registerIntegrator(), which returns an unsinged int
    to be used to access the variables for that integrator. The same sequence of registerIntegrator() calls will produce
    the same set of handles, so they can be used to read existing state values after loading data from a file.
    
    The state of current registered integrators is reset when a new IntegratorData is constructed. This is consistent
    with the above use-case, as the construction of a new IntegratorData means the construction of a new SystemData, and
    hence a new series of constructed Integrators, which will then re-register.
    
    \ingroup data_structs
*/
class IntegratorData
    {
    public:
        //! Constructs an empty list with no integrator variables
        IntegratorData() : m_num_registered(0) {}
       
        //! Constructs an IntegratorData from a given set of IntegratorVariables
        IntegratorData(const std::vector<IntegratorVariables>& variables)
            : m_num_registered(0)
            {
            m_integrator_variables = variables;
            }

        //! Destructor
        ~IntegratorData() {}
        
        //! Register an integrator (should occur during integrator construction)
        unsigned int registerIntegrator();
        
        //! Get the number of integrator variables
        /*! \return Number of integrator variables present
        */
        unsigned int getNumIntegrators() const
            {
            return (unsigned int)m_integrator_variables.size();
            }
            
        //! Load a number of integrator variables
        /*! \param n Number of variables to load
            When loading from a file, a given number of integrator variables must be preloaded without registering them.
            This method does that. After calling load(n), setIntegratorVariables() can be called for \a i up to \a n-1
        */
        void load(unsigned int n)
            {
            m_integrator_variables.resize(n);
            }
        
        //! Get a collection of integrator variables
        /*! \param i access integrator variables for integrator i
        */
        const IntegratorVariables& getIntegratorVariables(unsigned int i) const
            {
            assert(i < m_integrator_variables.size()); return m_integrator_variables[i];
            }
            
        //! Set a collection of integrator variables
        /*! \param i set integrator variables for integrator i
            \param v Variables to set
        */
        void setIntegratorVariables(unsigned int i, const IntegratorVariables& v)
            {
            assert(i < m_integrator_variables.size()); m_integrator_variables[i] = v;
            }
            
    private:
        unsigned int m_num_registered;                                  //!< Number of integrators that have registered
        std::vector<IntegratorVariables> m_integrator_variables;        //!< List of the integrator variables defined
        
    };

//! Exports IntegratorData to python
void export_IntegratorData();

#endif

