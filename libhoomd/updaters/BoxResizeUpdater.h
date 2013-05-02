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

/*! \file BoxResizeUpdater.h
    \brief Declares an updater that resizes the simulation box of the system
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include <boost/shared_ptr.hpp>

#include "Updater.h"
#include "Variant.h"

#ifndef __BOXRESIZEUPDATER_H__
#define __BOXRESIZEUPDATER_H__

//! Updates the simulation box over time
/*! This simple updater gets the box lengths from specified variants and sets those box sizes
    over time. As an option, particles can be rescaled with the box lengths or left where they are.

    \ingroup updaters
*/
class BoxResizeUpdater : public Updater
    {
    public:
        //! Constructor
        BoxResizeUpdater(boost::shared_ptr<SystemDefinition> sysdef,
                         boost::shared_ptr<Variant> Lx,
                         boost::shared_ptr<Variant> Ly,
                         boost::shared_ptr<Variant> Lz,
                         boost::shared_ptr<Variant> xy,
                         boost::shared_ptr<Variant> xz,
                         boost::shared_ptr<Variant> yz);

        //! Destructor
        virtual ~BoxResizeUpdater();
        
        //! Sets parameter flags
        void setParams(bool scale_particles);
        
        //! Take one timestep forward
        virtual void update(unsigned int timestep);
        
    private:
        boost::shared_ptr<Variant> m_Lx;    //!< Box Lx vs time
        boost::shared_ptr<Variant> m_Ly;    //!< Box Ly vs time
        boost::shared_ptr<Variant> m_Lz;    //!< Box Lz vs time
        boost::shared_ptr<Variant> m_xy;    //!< Box xy tilt factor vs time
        boost::shared_ptr<Variant> m_xz;    //!< Box xz tilt factor vs time
        boost::shared_ptr<Variant> m_yz;    //!< Box yz tilt factor vs time
        bool m_scale_particles;                //!< Set to true if particle positions are to be scaled as well
    };

//! Export the BoxResizeUpdater to python
void export_BoxResizeUpdater();

#endif

