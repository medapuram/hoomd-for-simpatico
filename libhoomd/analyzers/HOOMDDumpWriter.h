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

/*! \file HOOMDDumpWriter.h
    \brief Declares the HOOMDDumpWriter class
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include <string>

#include <boost/shared_ptr.hpp>

#include "Analyzer.h"

#ifndef __HOOMD_DUMP_WRITER_H__
#define __HOOMD_DUMP_WRITER_H__

//! Analyzer for writing out HOOMD  dump files
/*! HOOMDDumpWriter can be used to write out xml files containing various levels of information
    of the current time step of the simulation. At a minimum, the current time step and box
    dimensions are output. Optionally, particle positions, velocities and types can be included
    in the file.

    Usage:<br>
    Construct a HOOMDDumpWriter, attaching it to a ParticleData and specifying a base file name.
    Call analyze(timestep) to output a dump file with the state of the current time step
    of the simulation. It will create base_file.timestep.xml where timestep is a 0-padded
    10 digit number. The 0 padding is so files sorted "alphabetically" will be read in
    numberical order.

    To include positions, velocities and types, see: setOutputPosition() setOutputVelocity()
    and setOutputType(). Similarly, walls and bonds can be included with setOutputWall() and
    setOutputBond().

    Future versions will include the ability to dump forces on each particle to the file also.

    For information on the structure of the xml file format: see \ref page_dev_info
    Although, HOOMD's  user guide probably has a more up to date documentation on the format.
    \ingroup analyzers
*/
class HOOMDDumpWriter : public Analyzer
    {
    public:
        //! Construct the writer
        HOOMDDumpWriter(boost::shared_ptr<SystemDefinition> sysdef, std::string base_fname);

        //! Destructor
        ~HOOMDDumpWriter();

        //! Write out the data for the current timestep
        void analyze(unsigned int timestep);
        //! Enables/disables the writing of the particle positions
        void setOutputPosition(bool enable);
        //! Enables/disables the writing of particle images
        void setOutputImage(bool enable);
        //! Enables/disables the writing of particle velocities
        void setOutputVelocity(bool enable);
        //! Enables/disables the writing of particle masses
        void setOutputMass(bool enable);
        //! Enables/disables the writing of particle diameters
        void setOutputDiameter(bool enable);
        //! Enables/disables the writing of particle types
        void setOutputType(bool enable);
        //! Enables/disables the writing of bonds
        void setOutputBond(bool enable);
        //! Enables/disables the writing of angles
        void setOutputAngle(bool enable);
        //! Enables/disables the writing of dihedrals
        void setOutputDihedral(bool enable);
        //! Enables/disables the writing of impropers
        void setOutputImproper(bool enable);
        //! Enables/disables the writing of walls
        void setOutputWall(bool enable);
        //! Enables/disables the writing of acceleration
        void setOutputAccel(bool enable);
        //! Enables/disables the writing of body
        void setOutputBody(bool enable);
        //! Enables/disables the writing of charge
        void setOutputCharge(bool enable);
        //! Enables/disables the writing of orientation
        void setOutputOrientation(bool enable);
        //! Enables/disables the writing of moment_inertia
        void setOutputMomentInertia(bool enable);
        //! Sets the vizsigma value to write
        void setVizSigma(Scalar vizsigma)
            {
            m_vizsigma = vizsigma;
            m_vizsigma_set = true;
            }
        
        //! Writes a file at the current time step
        void writeFile(std::string fname, unsigned int timestep);
    private:
        std::string m_base_fname;   //!< String used to store the file name of the XML file
        bool m_output_position;     //!< true if the particle positions should be written
        bool m_output_image;        //!< true if the particle positions should be written
        bool m_output_velocity;     //!< true if the particle velocities should be written
        bool m_output_mass;         //!< true if the particle masses should be written
        bool m_output_diameter;     //!< true if the particle diameters should be written
        bool m_output_type;         //!< true if the particle types should be written
        bool m_output_bond;         //!< true if the bonds should be written
        bool m_output_angle;        //!< true if the angles should be written
        bool m_output_wall;         //!< true if the walls should be written
        bool m_output_dihedral;     //!< true if dihedrals should be written
        bool m_output_improper;     //!< true if impropers should be written
        bool m_output_accel;        //!< true if acceleration should be written
        bool m_output_body;         //!< true if body should be written
        bool m_output_charge;       //!< true if body should be written
        bool m_output_orientation;  //!< true if orientation should be written
        bool m_output_moment_inertia;  //!< true if moment_inertia should be written
        Scalar m_vizsigma;          //!< vizsigma value to write out to xml files
        bool m_vizsigma_set;        //!< true if vizsigma has been set
        };

//! Exports the HOOMDDumpWriter class to python
void export_HOOMDDumpWriter();

#endif

