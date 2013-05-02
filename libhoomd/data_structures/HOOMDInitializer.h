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

/*! \file HOOMDInitializer.h
    \brief Declares the HOOMDInitializer class
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include "ParticleData.h"
#include "WallData.h"
#include "BondData.h"
#include "AngleData.h"
#include "DihedralData.h"
#include "xmlParser.h"

#include <string>
#include <vector>
#include <map>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#ifndef __HOOMD_INITIALIZER_H__
#define __HOOMD_INITIALIZER_H__

//! Forward declarations
class ExecutionConfiguation;
class SnapshotSystemData;

//! Initializes particle data from a Hoomd input file
/*! The input XML file format is identical to the output XML file format that HOOMDDumpWriter writes.
    For more information on the XML file format design see \ref page_dev_info. Although, HOOMD's
    user guide probably has a more up to date documentation on the format.

    When HOOMDInitializer is instantiated, it reads in the XML file specified in the constructor
    and parses it into internal data structures. The initializer is then ready to be passed
    to ParticleData which will then make the needed calls to copy the data into its representation.

    HOOMD's XML file format and this class are designed to be very extensible. Parsers for inidividual
    XML nodes are written in separate functions and stored by name in the map \c m_parser_map. As the
    main parser loops through, it reads in xml nodes and fires of parsers from this map to parse each
    of them. Adding a new node to the file format parser is as simple as adding a new node parser function
    (like parsePositionNode()) and adding it to the map in the constructor.

    \ingroup data_structs
*/
class HOOMDInitializer 
    {
    public:
        //! Loads in the file and parses the data
        HOOMDInitializer(boost::shared_ptr<const ExecutionConfiguration> exec_conf,
                         const std::string &fname);

        //! Returns the timestep of the simulation
        virtual unsigned int getTimeStep() const;
        
        //! Sets the timestep of the simulation
        virtual void setTimeStep(unsigned int ts);

        //! initializes a snapshot with the particle data
        virtual boost::shared_ptr<SnapshotSystemData> getSnapshot() const;

        //! simple vec for storing particle data
        struct vec
            {
            //! Default construtor
            vec() : x(0.0), y(0.0), z(0.0)
                {
                }
            //! Constructs a vec with given components
            /*! \param xp x-component
                \param yp y-component
                \param zp z-component
            */
            vec(Scalar xp, Scalar yp, Scalar zp) : x(xp), y(yp), z(zp)
                {
                }
            Scalar x;   //!< x-component
            Scalar y;   //!< y-component
            Scalar z;   //!< z-component
            };
            
        //! simple integer vec for storing particle data
        struct vec_int
            {
            //! Default construtor
            vec_int() : x(0), y(0), z(0)
                {
                }
            //! Constructs a vec with given components
            /*! \param xp x-component
                \param yp y-component
                \param zp z-component
            */
            vec_int(int xp, int yp, int zp) : x(xp), y(yp), z(zp)
                {
                }
            int x;  //!< x-component
            int y;  //!< y-component
            int z;  //!< z-component
            };
            
        //! Access the read particle positions
        const std::vector< vec >& getPos() { return m_pos_array; }
        
        //! Access the read images
        const std::vector< vec_int >& getImage() { return m_image_array; }

    private:
        //! Helper function to read the input file
        void readFile(const std::string &fname);
        //! Helper function to parse the box node
        void parseBoxNode(const XMLNode& node);
        //! Helper function to parse the position node
        void parsePositionNode(const XMLNode& node);
        //! Helper function to parse the image node
        void parseImageNode(const XMLNode& node);
        //! Helper function to parse the velocity node
        void parseVelocityNode(const XMLNode& node);
        //! Helper function to parse the mass node
        void parseMassNode(const XMLNode& node);
        //! Helper function to parse diameter node
        void parseDiameterNode(const XMLNode& node);
        //! Helper function to parse the type node
        void parseTypeNode(const XMLNode& node);
        //! Helper function to parse the body node
        void parseBodyNode(const XMLNode& node);
        //! Helper function to parse the bonds node
        void parseBondNode(const XMLNode& node);
        //! Helper function to parse the angle node
        void parseAngleNode(const XMLNode& node);
        //! Helper function to parse the dihedral node
        void parseDihedralNode(const XMLNode& node);
        //! Helper function to parse the improper node
        void parseImproperNode(const XMLNode& node);
        //! Parse charge node
        void parseChargeNode(const XMLNode& node);
        //! Parse wall node
        void parseWallNode(const XMLNode& node);
        //! Parse orientation node
        void parseOrientationNode(const XMLNode& node);
        //! Parse moment inertia node
        void parseMomentInertiaNode(const XMLNode& node);
        
        //! Helper function for identifying the particle type id
        unsigned int getTypeId(const std::string& name);
        //! Helper function for identifying the bond type id
        unsigned int getBondTypeId(const std::string& name);
        //! Helper function for identifying the angle type id
        unsigned int getAngleTypeId(const std::string& name);
        //! Helper function for identifying the dihedral type id
        unsigned int getDihedralTypeId(const std::string& name);
        //! Helper function for identifying the improper type id
        unsigned int getImproperTypeId(const std::string& name);
        
        std::map< std::string, boost::function< void (const XMLNode&) > > m_parser_map; //!< Map for dispatching parsers based on node type
        
        BoxDim m_box;   //!< Simulation box read from the file
        bool m_box_read;    //!< Stores the box we read in
        
        unsigned int m_num_dimensions;              //!< number of spatial dimensions
        std::vector< vec > m_pos_array;             //!< positions of all particles loaded
        std::vector< vec_int > m_image_array;       //!< images of all particles loaded
        std::vector< vec > m_vel_array;             //!< velocities of all particles loaded
        std::vector< Scalar > m_mass_array;         //!< masses of all particles loaded
        std::vector< Scalar > m_diameter_array;     //!< diameters of all particles loaded
        std::vector< unsigned int > m_type_array;   //!< type values for all particles loaded
        std::vector< unsigned int > m_body_array;   //!< body values for all particles loaded
        std::vector< Scalar > m_charge_array;       //!< charge of the particles loaded
        std::vector< Wall > m_walls;                //!< walls loaded from the file
        std::vector< Bond > m_bonds;                //!< Bonds read in from the file
        std::vector< Angle > m_angles;              //!< Angle read in from the file
        std::vector< Dihedral > m_dihedrals;        //!< Dihedral read in from the file
        std::vector< Dihedral > m_impropers;        //!< Improper read in from the file
        unsigned int m_timestep;                    //!< The time stamp
        
        std::vector<std::string> m_type_mapping;          //!< The created mapping between particle types and ids
        std::vector<std::string> m_bond_type_mapping;     //!< The created mapping between bond types and ids
        std::vector<std::string> m_angle_type_mapping;    //!< The created mapping between angle types and ids
        std::vector<std::string> m_dihedral_type_mapping; //!< The created mapping between dihedral types and ids
        std::vector<std::string> m_improper_type_mapping; //!< The created mapping between improper types and ids
        
        std::vector<Scalar4> m_orientation;             //!< Orientation of each particle
        std::vector<InertiaTensor> m_moment_inertia;    //!< Inertia tensor for each particle

        boost::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< The execution configuration
    };

//! Exports HOOMDInitializer to python
void export_HOOMDInitializer();

#endif



