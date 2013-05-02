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

/*! \file RandomGenerator.h
    \brief Contains declarations for RandomGenerator and related classes.
 */

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include "ParticleData.h"
#include "BondData.h"

#include "SnapshotSystemData.h"

#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>

#ifndef __RANDOM_GENERATOR_H__
#define __RANDOM_GENERATOR_H__

//! Stores particles as they are generated in RandomGenerator
/*! \ingroup data_structs
    GeneratedParticles is a holding area where particles are stored when
    being generated by RandomGenerator and ParticleGenerator instances.

    It includes helper functions and data structures for placing particles
    that do not overlap. These helpers use the radius for each particle type
    as specified by setSeparationRadius(). Every particle type that will be
    generated must be specified before generation can begin.

    After all particles are placed in GeneratedParticles, RandomGenerator will
    then translate that data over to ParticleData in the initializer.
*/
class GeneratedParticles
    {
    public:
        //! Stores a single particle in GeneratedParticles
        struct particle
            {
            //! Default constructor
            particle() : x(0.0), y(0.0), z(0.0), ix(0), iy(0), iz(0), type("") {}
            Scalar x;   //!< X-coordinate
            Scalar y;   //!< Y-coordinate
            Scalar z;   //!< Z-coordinate
            int ix;     //!< Shift in x direction
            int iy;     //!< Shift in y direction
            int iz;     //!< Shift in z direction
            std::string type;   //!< Particle's type name
            unsigned int type_id;   //!< Particle's type id
            };
            
        //! Constructor
        GeneratedParticles(boost::shared_ptr<const ExecutionConfiguration> exec_conf, unsigned int n_particles, const BoxDim& box, const std::map< std::string, Scalar >& radii);
        //! Empty constructor
        /*! Included so that GeneratedParticles can be stored in a vector.
        */
        GeneratedParticles() { }
        
        //! Check if a particle can be placed while obeying the separation radii
        bool canPlace(const particle& p);
        
        //! Place a particle
        void place(const particle& p, unsigned int idx);
        
        //! Undo the placement of a particle
        void undoPlace(unsigned int idx);
        
        //! Get the box
        const BoxDim& getBox()
            {
            return m_box;
            }
            
        //! Add a bond
        void addBond(unsigned int a, unsigned int b, const std::string& type="");
        
    private:
        friend class RandomGenerator;
       
        boost::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< The execution configuration
        std::vector<particle> m_particles;                  //!< The generated particles
        BoxDim m_box;                                       //!< Box the particles are in
        std::vector< std::vector<unsigned int> > m_bins;    //!< Bins the particles are placed in for efficient distance checks
        int m_Mx;       //!< Number of bins in the x direction
        int m_My;       //!< Number of bins in the y direction
        int m_Mz;       //!< Number of bins in the z direction
        std::map< std::string, Scalar > m_radii;    //!< Separation radii accessed by particle type
        
        //! Structure representing a single bond
        struct bond
            {
            //! Default constructor
            bond() : tag_a(0), tag_b(0), type(""), type_id(0)
                {
                }
                
            //! Construct a bond between two particles
            /*! \param a tag of the first particle in the bond
                \param b tag of the second particle in the bond
                \param _type the type name of the bond
            */
            bond(unsigned int a, unsigned int b, const std::string& _type) : tag_a(a), tag_b(b), type(_type), type_id(0)
                {
                }
            unsigned int tag_a;     //!< First particle in the bond
            unsigned int tag_b;     //!< Second particle in the bond
            std::string type;       //!< Type of the bond
            unsigned int type_id;   //!< Type id of the bond (assigned in RandomGenerator::generate())
            };
            
        std::vector< bond > m_bonds;    //!< Bonds read in from the file
    };

//! Abstract interface for classes that generate particles
/*! \ingroup data_structs
    A ParticleGenerator is the workhorse that actually chooses where to place particles.
    A single ParticleGenerator should only place a small number of inter-related particles
    on each call to generateParticles() (i.e. a single polymer or a small cluster of particles).
    Larger systems are to be built from multiple calls of generateParticles() by RandomGenerator.
*/
class ParticleGenerator
    {
    public:
        //! Destructor
        virtual ~ParticleGenerator() {}
        
        //! Returns the number of particles that will be generated
        /*! Derived classes must implement this method so that RandomGenerator
            can properly allocate the space for the particles.
        
            The value returned by this function cannot vary from call to call
            for a particular instance of ParticleGenerator. Once instantiated,
            a ParticleGenerator must always generate the same number of particles
            each time it is called.
        */
        virtual unsigned int getNumToGenerate()=0;
        
        //! Actually generate the requested particles
        /*! \param particles Place generated particles here after a GeneratedParticles::canPlace() check
            \param rnd Random number generator to use
            \param start_idx Starting index to generate particles at
            Derived classes must implement this method. RandomGenerator will
            call it to generate the particles. Particles should be placed at indices
            \a start_idx, \a start_idx + 1, ... \a start_idx + getNumToGenerate()-1
        */
        virtual void generateParticles(GeneratedParticles& particles, boost::mt19937& rnd, unsigned int start_idx)=0;
    };

//! Generates random polymers
/*! \ingroup data_structs
    This ParticleGenerator can be used to generate systems of bead-spring polymers of any combination of
    partile types specified in an array.
*/
class PolymerParticleGenerator : public ParticleGenerator
    {
    public:
        //! Constructor
        PolymerParticleGenerator(boost::shared_ptr<const ExecutionConfiguration> exec_conf, Scalar bond_len, const std::vector<std::string>& types, const std::vector<unsigned int>& bond_a, const std::vector<unsigned int>& bond_b, const std::vector<string>& bond_type, unsigned int max_attempts);
        
        //! Returns the number of particles in each polymer
        virtual unsigned int getNumToGenerate()
            {
            return (unsigned int)m_types.size();
            }
            
        //! Generates a single polymer
        virtual void generateParticles(GeneratedParticles& particles, boost::mt19937& rnd, unsigned int start_idx);
        
    private:
        boost::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< Execution configuration for messaging
        Scalar m_bond_len;                  //!< Bond length
        std::vector<std::string> m_types;   //!< Particle types for each polymer bead
        std::vector<unsigned int> m_bond_a; //!< First particle in the bond pair
        std::vector<unsigned int> m_bond_b; //!< Second particle in the bond pair
        std::vector<string> m_bond_type;    //!< Type name of the bond
        unsigned int m_max_attempts;        //!< Number of attemps to make for each particle placement
        
        //! helper function to place particles recursively
        bool generateNextParticle(GeneratedParticles& particles, boost::mt19937& rnd, unsigned int i, unsigned int start_idx, const GeneratedParticles::particle& prev_particle);
        
    };


//! Generates a random particle system given a set of ParticleGenerator classes
/*! \ingroup data_structs
    RandomGenerator is the high level Initializer that brings all the pieces together to
    generate a random system of particles. The structure and types of the particles generated
    (i.e. a polymer system of A6B7A6 polymers) is determined by the ParticleGenerator classes
    that are added.

    Mixture systems can be created by adding more than one ParticleGenerator.
    class. Testing should be done to confirm this, but it is probably best to
    add the largest objects first so that smaller ones can be generated around them.

    By default, bonds are named "bond". This can be changed by calling setBondType().

    \b Usage:<br>
    Before the initializer can be passed to a ParticleData for initialization, the following
    steps must be performed.
     -# Contstruct a RandomGenerator (duh) with a given box size
     -# Set radii for all particle types to be generated
     -# Construct and add any number of ParticleGenerator instances to the RandomGenerator
     -# Call generate() to actually place the particles
*/
class RandomGenerator
    {
    public:
        //! Set the parameters
        RandomGenerator(boost::shared_ptr<const ExecutionConfiguration> exec_conf,
                        const BoxDim& box,
                        unsigned int seed);

        //! Empty Destructor
        virtual ~RandomGenerator() { }
        
        //! initializes a snapshot with the particle data
        virtual boost::shared_ptr<SnapshotSystemData> getSnapshot() const;

        //! Sets the separation radius for a particle
        void setSeparationRadius(string type, Scalar radius);
        
        //! Adds a generator
        void addGenerator(unsigned int repeat, boost::shared_ptr<ParticleGenerator> generator);
        
        //! Place the particles
        void generate();
       
    private:
        boost::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< The execution configuration
        BoxDim m_box;                                       //!< Precalculated box
        unsigned int m_seed;                                //!< Random seed to use
        GeneratedParticles m_data;                          //!< Actual particle data genreated
        std::map< std::string, Scalar > m_radii;            //!< Separation radii accessed by particle type
        std::vector< boost::shared_ptr<ParticleGenerator> > m_generators;   //!< Generators to place particles
        std::vector< unsigned int > m_generator_repeat;     //!< Repeat count for each generator
        std::vector<std::string> m_type_mapping;            //!< The created mapping between particle types and ids
        std::vector<std::string> m_bond_type_mapping;       //!< The created mapping between bond types and ids
        
        //! Helper function for identifying the particle type id
        unsigned int getTypeId(const std::string& name);
        //! Helper function for identifying the bond type id
        unsigned int getBondTypeId(const std::string& name);
    };

//! Exports RandomGenerator and related classes to python
void export_RandomGenerator();

#endif
