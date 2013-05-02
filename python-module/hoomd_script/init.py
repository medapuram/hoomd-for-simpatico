# -- start license --
# Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
# (HOOMD-blue) Open Source Software License Copyright 2008-2011 Ames Laboratory
# Iowa State University and The Regents of the University of Michigan All rights
# reserved.

# HOOMD-blue may contain modifications ("Contributions") provided, and to which
# copyright is held, by various Contributors who have granted The Regents of the
# University of Michigan the right to modify and/or distribute such Contributions.

# You may redistribute, use, and create derivate works of HOOMD-blue, in source
# and binary forms, provided you abide by the following conditions:

# * Redistributions of source code must retain the above copyright notice, this
# list of conditions, and the following disclaimer both in the code and
# prominently in any materials provided with the distribution.

# * Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions, and the following disclaimer in the documentation and/or
# other materials provided with the distribution.

# * All publications and presentations based on HOOMD-blue, including any reports
# or published results obtained, in whole or in part, with HOOMD-blue, will
# acknowledge its use according to the terms posted at the time of submission on:
# http://codeblue.umich.edu/hoomd-blue/citations.html

# * Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
# http://codeblue.umich.edu/hoomd-blue/

# * Apart from the above required attributions, neither the name of the copyright
# holder nor the names of HOOMD-blue's contributors may be used to endorse or
# promote products derived from this software without specific prior written
# permission.

# Disclaimer

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
# WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -- end license --

# Maintainer: joaander / All Developers are free to add commands for new features

from optparse import OptionParser;

import hoomd;

import math;
import sys;
import gc;
import os;

from hoomd_script import util;
from hoomd_script import globals;
from hoomd_script import data;

## \package hoomd_script.init
# \brief Data initialization commands
#
# Commands in the init package initialize the particle system. Initialization via
# any of the commands here must be done before any other command in hoomd_script can
# be run.
#
# \sa \ref page_quick_start

## Tests if the system has been initialized
#
# Returns True if a previous init.create* or init.read* command has completed successfully and initialized the system.
# Returns False otherwise.
def is_initialized():
    if globals.system is None:
        return False;
    else:
        return True;

## Resets all hoomd_script variables
#
# After calling init.reset() all global variables used in hoomd_script are cleared and all allocated
# memory is freed so the simulation can begin anew without needing to launch hoomd again.
#
# \note There is a very important memory management issue that must be kept in mind when using
# reset(). If you have saved a variable such as an integrator or a force for changing parameters, 
# that saved object \b must be deleted before the reset() command is called. If all objects are 
# not deleted, then a memory leak will result causing repeated runs of even a small simulation 
# to eventually run the system out of memory. reset() will throw an error if it detects that this 
# is the case.
#
# \note When using the python data access in hoomd scripts, iterators must also be deleted
# \code
# for p in sysdef.particles:
#   # do something
#
# del p
# init.reste()
# \endcode
#
# \b Example:
# \code
# init.create_random(N=1000, phi_p = 0.2)
# lj = pair.lj(r_cut=3.0)
# .... setup and run simulation
# del lj
# init.reset()
# init.create_random(N=2000, phi_p = 0.2)
# .... setup and run simulation
# \endcode
def reset():
    if not is_initialized():
        globals.msg.warning("Trying to reset an uninitialized system\n");
        return;

    # perform some reference counting magic to verify that the user has cleared all saved variables
    sysdef = globals.system_definition;
    globals.clear();
    
    gc.collect();
    count = sysdef.getPDataRefs()

    # note: the check should be against 2, getrefcount counts the temporary reference 
    # passed to it in the argument
    expected_count = 6
    if count != expected_count:
        globals.msg.warning("Not all saved variables were cleared before calling reset()\n");
        globals.msg.warning(str(count-expected_count) + " references to the particle data still exist somewhere\n");
        globals.msg.warning("Going to try and reset anyways, further errors (such as out of memory) may result\n");

    del sysdef
    gc.collect();
    gc.collect();

## Create an empty system
#
# \param N Number of particles to create
# \param box A 3-tuple of floats specifying the box lengths (Lx, Ly, Lz) (in distance units)
# \param n_particle_types Number of particle types to create
# \param n_bond_types Number of bond types to create
# \param n_angle_types Number of angle types to create
# \param n_dihedral_types Number of dihedral types to create
# \param n_improper_types Number of improper types to create
#
# \b Examples:
# \code
# system = init.create_empty(N=1000, box=(10, 10, 10)
# system = init.create_empty(N=64000, box=(20,20,20), n_particle_types=2)
# system = init.create_empty(N=64000, box=(20,20,20), n_bond_types=1, n_dihedral_types=2, n_improper_types=4)
# \endcode
#
# After init.create_empty returns, the requested number of particles will have been created <b>but all of them
# will have <i> DEFAULT VALUES</i> </b> and further initialization \b MUST be performed. See hoomd_script.data
# for full details on how such initialization can be performed.
#
# Specifically, all created particles will be:
# - At position 0,0,0
# - Have velocity 0,0,0
# - In box image 0,0,0
# - Have the type 'A'
# - Have charge 0
# - Have a mass of 1.0
#
# Furthermore, as a current limitation in the way create_empty works the defined particle type names are:
# 'A', 'B', 'C', ... up to \a n_particle_types types. An intuitive way for resetting these types to user-defined ones
# will come in a future enhancement.
# 
# \note The resulting empty system must have its particles fully initialized via python code, \b BEFORE
# any other hoomd_script commands are executed. As as example of what might go wrong, if the pair.lj command were to be
# run before the initial particle positions were set, \b all particles would have position 0,0,0 and the memory 
# initialized by the neighbor list would be so large that the memory allocation would fail.
#
# \sa hoomd_script.data
def create_empty(N, box, n_particle_types=1, n_bond_types=0, n_angle_types=0, n_dihedral_types=0, n_improper_types=0):
    util.print_status_line();
    
    # check if initialization has already occurred
    if is_initialized():
        globals.msg.error("Cannot initialize more than once\n");
        raise RuntimeError('Error initializing');
    
    my_exec_conf = _create_exec_conf();

    # create the empty system
    boxdim = hoomd.BoxDim(float(box[0]), float(box[1]), float(box[2]));
    globals.system_definition = hoomd.SystemDefinition(N,
                                                       boxdim,
                                                       n_particle_types,
                                                       n_bond_types,
                                                       n_angle_types,
                                                       n_dihedral_types,
                                                       n_improper_types,
                                                       my_exec_conf);
    
    # initialize the system
    globals.system = hoomd.System(globals.system_definition, 0);
    
    _perform_common_init_tasks();
    return data.system_data(globals.system_definition);

## Reads initial system state from an XML file
#
# \param filename File to read
# \param time_step (if specified) Time step number to use instead of the one stored in the XML file
#
# \b Examples:
# \code
# init.read_xml(filename="data.xml")
# init.read_xml(filename="directory/data.xml")
# init.read_xml(filename="restart.xml", time_step=0)
# system = init.read_xml(filename="data.xml")
# \endcode
#
# All particles, bonds, etc...  are read from the XML file given, 
# setting the initial condition of the simulation.
# After this command completes, the system is initialized allowing 
# other commands in hoomd_script to be run. For more details
# on the file format read by this command, see \ref page_xml_file_format.
#
# All values are read in native units, see \ref page_units for more information.
#
# If \a time_step is specified, its value will be used as the initial time 
# step of the simulation instead of the one read from the XML file.
#
# The result of init.read_xml can be saved in a variable and later used to read and/or change particle properties
# later in the script. See hoomd_script.data for more information.
#
# \sa dump.xml
def read_xml(filename, time_step = None):
    util.print_status_line();
    
    # initialize GPU/CPU execution configuration and MPI early
    my_exec_conf = _create_exec_conf();

    # check if initialization has already occured
    if is_initialized():
        globals.msg.error("Cannot initialize more than once\n");
        raise RuntimeError("Error creating random polymers");

    # read in the data
    initializer = hoomd.HOOMDInitializer(my_exec_conf,filename);
    snapshot = initializer.getSnapshot()

    my_domain_decomposition = _create_domain_decomposition(snapshot.global_box);

    if my_domain_decomposition is not None:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf, my_domain_decomposition);
    else:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf);

    # initialize the system
    if time_step is None:
        globals.system = hoomd.System(globals.system_definition, initializer.getTimeStep());
    else:
        globals.system = hoomd.System(globals.system_definition, time_step);
    
    _perform_common_init_tasks();
    return data.system_data(globals.system_definition);

## Reads initial system state from a binary file
#
# \param filename File to read
#
# \b Examples:
# \code
# init.read_bin(filename="data.bin.gz")
# init.read_bin(filename="directory/data.bin")
# system = init.read_bin(filename="data.bin.gz")
# \endcode
#
# All particles, bonds, etc...  are read from the binary file given, setting the initial condition of the simulation.
# Binary restart files also include state information needed to continue integrating time forward as if the previous job
# had never stopped. For more information see dump.bin.
#
# After this command completes, the system is initialized allowing other commands in hoomd_script to be run.
#
# The presence or lack of a .gz extension determines whether init.read_bin will attempt to decompress the %data
# before reading it.
#
# The result of init.read_bin can be saved in a variable and later used to read and/or change particle properties
# later in the script. See hoomd_script.data for more information.
#
# \sa dump.bin
def read_bin(filename):
    util.print_status_line();
    
    # initialize GPU/CPU execution configuration and MPI early
    my_exec_conf = _create_exec_conf();

    # check if initialization has already occurred
    if is_initialized():
        globals.msg.error("Cannot initialize more than once\n");
        raise RuntimeError('Error initializing');

    # read in the data
    initializer = hoomd.HOOMDBinaryInitializer(my_exec_conf,filename);
    snapshot = initializer.getSnapshot()

    my_domain_decomposition = _create_domain_decomposition(snapshot.global_box);

    if my_domain_decomposition is not None:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf, my_domain_decomposition);
    else:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf);

    # initialize the system
    if time_step is None:
        globals.system = hoomd.System(globals.system_definition, initializer.getTimeStep());
    else:
        globals.system = hoomd.System(globals.system_definition, time_step);
 
    _perform_common_init_tasks();
    return data.system_data(globals.system_definition);

## Generates N randomly positioned particles of the same type
#
# \param N Number of particles to create
# \param phi_p Packing fraction of particles in the simulation box (unitless)
# \param name Name of the particle type to create
# \param min_dist Minimum distance particles will be separated by (in distance units)
#
# \b Examples:
# \code
# init.create_random(N=2400, phi_p=0.20)
# init.create_random(N=2400, phi_p=0.40, min_dist=0.5)
# system = init.create_random(N=2400, phi_p=0.20)
# \endcode
#
# \a N particles are randomly placed in the simulation box. The 
# dimensions of the created box are such that the packing fraction
# of particles in the box is \a phi_p. The number density \e n
# is related to the packing fraction by \f$n = 6/\pi \cdot \phi_P\f$
# assuming the particles have a radius of 0.5.
# All particles are created with the same type, given by \a name.
# 
# The result of init.create_random can be saved in a variable and later used to read and/or change particle properties
# later in the script. See hoomd_script.data for more information.
#
def create_random(N, phi_p, name="A", min_dist=0.7):
    util.print_status_line();

    # initialize GPU/CPU execution configuration and MPI early
    my_exec_conf = _create_exec_conf();

    # check if initialization has already occured
    if is_initialized():
        globals.msg.error("Cannot initialize more than once\n");
        raise RuntimeError("Error initializing");

    # abuse the polymer generator to generate single particles
    
    # calculat the box size
    L = math.pow(math.pi/6.0*N / phi_p, 1.0/3.0);
    box = hoomd.BoxDim(L);
    
    # create the generator
    generator = hoomd.RandomGenerator(my_exec_conf, box, 12345);
    
    # build type list
    type_vector = hoomd.std_vector_string();
    type_vector.append(name);
    
    # empty bond lists for single particles
    bond_ab = hoomd.std_vector_uint();
    bond_type = hoomd.std_vector_string();
        
    # create the generator
    generator.addGenerator(int(N), hoomd.PolymerParticleGenerator(my_exec_conf, 1.0, type_vector, bond_ab, bond_ab, bond_type, 100));
    
    # set the separation radius
    generator.setSeparationRadius(name, min_dist/2.0);
        
    # generate the particles
    generator.generate();

    # initialize snapshot
    snapshot = generator.getSnapshot()
            
    my_domain_decomposition = _create_domain_decomposition(snapshot.global_box);
    if my_domain_decomposition is not None:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf, my_domain_decomposition);
    else:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf);

    # initialize the system
    globals.system = hoomd.System(globals.system_definition, 0);
    
    _perform_common_init_tasks();
    return data.system_data(globals.system_definition);

## Generates any number of randomly positioned polymers of configurable types
#
# \param box BoxDim specifying the simulation box to generate the polymers in
# \param polymers Specification for the different polymers to create (see below)
# \param separation Separation radii for different particle types (see below)
# \param seed Random seed to use
#
# Any number of polymers can be generated, of the same or different types, as 
# specified in the argument \a polymers. Parameters for each polymer include
# bond length, particle type list, bond list, and count.
#
# The syntax is best shown by example. The below line specifies that 600 block copolymers
# A6B7A6 with a %bond length of 1.2 be generated.
# \code
# polymer1 = dict(bond_len=1.2, type=['A']*6 + ['B']*7 + ['A']*6, 
#                 bond="linear", count=600)
# \endcode
# Here is an example for a second polymer, specifying just 100 polymers made of 5 B beads
# bonded in a branched pattern
# \code
# polymer2 = dict(bond_len=1.2, type=['B']*5, 
#                 bond=[(0, 1), (1,2), (1,3), (3,4)] , count=100)
# \endcode
# The \a polymers argument can be given a list of any number of polymer types specified
# as above. \a count randomly generated polymers of each type in the list will be
# generated in the system.
#
# In detail: 
# - \a bond_len defines the %bond length of the generated polymers. This should 
#   not necessarily be set to the equilibrium %bond length! The generator is dumb and doesn't know
#   that bonded particles can be placed closer together than the separation (see below). Thus
#   \a bond_len must be at a minimum set at twice the value of the largest separation radius. An 
#   error will be generated if this is not the case.
# - \a type is a python list of strings. Each string names a particle type in the order that
#   they will be created in generating the polymer.
# - \a %bond can be specified as "linear" in which case the generator connects all particles together
#   with bonds to form a linear chain. \a %bond can also be given a list if python tuples (see example
#   above). 
#   - Each tuple in the form of \c (a,b) specifies that particle \c a of the polymer be bonded to
#   particle \c b. These bonds are given the default type name of 'polymer' to be used when specifying parameters to 
#   bond forces such as bond.harmonic.
#   - A tuple with three elements (a,b,type) can be used as above, but with a custom name for the bond. For example,
#   a simple branched polymer with different bond types on each branch could be defined like so:
#\code
#bond=[(0,1), (1,2), (2,3,'branchA'), (3,4,'branchA), (2,5,'branchB'), (5,6,'branchB')]
#\endcode
# 
# \a separation \b must contain one entry for each particle type specified in \a polymers
# ('A' and 'B' in the examples above). The value given is the separation radius of each
# particle of that type. The generated polymer system will have no two overlapping 
# particles.
#
# \b Examples:
# \code
# init.create_random_polymers(box=hoomd.BoxDim(35), 
#                             polymers=[polymer1, polymer2], 
#                             separation=dict(A=0.35, B=0.35));
# 
# init.create_random_polymers(box=hoomd.BoxDim(31), 
#                             polymers=[polymer1], 
#                             separation=dict(A=0.35, B=0.35), seed=52);
# 
# init.create_random_polymers(box=hoomd.BoxDim(18,10,25), 
#                             polymers=[polymer2], 
#                             separation=dict(A=0.35, B=0.35), seed=12345);
# \endcode
#
# With all other parameters the same, create_random_polymers will always create the
# same system if \a seed is the same. Set a different \a seed (any integer) to create
# a different random system with the same parameters. Note that different versions
# of HOOMD \e may generate different systems even with the same seed due to programming
# changes.
#
# \note 1. For relatively dense systems (packing fraction 0.4 and higher) the simple random
# generation algorithm may fail to find room for all the particles and print an error message. 
# There are two methods to solve this. First, you can lower the separation radii allowing particles 
# to be placed closer together. Then setup integrate.nve with the \a limit option set to a 
# relatively small value. A few thousand time steps should relax the system so that the simulation can be
# continued without the limit or with a different integrator. For extremely troublesome systems,
# generate it at a very low density and shrink the box with the command update.box_resize
# to the desired final size.
#
# \note 2. The polymer generator always generates polymers as if there were linear chains. If you 
# provide a non-linear %bond topology, the bonds in the initial configuration will be stretched 
# significantly. This normally doesn't pose a problem for harmonic bonds (bond.harmonic) as
# the system will simply relax over a few time steps, but can cause the system to blow up with FENE 
# bonds (bond.fene). 
#
# \note 3. While the custom %bond list allows you to create ring shaped polymers, testing shows that
# such conformations have trouble relaxing and get stuck in tangled configurations. If you need 
# to generate a configuration of rings, you may need to write your own specialized initial configuration
# generator that writes HOOMD XML input files (see \ref page_xml_file_format). HOOMD's built-in polymer generator
# attempts to be as general as possible, but unfortunately cannot work in every possible case.
# 
# The result of init.create_random_polymers can be saved in a variable and later used to read and/or change particle
# properties later in the script. See hoomd_script.data for more information.
#
def create_random_polymers(box, polymers, separation, seed=1):
    util.print_status_line();
    
    # initialize GPU/CPU execution configuration and MPI early
    my_exec_conf = _create_exec_conf();

    # check if initialization has already occured
    if is_initialized():
        globals.msg.error("Cannot initialize more than once\n");
        raise RuntimeError("Error creating random polymers");

    if type(polymers) != type([]) or len(polymers) == 0:
        globals.msg.error("polymers specified incorrectly. See the hoomd_script documentation\n");
        raise RuntimeError("Error creating random polymers");
    
    if type(separation) != type(dict()) or len(separation) == 0:
        globals.msg.error("polymers specified incorrectly. See the hoomd_script documentation\n");
        raise RuntimeError("Error creating random polymers");
    
    # create the generator
    generator = hoomd.RandomGenerator(my_exec_conf,box, seed);
    
    # make a list of types used for an eventual check vs the types in separation for completeness
    types_used = [];
    
    # track the minimum bond length
    min_bond_len = None;
    
    # build the polymer generators
    for poly in polymers:
        type_list = [];
        # check that all fields are specified
        if not 'bond_len' in poly:
            globals.msg.error('Polymer specification missing bond_len\n');
            raise RuntimeError("Error creating random polymers");
        
        if min_bond_len is None:
            min_bond_len = poly['bond_len'];
        else:
            min_bond_len = min(min_bond_len, poly['bond_len']);
        
        if not 'type' in poly:
            globals.msg.error('Polymer specification missing type\n');
            raise RuntimeError("Error creating random polymers");
        if not 'count' in poly:
            globals.msg.error('Polymer specification missing count\n');
            raise RuntimeError("Error creating random polymers");
        if not 'bond' in poly:
            globals.msg.error('Polymer specification missing bond\n');
            raise RuntimeError("Error creating random polymers");
                
        # build type list
        type_vector = hoomd.std_vector_string();
        for t in poly['type']:
            type_vector.append(t);
            if not t in types_used:
                types_used.append(t);
        
        # build bond list
        bond_a = hoomd.std_vector_uint();
        bond_b = hoomd.std_vector_uint();
        bond_name = hoomd.std_vector_string();
        
        # if the bond setting is 'linear' create a default set of bonds
        if poly['bond'] == 'linear':
            for i in range(0,len(poly['type'])-1):
                bond_a.push_back(i);
                bond_b.push_back(i+1);
                bond_name.append('polymer')
        #if it is a list, parse the user custom bonds
        elif type(poly['bond']) == type([]):
            for t in poly['bond']:
                # a 2-tuple gets the default 'polymer' name for the bond
                if len(t) == 2:
                    a,b = t;
                    name = 'polymer';
                # and a 3-tuple specifies the name directly
                elif len(t) == 3:
                    a,b,name = t;
                else:
                    globals.msg.error('Custom bond ' + str(t) + ' must have either two or three elements\n');
                    raise RuntimeError("Error creating random polymers");
                                    
                bond_a.push_back(a);
                bond_b.push_back(b);
                bond_name.append(name);
        else:
            globals.msg.error('Unexpected argument value for polymer bond\n');
            raise RuntimeError("Error creating random polymers");
        
        # create the generator
        generator.addGenerator(int(poly['count']), hoomd.PolymerParticleGenerator(my_exec_conf, poly['bond_len'], type_vector, bond_a, bond_b, bond_name, 100));
        
        
    # check that all used types are in the separation list
    for t in types_used:
        if not t in separation:
            globals.msg.error("No separation radius specified for type " + str(t) + "\n");
            raise RuntimeError("Error creating random polymers");
            
    # set the separation radii, checking that it is within the minimum bond length
    for t,r in separation.items():
        generator.setSeparationRadius(t, r);
        if 2*r >= min_bond_len:
            globals.msg.error("Separation radius " + str(r) + " is too big for the minimum bond length of " + str(min_bond_len) + " specified\n");
            raise RuntimeError("Error creating random polymers");
        
    # generate the particles
    generator.generate();

    # copy over data to snapshot
    snapshot = generator.getSnapshot()
        
    my_domain_decomposition = _create_domain_decomposition(snapshot.global_box);

    if my_domain_decomposition is not None:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf, my_domain_decomposition);
    else:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf);

    # initialize the system
    globals.system = hoomd.System(globals.system_definition, 0);
    
    _perform_common_init_tasks();
    return data.system_data(globals.system_definition);

## Initializes the system from a snapshot
# 
# \param snapshot The snapshot to initialize the system from
#
# Snapshots temporarily store system %data. Snapshots contain the complete simulation state in a
# single object. They can be used to start or restart a simulation.
#
# Example use cases in which a simulation may be started from a snapshot include user code that generates initial
# particle positions.
#
# \note Snapshots do not yet have a python API, so they can only be generated by C++ plugins. A future version of 
#       HOOMD-blue will allow fast access to snapshot data in python.
#
# **Example:**
# \code
# snapshot = my_system_create_routine(.. parameters ..)
# system = init.read_snapshot(snapshot)
# \endcode
#
# \sa hoomd_script.data
def read_snapshot(snapshot):
    util.print_status_line();

    # initialize GPU/CPU execution configuration and MPI early
    my_exec_conf = _create_exec_conf();

    # check if initialization has already occured
    if is_initialized():
        globals.msg.error("Cannot initialize more than once\n");
        raise RuntimeError("Error creating random polymers");

    my_domain_decomposition = _create_domain_decomposition(snapshot.global_box);

    if my_domain_decomposition is not None:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf, my_domain_decomposition);
    else:
        globals.system_definition = hoomd.SystemDefinition(snapshot, my_exec_conf);

    # initialize the system
    globals.system = hoomd.System(globals.system_definition, time_step);
    
    _perform_common_init_tasks();
    return data.system_data(globals.system_definition);

 
## Performs common initialization tasks
#
# \internal
# Initialization tasks that are performed for every simulation are to
# be done here. For example, setting up communication, registering the
# SFCPackUpdater, initializing the log writer, etc...
def _perform_common_init_tasks():
    from hoomd_script import update;
    from hoomd_script import group;
    from hoomd_script import compute;

    # create the sorter, using the evil import __main__ trick to provide the user with a default variable
    import __main__;
    __main__.sorter = update.sort();
    
    # create the default compute.thermo on the all group
    util._disable_status_lines = True;
    all = group.all();
    compute._get_unique_thermo(group=all);
    util._disable_status_lines = False;

    # set up Communicator, and register it with the System 
    if hoomd.is_MPI_available():
        cpp_decomposition = globals.system_definition.getParticleData().getDomainDecomposition();
        if cpp_decomposition is not None:
            # create the c++ Communicator
            if not globals.exec_conf.isCUDAEnabled():
                cpp_communicator = hoomd.Communicator(globals.system_definition, cpp_decomposition)
            else:
                cpp_communicator = hoomd.CommunicatorGPU(globals.system_definition, cpp_decomposition)

            # set Communicator in C++ System
            globals.system.setCommunicator(cpp_communicator)

 
## Initializes the execution configuration
#
# \internal
def _create_exec_conf():
    # use a cached execution configuration if available
    if globals.exec_conf is not None:
        return globals.exec_conf

    mpi_available = hoomd.is_MPI_available();
    
    # set the openmp thread limits
    if globals.options.ncpu is not None:
        if globals.options.ncpu > hoomd.get_num_procs():
            globals.msg.warning("Requesting more CPU cores than there are available in the system\n");
        hoomd.set_num_threads(globals.options.ncpu);

    # if no command line options were specified, create a default ExecutionConfiguration
    if globals.options.mode is None:
        if mpi_available:
            if globals.options.nrank is not None:
                exec_conf = hoomd.ExecutionConfiguration(globals.options.min_cpu, globals.options.ignore_display, globals.msg, True, globals.options.nrank);
            else:
                exec_conf = hoomd.ExecutionConfiguration(globals.options.min_cpu, globals.options.ignore_display, globals.msg,True);
        else:
            exec_conf = hoomd.ExecutionConfiguration(globals.options.min_cpu, globals.options.ignore_display, globals.msg);
    else:
        # determine the GPU on which to execute
        if globals.options.gpu is not None:
            gpu_id = int(globals.options.gpu);
        else:
            gpu_id = -1;
        
        # create the specified configuration
        if globals.options.mode == "cpu":
            if mpi_available:
                if globals.options.nrank is not None:
                    exec_conf = hoomd.ExecutionConfiguration(hoomd.ExecutionConfiguration.executionMode.CPU, gpu_id, globals.options.min_cpu, globals.options.ignore_display, globals.msg, True, globals.options.nrank); 
                else:
                    exec_conf = hoomd.ExecutionConfiguration(hoomd.ExecutionConfiguration.executionMode.CPU, gpu_id, globals.options.min_cpu, globals.options.ignore_display, globals.msg, True);
            else:
                exec_conf = hoomd.ExecutionConfiguration(hoomd.ExecutionConfiguration.executionMode.CPU, gpu_id, globals.options.min_cpu, globals.options.ignore_display, globals.msg);
        elif globals.options.mode == "gpu":
            if mpi_available:
                if globals.options.nrank is not None:
                    exec_conf = hoomd.ExecutionConfiguration(hoomd.ExecutionConfiguration.executionMode.GPU, gpu_id, globals.options.min_cpu, globals.options.ignore_display, globals.msg, True, globals.options.nrank);
                else:
                    exec_conf = hoomd.ExecutionConfiguration(hoomd.ExecutionConfiguration.executionMode.GPU, gpu_id, globals.options.min_cpu, globals.options.ignore_display, globals.msg, True);
            else:
                exec_conf = hoomd.ExecutionConfiguration(hoomd.ExecutionConfiguration.executionMode.GPU, gpu_id, globals.options.min_cpu, globals.options.ignore_display, globals.msg);
        else:
            raise RuntimeError("Error initializing");
    
    # if gpu_error_checking is set, enable it on the GPU
    if globals.options.gpu_error_checking:
       exec_conf.setCUDAErrorChecking(True);
    
    globals.exec_conf = exec_conf;

    return exec_conf;

## Create a DomainDecomposition object
# \internal 
def _create_domain_decomposition(box):
        if not hoomd.is_MPI_available():
            return None

        # default values for arguents
        nx = ny = nz = 0
        linear = False

        if globals.options.nx is not None:
            nx = globals.options.nx
        if globals.options.ny is not None:
            ny = globals.options.ny
        if globals.options.nz is not None:
            nz = globals.options.nz
        if globals.options.linear is not None:
            linear = globals.options.linear

        if linear is True:
            # set up linear decomposition
            nz = globals.exec_conf.getNRanks()
  
        # if we are only running on one processor, we use optimized code paths
        # for single-GPU execution
        if globals.exec_conf.getNRanks() == 1:
            return None

        # initialize domain decomposition
        return hoomd.DomainDecomposition(globals.exec_conf, box.getL(), nx, ny, nz);

