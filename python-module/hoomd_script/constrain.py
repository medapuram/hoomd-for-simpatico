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

## \package hoomd_script.constrain
# \brief Commands that create constraint forces on particles
#
# Constraint forces %constrain a given set of particle to a given surface, to have some relative orientation,
# or impose some other type of constraint. For example, a group of particles can be constrained to the surface of a
# sphere with constrain.sphere.
#
# As with other force commands in hoomd_script, multiple constrain commands can be issued to specify multiple
# constraints, which are additively applied. Note, however, that not all constraints specified in this manner will
# be valid if two separate constrain commands operate on the same particles.
#
# The degrees of freedom removed from the system by constraints are correctly taken into account when computing the
# temperature for thermostatting and logging.
#

from hoomd_script import globals;
from hoomd_script import force;
import hoomd;
from hoomd_script import util;
from hoomd_script import init;
from hoomd_script import data;

## \internal
# \brief Base class for constraint forces
#
# A constraint_force in hoomd_script reflects a ForceConstraint in c++. It is responsible
# for all high-level management that happens behind the scenes for hoomd_script
# writers. 1) The instance of the c++ constraint force itself is tracked and added to the
# System 2) methods are provided for disabling the force from being added to the
# net force on each particle
class _constraint_force:
    ## \internal
    # \brief Constructs the constraint force
    #
    # \param name name of the constraint force instance 
    #
    # Initializes the cpp_force to None.
    # If specified, assigns a name to the instance
    # Assigns a name to the force in force_name;
    def __init__(self):
        # check if initialization has occured
        if not init.is_initialized():
            globals.msg.error("Cannot create force before initialization\n");
            raise RuntimeError('Error creating constraint force');
        
        self.cpp_force = None;

        # increment the id counter
        id = _constraint_force.cur_id;
        _constraint_force.cur_id += 1;
        
        self.force_name = "constraint_force%d" % (id);
        self.enabled = True;
        globals.constraint_forces.append(self);
        
        # create force data iterator
        self.forces = data.force_data(self);

    ## \var enabled
    # \internal
    # \brief True if the force is enabled

    ## \var cpp_force
    # \internal
    # \brief Stores the C++ side ForceCompute managed by this class
    
    ## \var force_name
    # \internal
    # \brief The Force's name as it is assigned to the System

    ## \internal
    # \brief Checks that proper initialization has completed
    def check_initialization(self):
        # check that we have been initialized properly
        if self.cpp_force is None:
            globals.msg.error('Bug in hoomd_script: cpp_force not set, please report\n');
            raise RuntimeError();
        

    ## Disables the force
    #
    # \b Examples:
    # \code
    # force.disable()
    # \endcode
    #
    # Executing the disable command will remove the force from the simulation.
    # Any run() command executed after disabling a force will not calculate or 
    # use the force during the simulation. A disabled force can be re-enabled
    # with enable()
    #
    # To use this command, you must have saved the force in a variable, as 
    # shown in this example:
    # \code
    # force = constrain.some_force()
    # # ... later in the script
    # force.disable()
    # \endcode
    def disable(self):
        util.print_status_line();
        self.check_initialization();
            
        # check if we are already disabled
        if not self.enabled:
            globals.msg.warning("Ignoring command to disable a force that is already disabled");
            return;
        
        self.enabled = False;
        
        # remove the compute from the system
        globals.system.removeCompute(self.force_name);

    ## Benchmarks the force computation
    # \param n Number of iterations to average the benchmark over
    #
    # \b Examples:
    # \code
    # t = force.benchmark(n = 100)
    # \endcode
    #
    # The value returned by benchmark() is the average time to perform the force 
    # computation, in milliseconds. The benchmark is performed by taking the current
    # positions of all particles in the simulation and repeatedly calculating the forces
    # on them. Thus, you can benchmark different situations as you need to by simply 
    # running a simulation to achieve the desired state before running benchmark().
    #
    # \note
    # There is, however, one subtle side effect. If the benchmark() command is run 
    # directly after the particle data is initialized with an init command, then the 
    # results of the benchmark will not be typical of the time needed during the actual
    # simulation. Particles are not reordered to improve cache performance until at least
    # one time step is performed. Executing run(1) before the benchmark will solve this problem.
    #
    # To use this command, you must have saved the force in a variable, as
    # shown in this example:
    # \code
    # force = pair.some_force()
    # # ... later in the script
    # t = force.benchmark(n = 100)
    # \endcode
    def benchmark(self, n):
        self.check_initialization();
        
        # run the benchmark
        return self.cpp_force.benchmark(int(n))

    ## Enables the force
    #
    # \b Examples:
    # \code
    # force.enable()
    # \endcode
    #
    # See disable() for a detailed description.
    def enable(self):
        util.print_status_line();
        self.check_initialization();
            
        # check if we are already disabled
        if self.enabled:
            globals.msg.warning("Ignoring command to enable a force that is already enabled");
            return;
        
        # add the compute back to the system
        globals.system.addCompute(self.cpp_force, self.force_name);
        
        self.enabled = True;
        
# set default counter
_constraint_force.cur_id = 0;


## Constrain particles to the surface of a sphere
#
# The command constrain.sphere specifies that forces will be applied to all particles in the given group to constrain
# them to a sphere.
#
class sphere(_constraint_force):
    ## Specify the %sphere constraint %force
    #
    # \param group Group on which to apply the constraint
    # \param P (x,y,z) tuple indicating the position of the center of the sphere (in distance units)
    # \param r Radius of the sphere (in distance units)
    #
    # \b Examples:
    # \code
    # constrain.sphere(group=groupA, P=(0,10,2), r=10)
    # \endcode
    def __init__(self, group, P, r):
        util.print_status_line();
        
        # initialize the base class
        _constraint_force.__init__(self);
        
        # create the c++ mirror class
        P = hoomd.make_scalar3(P[0], P[1], P[2]);
        if not globals.exec_conf.isCUDAEnabled():
            self.cpp_force = hoomd.ConstraintSphere(globals.system_definition, group.cpp_group, P, r);
        else:
            self.cpp_force = hoomd.ConstraintSphereGPU(globals.system_definition, group.cpp_group, P, r);

        globals.system.addCompute(self.cpp_force, self.force_name);


