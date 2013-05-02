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

/*! \file Profiler.h
    \brief Declares the Profiler class
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include "ClockSource.h"
#include "ExecutionConfiguration.h"

#ifdef ENABLE_CUDA
#include <cuda_runtime.h>
#include <boost/bind.hpp>
#endif

#include <string>
#include <stack>
#include <map>
#include <iostream>
#include <cassert>

#ifdef VTRACE
#include <vt_user.h>
#endif

#ifndef __PROFILER_H__
#define __PROFILER_H__

// forward declarations
class ProfileDataElem;
class Profiler;

/*! \ingroup hoomd_lib
    @{
*/

/*! \defgroup utils Utility classes
    \brief General purpose utility classes
*/

/*! @}
*/

//! Internal class for storing profile data
/*! This is a simple utility class, so it is fully public. It is really only designed to be used in
    concert with the Profiler class.
    \ingroup utils
*/
class ProfileDataElem
    {
    public:
        //! Constructs an element with zeroed counters
        ProfileDataElem() : m_start_time(0), m_elapsed_time(0), m_flop_count(0), m_mem_byte_count(0)
            {}
            
        //! Returns the total elapsed time of this nodes children
        int64_t getChildElapsedTime() const;
        //! Returns the total flop count of this node + children
        int64_t getTotalFlopCount() const;
        //! Returns the total memory byte count of this node + children
        int64_t getTotalMemByteCount() const;
        
        //! Output helper function
        void output(std::ostream &o, const std::string &name, int tab_level, int64_t total_time, int name_width) const;
        //! Another output helper function
        void output_line(std::ostream &o,
                         const std::string &name,
                         double sec,
                         double perc,
                         double flops,
                         double bytes,
                         unsigned int name_width) const;
        
        std::map<std::string, ProfileDataElem> m_children; //!< Child nodes of this profile
        
        int64_t m_start_time;   //!< The start time of the most recent timed event
        int64_t m_elapsed_time; //!< A running total of elapsed running time
        int64_t m_flop_count;   //!< A running total of floating point operations
        int64_t m_mem_byte_count;   //!< A running total of memory bytes transferred

#ifdef VTRACE
        unsigned int m_vt_id;   //!< Vampir Trace identifier
#endif
    };



//! A class for doing coarse-level profiling of code
/*! Stores and organizes a tree of profiles that can be created with a simple push/pop
    type interface. Any number of root profiles can be created via the default constructor
    Profiler::Profiler(). Sub-profiles are created with the push() member. They take a time sample on creation
    and take a second time sample when they are pop() ed. One can perform a guesstimate
    on memory bandwidth and FLOPS by calling the appropriate pop() member that takes
    a number of operations executed as a parameter.

    Pushing and popping the same profile tree over and over is important to do,
    the system will tally up total time in each slot. Pushing and popping different
    names on each pass will generate a jumbled mess.

    There are verions of push() and pop() that take in a reference to an ExecutionConfiguration.
    These methods automatically syncrhonize with the asynchronous GPU execution stream in order
    to provide accurate timing information.

    These profiles can of course be output via normal ostream operators.
    \ingroup utils
    */
class Profiler
    {
    public:
        //! Constructs an empty profiler and starts its timer ticking
        Profiler(const std::string& name = "Profile");
        //! Pushes a new sub-category into the current category
        void push(const std::string& name);
        //! Pops back up to the next super-category
        void pop(uint64_t flop_count = 0, uint64_t byte_count = 0);
        
        //! Pushes a new sub-category into the current category & syncs the GPUs
        void push(boost::shared_ptr<const ExecutionConfiguration> exec_conf, const std::string& name);
        //! Pops back up to the next super-category & syncs the GPUs
        void pop(boost::shared_ptr<const ExecutionConfiguration> exec_conf, uint64_t flop_count = 0, uint64_t byte_count = 0);
        
    private:
        ClockSource m_clk;  //!< Clock to provide timing information
        std::string m_name; //!< The name of this profile
        ProfileDataElem m_root; //!< The root profile element
        std::stack<ProfileDataElem *> m_stack;  //!< A stack of data elements for the push/pop structure
        
        //! Output helper function
        void output(std::ostream &o);
        
        //! friend operator to enable stream output
        friend std::ostream& operator<<(std::ostream &o, Profiler& prof);
    };

//! Exports the Profiler class to python
void export_Profiler();

//! Output operator for Profiler
std::ostream& operator<<(std::ostream &o, Profiler& prof);

/////////////////////////////////////
// Profiler inlines

inline void Profiler::push(boost::shared_ptr<const ExecutionConfiguration> exec_conf, const std::string& name)
    {
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        cudaThreadSynchronize();
#endif
    push(name);
   }

inline void Profiler::pop(boost::shared_ptr<const ExecutionConfiguration> exec_conf, uint64_t flop_count, uint64_t byte_count)
    {
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        cudaThreadSynchronize();
#endif
    pop(flop_count, byte_count);
    }

inline void Profiler::push(const std::string& name)
    {
    // sanity checks
    assert(!m_stack.empty());
    
    // pushing a new record on to the stack involves taking a time sample
    int64_t t = m_clk.getTime();
    
    ProfileDataElem *cur = m_stack.top();
    
    // then creating (or accessing) the named sample and setting the start time
    cur->m_children[name].m_start_time = t;
    
    // and updating the stack
    m_stack.push(&cur->m_children[name]);

#ifdef VTRACE
    // log VampirTrace event
    unsigned int id;
    id = VT_USER_DEF(name.c_str());
    VT_USER_START_ID(id);
    cur->m_children[name].m_vt_id = id;
#endif
    }

inline void Profiler::pop(uint64_t flop_count, uint64_t byte_count)
    {
    // sanity checks
    assert(!m_stack.empty());
    assert(!(m_stack.top() == &m_root));
    
    // popping up a level in the profile stack involves taking a time sample
    int64_t t = m_clk.getTime();
    
    // then increasing the elapsed time for the current item
    ProfileDataElem *cur = m_stack.top();
#ifdef VTRACE
    VT_USER_END_ID(cur->m_vt_id);
#endif
    cur->m_elapsed_time += t - cur->m_start_time;
    
    // and increasing the flop and mem counters
    cur->m_flop_count += flop_count;
    cur->m_mem_byte_count += byte_count;
    
    // and finally popping the stack so that the next pop will access the correct element
    m_stack.pop();
    }

#endif

