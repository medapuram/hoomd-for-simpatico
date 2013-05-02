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

/*! \file IMDInterface.h
    \brief Declares the IMDInterface class
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include <boost/shared_ptr.hpp>

#include "Analyzer.h"
#include "ConstForceCompute.h"

#ifndef __IMD_INTERFACE_H__
#define __IMD_INTERFACE_H__

//! Iterfaces with VMD through the IMD communcations port
/*! analyze() can be called very often. When not connected to
    VMD, it will do nothing. After a connection has been established,
    which can only happen during a call to analyze(), further calls will
    transmit particle positions to VMD.

    In its current implementation, only a barebones set of commands are
    supported. The sending of any command that is not understood will
    result in the socket closing the connection.
    \ingroup analyzers
*/
class IMDInterface : public Analyzer
    {
    public:
        //! Constructor
        IMDInterface(boost::shared_ptr<SystemDefinition> sysdef,
                     int port = 54321,
                     bool pause = false,
                     unsigned int rate=1,
                     boost::shared_ptr<ConstForceCompute> force = boost::shared_ptr<ConstForceCompute>(),
                     float force_scale=1.0);
        
        //! Destructor
        ~IMDInterface();
        
        //! Handle connection requests and send current positions if connected
        void analyze(unsigned int timestep);
    private:
        void *m_listen_sock;    //!< Socket we are listening on
        void *m_connected_sock; //!< Socket to transmit/receive data
        float *m_tmp_coords;    //!< Temporary holding location for coordinate data
        
        bool m_active;          //!< True if we have received a go command
        bool m_paused;          //!< True if we are paused
        unsigned int m_trate;   //!< Transmission rate
        unsigned int m_count;   //!< Count the number of times analyze() is called (used with trate)

        bool m_is_initialized;  //!< True if the interface has been initialized
        int m_port;             //!< Port to listen on
        
        boost::shared_ptr<ConstForceCompute> m_force;   //!< Force for applying IMD forces
        float m_force_scale;                            //!< Factor by which to scale all IMD forces
        
        //! Helper function that reads message headers and dispatches them to the relevant process functions
        void dispatch();
        //! Helper function to determine of messages are still available
        bool messagesAvailable();
        //! Process the IMD_DISCONNECT message
        void processIMD_DISCONNECT();
        //! Process the IMD_GO message
        void processIMD_GO();
        //! Process the IMD_KILL message
        void processIMD_KILL();
        //! Process the IMD_MDCOMM message
        void processIMD_MDCOMM(unsigned int n);
        //! Process the IMD_TRATE message
        void processIMD_TRATE(int rate);
        //! Process the IMD_PAUSE message
        void processIMD_PAUSE();
        //! Process the IMD_IOERROR message
        void processIMD_IOERROR();
        //! Process a dead connection
        void processDeadConnection();
        
        //! Helper function to establish a connection
        void establishConnectionAttempt();
        //! Helper function to send current data to VMD
        void sendCoords(unsigned int timestep);

        //! Initialize socket and internal state variables for communication
        void initConnection();
    };

//! Exports the IMDInterface class to python
void export_IMDInterface();

#endif

