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


#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4103 4244 )
#endif

#include <math.h>
#include "HOOMDDumpWriter.h"
#include "HOOMDInitializer.h"
#include "BondData.h"
#include "AngleData.h"

#include <iostream>
#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
using namespace boost::filesystem;
#include <boost/shared_ptr.hpp>
using namespace boost;

#include <fstream>
using namespace std;

//! Name the unit test module
#define BOOST_TEST_MODULE XMLReaderWriterTest
#include "boost_utf_configure.h"

/*! \file xml_reader_writer_test.cc
    \brief Unit tests for HOOMDDumpWriter and HOOMDumpReader
    \ingroup unit_tests
*/

//! Performs low level tests of HOOMDDumpWriter
BOOST_AUTO_TEST_CASE( HOOMDDumpWriterBasicTests )
    {
    InertiaTensor I;
    
    // start by creating a single particle system: see it the correct file is written
    BoxDim box(Scalar(35), Scalar(55), Scalar(125));
    
    // set some tilt factors
    box.setTiltFactors(Scalar(1.0),Scalar(0.5),Scalar(0.25));
    int n_types = 5;
    int n_bond_types = 2;
    int n_angle_types = 1;
    int n_dihedral_types = 1;
    int n_improper_types = 1;
    
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(4, box, n_types, n_bond_types, n_angle_types, n_dihedral_types, n_improper_types));
    shared_ptr<ParticleData> pdata = sysdef->getParticleData();
    
    // set recognizable values for the particle
    {
    ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(pdata->getAccelerations(), access_location::host, access_mode::readwrite);
    ArrayHandle<int3> h_image(pdata->getImages(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_body(pdata->getBodies(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_rtag(pdata->getRTags(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_charge(pdata->getCharges(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_diameter(pdata->getDiameters(), access_location::host, access_mode::readwrite);


    h_pos.data[0].x = Scalar(1.5);
    h_pos.data[0].y = Scalar(2.5);
    h_pos.data[0].z = Scalar(-5.5);
    
    h_image.data[0].x = -1;
    h_image.data[0].y = -5;
    h_image.data[0].z = 6;
    
    h_vel.data[0].x = Scalar(-1.5);
    h_vel.data[0].y = Scalar(-10.5);
    h_vel.data[0].z = Scalar(56.5);
    
    h_vel.data[0].w = Scalar(1.5); //mass
    
    h_diameter.data[0] = Scalar(3.5);
    
    h_pos.data[0].w = __int_as_scalar(3); //type

    h_body.data[0] = NO_BODY;
    
    I.set(0, 1, 2, 3, 4, 5);
    pdata->setInertiaTensor(0, I);

    h_pos.data[1].x = Scalar(1.5);
    h_pos.data[1].y = Scalar(2.5);
    h_pos.data[1].z = Scalar(-3.5);
    
    h_image.data[1].x = 10;
    h_image.data[1].y = 500;
    h_image.data[1].z = 900;
    
    h_vel.data[1].x = Scalar(-1.5);
    h_vel.data[1].y = Scalar(-10.5);
    h_vel.data[1].z = Scalar(5.5);
    
    h_vel.data[1].w = Scalar(2.5); /// mass
    
    h_diameter.data[1] = Scalar(4.5);
    
    h_pos.data[1].w = __int_as_scalar(0);
    
    h_body.data[1] = 1;

    I.set(5, 4, 3, 2, 1, 0);
    pdata->setInertiaTensor(1, I);

    h_pos.data[2].x = Scalar(-1.5);
    h_pos.data[2].y = Scalar(2.5);
    h_pos.data[2].z = Scalar(3.5);
    
    h_image.data[2].x = 10;
    h_image.data[2].y = 500;
    h_image.data[2].z = 900;
    
    h_vel.data[2].x = Scalar(-1.5);
    h_vel.data[2].y = Scalar(-10.5);
    h_vel.data[2].z = Scalar(5.5);
    
    h_vel.data[2].w = Scalar(2.5);
    
    h_diameter.data[2] = Scalar(4.5);
    
    h_pos.data[2].w = __int_as_scalar(1);
    
    h_body.data[2] = 1;
    
    I.set(1, 11, 21, 31, 41, 51);
    pdata->setInertiaTensor(2, I);

    h_pos.data[3].x = Scalar(-1.5);
    h_pos.data[3].y = Scalar(2.5);
    h_pos.data[3].z = Scalar(3.5);
    
    h_image.data[3].x = 105;
    h_image.data[3].y = 5005;
    h_image.data[3].z = 9005;
    
    h_vel.data[3].x = Scalar(-1.5);
    h_vel.data[3].y = Scalar(-10.5);
    h_vel.data[3].z = Scalar(5.5);
    
    h_vel.data[3].w = Scalar(2.5);
    
    h_diameter.data[3] = Scalar(4.5);
    
    h_pos.data[3].w = __int_as_scalar(1);
    
    h_body.data[3] = 0;
    
    I.set(51, 41, 31, 21, 11, 1);
    pdata->setInertiaTensor(3, I);
    
    }
    
    // add a couple walls for fun
    sysdef->getWallData()->addWall(Wall(1,0,0, 0,1,0));
    sysdef->getWallData()->addWall(Wall(0,1,0, 0,0,1));
    sysdef->getWallData()->addWall(Wall(0,0,1, 1,0,0));
    
    // add a few bonds too
    sysdef->getBondData()->addBond(Bond(0, 0, 1));
    sysdef->getBondData()->addBond(Bond(1, 1, 0));
    
    // and angles as well
    sysdef->getAngleData()->addAngle(Angle(0, 0, 1, 2));
    sysdef->getAngleData()->addAngle(Angle(0, 1, 2, 0));
    
    // and a dihedral
    sysdef->getDihedralData()->addDihedral(Dihedral(0, 0, 1, 2, 3));
    
    // and an improper
    sysdef->getImproperData()->addDihedral(Dihedral(0, 3, 2, 1, 0));
    
    // create the writer
    shared_ptr<HOOMDDumpWriter> writer(new HOOMDDumpWriter(sysdef, "test"));
    
    writer->setOutputPosition(false);
    
    // first test
        {
        // make sure the first output file is deleted
        remove_all("test.0000000000.xml");
        BOOST_REQUIRE(!exists("test.0000000000.xml"));
        
        // write the first output
        writer->analyze(0);
        
        // make sure the file was created
        BOOST_REQUIRE(exists("test.0000000000.xml"));
        
        // check the output line by line
        ifstream f("test.0000000000.xml");
        string line;
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<hoomd_xml version=\"1.5\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "<configuration time_step=\"0\" dimensions=\"3\" natoms=\"4\" >");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "<box lx=\"35\" ly=\"55\" lz=\"125\" xy=\"1\" xz=\"0.5\" yz=\"0.25\"/>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "</configuration>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "</hoomd_xml>");
        BOOST_REQUIRE(!f.bad());
        f.close();
        }
        
    // second test: test position
        {
        writer->setOutputPosition(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000010.xml");
        BOOST_REQUIRE(!exists("test.0000000010.xml"));
        
        // write the file
        writer->analyze(10);
        
        // make sure the file was created
        BOOST_REQUIRE(exists("test.0000000010.xml"));
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000010.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<position num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1.5 2.5 -5.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1.5 2.5 -3.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1.5 2.5 3.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1.5 2.5 3.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</position>");
        
        getline(f, line); // </configuration
        getline(f, line); // </HOOMD_xml
        f.close();
        }
        
    // third test: test velocity
        {
        writer->setOutputPosition(false);
        writer->setOutputVelocity(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000020.xml");
        BOOST_REQUIRE(!exists("test.0000000020.xml"));
        
        // write the file
        writer->analyze(20);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000020.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<velocity num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1.5 -10.5 56.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1.5 -10.5 5.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1.5 -10.5 5.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1.5 -10.5 5.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</velocity>");
        f.close();
        }
        
    // fourth test: the type array
        {
        writer->setOutputVelocity(false);
        writer->setOutputType(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000030.xml");
        BOOST_REQUIRE(!exists("test.0000000030.xml"));
        
        // write the file
        writer->analyze(30);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000030.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<type num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "D");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "A");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "B");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "B");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</type>");
        f.close();
        }
        
    // fifth test: the wall array
        {
        writer->setOutputType(false);
        writer->setOutputWall(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000040.xml");
        BOOST_REQUIRE(!exists("test.0000000040.xml"));
        
        // write the file
        writer->analyze(40);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000040.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<wall>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<coord ox=\"1\" oy=\"0\" oz=\"0\" nx=\"0\" ny=\"1\" nz=\"0\" />");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<coord ox=\"0\" oy=\"1\" oz=\"0\" nx=\"0\" ny=\"0\" nz=\"1\" />");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<coord ox=\"0\" oy=\"0\" oz=\"1\" nx=\"1\" ny=\"0\" nz=\"0\" />");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</wall>");
        f.close();
        }
        
    // sixth test: the bond array
        {
        writer->setOutputWall(false);
        writer->setOutputBond(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000050.xml");
        BOOST_REQUIRE(!exists("test.0000000050.xml"));
        
        // write the file
        writer->analyze(50);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000050.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<bond num=\"2\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "bondA 0 1");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "bondB 1 0");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</bond>");
        f.close();
        }
        
    // seventh test: the angle array
        {
        writer->setOutputBond(false);
        writer->setOutputAngle(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000060.xml");
        BOOST_REQUIRE(!exists("test.0000000060.xml"));
        
        // write the file
        writer->analyze(60);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000060.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<angle num=\"2\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "angleA 0 1 2");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "angleA 1 2 0");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</angle>");
        f.close();
        }
        
    // eighth test: test image
        {
        writer->setOutputAngle(false);
        writer->setOutputImage(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000070.xml");
        BOOST_REQUIRE(!exists("test.0000000070.xml"));
        
        // write the file
        writer->analyze(70);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000070.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<image num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1 -5 6");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "10 500 900");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "10 500 900");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "105 5005 9005");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</image>");
        f.close();
        }
        
    // nineth test: test mass
        {
        writer->setOutputImage(false);
        writer->setOutputMass(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000080.xml");
        BOOST_REQUIRE(!exists("test.0000000080.xml"));
        
        // write the file
        writer->analyze(80);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000080.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<mass num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "2.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "2.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "2.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</mass>");
        f.close();
        }
        
    // tenth test: test diameter
        {
        writer->setOutputMass(false);
        writer->setOutputDiameter(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000090.xml");
        BOOST_REQUIRE(!exists("test.0000000090.xml"));
        
        // write the file
        writer->analyze(90);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000090.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<diameter num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "3.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "4.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "4.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "4.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</diameter>");
        f.close();
        }
        
    // eleventh test: the dihedral array
        {
        writer->setOutputDiameter(false);
        writer->setOutputDihedral(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000100.xml");
        BOOST_REQUIRE(!exists("test.0000000100.xml"));
        
        // write the file
        writer->analyze(100);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000100.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<dihedral num=\"1\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "dihedralA 0 1 2 3");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</dihedral>");
        f.close();
        }
        
        
    // twelfth test: the improper array
        {
        writer->setOutputDihedral(false);
        writer->setOutputImproper(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000110.xml");
        BOOST_REQUIRE(!exists("test.0000000110.xml"));
        
        // write the file
        writer->analyze(110);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000110.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<improper num=\"1\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "dihedralA 3 2 1 0");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</improper>");
        f.close();
        }
    
    // thirteenth test: the body array
        {
        writer->setOutputImproper(false);
        writer->setOutputBody(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000120.xml");
        BOOST_REQUIRE(!exists("test.0000000120.xml"));
        
        // write the file
        writer->analyze(120);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000120.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<body num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-1");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "0");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</body>");
        f.close();
        }
        
    // fourteenth test: the moment_inertia array
        {
        writer->setOutputBody(false);
        writer->setOutputMomentInertia(true);
        
        // make sure the first output file is deleted
        remove_all("test.0000000130.xml");
        BOOST_REQUIRE(!exists("test.0000000130.xml"));
        
        // write the file
        writer->analyze(130);
        
        // assume that the first lines tested in the first case are still OK and skip them
        ifstream f("test.0000000130.xml");
        string line;
        getline(f, line); // <?xml
        getline(f, line); // <HOOMD_xml
        getline(f, line); // <Configuration
        getline(f, line); // <Box
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<moment_inertia num=\"4\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "0 1 2 3 4 5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "5 4 3 2 1 0");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1 11 21 31 41 51");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "51 41 31 21 11 1");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</moment_inertia>");
        f.close();
        }

    remove_all("test.0000000000.xml");
    remove_all("test.0000000010.xml");
    remove_all("test.0000000020.xml");
    remove_all("test.0000000030.xml");
    remove_all("test.0000000040.xml");
    remove_all("test.0000000050.xml");
    remove_all("test.0000000060.xml");
    remove_all("test.0000000070.xml");
    remove_all("test.0000000080.xml");
    remove_all("test.0000000090.xml");
    remove_all("test.0000000100.xml");
    remove_all("test.0000000110.xml");
    remove_all("test.0000000120.xml");
    remove_all("test.0000000130.xml");
    }

//! Tests the ability of HOOMDDumpWriter to handle tagged and reordered particles
BOOST_AUTO_TEST_CASE( HOOMDDumpWriter_tag_test )
    {
    // start by creating a single particle system: see it the correct file is written
    BoxDim box(Scalar(100.5), Scalar(120.5), Scalar(130.5));
    int n_types = 10;
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(6, box, n_types));
    shared_ptr<ParticleData> pdata = sysdef->getParticleData();
    
    // this is the shuffle order of the particles
    unsigned int tags[6] = { 5, 2, 3, 1, 0, 4 };
    unsigned int rtags[6] = { 4, 3, 1, 2, 5, 0 };
    
    {
    // set recognizable values for the particle
    ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(pdata->getVelocities(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_accel(pdata->getAccelerations(), access_location::host, access_mode::readwrite);
    ArrayHandle<int3> h_image(pdata->getImages(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_body(pdata->getBodies(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::readwrite);
    ArrayHandle<unsigned int> h_rtag(pdata->getRTags(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_charge(pdata->getCharges(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar> h_diameter(pdata->getDiameters(), access_location::host, access_mode::readwrite);

    for (int i = 0; i < 6; i++)
        {
        h_tag.data[i] = tags[i];
        unsigned int tag = tags[i];
        
        h_pos.data[i].x = Scalar(tag)+Scalar(0.5);
        h_pos.data[i].y = Scalar(tag)+Scalar(1.5);
        h_pos.data[i].z = Scalar(tag)+Scalar(2.5);
        
        h_image.data[i].x = tag - 10;
        h_image.data[i].y = tag - 11;
        h_image.data[i].z = tag + 50;
        
        h_vel.data[i].x = Scalar(tag)*Scalar(10.0);
        h_vel.data[i].y = Scalar(tag)*Scalar(11.0);
        h_vel.data[i].z = Scalar(tag)*Scalar(12.0);
        
        h_pos.data[i].w =__int_as_scalar(tag + 2);
        h_rtag.data[i] = rtags[i];
        }
    }
    
    // create the writer
    shared_ptr<HOOMDDumpWriter> writer(new HOOMDDumpWriter(sysdef, "test"));
    
    // write the file with all outputs enabled
    writer->setOutputPosition(true);
    writer->setOutputVelocity(true);
    writer->setOutputType(true);
    writer->setOutputImage(true);
    
    // now the big mess: check the file line by line
        {
        // make sure the first output file is deleted
        remove_all("test.0000000100.xml");
        BOOST_REQUIRE(!exists("test.0000000100.xml"));
        
        // write the first output
        writer->analyze(100);
        
        // make sure the file was created
        BOOST_REQUIRE(exists("test.0000000100.xml"));
        
        // check the output line by line
        ifstream f("test.0000000100.xml");
        string line;
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<hoomd_xml version=\"1.5\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "<configuration time_step=\"100\" dimensions=\"3\" natoms=\"6\" >");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "<box lx=\"100.5\" ly=\"120.5\" lz=\"130.5\" xy=\"0\" xz=\"0\" yz=\"0\"/>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<position num=\"6\">");
        BOOST_REQUIRE(!f.bad());
        
        // check all the positions
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "0.5 1.5 2.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "1.5 2.5 3.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "2.5 3.5 4.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "3.5 4.5 5.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "4.5 5.5 6.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "5.5 6.5 7.5");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "</position>");
        BOOST_REQUIRE(!f.bad());
        
        // check all the images
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<image num=\"6\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-10 -11 50");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-9 -10 51");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-8 -9 52");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-7 -8 53");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-6 -7 54");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "-5 -6 55");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "</image>");
        BOOST_REQUIRE(!f.bad());
        
        // check all velocities
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<velocity num=\"6\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "0 0 0");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "10 11 12");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "20 22 24");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "30 33 36");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "40 44 48");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "50 55 60");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</velocity>");
        
        // check all types
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "<type num=\"6\">");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "C");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "D");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "E");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "F");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "G");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "H");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line, "</type>");
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "</configuration>");
        BOOST_REQUIRE(!f.bad());
        
        getline(f, line);
        BOOST_CHECK_EQUAL(line,  "</hoomd_xml>");
        BOOST_REQUIRE(!f.bad());
        f.close();
        remove_all("test.0000000100.xml");
        }
    }

//! Test basic functionality of HOOMDInitializer
BOOST_AUTO_TEST_CASE( HOOMDInitializer_basic_tests )
    {
    // create a test input file
    ofstream f("test_input.xml");
    f << "<?xml version =\"1.0\" encoding =\"UTF-8\" ?>\n\
<hoomd_xml version=\"1.3\">\n\
<configuration time_step=\"150000000\" dimensions=\"2\">\n\
<box lx=\"20.05\" ly= \"32.12345\" lz=\"45.098\" xy=\".12\" xz=\".23\" yz=\".34\"/>\n\
<position >\n\
1.4 2.567890 3.45\n\
2.4 3.567890 4.45\n\
3.4 4.567890 5.45\n\
4.4 5.567890 6.45\n\
5.4 6.567890 7.45\n\
6.4 7.567890 8.45\n\
</position>\n\
<image>\n\
10 20 30\n\
11 21 31\n\
12 22 32\n\
13 23 33\n\
14 24 34\n\
15 25 35\n\
</image>\n\
<velocity>\n\
10.12 12.1567 1.056\n\
20.12 22.1567 2.056\n\
30.12 32.1567 3.056\n\
40.12 42.1567 4.056\n\
50.12 52.1567 5.056\n\
60.12 62.1567 6.056\n\
</velocity>\n\
<mass>\n\
1.0\n\
2.0\n\
3.0\n\
4.0\n\
5.0\n\
6.0\n\
</mass>\n\
<diameter>\n\
7.0\n\
8.0\n\
9.0\n\
10.0\n\
11.0\n\
12.0\n\
</diameter>\n\
<body>\n\
-1\n\
0\n\
1\n\
2\n\
3\n\
4\n\
</body>\n\
<type>\n\
5\n\
4\n\
3\n\
2\n\
1\n\
0\n\
</type>\n\
<charge>\n\
0.0\n\
10.0\n\
20.0\n\
30.0\n\
40.0\n\
50.0\n\
</charge>\n\
<moment_inertia>\n\
0 1 2 3 4 5\n\
10 11 12 13 14 15\n\
20 21 22 23 24 25\n\
30 31 32 33 34 35\n\
40 41 42 43 44 45\n\
50 51 52 53 54 55\n\
</moment_inertia>\n\
<wall>\n\
<coord ox=\"1.0\" oy=\"2.0\" oz=\"3.0\" nx=\"4.0\" ny=\"5.0\" nz=\"6.0\"/>\n\
<coord ox=\"7.0\" oy=\"8.0\" oz=\"9.0\" nx=\"10.0\" ny=\"11.0\" nz=\"-12.0\"/>\n\
</wall>\n\
<bond>\n\
bond_a 0 1\n\
bond_b 1 2\n\
bond_a 2 3\n\
bond_c 3 4\n\
</bond>\n\
<angle>\n\
angle_a 0 1 2\n\
angle_b 1 2 3\n\
angle_a 2 3 4\n\
</angle>\n\
<dihedral>\n\
di_a 0 1 2 3\n\
di_b 1 2 3 4\n\
</dihedral>\n\
<improper>\n\
im_a 3 2 1 0\n\
im_b 5 4 3 2\n\
</improper>\n\
</configuration>\n\
</hoomd_xml>" << endl;
    f.close();
    
    // now that we have created a test file, load it up into a pdata
    boost::shared_ptr<ExecutionConfiguration> exec_conf(new ExecutionConfiguration(ExecutionConfiguration::CPU));
    HOOMDInitializer init(exec_conf,"test_input.xml");
    boost::shared_ptr<SnapshotSystemData> snapshot;
    snapshot = init.getSnapshot();
    shared_ptr<SystemDefinition> sysdef(new SystemDefinition(snapshot));
    shared_ptr<ParticleData> pdata = sysdef->getParticleData();
    
    // verify all parameters
    BOOST_CHECK_EQUAL(init.getTimeStep(), (unsigned int)150000000);
    BOOST_CHECK_EQUAL(sysdef->getNDimensions(), (unsigned int)2);
    BOOST_CHECK_EQUAL(pdata->getN(), (unsigned int)6);
    BOOST_CHECK_EQUAL(pdata->getNTypes(), (unsigned int)6);
    MY_BOOST_CHECK_CLOSE(pdata->getGlobalBox().getL().x, 20.05, tol);
    MY_BOOST_CHECK_CLOSE(pdata->getGlobalBox().getL().y, 32.12345, tol);
    MY_BOOST_CHECK_CLOSE(pdata->getGlobalBox().getL().z, 45.098, tol);
    MY_BOOST_CHECK_CLOSE(pdata->getGlobalBox().getTiltFactorXY(), 0.12, tol);
    MY_BOOST_CHECK_CLOSE(pdata->getGlobalBox().getTiltFactorXZ(), 0.23, tol);
    MY_BOOST_CHECK_CLOSE(pdata->getGlobalBox().getTiltFactorYZ(), 0.34, tol);
    
    {
    ArrayHandle<Scalar4> h_pos(pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_vel(pdata->getVelocities(), access_location::host, access_mode::read);
    ArrayHandle<Scalar3> h_accel(pdata->getAccelerations(), access_location::host, access_mode::read);
    ArrayHandle<int3> h_image(pdata->getImages(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_body(pdata->getBodies(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_tag(pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_rtag(pdata->getRTags(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_charge(pdata->getCharges(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_diameter(pdata->getDiameters(), access_location::host, access_mode::read);

    for (int i = 0; i < 6; i++)
        {
        MY_BOOST_CHECK_CLOSE(h_pos.data[i].x, Scalar(i) + Scalar(1.4), tol);
        MY_BOOST_CHECK_CLOSE(h_pos.data[i].y, Scalar(i) + Scalar(2.567890), tol);
        MY_BOOST_CHECK_CLOSE(h_pos.data[i].z, Scalar(i) + Scalar(3.45), tol);
        
        BOOST_CHECK_EQUAL(h_image.data[i].x, 10 + i);
        BOOST_CHECK_EQUAL(h_image.data[i].y, 20 + i);
        BOOST_CHECK_EQUAL(h_image.data[i].z, 30 + i);
        
        MY_BOOST_CHECK_CLOSE(h_vel.data[i].x, Scalar(i+1)*Scalar(10.0) + Scalar(0.12), tol);
        MY_BOOST_CHECK_CLOSE(h_vel.data[i].y, Scalar(i+1)*Scalar(10.0) + Scalar(2.1567), tol);
        MY_BOOST_CHECK_CLOSE(h_vel.data[i].z, Scalar(i+1) + Scalar(0.056), tol);
        
        MY_BOOST_CHECK_CLOSE(h_vel.data[i].w, Scalar(i+1), tol); // mass
        
        MY_BOOST_CHECK_CLOSE(h_diameter.data[i], Scalar(i+7), tol);
        
        MY_BOOST_CHECK_CLOSE(h_charge.data[i], Scalar(i)*Scalar(10.0), tol);
        
        BOOST_CHECK_EQUAL(h_body.data[i], (unsigned int)(i-1));
        
        // checking that the type is correct becomes tricky because types are identified by
        // string
        ostringstream type_name;
        type_name << 5-i;   // the expected type is the integer 5-i
        BOOST_CHECK_EQUAL((unsigned int)__scalar_as_int(h_pos.data[i].w), pdata->getTypeByName(type_name.str()));
        BOOST_CHECK_EQUAL(h_tag.data[i], (unsigned int)i);
        BOOST_CHECK_EQUAL(h_rtag.data[i], (unsigned int)i);

        // check the moment_inertia values
        InertiaTensor I;
        I = pdata->getInertiaTensor(i);
        for (unsigned int c = 0; c < 6; c++)
            {
            MY_BOOST_CHECK_CLOSE(I.components[c], i*10 + c, tol);
            }
        }
    }
    
    // check the walls
    BOOST_REQUIRE_EQUAL(sysdef->getWallData()->getNumWalls(), (unsigned int)2);
    Wall wall1 = sysdef->getWallData()->getWall(0);
    MY_BOOST_CHECK_CLOSE(wall1.origin_x, 1.0, tol);
    MY_BOOST_CHECK_CLOSE(wall1.origin_y, 2.0, tol);
    MY_BOOST_CHECK_CLOSE(wall1.origin_z, 3.0, tol);
    // normals are made unit length when loaded, so these values differ from the ones in the file
    MY_BOOST_CHECK_CLOSE(wall1.normal_x, 0.455842306, tol);
    MY_BOOST_CHECK_CLOSE(wall1.normal_y, 0.569802882, tol);
    MY_BOOST_CHECK_CLOSE(wall1.normal_z, 0.683763459, tol);
    
    Wall wall2 = sysdef->getWallData()->getWall(1);
    MY_BOOST_CHECK_CLOSE(wall2.origin_x, 7.0, tol);
    MY_BOOST_CHECK_CLOSE(wall2.origin_y, 8.0, tol);
    MY_BOOST_CHECK_CLOSE(wall2.origin_z, 9.0, tol);
    // normals are made unit length when loaded, so these values differ from the ones in the file
    MY_BOOST_CHECK_CLOSE(wall2.normal_x, 0.523423923, tol);
    MY_BOOST_CHECK_CLOSE(wall2.normal_y, 0.575766315, tol);
    MY_BOOST_CHECK_CLOSE(wall2.normal_z, -0.628108707, tol);
    
    // check the bonds
    boost::shared_ptr<BondData> bond_data = sysdef->getBondData();
    
    // 4 bonds should have been read in
    BOOST_REQUIRE_EQUAL(bond_data->getNumBonds(), (unsigned int)4);
    
    // check that the types have been named properly
    BOOST_REQUIRE_EQUAL(bond_data->getNBondTypes(), (unsigned int)3);
    BOOST_CHECK_EQUAL(bond_data->getTypeByName("bond_a"), (unsigned int)0);
    BOOST_CHECK_EQUAL(bond_data->getTypeByName("bond_b"), (unsigned int)1);
    BOOST_CHECK_EQUAL(bond_data->getTypeByName("bond_c"), (unsigned int)2);
    
    BOOST_CHECK_EQUAL(bond_data->getNameByType(0), string("bond_a"));
    BOOST_CHECK_EQUAL(bond_data->getNameByType(1), string("bond_b"));
    BOOST_CHECK_EQUAL(bond_data->getNameByType(2), string("bond_c"));
    
    // verify each bond
    Bond b = bond_data->getBond(0);
    BOOST_CHECK_EQUAL(b.a, (unsigned int)0);
    BOOST_CHECK_EQUAL(b.b, (unsigned int)1);
    BOOST_CHECK_EQUAL(b.type, (unsigned int)0);
    
    b = bond_data->getBond(1);
    BOOST_CHECK_EQUAL(b.a, (unsigned int)1);
    BOOST_CHECK_EQUAL(b.b, (unsigned int)2);
    BOOST_CHECK_EQUAL(b.type, (unsigned int)1);
    
    b = bond_data->getBond(2);
    BOOST_CHECK_EQUAL(b.a, (unsigned int)2);
    BOOST_CHECK_EQUAL(b.b, (unsigned int)3);
    BOOST_CHECK_EQUAL(b.type, (unsigned int)0);
    
    b = bond_data->getBond(3);
    BOOST_CHECK_EQUAL(b.a, (unsigned int)3);
    BOOST_CHECK_EQUAL(b.b, (unsigned int)4);
    BOOST_CHECK_EQUAL(b.type, (unsigned int)2);
    
    // check the angles
    boost::shared_ptr<AngleData> angle_data = sysdef->getAngleData();
    
    // 3 angles should have been read in
    BOOST_REQUIRE_EQUAL(angle_data->getNumAngles(), (unsigned int)3);
    
    // check that the types have been named properly
    BOOST_REQUIRE_EQUAL(angle_data->getNAngleTypes(), (unsigned int)2);
    BOOST_CHECK_EQUAL(angle_data->getTypeByName("angle_a"), (unsigned int)0);
    BOOST_CHECK_EQUAL(angle_data->getTypeByName("angle_b"), (unsigned int)1);
    
    BOOST_CHECK_EQUAL(angle_data->getNameByType(0), string("angle_a"));
    BOOST_CHECK_EQUAL(angle_data->getNameByType(1), string("angle_b"));
    
    // verify each angle
    Angle a = angle_data->getAngle(0);
    BOOST_CHECK_EQUAL(a.a, (unsigned int)0);
    BOOST_CHECK_EQUAL(a.b, (unsigned int)1);
    BOOST_CHECK_EQUAL(a.c, (unsigned int)2);
    BOOST_CHECK_EQUAL(a.type, (unsigned int)0);
    
    a = angle_data->getAngle(1);
    BOOST_CHECK_EQUAL(a.a, (unsigned int)1);
    BOOST_CHECK_EQUAL(a.b, (unsigned int)2);
    BOOST_CHECK_EQUAL(a.c, (unsigned int)3);
    BOOST_CHECK_EQUAL(a.type, (unsigned int)1);
    
    a = angle_data->getAngle(2);
    BOOST_CHECK_EQUAL(a.a, (unsigned int)2);
    BOOST_CHECK_EQUAL(a.b, (unsigned int)3);
    BOOST_CHECK_EQUAL(a.c, (unsigned int)4);
    BOOST_CHECK_EQUAL(a.type, (unsigned int)0);
    
    // check the dihedrals
    boost::shared_ptr<DihedralData> dihedral_data = sysdef->getDihedralData();
    
    // 2 dihedrals should have been read in
    BOOST_REQUIRE_EQUAL(dihedral_data->getNumDihedrals(), (unsigned int)2);
    
    // check that the types have been named properly
    BOOST_REQUIRE_EQUAL(dihedral_data->getNDihedralTypes(), (unsigned int)2);
    BOOST_CHECK_EQUAL(dihedral_data->getTypeByName("di_a"), (unsigned int)0);
    BOOST_CHECK_EQUAL(dihedral_data->getTypeByName("di_b"), (unsigned int)1);
    
    BOOST_CHECK_EQUAL(dihedral_data->getNameByType(0), string("di_a"));
    BOOST_CHECK_EQUAL(dihedral_data->getNameByType(1), string("di_b"));
    
    // verify each dihedral
    Dihedral d = dihedral_data->getDihedral(0);
    BOOST_CHECK_EQUAL(d.a, (unsigned int)0);
    BOOST_CHECK_EQUAL(d.b, (unsigned int)1);
    BOOST_CHECK_EQUAL(d.c, (unsigned int)2);
    BOOST_CHECK_EQUAL(d.d, (unsigned int)3);
    BOOST_CHECK_EQUAL(d.type, (unsigned int)0);
    
    d = dihedral_data->getDihedral(1);
    BOOST_CHECK_EQUAL(d.a, (unsigned int)1);
    BOOST_CHECK_EQUAL(d.b, (unsigned int)2);
    BOOST_CHECK_EQUAL(d.c, (unsigned int)3);
    BOOST_CHECK_EQUAL(d.d, (unsigned int)4);
    BOOST_CHECK_EQUAL(d.type, (unsigned int)1);
    
    
    // check the impropers
    boost::shared_ptr<DihedralData> improper_data = sysdef->getImproperData();
    
    // 2 dihedrals should have been read in
    BOOST_REQUIRE_EQUAL(improper_data->getNumDihedrals(), (unsigned int)2);
    
    // check that the types have been named properly
    BOOST_REQUIRE_EQUAL(improper_data->getNDihedralTypes(), (unsigned int)2);
    BOOST_CHECK_EQUAL(improper_data->getTypeByName("im_a"), (unsigned int)0);
    BOOST_CHECK_EQUAL(improper_data->getTypeByName("im_b"), (unsigned int)1);
    
    BOOST_CHECK_EQUAL(improper_data->getNameByType(0), string("im_a"));
    BOOST_CHECK_EQUAL(improper_data->getNameByType(1), string("im_b"));
    
    // verify each dihedral
    d = improper_data->getDihedral(0);
    BOOST_CHECK_EQUAL(d.a, (unsigned int)3);
    BOOST_CHECK_EQUAL(d.b, (unsigned int)2);
    BOOST_CHECK_EQUAL(d.c, (unsigned int)1);
    BOOST_CHECK_EQUAL(d.d, (unsigned int)0);
    BOOST_CHECK_EQUAL(d.type, (unsigned int)0);
    
    d = improper_data->getDihedral(1);
    BOOST_CHECK_EQUAL(d.a, (unsigned int)5);
    BOOST_CHECK_EQUAL(d.b, (unsigned int)4);
    BOOST_CHECK_EQUAL(d.c, (unsigned int)3);
    BOOST_CHECK_EQUAL(d.d, (unsigned int)2);
    BOOST_CHECK_EQUAL(d.type, (unsigned int)1);
    
    // clean up after ourselves
    remove_all("test_input.xml");
    }

#ifdef WIN32
#pragma warning( pop )
#endif

