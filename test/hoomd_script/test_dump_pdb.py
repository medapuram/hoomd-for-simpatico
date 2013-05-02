# -*- coding: iso-8859-1 -*-
# Maintainer: joaander

from hoomd_script import *
import unittest
import os

# unit tests for dump.pdb
class dmp_pdb_tests (unittest.TestCase):
    def setUp(self):
        print
        init.create_random(N=100, phi_p=0.05);
        import __main__;
        __main__.sorter.set_params(grid=8)

    # tests basic creation of the dump
    def test(self):
        dump.pdb(filename="dump_pdb", period=100);
        run(101)
        os.remove("dump_pdb.0000000000.pdb")
    
    # tests variable periods
    def test_variable(self):
        dump.pdb(filename="dump_pdb", period=lambda n: n*100);
        run(101);
        os.remove("dump_pdb.0000000000.pdb")
    
    def tearDown(self):
        init.reset();

if __name__ == '__main__':
    unittest.main(argv = ['test.py', '-v'])

