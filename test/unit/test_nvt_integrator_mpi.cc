//! name the boost unit test module
#define BOOST_TEST_MODULE NVTUpdaterTestsMPI
#include "MPITestSetup.h"

#include "HOOMDMath.h"
#include "ExecutionConfiguration.h"
#include "System.h"
#include "TwoStepNVT.h"
#include "IntegratorTwoStep.h"
#include "AllPairPotentials.h"
#include "NeighborListBinned.h"
#include "RandomGenerator.h"

#include <boost/python.hpp>
#include <boost/mpi.hpp>
#include <boost/shared_ptr.hpp>

#include <math.h>

#include "Communicator.h"
#include "DomainDecomposition.h"

#ifdef ENABLE_CUDA
#include "CellListGPU.h"
#include "NeighborListGPUBinned.h"
#include "TwoStepNVTGPU.h"
#include "ComputeThermoGPU.h"
#include "CommunicatorGPU.h"
#endif

using namespace boost;

void test_nvt_integrator_mpi(boost::shared_ptr<ExecutionConfiguration> exec_conf)
{
    // initialize two random particle system
    Scalar phi_p = 0.2;
    unsigned int N = 20000;
    Scalar L = pow(M_PI/6.0/phi_p*Scalar(N),1.0/3.0);
    BoxDim box_g(L);
    RandomGenerator rand_init(exec_conf, box_g, 12345);
    std::vector<string> types;
    types.push_back("A");
    std::vector<uint> bonds;
    std::vector<string> bond_types;
    rand_init.addGenerator((int)N, boost::shared_ptr<PolymerParticleGenerator>(new PolymerParticleGenerator(1.0, types, bonds, bonds, bond_types, 100)));
    rand_init.setSeparationRadius("A", .4);

    rand_init.generate();

    boost::shared_ptr<SnapshotSystemData> snap;
    snap = rand_init.getSnapshot();

    boost::shared_ptr<DomainDecomposition> decomposition(new DomainDecomposition(exec_conf,snap->global_box.getL(), 0));

    shared_ptr<SystemDefinition> sysdef_1(new SystemDefinition(snap, exec_conf,decomposition));
    shared_ptr<SystemDefinition> sysdef_2(new SystemDefinition(snap, exec_conf));

    shared_ptr<ParticleData> pdata_1 = sysdef_1->getParticleData();
    shared_ptr<ParticleData> pdata_2 = sysdef_2->getParticleData();
 
    boost::shared_ptr<Communicator> comm;
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        comm = shared_ptr<Communicator>(new CommunicatorGPU(sysdef_1, decomposition));
    else
#endif
        comm = boost::shared_ptr<Communicator>(new Communicator(sysdef_1,decomposition));

    shared_ptr<ParticleSelector> selector_all_1(new ParticleSelectorTag(sysdef_1, 0, pdata_1->getNGlobal()-1));
    shared_ptr<ParticleGroup> group_all_1(new ParticleGroup(sysdef_1, selector_all_1));

    shared_ptr<ParticleSelector> selector_all_2(new ParticleSelectorTag(sysdef_2, 0, pdata_2->getNGlobal()-1));
    shared_ptr<ParticleGroup> group_all_2(new ParticleGroup(sysdef_2, selector_all_2));

    Scalar r_cut = Scalar(3.0);
    Scalar r_buff = Scalar(0.8);
    shared_ptr<NeighborList> nlist_1(new NeighborListBinned(sysdef_1, r_cut, r_buff));
    shared_ptr<NeighborList> nlist_2(new NeighborListBinned(sysdef_2, r_cut, r_buff));

    nlist_1->setStorageMode(NeighborList::full);
    nlist_2->setStorageMode(NeighborList::full);

    nlist_1->setCommunicator(comm);

    shared_ptr<PotentialPairLJ> fc_1 = shared_ptr<PotentialPairLJ>(new PotentialPairLJ(sysdef_1, nlist_1));
    shared_ptr<PotentialPairLJ> fc_2 = shared_ptr<PotentialPairLJ>(new PotentialPairLJ(sysdef_2, nlist_2));

    fc_1->setRcut(0, 0, r_cut);
    fc_2->setRcut(0, 0, r_cut);

    // setup some values for alpha and sigma
    Scalar epsilon = Scalar(1.0);
    Scalar sigma = Scalar(1.0);
    Scalar alpha = Scalar(1.0);
    Scalar lj1 = Scalar(4.0) * epsilon * pow(sigma,Scalar(12.0));
    Scalar lj2 = alpha * Scalar(4.0) * epsilon * pow(sigma,Scalar(6.0));
    fc_1->setParams(0,0,make_scalar2(lj1,lj2));
    fc_2->setParams(0,0,make_scalar2(lj1,lj2));

    Scalar deltaT = Scalar(0.005);
    Scalar Q = Scalar(2.0);
    Scalar T = Scalar(1.5/3.0);
    Scalar tau = sqrt(Q / (Scalar(3.0) * T));
    shared_ptr<VariantConst> T_variant_1(new VariantConst(T));
    shared_ptr<VariantConst> T_variant_2(new VariantConst(T));
    shared_ptr<IntegratorTwoStep> nvt_1(new IntegratorTwoStep(sysdef_1, deltaT));
    shared_ptr<IntegratorTwoStep> nvt_2(new IntegratorTwoStep(sysdef_2, deltaT));
    shared_ptr<ComputeThermo> thermo_1 = shared_ptr<ComputeThermo>(new ComputeThermo(sysdef_1,group_all_1));
    thermo_1->setCommunicator(comm);
    shared_ptr<ComputeThermo> thermo_2 = shared_ptr<ComputeThermo>(new ComputeThermo(sysdef_2,group_all_2));

    shared_ptr<TwoStepNVT> two_step_nvt_1;
    shared_ptr<TwoStepNVT> two_step_nvt_2;
#ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        {
        two_step_nvt_1 = boost::shared_ptr<TwoStepNVT>(new TwoStepNVTGPU(sysdef_1, group_all_1, thermo_1, tau, T_variant_1));
        two_step_nvt_2 = boost::shared_ptr<TwoStepNVT>(new TwoStepNVTGPU(sysdef_2, group_all_2, thermo_2, tau, T_variant_2));
        }
    else
#endif
    two_step_nvt_1 = boost::shared_ptr<TwoStepNVT>(new TwoStepNVT(sysdef_1, group_all_1, thermo_1, tau, T_variant_1));
    two_step_nvt_2 = boost::shared_ptr<TwoStepNVT>(new TwoStepNVT(sysdef_2, group_all_2, thermo_2, tau, T_variant_2));

    nvt_1->addIntegrationMethod(two_step_nvt_1);
    nvt_1->addForceCompute(fc_1);
    nvt_1->setCommunicator(comm);

    nvt_2->addIntegrationMethod(two_step_nvt_2);
    nvt_2->addForceCompute(fc_2);

    unsigned int ndof = nvt_1->getNDOF(group_all_1);
    BOOST_CHECK_EQUAL(ndof, nvt_2->getNDOF(group_all_2));
    thermo_1->setNDOF(ndof);
    thermo_2->setNDOF(ndof);

    nvt_1->prepRun(0);
    nvt_2->prepRun(0);

    SnapshotParticleData snap_1(N);
    SnapshotParticleData snap_2(N);
    for (int i=0; i< 100; i++)
        {
        // compare temperatures
        if (exec_conf->getRank() == 0)
            BOOST_CHECK_CLOSE(thermo_1->getTemperature(), thermo_2->getTemperature(), tol_small);
      
//       if (world->rank() ==0)
//           std::cout << "step " << i << std::endl;
        Scalar rough_tol = 15.0;

        // in the first five steps, compare all accelerations and velocities
        // beyond this number of steps, trajectories will generally diverge, since they are chaotic
        // compare the snapshot of the parallel simulation
        if (i < 5)
            {
            pdata_1->takeSnapshot(snap_1);
            // ... against the serial simulation
            pdata_2->takeSnapshot(snap_2);

            if (exec_conf->getRank() == 0)
                {
                // check position, velocity and acceleration
                for (unsigned int j = 0; j < N; j++)
                    {
                    // we do not check positions (or we would need to pull back vectors over the boundaries)
                    //MY_BOOST_CHECK_CLOSE(snap_1.pos[j].x, snap_2.pos[j].x, rough_tol);
                    //MY_BOOST_CHECK_CLOSE(snap_1.pos[j].y, snap_2.pos[j].y, rough_tol);
                    //MY_BOOST_CHECK_CLOSE(snap_1.pos[j].z, snap_2.pos[j].z, rough_tol);

                    MY_BOOST_CHECK_CLOSE(snap_1.vel[j].x, snap_2.vel[j].x, rough_tol);
                    MY_BOOST_CHECK_CLOSE(snap_1.vel[j].y, snap_2.vel[j].y, rough_tol);
                    MY_BOOST_CHECK_CLOSE(snap_1.vel[j].z, snap_2.vel[j].z, rough_tol);
     
                    MY_BOOST_CHECK_CLOSE(snap_1.accel[j].x, snap_2.accel[j].x, rough_tol);
                    MY_BOOST_CHECK_CLOSE(snap_1.accel[j].y, snap_2.accel[j].y, rough_tol);
                    MY_BOOST_CHECK_CLOSE(snap_1.accel[j].z, snap_2.accel[j].z, rough_tol);
                    }
                } 
            }
        nvt_1->update(i);
        nvt_2->update(i);
        }

}

//! Tests MPI domain decomposition with NVT integrator
BOOST_AUTO_TEST_CASE( DomainDecomposition_NVT_test )
    {
    test_nvt_integrator_mpi(exec_conf_cpu);
    }

#ifdef ENABLE_CUDA
//! Tests MPI domain decomposition with NVT integrator on the GPU
BOOST_AUTO_TEST_CASE( DomainDecomposition_NVT_test_GPU )
    {
    test_nvt_integrator_mpi(exec_conf_gpu);
    }
#endif
