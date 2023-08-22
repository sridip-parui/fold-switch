#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, master_runner
from meld import comm, vault
from meld import system
from meld import parse
import meld.system.montecarlo as mc
from meld.system.restraints import LinearRamp, ConstantRamp
from collections import namedtuple
import glob as glob
from typing import List, Optional, TextIO, NewType
from collections import namedtuple
from meld.system import System
from meld.system.restraints import (
    RestraintGroup,
    Restraint,
    TorsionRestraint,
    DistanceRestraint,
    RdcRestraint,
    TimeRamp,
    ConstantRamp,
    RestraintScaler,
)
from meld.system.patchers import RdcAlignmentPatcher

#Inputs for simulation
N_REPLICAS = 32
N_STEPS = 50000
BLOCK_SIZE = 100
#Inputs for restraints
contact = [ 'contact_model1.dat', 'contact_model2.dat' ]  # files for contacts, here model1, model2 correspond to alpha and beta respectively
rest_activation = 0.98 # fraction of restraints to be activated (Please note force constant also depends on frac_activation)
k_model1 = 1700.0 # force constant for alpha (kb*gamma gamma=Eb/Ea)
k_model2 = 250.0 # force constant for beta

def file_length(filename):
        with open(filename) as f:
                for i, l in enumerate(f):
                        pass
        return i + 1

def get_CA_restraints_CTD(filename, s, ramp, force_constant, scaler, frac_activation):
    r1 = 0.0
    r2 = 0.0
    r3 = 0.8
    r4 = 1.0
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    rest_group = []
    dists = []
    for line in lines:
            cols = line.split()
            i = int(cols[0])
            j = int(cols[1])
            k = force_constant
            d = DistanceRestraint(s, scaler, ramp, i, "CA", j, "CA", r1,
                 r2, r3, r4, k)

            rest_group.append(d) 

    if file_length(contact[0])<file_length(contact[1]):
        active=int(file_length(contact[0]) * frac_activation)
    else:
        active=int(file_length(contact[1]) * frac_activation)

    print ('CA restraints', len(rest_group), 'enforce', active)
    dists.append(s.restraints.create_restraint_group(rest_group, active))
    return dists



def get_CA_restraints_NTD(filename, s, scaler):
    dists = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        cols = line.split()
        i = int(cols[0])
        j = int(cols[1])
        r = (float(cols[2])/10.0)
#        r2= ((float(cols[2])/10.0) - 0.1)
#        r3= ((float(cols[2])/10.0) + 0.1)
#        r4= r3+0.2

        rest = s.restraints.create_restraint('distance', scaler, LinearRamp(0, 100, 0, 1),
                                             r1=0.0, r2=r-0.10, r3=r+0.10, r4=r+0.30, k=250,
                                             atom_1_res_index=i, atom_2_res_index=j,
                                             atom_1_name='CA', atom_2_name='CA')
        dists.append(s.restraints.create_restraint_group([rest], 1))
    return dists

def gen_state_templates(index,templates):
    n_templates = len(templates)
    print (index,n_templates,((n_templates-1)-(index%n_templates)))
    a = system.ProteinMoleculeFromPdbFile(templates[(n_templates-1)-(index%n_templates)])
    # Note that it does not matter which ff we use here to build as that info is not passed on
    b = system.SystemBuilder(forcefield="ff14sbside")
    c = b.build_system_from_molecules([a])
    pos = c._coordinates
    c._box_vectors=np.array([0.,0.,0.])
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos,vel,alpha,energy,c._box_vectors)


def setup_system():
    # load the sequence
    templates = glob.glob('TEMPLATES/*.pdb')
    top_files = glob.glob('TEMPLATES/*.top')
    crd_files = glob.glob('TEMPLATES/*.crd')
    # build the system
    s = system.builder.load_amber_system(top_files[0], crd_files[0])
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 0.8, 300., 500.)
    n_res = s.residue_numbers[-1]


# CA restraints for CTD
    distance_scaler_CTD = s.restraints.create_scaler('nonlinear', alpha_min=0.2, alpha_max=0.8, factor=4.0)
    dis_CTD = get_CA_restraints_CTD(contact[0], s, ramp=LinearRamp(0, 100, 0, 1), force_constant=k_model1,
                        scaler=distance_scaler_CTD, frac_activation=rest_activation)
    dis_CTD.extend(get_CA_restraints_CTD(contact[1], s, ramp=LinearRamp(0, 100, 0, 1), force_constant=k_model2,
                             scaler=distance_scaler_CTD, frac_activation=rest_activation))
    print ('dis_group', len(dis_CTD))
    s.restraints.add_selectively_active_collection(dis_CTD,1)
   
# CA restraints for NTD
    distance_scaler_NTD = s.restraints.create_scaler('constant')
    dis_NTD = get_CA_restraints_NTD(
        'contact_NTD.dat', s, distance_scaler_NTD)
    n_distances_keep = int(len(dis_NTD) * 1.0)
    print ('NTD', len(dis_NTD), 'enforce', n_distances_keep)
    s.restraints.add_selectively_active_collection(
        dis_NTD, n_distances_keep)

    # setup mcmc at startup
    movers = []
    n_atoms = s.n_atoms
    for i in range(1, n_res + 1):
        n = s.index_of_atom(i, 'N') - 1
        ca = s.index_of_atom(i, 'CA') - 1
        c = s.index_of_atom(i, 'C') - 1

        mover = mc.DoubleTorsionMover(n, ca, list(range(ca, n_atoms)),
                                      ca, c, list(range(c, n_atoms)))

        movers.append((mover, 1))

    sched = mc.MonteCarloScheduler(movers, n_res * 60)

    # create the options
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck'
    options.use_big_timestep = False
    options.use_bigger_timestep = True
    options.cutoff = 1.8

    options.use_amap = False
    options.amap_alpha_bias = 0.0
    options.amap_beta_bias = 0.0
    options.timesteps = 11111
    options.minimize_steps = 20000
    #print ("Heads up! using MC minimizer!")
    print ("!Do not use MC Minimizer!")
    #options.min_mc = sched
    options.min_mc = None
    options.run_mc = None

    # create a store
    store = vault.DataStore(s.n_atoms, N_REPLICAS,
                            s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=100)
    policy = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(
        n_replicas=N_REPLICAS, adaptation_policy=policy)

    remd_runner = master_runner.MasterReplicaExchangeRunner(
        N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
#    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    states = [gen_state_templates(i,templates) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()
