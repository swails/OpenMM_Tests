#!/usr/bin/env python
from __future__ import division

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

from chemistry.amber.openmmloader import OpenMMAmberParm as AmberParm

parm = AmberParm('4LYT.solv10.parm7', '4LYT.solv10.equil.rst7')

system = parm.createSystem(nonbondedCutoff=8.0*angstroms, nonbondedMethod=PME,
                           rigidWater=True)
integrator = VerletIntegrator(1.0e-10*picoseconds)

context = Context(system, integrator)
context.setPositions(parm.positions)

state = context.getState(getEnergy=True)
