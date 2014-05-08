"""
Tests various options for HEWL systems and compares their results to what Amber
computes.
"""
from copy import copy

from utils import get_fn, colorize_error, red, green, colorize_list

# Now import the ParmEd and OpenMM utilities
from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
try:
    from chemistry.amber.openmmloader import OpenMMAmberParm as AmberParm
except ImportError:
    raise ImportError('Could not import ParmEd modules. Make sure you have '
                      'installed ParmEd!')
try:
    from scipy.io.netcdf import netcdf_file
except ImportError:
    raise ImportError('Could not import NetCDF functionality from scipy!')

class TestPME(object):
    """
    Tests the PME implementation on HEWL, a large-ish protein solvated in water
    """
    def __init__(self):
        parmname = get_fn('4LYT.solv10.parm7')
        rst7name = get_fn('4LYT.solv10.equil.rst7')
        self.parm = AmberParm(parmname, rst7name)
        self.system = self.parm.createSystem(nonbondedCutoff=8.0*u.angstrom,
                                             nonbondedMethod=app.PME)
        # Also load objects using the app layer
        self.parmapp = app.AmberPrmtopFile(parmname)
        self.crdapp = crd = app.AmberInpcrdFile(rst7name, loadVelocities=True,
                                                loadBoxVectors=True)
        self.systemapp = self.parmapp.createSystem(
                                nonbondedCutoff=8.0*u.angstroms,
                                nonbondedMethod=app.PME
        )
        self.systemapp.setDefaultPeriodicBoxVectors(*crd.getBoxVectors())
    
    def run(self, platform, precision=None, devices=None,
            use_dispersion_correction=True):
        """
        Runs the test on the given platform with the given precision model.

        Parameters
        ----------
        platform : str
            Name of the OpenMM platform to use (CPU, Reference, CUDA, or OpenCL)

        precision : str
            Precision model to use for CUDA or OpenCL (single, double, or mixed)

        devices : int or tuple of ints
            Which GPUs to run on (a tuple will run in parallel)
        """
        if not platform in ('CPU', 'Reference', 'CUDA', 'OpenCL'):
            raise ValueError('Platform %s not recognized.' % platform)
        # Define our platform and our platform properties
        plat = mm.Platform.getPlatformByName(platform)
        properties = None
        if platform == 'CUDA':
            if not precision in ('mixed', 'double', 'single'):
                raise ValueError('You must set the precision to single, '
                                 'double, or mixed for the CUDA platform.')
            properties = dict(CudaPrecision=precision)
            if devices is not None:
                properties['CudaDeviceIndex'] = str(devices)
        elif platform == 'OpenCL':
            if not precision in ('mixed', 'double', 'single'):
                raise ValueError('You must set the precision to single, '
                                 'double, or mixed for the CUDA platform.')
            properties = dict(OpenCLPrecision=precision)
            if devices is not None:
                properties['OpenCLDeviceIndex'] = str(devices)

        # Create a new system with no charges so we can compare vdW and EEL
        # energies to Amber
        parmcopy = copy(self.parm)
        for i in range(len(parmcopy.parm_data['CHARGE'])):
            parmcopy.parm_data['CHARGE'][i] = 0.0
        system = parmcopy.createSystem(nonbondedCutoff=8.0*u.angstroms,
                                       nonbondedMethod=app.PME)
        system.setDefaultPeriodicBoxVectors(*self.parm.box_vectors)

        # Test serialization
        xmlsys = mm.XmlSerializer.deserialize(
                        mm.XmlSerializer.serialize(self.system)
        )

        # Loop through all systems and turn on or off the dispersion correction
        print 'Trying to set PME parameters...',
        succeeded = None
        for sysmod in (self.system, self.systemapp, system, xmlsys):
            for force in sysmod.getForces():
                if isinstance(force, mm.NonbondedForce):
                    force.setUseDispersionCorrection(use_dispersion_correction)
                    # See if we can set the PME parameters
                    try:
                        force.setPMEParameters(3.285326106/u.nanometers,
                                               60, 64, 60)
                        succeeded = True
                    except AttributeError:
                        # This version of OpenMM does not support setting PME
                        # parameters
                        succeeded = False
        if succeeded:
            print 'Changed.'
        elif succeeded is None:
            print 'No NonbondedForce detected.'
        else:
            print 'OpenMM is too old. Could not change PME parameters.'
        # Define some integrators
        dummyint1 = mm.VerletIntegrator(1.0e-6*u.picoseconds)
        dummyint2 = mm.VerletIntegrator(1.0e-6*u.picoseconds)
        dummyint3 = mm.VerletIntegrator(1.0e-6*u.picoseconds)
        dummyint4 = mm.VerletIntegrator(1.0e-6*u.picoseconds)
        # Define the contexts
        if properties is None:
            context1 = mm.Context(self.system, dummyint1, plat)
            context2 = mm.Context(system, dummyint2, plat)
            context3 = mm.Context(self.systemapp, dummyint3, plat)
            context4 = mm.Context(xmlsys, dummyint4, plat)
        else:
            context1 = mm.Context(self.system, dummyint1, plat, properties)
            context2 = mm.Context(system, dummyint2, plat, properties)
            context3 = mm.Context(self.systemapp, dummyint3, plat,
                                  properties)
            context4 = mm.Context(xmlsys, dummyint4, plat, properties)
        # Set the context positions
        context1.setPositions(self.parm.positions)
        context2.setPositions(self.parm.positions)
        context3.setPositions(self.crdapp.getPositions())
        context4.setPositions(self.parm.positions)
        # Get the energies
        eunit = u.kilocalories_per_mole
        state = context1.getState(getEnergy=True, getForces=True,
                                  enforcePeriodicBox=True)
        funit = eunit / u.angstrom
        forces = state.getForces().value_in_unit(funit)
        tote = state.getPotentialEnergy().value_in_unit(eunit)
        state = context3.getState(getEnergy=True, enforcePeriodicBox=True)
        toteapp = state.getPotentialEnergy().value_in_unit(eunit)
        state = context4.getState(getEnergy=True, enforcePeriodicBox=True)
        xmltote = state.getPotentialEnergy().value_in_unit(eunit)
        # Now get the decomposed energies from both the system and the
        # deserialized system to check that serialization and deserialization
        # behave as expected with these force objects and force groups
        state = context1.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.BOND_FORCE_GROUP)
        bonde = state.getPotentialEnergy().value_in_unit(eunit)
        state = context4.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.BOND_FORCE_GROUP)
        xmlbonde = state.getPotentialEnergy().value_in_unit(eunit)

        state = context1.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.ANGLE_FORCE_GROUP)
        anglee = state.getPotentialEnergy().value_in_unit(eunit)
        state = context4.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.ANGLE_FORCE_GROUP)
        xmlanglee = state.getPotentialEnergy().value_in_unit(eunit)

        state = context1.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.DIHEDRAL_FORCE_GROUP)
        dihede = state.getPotentialEnergy().value_in_unit(eunit)
        state = context4.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.DIHEDRAL_FORCE_GROUP)
        xmldihede = state.getPotentialEnergy().value_in_unit(eunit)

        state = context1.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.NONBONDED_FORCE_GROUP)
        nonbe = state.getPotentialEnergy().value_in_unit(eunit)
        state = context4.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.NONBONDED_FORCE_GROUP)
        xmlnonbe = state.getPotentialEnergy().value_in_unit(eunit)

        state = context2.getState(getEnergy=True, enforcePeriodicBox=True,
                                  groups=2**self.parm.NONBONDED_FORCE_GROUP)
        vdwe = state.getPotentialEnergy().value_in_unit(eunit)

        eele = nonbe - vdwe
        # Now get the sander forces and compare them
        traj = netcdf_file(get_fn('sander_pme.nc'), 'r')
        sander_forces = traj.variables['forces'][0]
        maxdif = [abs(ofrc-sfrc)
                  for ofrc, sfrc in zip(forces[0], sander_forces[0])]
        maxrel = [abs(ofrc-sfrc)/sfrc
                  for ofrc, sfrc in zip(forces[0], sander_forces[0])]
        avgdif = [0, 0, 0]
        avgrel = [0, 0, 0]
        n = 0
        for ofrcs, sfrcs in zip(forces, sander_forces):
            for i, sfrc in enumerate(sfrcs):
                ofrc = ofrcs[i]
                dif = abs(ofrc-sfrc)
                rel = dif/sfrc
                maxdif[i] = max(maxdif[i], dif)
                maxrel[i] = max(maxrel[i], rel)
                avgdif[i] += dif
                avgrel[i] += rel
            n += 1
        avgdif = [x/n for x in avgdif]
        avgrel = [x/n for x in avgrel]
        # The sander energies are:
# Etot   =    -69285.4160  EKtot   =         0.0000  EPtot      =    -69285.4160
# BOND   =       404.9439  ANGLE   =      1003.4499  DIHED      =      2231.7367
# 1-4 NB =       440.7084  1-4 EEL =      3818.2959  VDWAALS    =      8271.5191
# EELEC  =    -85456.0701  EHBOND  =         0.0000  RESTRAINT  =         0.0000
        sander = dict(bond=404.9439, angle=1003.4499, dihedral=2231.7367,
                      vdw=8271.5191+440.7084, eel=-85456.0701+3818.2959,
                      total=-69285.4160)
        if not use_dispersion_correction:
            # Without the long-range dispersion correction, VDWAALS = 8943.8420
            sander['total'] -= sander['vdw']
            sander['vdw'] = 8943.8420 + 440.7084
            sander['total'] += sander['vdw']
        bonddif = bonde - sander['bond']
        angledif = anglee - sander['angle']
        diheddif = dihede - sander['dihedral']
        vdwdif = vdwe - sander['vdw'] # includes 1-4 also
        eeldif = eele - sander['eel'] # Includes 1-4 also
        totaldif = tote - sander['total']
        appdif = tote - toteapp
        print 'Energy differences compared to sander/Amber (kcal/mol)'
        print '             Absolute     Relative    sander'
        print '------------------------------------------------------'
        print 'Bond     =', colorize_error(bonddif), \
              colorize_error(bonddif/sander['bond'], 1e-6), \
              '%12.4f'%sander['bond']
        print 'Angle    =', colorize_error(angledif), \
              colorize_error(angledif/sander['angle'], 1e-6), \
              '%12.4f'%sander['angle']
        print 'Dihedral =', colorize_error(diheddif), \
              colorize_error(diheddif/sander['dihedral'], 1e-6), \
              '%12.4f'%sander['dihedral']
        if use_dispersion_correction:
            # The dispersion correction in Amber neglects the repulsive part of
            # the correction, but OpenMM does not. Therefore, when we are using
            # the dispersion correction we should allow for a slightly larger
            # energy difference.
            print 'vdWaals  =', colorize_error(vdwdif, 1.0), \
                   colorize_error(vdwdif/sander['vdw'], 1e-4), \
                   '%12.4f'%sander['vdw']
        else:
            print 'vdWaals  =', colorize_error(vdwdif, 1e-2), \
                   colorize_error(vdwdif/sander['vdw'], 1e-6), \
                   '%12.4f'%sander['vdw']
        print 'Elec     =', colorize_error(eeldif, 4e0), \
               colorize_error(eeldif/sander['eel'], 1e-3), '%12.4f'%sander['eel']
        print 'Total    =', colorize_error(totaldif, 4e0), \
               colorize_error(totaldif/sander['total'], 1e-3), \
               '%12.4f'%sander['total']
        print ''
        print 'Difference b/w ParmEd and OpenMM App layer'
        print '------------------------------------------'
        print 'Total    =', colorize_error(appdif, tolerance=5e-5)
        print ''
        print 'Difference b/w sander and OpenMM forces'
        print '---------------------------------------'
        print 'Maximum deviation = [%12s, %12s, %12s]' % colorize_list(maxdif,2e0)
        print 'Maximum rel. dev. = [%12s, %12s, %12s]' % colorize_list(maxrel,2e0)
        print 'Average deviation = [%12s, %12s, %12s]' % colorize_list(avgdif,1e-1)
        print 'Average rel. dev. = [%12s, %12s, %12s]' % colorize_list(avgrel,5e-1)

        # Now test serialization
        CUTOFF = 1e-5
        CUTOFFNB = 1e-2
        print ''
        print 'Serialization tests'
        print '-------------------'
        print 'Bond........',
        if abs(xmlbonde - bonde) < CUTOFF:
            print green('OK')
        else:
            dif = xmlbonde - bonde
            print red('off by %.4e (%f%%)' % (dif, 100*dif/(max(bonde,xmlbonde))))
        print 'Angle.......',
        if abs(xmlanglee - anglee) < CUTOFF:
            print green('OK')
        else:
            dif = xmlanglee - anglee
            print red('off by %.4e (%f%%)' % (dif,100*dif/(max(anglee,xmlanglee))))
        print 'Dihedral....',
        if abs(xmldihede - dihede) < CUTOFF:
            print green('OK')
        else:
            dif = xmldihede - dihede
            print red('off by %.4e (%f%%)' % (dif,100*dif/(max(dihede,xmldihede))))
        print 'Nonbonded...',
        if abs(xmlnonbe - nonbe) < CUTOFFNB:
            print green('OK')
        else:
            dif = xmlnonbe - nonbe
            print red('off by %.4e (%f%%)' % (dif,100*dif/(max(nonbe,xmlnonbe))))

# Log of available tests
tests = (TestPME,)
