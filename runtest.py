#!/usr/bin/env python

"""
This is the driver program that actually runs the tests.
"""
# Register the test directory
from argparse import ArgumentParser
import os
import sys
import testcases
import testcases.utils as utils
import warnings
fpath, fname = os.path.split(sys.argv[0])
utils.register_test_directory(os.path.join(fpath, 'test_files'))

parser = ArgumentParser()
group = parser.add_argument_group('Available test cases',
                'All available test cases that can be run with OpenMM.')
group.add_argument('--all-tests', dest='all_tests', default=False,
                   action='store_true', help='Run all available tests')
group.add_argument('--amber-hewl-pme', dest='amber_hewl_pme', default=False,
                   action='store_true', help='''Run the HEWL tests from the
                   Amber starting files in explicit solvent.''')
group.add_argument('--skip-slow' dest='skip_slow', default=False,
                   action='store_true', help='''Skip platforms for tests that
                   are known to be very slow (CustomGBForce on Reference or CPU
                   platforms on big systems, for instance)''')
group = parser.add_argument_group('Platforms', '''Options controlling which
                                  platforms get tested.''')
group.add_argument('--list-platforms', dest='list_platforms', default=False,
                   action='store_true', help='''List all available platforms on
                   this machine and quit.''')
group.add_argument('--all-platforms', dest='all_platforms', default=False,
                   action='store_true', help='Test all platforms')
group.add_argument('--reference', dest='reference', default=False,
                   action='store_true', help='Test the Reference platform')
group.add_argument('--cpu', dest='cpu', default=False,
                   action='store_true', help='Test the CPU platform')
group.add_argument('--cuda', dest='cuda', default=False,
                   action='store_true', help='Test the CUDA platform')
group.add_argument('--opencl', dest='opencl', default=False,
                   action='store_true', help='Test the OpenCL platform')

opt = parser.parse_args()

platmap = dict(CPU=opt.cpu, CUDA=opt.cuda, Reference=opt.reference,
               OpenCL=opt.opencl)

if opt.list_platforms:
    print 'Available platforms are:'
    print '\t' + '\n\t'.join(testcases.available_platforms)
    sys.exit(0)

if opt.cuda and not 'CUDA' in testcases.available_platforms:
    warnings.warn('CUDA platform is not available on this machine')
    opt.cuda = False
if opt.opencl and not 'OpenCL' in testcases.available_platforms:
    warnings.warn('OpenCL platform is not available on this machine')
    opt.opencl = False
if opt.cpu and not 'OpenCL' in testcases.available_platforms:
    warnings.warn('CPU platform is not available in this OpenMM release.')
    opt.cpu = False

# Make a list of the platforms we want to test
testplatforms = []
for plat in testcases.available_platforms:
    if platmap[plat] or opt.all_platforms:
        testplatforms.append(plat)
if not testplatforms:
    sys.exit('No available platforms selected for testing')
print 'Testing the following platforms:'
print '\t' + '\n\t'.join(testplatforms) + '\n'

def runtest(test, parallel=None, **kwargs):
    """
    Runs an arbitrary test case through all testing platforms and all precision
    models
    """
    global testplatforms
    for plat in testplatforms:
        print '\n\tTesting %s platform' % plat
        print '\t----------------' + '-'*len(plat)
        if plat in ('CUDA', 'OpenCL'):
            for precision in ('double', 'mixed', 'single'):
                print '\t%s precision' % precision
                test.run(plat, precision, devices=parallel, **kwargs)
        else:
            test.run(plat, **kwargs)

# Now go through all of the tests
if opt.all_tests or opt.amber_hewl_pme:
    print '\nTesting the HEWL PME test case (with long-range correction)'
    print '-----------------------------------------------------------'
    test = testcases.testhewl.TestPME()
    runtest(test, use_dispersion_correction=True)
    print '\nTesting the HEWL PME test case (without long-range correction)'
    print '--------------------------------------------------------------'
    runtest(test, use_dispersion_correction=False)
    print '='*80
