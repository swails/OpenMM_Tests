"""
Various utilities for running the tests
"""

import os
TEST_FILE_DIRECTORY = None

TOTAL_FAILED = 0
TOTAL_PASSED = 0

TOTAL_TEST_SUMMARY = dict()

def register_test_directory(folder):
    """ Registers the directory where the test files are kept """
    global TEST_FILE_DIRECTORY
    TEST_FILE_DIRECTORY = os.path.abspath(folder)

def get_fn(fname, status='OLD'):
    """ Returns the name of the file in the test file directory
    
    Parameters
    ----------
    fname : str
        Name of the file to find
    status : 'OLD' or 'NEW' (case-insensitive)
        Whether the file should exist ('OLD') or 

    Returns
    -------
    string
        Name of the desired file with the full path

    Raises
    ------
    If "register_test_directory" has not yet been called, RuntimeError is raised
    If status is not OLD or NEW, ValueError is raised
    If status is OLD and the file does not exist, IOError is raised
    """
    global TEST_FILE_DIRECTORY
    if TEST_FILE_DIRECTORY is None:
        raise RuntimeError('The test file directory was not registered!')
    status = status.upper()
    if status not in ('OLD', 'NEW'):
        raise ValueError('status must be OLD or NEW')
    fname = os.path.join(TEST_FILE_DIRECTORY, fname)
    if not os.path.exists(fname) and status == 'OLD':
        raise IOError('Could not find %s' % fname)
    return fname

def red(text, count=True):
    """ Return text red """
    global TOTAL_FAILED
    TOTAL_FAILED += int(bool(count))
    return '\033[91m%s\033[0m' % text

def green(text, count=True):
    """ Return text green """
    global TOTAL_PASSED
    TOTAL_PASSED += int(bool(count))
    return '\033[92m%s\033[0m' % text

def colorize_error(error, tolerance=0.001):
    """
    Print out the error in either green or red (agrees within tolerance or
    differs within tolerance, respectively)

    Parameters
    ----------
    error : float
        The error to print
    tolerance : float
        The cutoff by which we declare an error is "ok" (green) or "not" (red)

    Notes
    -----
    The absolute value of the error is compared to the tolerance. Tolerance is
    inclusive (i.e., if error == tolerance, it is printed as green)
    """
    if abs(error) <= tolerance:
        return green('%12.4e' % error)
    return red('%12.4e' % error)

def colorize_list(numbers, tolerance=0.001):
    """
    Same as colorize_error, but for a list of numbers

    Parameters
    ----------
    error : float
        The error to print
    tolerance : float
        The cutoff by which we declare an error is "ok" (green) or "not" (red)

    Returns
    -------
    tuple of colorized strings for each number

    Notes
    -----
    The absolute value of the error is compared to the tolerance. Tolerance is
    inclusive (i.e., if error == tolerance, it is printed as green)
    """
    retval = []
    for number in numbers:
        if abs(number) <= tolerance:
            retval.append(green('%12.4e' % number))
        else:
            retval.append(red('%12.4e' % number))
    return tuple(retval)

def summarize(testname):
    """
    Summarizes the number of times we failed a comparison (RED) and the number
    of times we passed a comparison (GREEN). Reset the tallies so we can use
    this again
    """
    global TOTAL_FAILED, TOTAL_PASSED, TOTAL_TEST_SUMMARY
    failmsg = red('%d total failures' % TOTAL_FAILED, count=False)
    passmsg = green('%d total successes' % TOTAL_PASSED, count=False)
    TOTAL_TEST_SUMMARY[testname] = '%s\n%s' % (failmsg, passmsg)
    print failmsg
    print passmsg
    TOTAL_FAILED = TOTAL_PASSED = 0
