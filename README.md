OpenMM_Tests
============

A series of tests for OpenMM that tests OpenMM's comparison with other MD programs in terms of accuracy.


To run
======

To run the tests, select the tests you want to run and the platforms you want to test. The GPU platforms will have all precision models tested by default.

For a full listing of help from the driver script, use the `--help` flag.

As an example:

```
./runtest.py --all-platforms --all-tests
```

will run all of the tests on all of the platforms.
