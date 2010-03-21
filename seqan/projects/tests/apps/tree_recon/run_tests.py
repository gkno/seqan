#!/usr/bin/env python
"""Execute the tests for the tree_recomb program.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import os.path
import sys

import app_tests
 
# Path of the binary under test, relative to the checkout.
BINARY = 'projects/library/cmake/apps/tree_recon'


def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for tree_recomb'
    print '=============================='
    print
    
    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'projects/tests/apps/tree_recon')  # tests dir

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []
    for i in [1, 2, 3]:
        conf = app_tests.TestConf(
            program=os.path.join(ph.binary_base_path, BINARY),
            args=['-m', ph.inFile('example%d.dist' % i),
                  '-o', ph.outFile('example%d.out' % i)],
            to_diff=[(ph.inFile('example%d.out' % i),
                      ph.outFile('example%d.out' % i))])
        conf_list.append(conf)
    for i in [1, 2, 3]:
        for b in ['nj', 'min', 'max', 'avg', 'wavg']:
            conf = app_tests.TestConf(
                program=os.path.join(ph.binary_base_path, BINARY),
                args=['-b', b,
                      '-m', ph.inFile('example%d.dist' % i),
                      '-o', ph.outFile('example%d.%s.out' % (i, b))],
                to_diff=[(ph.inFile('example%d.%s.out' % (i, b)),
                          ph.outFile('example%d.%s.out' % (i, b)))])
            conf_list.append(conf)
    for i in [1, 2, 3]:
        for f in ['dot', 'newick']:
            conf = app_tests.TestConf(
                program=os.path.join(ph.binary_base_path, BINARY),
                args=['-f', f,
                      '-m', ph.inFile('example%d.dist' % i),
                      '-o', ph.outFile('example%d.%s.out' % (i, f))],
                to_diff=[(ph.inFile('example%d.%s.out' % (i, f)),
                          ph.outFile('example%d.%s.out' % (i, f)))])
            conf_list.append(conf)
    
    # Execute the tests.
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['tree_recomb'] + conf.args),
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'
    print
    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    print >>sys.stderr, sys.argv
    if len(sys.argv) != 3:
        print >>sys.stderr, 'ERROR: Invalid arguments!'
        print >>sys.stderr, 'Usage: run_tests SOURCE_ROOT_PATH BINARY_ROOT_PATH'
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
