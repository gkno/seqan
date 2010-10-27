#!/usr/bin/env python
"""Execute the tests for the sak program.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os.path
import sys

import seqan.app_tests as app_tests

def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for sak'
    print '======================'
    print

    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'projects/tests/apps/sak')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'projects/library/cmake/apps', 'sak')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # ============================================================
    # Run on DNA (Adenoviruses).
    # ============================================================

    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-o', ph.outFile('adeno.all.out')],
        to_diff=[(ph.inFile('adeno.all.out'),
                  ph.outFile('adeno.all.out'))])
    conf_list.append(conf)
    
    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-s', '1',
              '-o', ph.outFile('adeno.seq1.out')],
        to_diff=[(ph.inFile('adeno.seq1.out'),
                  ph.outFile('adeno.seq1.out'))])
    conf_list.append(conf)
    
    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-ss', '1', '2',
              '-o', ph.outFile('adeno.seq1-2.out')],
        to_diff=[(ph.inFile('adeno.seq1-2.out'),
                  ph.outFile('adeno.seq1-2.out'))])
    conf_list.append(conf)
    
    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-s', '3',
              '-o', ph.outFile('adeno.seq3.out')],
        to_diff=[(ph.inFile('adeno.seq3.out'),
                  ph.outFile('adeno.seq3.out'))])
    conf_list.append(conf)
    
    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-sn', 'gi|9626621',
              '-o', ph.outFile('adeno.sn.out')],
        to_diff=[(ph.inFile('adeno.sn.out'),
                  ph.outFile('adeno.sn.out'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-s', '1',
              '-i', '5', '25',
              '-o', ph.outFile('adeno.s1i5-25.out')],
        to_diff=[(ph.inFile('adeno.s1i5-25.out'),
                  ph.outFile('adeno.s1i5-25.out'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-ss', '1', '2',
              '-i', '5', '25',
              '-o', ph.outFile('adeno.s1-2i5-25.out')],
        to_diff=[(ph.inFile('adeno.s1-2i5-25.out'),
                  ph.outFile('adeno.s1-2i5-25.out'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        args=[ph.inFile('adeno.fasta'),
              '-s', '1',
              '-rc',
              '-o', ph.outFile('adeno.s1rc.out')],
        to_diff=[(ph.inFile('adeno.s1rc.out'),
                  ph.outFile('adeno.s1rc.out'))])
    conf_list.append(conf)

    # Execute the tests.
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['sak'] + conf.args),
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
    sys.exit(app_tests.main(main))
