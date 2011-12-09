#!/usr/bin/env python
"""Execute the tests for micro_razers.

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

    print 'Executing test for micro_razers'
    print '========================='
    print
    
    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'core/apps/micro_razers/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'core/apps/micro_razers', 'micro_razers')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # ============================================================
    # First Section.
    # ============================================================

    # Run with default options.
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('se-adeno-reads36_1_default.stdout'),
        args=[ph.inFile('adeno-genome.fa'),
              ph.inFile('adeno-reads36_1.fa'),
              '-o', ph.outFile('se-adeno-reads36_1_default.out' )],
        to_diff=[(ph.inFile('se-adeno-reads36_1_default.out' ),
                  ph.outFile('se-adeno-reads36_1_default.out' )),
                 (ph.inFile('se-adeno-reads36_1_default.stdout' ),
                  ph.outFile('se-adeno-reads36_1_default.stdout' ))])
    conf_list.append(conf)
    
    # Run with different seed lengths
    for sl in range(14,21):
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile('se-adeno-reads36_1_sl%d.stdout' % sl),
            args=['-sL', str(sl), 
                  ph.inFile('adeno-genome.fa'),
                  ph.inFile('adeno-reads36_1.fa'),
                  '-o', ph.outFile('se-adeno-reads36_1_sl%d.out' % sl)],
            to_diff=[(ph.inFile('se-adeno-reads36_1_sl%d.out' % sl),
                      ph.outFile('se-adeno-reads36_1_sl%d.out' % sl)),
                     (ph.inFile('se-adeno-reads36_1_sl%d.stdout' % sl),
                      ph.outFile('se-adeno-reads36_1_sl%d.stdout' % sl))])
        conf_list.append(conf)
    
        # allow error in seed
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile('se-adeno-reads36_1_sl%d_se.stdout' % sl),
            args=['-sL', str(sl), '-sE',
                  ph.inFile('adeno-genome.fa'),
                  ph.inFile('adeno-reads36_1.fa'),
                  '-o', ph.outFile('se-adeno-reads36_1_sl%d_se.out' % sl)],
            to_diff=[(ph.inFile('se-adeno-reads36_1_sl%d_se.out' % sl),
                      ph.outFile('se-adeno-reads36_1_sl%d_se.out' % sl)),
                     (ph.inFile('se-adeno-reads36_1_sl%d_se.stdout' % sl),
                      ph.outFile('se-adeno-reads36_1_sl%d_se.stdout' % sl))])
        conf_list.append(conf)


    # change maxhits parameter 
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('se-adeno-reads36_1_sl18_m20_pa.stdout' ),
        args=['-sL', str(18), '-m', str(20), '-pa',
              ph.inFile('adeno-genome.fa'),
              ph.inFile('adeno-reads36_1.fa'),
              '-o', ph.outFile('se-adeno-reads36_1_sl18_m20_pa.out' )],
        to_diff=[(ph.inFile('se-adeno-reads36_1_sl18_m20_pa.out' ),
                  ph.outFile('se-adeno-reads36_1_sl18_m20_pa.out' )),
                 (ph.inFile('se-adeno-reads36_1_sl18_m20_pa.stdout' ),
                  ph.outFile('se-adeno-reads36_1_sl18_m20_pa.stdout' ))])
    conf_list.append(conf)

    # ============================================================
    # Execute the tests.
    # ============================================================
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['micro_razers'] + conf.args),
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'

    # Cleanup.
    ph.deleteTempDir()

    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    sys.exit(app_tests.main(main))
