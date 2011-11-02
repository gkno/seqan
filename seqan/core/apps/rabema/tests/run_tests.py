#!/usr/bin/env python
"""Execute the tests for rabema.

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

    print 'Executing test for rabema'
    print '========================='
    print
    
    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'core/apps/rabema/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'core/apps/rabema', 'rabema')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # ============================================================
    # Build Gold Standard
    # ============================================================

    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('gold-adeno-hamming-08.stdout'),
        args=['build_standard', '-d', 'hamming', '-e', '8',
              '-o', ph.outFile('gold-adeno-hamming-08.wit'),
              ph.inFile('adeno-genome.fa'),
              ph.inFile('gold-adeno-hamming-08.sam')],
        to_diff=[(ph.inFile('gold-adeno-hamming-08.stdout'),
                  ph.outFile('gold-adeno-hamming-08.stdout')),
                 (ph.inFile('gold-adeno-hamming-08.wit'),
                  ph.outFile('gold-adeno-hamming-08.wit'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('gold-adeno-edit-08.stdout'),
        args=['build_standard', '-d', 'edit', '-e', '8',
              '-o', ph.outFile('gold-adeno-edit-08.wit'),
              ph.inFile('adeno-genome.fa'),
              ph.inFile('gold-adeno-edit-08.sam')],
        to_diff=[(ph.inFile('gold-adeno-edit-08.stdout'),
                  ph.outFile('gold-adeno-edit-08.stdout')),
                 (ph.inFile('gold-adeno-edit-08.wit'),
                  ph.outFile('gold-adeno-edit-08.wit'))])
    conf_list.append(conf)

    # ============================================================
    # Compare.
    # ============================================================

    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('razers2-adeno-hamming-08.stdout'),
        args=['compare', '-d', 'hamming', '-e', '8',
              ph.inFile('adeno-genome.fa'),
              ph.inFile('razers2-adeno-hamming-08.sam'),
              ph.inFile('gold-adeno-hamming-08.wit')],
        to_diff=[(ph.inFile('razers2-adeno-hamming-08.stdout'),
                  ph.outFile('razers2-adeno-hamming-08.stdout'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('razers2-adeno-hamming-04.stdout'),
        args=['compare', '-d', 'hamming', '-e', '8',
              ph.inFile('adeno-genome.fa'),
              ph.inFile('razers2-adeno-hamming-04.sam'),
              ph.inFile('gold-adeno-hamming-08.wit')],
        to_diff=[(ph.inFile('razers2-adeno-hamming-04.stdout'),
                  ph.outFile('razers2-adeno-hamming-04.stdout'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('razers2-adeno-edit-08.stdout'),
        args=['compare', '-d', 'edit', '-e', '8',
              ph.inFile('adeno-genome.fa'),
              ph.inFile('razers2-adeno-edit-08.sam'),
              ph.inFile('gold-adeno-edit-08.wit')],
        to_diff=[(ph.inFile('razers2-adeno-edit-08.stdout'),
                  ph.outFile('razers2-adeno-edit-08.stdout'))])
    conf_list.append(conf)

    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('razers2-adeno-edit-04.stdout'),
        args=['compare', '-d', 'edit', '-e', '8',
              ph.inFile('adeno-genome.fa'),
              ph.inFile('razers2-adeno-edit-04.sam'),
              ph.inFile('gold-adeno-edit-08.wit')],
        to_diff=[(ph.inFile('razers2-adeno-edit-04.stdout'),
                  ph.outFile('razers2-adeno-edit-04.stdout'))])

    conf_list.append(conf)

    # Execute the tests.
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['rabema'] + conf.args),
        if res:
             print 'OK'
        else:
            failures += 1
            print 'FAILED'

    print '=============================='
    print '     total tests: %d' % len(conf_list)
    print '    failed tests: %d' % failures
    print 'successful tests: %d' % (len(conf_list) - failures)
    print '=============================='
    # Compute and return return code.
    return failures != 0


if __name__ == '__main__':
    sys.exit(app_tests.main(main))