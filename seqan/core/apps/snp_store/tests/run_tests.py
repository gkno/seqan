#!/usr/bin/env python
"""Execute the tests for snp_store.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os.path
import sys

# Automagically add util/py_lib to PYTHONPATH environment variable.
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..',
                                    '..', '..', 'util', 'py_lib'))
sys.path.insert(0, path)

import seqan.app_tests as app_tests

def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for snp_store'
    print '========================='
    print
    
    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'core/apps/snp_store/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'core/apps/snp_store', 'snp_store')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    # ============================================================
    # First Section.
    # ============================================================

    # App TestConf objects to conf_list, just like this for each
    # test you want to run.
    # default
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('snp_store_default.stdout'),
        args=[ph.inFile('human-chr22-inf2.fa'),
              ph.inFile('human-reads2.gff'),
              '-id', ph.outFile('indels_default.out'),
              '-o', ph.outFile('snps_default.out')],
        to_diff=[(ph.inFile('snp_store_default.stdout'),
                  ph.outFile('snp_store_default.stdout')),
#                 (ph.inFile('snps_default.out'),
#                  ph.outFile('snps_default.out')),
                 (ph.inFile('indels_default.out'),
                  ph.outFile('indels_default.out'))])
    conf_list.append(conf)
    
    # test 2
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('snp_store_realign.stdout'),
        args=[ph.inFile('human-chr22-inf2.fa'),
              ph.inFile('human-reads2.sam'),
              '-re', '-if', str(1),
              '-id', ph.outFile('indels_realign.out'),
              '-o', ph.outFile('snps_realign.out')],
        to_diff=[(ph.inFile('snp_store_realign.stdout'),
                 ph.outFile('snp_store_realign.stdout')),
#                 (ph.inFile('snps_realign.out'),
#                  ph.outFile('snps_realign.out')),
                 (ph.inFile('indels_realign.out'),
                  ph.outFile('indels_realign.out'))])
    conf_list.append(conf)

    # test 3
    conf = app_tests.TestConf(
        program=path_to_program,
        redir_stdout=ph.outFile('snp_store_realign_m0mp1oa.stdout'),
        args=[ph.inFile('human-chr22-inf2.fa'),
              ph.inFile('human-reads2.sam'),
              '-re', '-if', str(1), '-oa', '-hq', '-it', str(1), '-mp', str(1), '-m', str(1),
              '-id', ph.outFile('indels_realign_m0mp1oa.out'),
              '-o', ph.outFile('snps_realign_m0mp1oa.out')],
        to_diff=[(ph.inFile('snp_store_realign_m0mp1oa.stdout'),
                  ph.outFile('snp_store_realign_m0mp1oa.stdout')),
#                 (ph.inFile('snps_realign_m0mp1oa.out'),
#                  ph.outFile('snps_realign_m0mp1oa.out')),
                 (ph.inFile('indels_realign_m0mp1oa.out'),
                  ph.outFile('indels_realign_m0mp1oa.out'))])
    conf_list.append(conf)



    # ============================================================
    # Execute the tests.
    # ============================================================
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['snp_store'] + conf.args),
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
