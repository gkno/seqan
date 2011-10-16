#!/usr/bin/env python
"""Execute the tests for the pair_align program.

The golden test outputs are generated by the script generate_outputs.sh.

You have to give the root paths to the source and the binaries as arguments to
the program.  These are the paths to the directory that contains the 'projects'
directory.

Usage:  run_tests.py SOURCE_ROOT_PATH BINARY_ROOT_PATH
"""
import logging
import os.path
import sys
import zipfile

import seqan.app_tests as app_tests

def unzip_file_into_dir(file, dir):
    if not os.path.exists(dir):
        os.mkdir(dir, 0777)
    zfobj = zipfile.ZipFile(file)
    for name in zfobj.namelist():
        if name.endswith('/'):
            if not os.path.exists(os.path.join(dir, name)):
                os.mkdir(os.path.join(dir, name))
        else:
            outfile = open(os.path.join(dir, name), 'wb')
            outfile.write(zfobj.read(name))
            outfile.close()

def main(source_base, binary_base):
    """Main entry point of the script."""

    print 'Executing test for dfi'
    print '======================'
    print

    ph = app_tests.TestPathHelper(
        source_base, binary_base,
        'core/apps/dfi/tests')  # tests dir

    # ============================================================
    # Auto-detect the binary path.
    # ============================================================

    path_to_program = app_tests.autolocateBinary(
      binary_base, 'core/apps/dfi', 'dfi')

    # ============================================================
    # Built TestConf list.
    # ============================================================

    # Build list with TestConf objects, analoguely to how the output
    # was generated in generate_outputs.sh.
    conf_list = []

    paramsFile = open(ph.inFile("params.txt"), "r")

    # unzip datasets into test directory
    unzip_file_into_dir(ph.inFile("datasets.zip"), ph.inFile(""))
    unzip_file_into_dir(ph.inFile("results.zip"), ph.inFile(""))

    # We run the following for all read lengths we have reads for.
    for paramLine in paramsFile.readlines():
        params = paramLine.split();
        # Run with options parsed from params file.
        conf = app_tests.TestConf(
            program=path_to_program,
            redir_stdout=ph.outFile(params[0]),
            args=[ph.inFile(params[1]),
                  ph.inFile(params[2])] + params[3:],
            to_diff=[(ph.inFile(params[0]),
                      ph.outFile(params[0]))])
        conf_list.append(conf)

    # Execute the tests.
    failures = 0
    for conf in conf_list:
        res = app_tests.runTest(conf)
        # Output to the user.
        print ' '.join(['dfi'] + conf.args),
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
