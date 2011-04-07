#!/usr/bin/env python
"""SeqAn code generation from templates / skelletons.

This module contains code to help the creation of modules, tests, apps etc.
It can be called directly or imported and the main() function can be called.

It will perform the following replacements:

  %(AUTHOR)s  will be replaced by the author's name, either given on command
              line or taken from environment variable SEQAN_AUTHOR.

  %(NAME)s    will be replaced by the name of the generated code.
  %(TITLE)s   will be replaced by the name of the generated, but centered in
              74 characters, to be used in the file header comment.

  %(YEAR)d    will be replaced by the current year.
  %(DATE)s    will be replaced by the current date.
  %(TIME)s    will be replaced by the current time.

  %(HEADER_GUARD)s  will be replaced by the UPPER_CASE_PATH_H_ to the file.

  %(CMAKE_PROJECT_PATH)s  will be replaced by lower_case_path to the target
                          directory.

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import datetime
import optparse
import os
import os.path
import sys
import string

import paths

# Length of the header comment.
HEADER_CENTER_WIDTH = 74

# Fallback for author string if neither given on command line or environment
# Variable SEQAN_AUTHOR.
DEFAULT_AUTHOR = 'Your Name <your.email@example.net>'

# Program usage string for command line parser.
USAGE = """
Usage: %prog [options] create [module|test|app|demo|repository] NAME LOCATION
       %prog [options] info [module|test|app|demo|repository]
""".strip()

# Program description, used for command line parser.  Will be wrapped by, though.
DESCRIPTION = """
SeqAn code generator.

The create command uses the template to code of the given type with the given
NAME in the given LOCATION.  The info command displays information on the
template.  The LOCATION is the repository name and could be "core", "extras" or
"sandbox/fub_students".
""".strip()

def createDirectory(path, dry_run=False):
    print 'mkdir(%s)' % path
    print
    if not dry_run:
        os.mkdir(path)

def configureFile(target_file, source_file, replacements, dry_run):
    print 'Configuring file.'
    print '  Source:', source_file
    print '  Target:', target_file
    print
    if os.path.exists(target_file):
        msg = 'Target file already exists.  Move it away and call the script again.'
        print >>sys.stderr, msg
        return 1

    with open(source_file, 'rb') as f:
        contents = f.read()
    target_contents = contents % replacements
    if dry_run:
        print 'The contents of the target file are:'
        print '-' * 78
        print target_contents
        print '-' * 78
    else:
        with open(target_file, 'wb') as f:
            f.write(target_contents)
    return 0

def _pathToIdentifier(relative_path):
    result = relative_path.replace('/', '_')
    result = result.replace('-', '_')
    result = result.replace('.', '_')
    result = result.replace(' ', '_')
    return result

def buildReplacements(type_, name, location, target_file, options):
    result = {}
    result['AUTHOR'] = options.author
    result['YEAR'] = datetime.date.today().year
    result['TIME'] = datetime.datetime.now().strftime('%H:%M')
    result['DATE'] = datetime.date.today().strftime('%Y-%m-%d')
    result['NAME'] = name
    result['TITLE'] = name.center(HEADER_CENTER_WIDTH).rstrip()
    path = os.path.relpath(target_file, paths.repositoryRoot())
    guard = _pathToIdentifier(path).upper()
    result['HEADER_GUARD'] = guard + '_'
    path = os.path.relpath(os.path.dirname(target_file),
                           paths.repositoryRoot())
    cmake_project_name = _pathToIdentifier(path)
    result['CMAKE_PROJECT_NAME'] = cmake_project_name
    if type_ == 'repository':
        result['REPOSITORY_PSEUDO_TARGET_NAME'] = string.capwords(name.replace('/', ' ')).replace(' ', '')
    return result

def _checkTargetPaths(target_path):
    """Check that the path does not exist but its parent does."""
    # Check that the given path does not exist yet.
    if os.path.exists(target_path):
        msg = 'The path %s already exists. Move it and call this script again.'
        print >>sys.stderr, msg % target_path
        return False
    # Check that the parent path already exists.
    if not os.path.exists(os.path.dirname(target_path)):
        msg = 'The parent of the target path does not exist yet: %s'
        print >>sys.stderr, msg % os.path.dirname(target_path)
        print >>sys.stderr, 'Please create it can call this script again.'
        return False
    return True

def createModule(name, location, options):
    include_path = paths.pathToInclude(location)
    seqan_path = os.path.join(include_path, 'seqan')
    module_path = os.path.join(seqan_path, name)
    header_path = os.path.join(seqan_path, '%s.h' % name)
    print 'Creating module in %s' % module_path
    if not options.cmakelists_only and not _checkTargetPaths(module_path):
        return 1
    if not options.cmakelists_only and not _checkTargetPaths(header_path):
        return 1
    print '  Module path is: %s' % module_path
    print '  Module header path is: %s' % header_path
    print ''
    if not options.cmakelists_only:
        # Create directory.
        createDirectory(module_path, options.dry_run)
        # Copy over module header.
        source_file = paths.pathToTemplate('module_template', 'module.h')
        target_file = header_path
        replacements = buildReplacements('module', name, seqan_path, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
        # Copy over header inside module.
        source_file = paths.pathToTemplate('module_template', 'header.h')
        target_file = os.path.join(module_path, '%s_base.h' % name)
        replacements = buildReplacements('module', name, seqan_path, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
        # Copy over INFO file for app and perform replacements.
        source_file = paths.pathToTemplate('app_template', 'INFO')
        target_file = os.path.join(target_path, 'INFO')
        replacements = buildReplacements('app', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
    return 0

def createTest(name, location, options):
    target_path = paths.pathToTest(location, name)
    print 'Creating test in %s' % target_path
    if not options.cmakelists_only and not _checkTargetPaths(target_path):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if not options.cmakelists_only:
        # Create directory.
        createDirectory(target_path, options.dry_run)
        # Copy over .cpp file for test and perform replacements.
        source_file = paths.pathToTemplate('test_template', 'test.cpp')
        target_file = os.path.join(target_path, 'test_%s.cpp' % name)
        replacements = buildReplacements('test', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
        # Copy over .h file for test and perform replacements.
        source_file = paths.pathToTemplate('test_template', 'test.h')
        target_file = os.path.join(target_path, 'test_%s.h' % name)
        replacements = buildReplacements('test', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
    # Copy over CMakeLists.txt file for test and perform replacements.
    source_file = paths.pathToTemplate('test_template', 'CMakeLists.txt')
    target_file = os.path.join(target_path, 'CMakeLists.txt')
    replacements = buildReplacements('test', name, location, target_file, options)
    res = configureFile(target_file, source_file, replacements, options.dry_run)
    if res: return res
    return 0

def createApp(name, location, options):
    target_path = paths.pathToApp(location, name)
    print 'Creating app in %s' % target_path
    if not options.cmakelists_only and not _checkTargetPaths(target_path):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if not options.cmakelists_only:
        # Create directory.
        createDirectory(target_path, options.dry_run)
        # Copy over .cpp file for app and perform replacements.
        source_file = paths.pathToTemplate('app_template', 'app.cpp')
        target_file = os.path.join(target_path, '%s.cpp' % name)
        replacements = buildReplacements('app', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
        # Copy over .h file for app and perform replacements.
        source_file = paths.pathToTemplate('app_template', 'app.h')
        target_file = os.path.join(target_path, '%s.h' % name)
        replacements = buildReplacements('app', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
        # Copy over INFO file for app and perform replacements.
        source_file = paths.pathToTemplate('app_template', 'INFO')
        target_file = os.path.join(target_path, 'INFO')
        replacements = buildReplacements('app', name, location, target_file, options)
        res = configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
    # Copy over CMakeLists.txt file for app and perform replacements.
    source_file = paths.pathToTemplate('app_template', 'CMakeLists.txt')
    target_file = os.path.join(target_path, 'CMakeLists.txt')
    replacements = buildReplacements('app', name, location, target_file, options)
    res = configureFile(target_file, source_file, replacements, options.dry_run)
    if res: return res
    return 0

def createDemo(name, location, options):
    target_path = paths.pathToDemo(location, name)
    print 'Creating app in %s' % target_path
    if not options.cmakelists_only and not _checkTargetPaths(target_path):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if not options.cmakelists_only:
        # Copy over .cpp file for app and perform replacements.
        source_file = paths.pathToTemplate('demo_template', 'demo.cpp')
        target_file = os.path.join(target_path)
        replacements = buildReplacements('app', name, location, target_file, options)
        configureFile(target_file, source_file, replacements, options.dry_run)
        if res: return res
    return 0

def createRepository(location, options):
    print 'Creating module %s' % location
    target_path = paths.pathToRepository(location)
    if not _checkTargetPaths(target_path):
        return 1
    print '  Target path is: %s' % target_path
    print ''
    if not options.cmakelists_only:
        # Create directories.
        createDirectory(target_path, options.dry_run)
        createDirectory(os.path.join(target_path, 'apps'), options.dry_run)
        createDirectory(os.path.join(target_path, 'demos'), options.dry_run)
        createDirectory(os.path.join(target_path, 'include'), options.dry_run)
        createDirectory(os.path.join(target_path, 'include', 'seqan'), options.dry_run)
        createDirectory(os.path.join(target_path, 'tests'), options.dry_run)
    # Copy over file ${REPOSITORY}/CMakeLists.txt.
    target_file = os.path.join(target_path, 'CMakeLists.txt')
    source_file = paths.pathToTemplate('repository_template', 'CMakeLists.txt')
    replacements = buildReplacements('repository', location, target_path, target_file, options)
    configureFile(target_file, source_file, replacements, options.dry_run)
    # Copy over file ${REPOSITORY}/apps/CMakeLists.txt.
    target_file = os.path.join(target_path, 'apps', 'apps_CMakeLists.txt')
    source_file = paths.pathToTemplate('repository_template', 'CMakeLists.txt')
    replacements = buildReplacements('repository', location, target_path, target_file, options)
    configureFile(target_file, source_file, replacements, options.dry_run)
    # Copy over file ${REPOSITORY}/tests/CMakeLists.txt.
    target_file = os.path.join(target_path, 'tests', 'tests_CMakeLists.txt')
    source_file = paths.pathToTemplate('repository_template', 'CMakeLists.txt')
    replacements = buildReplacements('repository', location, target_path, target_file, options)
    configureFile(target_file, source_file, replacements, options.dry_run)
    # Copy over file ${REPOSITORY}/demos/CMakeLists.txt.
    target_file = os.path.join(target_path, 'demos', 'demos_CMakeLists.txt')
    source_file = paths.pathToTemplate('repository_template', 'CMakeLists.txt')
    replacements = buildReplacements('repository', location, target_path, target_file, options)
    configureFile(target_file, source_file, replacements, options.dry_run)
    return 0

def main():
    # Parse arguments.
    parser = optparse.OptionParser(usage=USAGE, description=DESCRIPTION)
    parser.add_option('-s', '--skel-root', dest='skel_root',
                      help=('Set path to the directory where the skelletons '
                            'live in.  Taken from environment variable '
                            'SEQAN_SKELS if available.'),
                      default=os.environ.get('SEQAN_SKELS',
                                             paths.pathToSkelletons()))
    parser.add_option('-a', '--author', dest='author',
                      help=('Set author to use.  Should have the format USER '
                            '<EMAIL>.  Taken from environment variable '
                            'SEQAN_AUTHOR if it exists.'),
                      default=os.environ.get('SEQAN_AUTHOR', DEFAULT_AUTHOR))
    parser.add_option('-d', '--dry-run', dest='dry_run', action='store_true',
                      help='Do not change anything, just simulate.',
                      default=False)
    parser.add_option('-c', '--cmakelists-only', dest='cmakelists_only',
                      action='store_true',
                      help='Only create CMakeLists.txt files',
                      default=False)
    options, args = parser.parse_args()
    if not args:
        parser.print_help(file=sys.stderr)
        return 1
    if args[0] == 'create':
        if len(args) < 3:
            print >>sys.stderr, 'Invalid argument count!'
            return 1
        if args[1] not in ['module', 'test', 'app', 'demo', 'repository']:
            print >>sys.stderr, 'Invalid template "%s".' % args[1]
            return 1
        if args[1] == 'repository':
            if len(args) != 3:
                print >>sys.stderr, 'Invalid argument count!'
                return 1
            return createRepository(args[2], options)
        elif len(args) != 4:
            print >>sys.stderr, 'Invalid argument count!'
            return 1
        create_methods = {
            'module' : createModule,
            'test': createTest,
            'app': createApp,
            'demo': createDemo,
            }
        return create_methods[args[1]](args[2], args[3], options)
    elif args[0] == 'info':
        print >>sys.stderr, 'Command "info" is not implemented yet.'
        return 1
    else:
        print >>sys.stderr, 'Unknown command "%s".' % args[0]
        parser.print_help(file=sys.stderr)
        return 1

if __name__ == '__main__':
   sys.exit(main())

