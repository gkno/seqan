#!/usr/bin/env python
"""pyclangcheck driver code

This code is the driver code for the pyclangcheck tool.

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

from __future__ import with_statement

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import optparse
import os
import pickle
import sys

import clang.cindex as ci


def _hasFileLocation(node):
    """Return True if node has a file lcoation."""
    if not hasattr(node, 'location'):
        return False
    if not hasattr(node.location, 'file'):
        return False
    if not node.location.file:
        return False
    if not hasattr(node.location.file, 'name'):
        return False
    if not node.location.file.name:
        return False
    return True


class CollectCompoundStatementNodeVisitor(object):
    def __init__(self, options):
        self.options = options
        self.stack = []
        self.ranges = []
    
    def enterNode(self, node):
        self.stack.append(node)
        ## print '  ' * len(self.stack), node.kind,
        num_children = len([x for x in node.get_children()])
        ## if _hasFileLocation(node):
        ##     print node.location.file.name, '%d-%d' % (node.extent.start.line, node.extent.end.line)
        ## else:
        ##     print
        # Only add range for statements that are no compound statements.  Add
        # for empty compounds.
        if not node.kind.is_statement():
            ## print 'skipping, no statement'
            return
        if node.kind == ci.CursorKind.COMPOUND_STMT and num_children > 0:
            ## print 'skipping, non-empty compound statement', num_children
            return
        # Only add if has file location.
        if _hasFileLocation(node):
            self.ranges.append((node.location.file.name, node.extent.start.line, node.extent.end.line))

    def exitNode(self, node):
        self.stack.pop()


class VisitAllowedRule(object):
    def __init__(self, options):
        self.options = options
        self.include_dirs = [os.path.abspath(x) for x in options.include_dirs]
        self.cache = {}

    def visitAllowed(self, node):
        """Return True if visiting is allowed."""
        # TODO(holtgrew): For this application, stopping at compound statements would be enough.
        print node._kind_id
        if node.kind == ci.CursorKind.TRANSLATION_UNIT:
            return True
        if not _hasFileLocation(node):
            return False
        if self.cache.has_key(node.location.file.name):
            return self.cache[node.location.file.name]
        # Check whether node's location is below the include directories.
        filename = os.path.abspath(node.location.file.name)
        result = False
        for x in self.include_dirs:
            if filename.startswith(x):
                # print filename, x
                result = True
                break
        self.cache[node.location.file.name] = result
        return result


class AstTraverser(object):
    def __init__(self, node_visitor, options):
        self.node_visitor = node_visitor
        self.options = options
        self.visit_allowed_rule = VisitAllowedRule(options)

    def _recurse(self, node):
        if not self.visit_allowed_rule.visitAllowed(node):
            return False  # We did not visit this node.
        self.node_visitor.enterNode(node)
        for c in node.get_children():
            self._recurse(c)
        self.node_visitor.exitNode(node)
        return True

    def run(self, filename):
        index = ci.Index.create()
        args = ['-I%s' % s for s in self.options.include_dirs]
        # print args
        tu = index.parse(filename, args=args)
        print 'Translation unit: %s.' % tu.spelling
        return self._recurse(tu.cursor)
    
    @classmethod
    def visitFile(klass, filename, node_visitor, options):
        traverser = AstTraverser(node_visitor, options)
        res = traverser.run(filename)
        return res == True


def main():
    # Parse command line arguments.
    parser = optparse.OptionParser("USAGE: %prog [options] -s file.cpp")
    parser.add_option('-I', '--include-dir', dest='include_dirs', default=[],
                      type='string', help='Specify include directories',
                      action='append')
    parser.add_option('-s', '--src-file', dest='source_files', default=[],
                      type='string', help='Specify compilation units.',
                      action='append')
    parser.add_option('-l', '--location-file', dest='location_file',
                      default='locations.dat', type='string',
                      help='Path to file with compound statement locations.')
    parser.add_option('-g', '--gcov-file', dest='gcov_files', default=[],
                      type='string', help='Specify gcov files to process.',
                      action='append')
    options, args = parser.parse_args()
    if len(args) != 0:
        parser.error('Incorrect number of arguments!')
        return 1

    options.include_dirs += [os.path.abspath(os.path.dirname(s)) for s in options.source_files]

    if not options.source_files and not options.gcov_files:
        parser.error('Neither source nor gcov file given!')
        return 1

    if options.source_files:
        # If any source file is given, all given source files are parsed and all
        # lines with compound statements in all included files are written to
        # the location file.
        print >>sys.stderr, 'Building Locations'
        print >>sys.stderr, '=================='

        # Fire off AST traversal.
        print >>sys.stderr, 'AST Traversal'
        node_visitor = CollectCompoundStatementNodeVisitor(options)
        for src in options.source_files:
            print >>sys.stderr, '  Compilation Unit', src
            AstTraverser.visitFile(src, node_visitor, options)

        # Convert locations into points.
        locations = {}
        for filename, start, stop in node_visitor.ranges:
            filename = os.path.abspath(filename)
            for i in range(start, stop + 1):
                locations.setdefault(filename, set()).add(i)

        # Write out the source locations.
        print >>sys.stderr, 'Writing out locations to', options.location_file
        with open(options.location_file, 'wb') as f:
            pickle.dump(locations, f)

    if options.gcov_files:
        # If no source files and gcov files are given then
        print >>sys.stderr, 'Updating gcov Results'
        print >>sys.stderr, '====================='

        if not options.source_files:
            print >>sys.stderr, 'Loading locations from', options.location_file
            with open(options.location_file, 'rb') as f:
                locations = pickle.load(f)

        for filename in options.gcov_files:
            filename = os.path.abspath(filename)
            print >>sys.stderr, 'Processing', filename
            with open(filename, 'rb') as f:
                lines = f.readlines()
            pos0 = lines[0].find(':')
            pos1 = lines[0].find(':', pos0 + 1)
            source = None
            result = []
            skip = False
            for i, line in enumerate(lines):
                coverage = line[:pos0]
                lineno = int(line[pos0 + 1:pos1].strip())
                slineno = line[pos0 + 1:pos1]
                txt = line[pos1 + 1:]
                if txt.startswith('Source:'):
                    source = os.path.abspath(txt[len('Source:'):].strip())
                    if not locations.has_key(source):
                        print >>sys.stderr, '  Skipping.'
                        skip = True
                        break
                if not source or lineno == 0:
                    result.append(line)
                    continue  # Proceed only if in file.
                if lineno in locations[source] and coverage.strip() == '-':
                    coverage = ('%%%ds' % pos0) % '#####'
                result.append(':'.join([coverage, slineno, txt]))
            # Write back file if not skipped.
            if skip:
                continue
            print ''.join(result)

if __name__ == '__main__':
    sys.exit(main())
