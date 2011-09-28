#!/usr/bin/env python
"""pyclangcheck driver code

This code is the driver code for the pyclangcheck tool.

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import optparse
import os
import os.path
import sys

import rules
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


class CollectViolationsVisitor(object):
    """Visitor for AST nodes that collects rule violations."""
    
    def __init__(self, options, rules):
        self.options = options
        self.rules = []
        self.stack = []
        self.violations = []
    
    def enterNode(self, node):
        """Called when a node is entered ("pre-order" traversal)."""
        self.stack.append(node)
        if self.options.verbosity >= 2:
            print ' ' * len(self.stack), 'Entering', node.kind, node.spelling
        violations = []
        for rule in self.rules:
            violations += rule.check(node)
        for v in violations:
            if self.options.verbosity >= 2:
                print v
            self.violations[v.key()] = v

    def exitNode(self, node):
        """Called when a node is left ("post-order" traversa)."""
        self.stack.pop()


class VisitAllowedRule(object):
    """Decides whether a AST node and its children is visited."""
    
    def __init__(self, options):
        self.options = options
        self.include_dirs = [os.path.abspath(x) for x in options.include_dirs]
        self.cache = {}

    def visitAllowed(self, node):
        """Return True if visiting is allowed."""
        # would be enough.  Visit if translation unit.
        if node.kind == ci.CursorKind.TRANSLATION_UNIT:
            return True
        # Don't visit if it has no location (built-in).
        if not _hasFileLocation(node):
            return False
        # Try to hit cache.
        if self.cache.has_key(node.location.file.name):
            return self.cache[node.location.file.name]
        # Check whether node's location is below the include directories.  It is
        # only visited if this is the case.
        filename = os.path.abspath(node.location.file.name)
        result = False
        for x in self.include_dirs:
            if filename.startswith(x):
                # print filename, x
                result = True
                break
        self.cache[node.location.file.name] = result  # Save in cache.
        return result


class AstTraverser(object):
    """Traverses AST tree and applies given visitor object."""
    
    def __init__(self, node_visitor, options):
        self.node_visitor = node_visitor
        self.options = options
        self.visit_allowed_rule = VisitAllowedRule(options)

    def _recurse(self, node):
        """Recursion helper."""
        if not self.visit_allowed_rule.visitAllowed(node):
            return False  # We did not visit this node.
        self.node_visitor.enterNode(node)
        for c in node.get_children():
            self._recurse(c)
        self.node_visitor.exitNode(node)
        return True

    def run(self, filename):
        """Main entry point."""
        index = ci.Index.create()
        args = ['-I%s' % s for s in self.options.include_dirs]
        # print args
        tu = index.parse(filename, args=args)
        if self.options.verbosity >= 1:
            print 'Translation unit: %s.' % tu.spelling
        return self._recurse(tu.cursor)
    
    @classmethod
    def visitFile(klass, filename, node_visitor, options):
        """Don't instantiate AstTraverser yourself, use this function."""
        if options.verbosity >= 1:
            print >>sys.stderr, 'Checking', filename
        traverser = AstTraverser(node_visitor, options)
        res = traverser.run(filename)
        return res != True


def main():
    # ========================================================================
    # Parse command line arguments.
    # ========================================================================
    parser = optparse.OptionParser("USAGE: %prog [options] file.cpp")
    parser.add_option('-s', '--source-file', dest='source_files', default=[],
                      type='string', help='Specify source (.cpp) files.',
                      action='append')
    parser.add_option('-I', '--include-dir', dest='include_dirs', default=[],
                      type='string', help='Specify include directories',
                      action='append')
    parser.add_option('-q', '--quiet', dest='verbosity', default=1,
                      action='store_const', const=0, help='Fewer message.')
    parser.add_option('-v', '--verbose', dest='verbosity', default=1,
                      action='store_const', const=2, help='More messages.')
    options, args = parser.parse_args()

    if len(args) != 0:
        parser.error('Incorrect number of arguments!')
        return 1

    # Recursion Rule: Only check symbols within the include directories.
    recurse_rules = []
    recurse_rules.append(rules.InIncludeDirsRule(options.include_dirs + ['.']))
    # Define symbol naming rules.
    R = rules.GenericSymbolNameRule
    r = rules
    ck = ci.CursorKind
    check_rules = [
        # TODO(holtgrew): There currently is no struct template node type, only struct, which is also used for templates.
        R(ck.STRUCT_DECL                          , r.RE_STRUCT       , r.RULE_NAMING_STRUCT                ),
        R(ck.UNION_DECL                           , r.RE_TYPE         , r.RULE_NAMING_UNION                 ),
        R(ck.CLASS_DECL                           , r.RE_TYPE         , r.RULE_NAMING_CLASS                 ),
        R(ck.ENUM_DECL                            , r.RE_TYPE         , r.RULE_NAMING_ENUM                  ),
        # TODO(holtgrew): Analyze variable type to enforce constant syntax for constants.
        R(ck.FIELD_DECL                           , r.RE_VARIABLE     , r.RULE_NAMING_FIELD                 ),
        R(ck.ENUM_CONSTANT_DECL                   , r.RE_CONSTANT     , r.RULE_NAMING_ENUM_CONSTANT         ),
        R(ck.FUNCTION_DECL                        , r.RE_FUNCTION     , r.RULE_NAMING_FUNCTION              ),
        R(ck.VAR_DECL                             , r.RE_VARIABLE     , r.RULE_NAMING_VARIABLE              ),
        R(ck.PARM_DECL                            , r.RE_VARIABLE     , r.RULE_NAMING_PARAMETER             ),
        # R(ck.TYPEDEF_DECL                         , r.RE_FUNCTION     , r.RULE_NAMING_TYPEDEF               ),
        R(ck.CXX_METHOD                           , r.RE_FUNCTION     , r.RULE_NAMING_CXX_METHOD            ),
        R(ck.TEMPLATE_TYPE_PARAMETER              , r.RE_TYPE         , r.RULE_NAMING_TPL_TYPE_PARAMETER    ),
        R(ck.TEMPLATE_NON_TYPE_PARAMETER          , r.RE_CONSTANT     , r.RULE_NAMING_TPL_NON_TYPE_PARAMETER),
        R(ck.TEMPLATE_TEMPLATE_PARAMTER           , r.RE_TYPE         , r.RULE_NAMING_TPL_TPL_PARAMETER     ),
        R(ck.FUNCTION_TEMPLATE                    , r.RE_FUNCTION     , r.RULE_NAMING_FUNCTION_TPL          ),
        R(ck.CLASS_TEMPLATE                       , r.RE_TYPE_TEMPLATE, r.RULE_NAMING_CLASS_TPL             ),
        R(ck.CLASS_TEMPLATE_PARTIAL_SPECIALIZATION, r.RE_TYPE_TEMPLATE, r.RULE_NAMING_CLASS_TPL_SPEC        ),
    ]

    # Fire off AST traversal.
    node_visitor = CollectViolationsVisitor(options, check_rules)
    for filename in options.source_files:
        res = AstTraverser.visitFile(filename, node_visitor, options)
        if res:
            break
    print node_visitor.violations
    return len(node_visitor.violations) > 0


if __name__ == '__main__':
    sys.exit(main())
