#!/usr/bin/env python
"""pyclangcheck driver code

This code is the driver code for the pyclangcheck tool.

Copyright: (c) 2010, Knut Reinert, FU Berlin
License:   3-clause BSD (see LICENSE)
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import optparse
import sys

import rules
import clang.cindex as ci


class RulesDispatcher(object):
    def __init__(self, rules, options):
        self.rules = rules
        self.options = options
        self.violations = {}
    
    def enterNode(self, node):
        pass

    def exitNode(self, node):
        pass
    
    def check(self, node):
        violations = []
        for x in self.rules:
            if x.allowVisit(node):
                violations += x.check(node)
        for v in violations:
            print v
            self.violations[v.key()] = v


class AstVisitor(object):
    def __init__(self, dispatcher, recurse_rules, check_rules, options):
        self.dispatcher = dispatcher
        self.recurse_rules = recurse_rules
        self.check_rules = check_rules
        self.options = options

    def _recurseAllowed(self, node):
        for r in self.recurse_rules:
            if not r.allowRecurse(node):
                # print 'NOT recursing into', ci.Cursor_displayname(node), node.kind
                return False  # visiting not allowed
        return True
    
    def _visitAllowed(self, node):
        # if rules._hasFileLocation(node):
        #     print node.location.file.name
        for r in self.recurse_rules:
            if not r.allowVisit(node):
                # print 'NOT visiting', ci.Cursor_displayname(node), node.kind
                return False  # visiting not allowed
        return True
    
    def _recurse(self, node):
        if not self._visitAllowed(node):
            return False  # We did not visit this node.
        self.dispatcher.enterNode(node)
        self.dispatcher.check(node)
        if self._recurseAllowed(node):
            for c in node.get_children():
                self._recurse(c)
        self.dispatcher.exitNode(node)
        return True
    
    def run(self, filename):
        index = ci.Index.create()
        args = ['-I%s' % s for s in self.options.include_dirs]
        # print args
        tu = index.parse(filename, args=args)
        print 'Translation unit: %s.' % tu.spelling
        return self._recurse(tu.cursor)
    
    @classmethod
    def visitFile(klass, filename, options, recurse_rules, check_rules):
        dispatcher = RulesDispatcher(check_rules, options)
        visitor = AstVisitor(dispatcher, recurse_rules, check_rules, options)
        res = visitor.run(filename)
        for k in sorted(dispatcher.violations.keys()):
            print dispatcher.violations[k]
        return res == True


def main():
    # Parse command line arguments.
    parser = optparse.OptionParser("USAGE: %prog [options] file.cpp")
    parser.add_option('-I', '--include-dir', dest='include_dirs', default=[],
                      type='string', help='Specify include directories',
                      action='append')
    options, args = parser.parse_args()
    if len(args) != 1:
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
    return AstVisitor.visitFile(args[0], options, recurse_rules, check_rules)

if __name__ == '__main__':
    sys.exit(main())
