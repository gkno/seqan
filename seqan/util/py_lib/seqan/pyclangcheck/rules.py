#!/usr/bin/env python

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os.path
import re

import clang.cindex as ci

RULE_NAMING_CONSTANT = 'naming.constant'
RULE_NAMING_STRUCT = 'naming.struct'
RULE_NAMING_UNION = 'naming.union'
RULE_NAMING_CLASS = 'naming.class'
RULE_NAMING_ENUM = 'naming.enum'
RULE_NAMING_FIELD = 'naming.field'
RULE_NAMING_ENUM_CONSTANT = 'naming.enum_constant'
RULE_NAMING_VARIABLE = 'naming.variable'
RULE_NAMING_FUNCTION = 'naming.function'
RULE_NAMING_PARAMETER = 'naming.parameter'
RULE_NAMING_VARIABLE = 'naming.variable'
RULE_NAMING_CXX_METHOD = 'naming.method'
RULE_NAMING_TYPEDEF = 'naming.typedef'
RULE_NAMING_TPL_NON_TYPE_PARAMETER = 'naming.tpl_nontype_param'
RULE_NAMING_TPL_TYPE_PARAMETER = 'naming.tpl_type_param'
RULE_NAMING_TPL_TPL_PARAMETER = 'naming.tpl_tpl_param'
RULE_NAMING_FUNCTION_TPL = 'naming.function_tpl'
RULE_NAMING_CLASS_TPL = 'naming.class_tpl'
RULE_NAMING_CLASS_TPL_SPEC = 'naming.class_tpl_spec'

RE_CONSTANT = r'^[A-Z]([_A-Z0-9])*_*$'
RE_VARIABLE = r'^_?[a-z]([_a-zA-Z0-9])*_*$'
RE_FUNCTION = r'operator.*|^_?[a-z]([_a-zA-Z0-9])*\(.*\)_*$'
RE_TYPE = r'^[A-Z]([a-zA-Z0-9])*_*$'
RE_TYPE_TEMPLATE = r'^[A-Z]([a-zA-Z0-9])*_*<.*>$'
RE_STRUCT = r'^[A-Z]([a-zA-Z0-9])*_*(<.*>)?$'

RULE_TEXTS = {
    RULE_NAMING_CONSTANT: 'Constant names must be all upper-case, separated by underscores.',
    RULE_NAMING_STRUCT: 'Struct names must be all camel case, starting with upper case.',
    RULE_NAMING_UNION: 'Union names must be all camel case, starting with upper case.',
    RULE_NAMING_CLASS: 'Class names must be all camel case, starting with upper case.',
    RULE_NAMING_ENUM: 'Enum names must be all camel case, starting with upper case.',
    RULE_NAMING_FIELD: 'Field names must be camel case case, starting with lower case.',
    RULE_NAMING_ENUM_CONSTANT: 'Enum constant names must be all upper-case, separated by underscores.',
    RULE_NAMING_VARIABLE: 'Variable names must be camel case case, starting with lower case.',
    RULE_NAMING_FUNCTION: 'Function names must be camel case case, starting with lower case.',
    RULE_NAMING_PARAMETER: 'Parameter names must be camel case case, starting with lower case.',
    RULE_NAMING_CXX_METHOD: 'Method names must be camel case case, starting with lower case.',
    RULE_NAMING_TYPEDEF: 'Typedef names must be all camel case, starting with upper case.',
    RULE_NAMING_TPL_NON_TYPE_PARAMETER: 'Template non-type parameters must be all upper-case, separated by underscores.',
    RULE_NAMING_TPL_TYPE_PARAMETER: 'Template type parameter names must be all camel case, starting with upper case.',
    RULE_NAMING_TPL_TPL_PARAMETER: 'Template template parameter names must be all camel case, starting with upper case.',
    RULE_NAMING_FUNCTION_TPL: 'Function template names must be camel case case, starting with lower case.',
    RULE_NAMING_CLASS_TPL: 'Class template names must be all camel case, starting with upper case.',
    RULE_NAMING_CLASS_TPL_SPEC: 'Partial specialization names must be all camel case, starting with upper case.',
}


class RuleViolation(object):
    def __init__(self, rule_id, violator, file, line, column):
        self.rule_id = rule_id
        self.violator = violator
        self.file = file
        self.line = line
        self.column = column
    
    def key(self):
        return (self.rule_id, self.line, self.column, self.violator)
    
    def __str__(self):
        msg = '%s [%s:%d/%d] "%s": %s'
        return msg % (self.rule_id, self.file, self.line, self.column,
                      self.violator, RULE_TEXTS[self.rule_id])


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


class AllYesRule(object):
    """A rule that allows all visiting and checks always evaluate to True."""
    def allowVisit(self, node):
        return True

    def allowRecurse(self, node):
        return True
    
    def check(self, node):
        return []


class GenericSymbolNameRule(AllYesRule):
    def __init__(self, kind, regular_ex, rule_name):
        self.kind = kind
        self.regular_ex = regular_ex
        self.rule_name = rule_name

    def allowVisit(self, node):
        if not _hasFileLocation(node):
            return False
        displayname = ci.Cursor_displayname(node)
        if not displayname:
            return False  # Ignore empty symbols.
        # print 'allow visit template type?', displayname, node.kind
        if node.kind == self.kind:
            return True
        return False
    
    def check(self, node):
        displayname = ci.Cursor_displayname(node)
        print 'checking', displayname
        import pdb; pdb.set_trace()
        if not re.match(self.regular_ex, displayname):
            v = RuleViolation(
                self.rule_name, displayname, node.location.file.name,
                node.location.line, node.location.column)
            return [v]
        return []


class InIncludeDirsRule(AllYesRule):
    """Rule to block visiting and recursion outside include dirs."""
    
    def __init__(self, include_dirs):
        self.include_dirs = [os.path.abspath(x) for x in include_dirs]
        self.cache = {}
    
    def allowVisit(self, node):
        """Return True if visiting is allowed."""
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
    
    def allowRecurse(self, node):
        """Return True if we want to recurse below node."""
        return self.allowVisit(node)
