#!/usr/bin/env python
"""Code related to violations and suppressions."""

from __future__ import with_statement

import os
import os.path
import sys

import rules


class RuleViolation(object):
    def __init__(self, rule_id, violator, file, line, column):
        self.rule_id = rule_id
        self.violator = violator
        self.file = file
        self.line = line
        self.column = column
    
    def key(self):
        return (self.file, self.line, self.column, self.rule_id, self.violator)
    
    def __str__(self):
        msg = '[%s:%d/%d] %s "%s": %s'
        return msg % ('/'.join(self.file.split('/')[-2:]), self.line, self.column,
                      self.rule_id, self.violator, rules.RULE_TEXTS[self.rule_id])


class NolintManager(object):
    """Manage the lines ending in '//nolint'."""

    def __init__(self):
        self.locations = {}

    def hasNolint(self, filename, lineno):
        filename = os.path.abspath(filename)
        # Ensure that the nolint lines are registered in self.locations[filename].
        if not self.locations.has_key(filename):
            line_set = set()
            with open(filename, 'rb') as f:
                line_no = 0
                for line in f:
                    line_no += 1
                    if line.strip().endswith('// nolint'):
                        ## print 'nolint', filename, line_no
                        line_set.add(line_no)
            self.locations[filename] = line_set
        # Query self.locations[filename].
        return lineno in self.locations[filename]
