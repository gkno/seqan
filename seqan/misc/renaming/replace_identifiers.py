#!/usr/bin/env python

import sys
import re

from helpers import *

PROGRAM_USAGE = """
SeqAn script to replace invalid identifiers (previously collected) in the SeqAn
codebase.

USAGE: replace_identifiers.py BASE_PATH REPLACEMENTS

BASE_PATH is the root path of all the folders to be searched.
REPLACEMENTS is a file of `key:value' pairs which contain the invalid
    identifier and the replacement string.
""".strip()

def replace_all(text, subst):
    """
    Perform the substitutions given by the dictionary ``subst`` on ``text``.
    """
    for old in subst.keys():
        text = old.sub(subst[old], text)

    return text


def validate_file(file, subst):
    """
    Perform the substitutions given by the dictionary ``subst`` on ``file``.
    """
    code = ''
    with open(file, 'r') as f:
        code = f.read()

    old_len = len(code)

    replaced = replace_all(code, subst)
    assert old_len == len(replaced)

    open(file, 'w').write(replaced)


def build_subst_table(file):
    """
    Read the substitutions defined in ``file`` and build a substitution table.
    """
    f = open(file, 'r')
    table = {}

    for line in f:
        old, new = line.rstrip('\r\n').split(':')
        table[re.compile(r'\b%s\b' % old.strip())] = new.strip()

    return table


def main():
    if len(sys.argv) != 3:
        print >>sys.stderr, 'ERROR: Invalid number of arguments.'
        print >>sys.stderr, PROGRAM_USAGE
        return 1

    project_path = sys.argv[1]
    replacements_file = sys.argv[2]

    substitutions = build_subst_table(replacements_file)

    for file in all_files(project_path):
        validate_file(file, substitutions)

    return 0


if __name__ == '__main__':
    sys.exit(main())
