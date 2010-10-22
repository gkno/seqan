#!/usr/bin/env python

import sys
import re

from helpers import *

def replace_all(text, subst):
    for old in subst.keys():
        text = old.sub(subst[old], text)

    return text


def validate_file(file, subst):
    code = ''
    with open(file, 'r') as f:
        code = f.read()

    old_len = len(code)

    replaced = replace_all(code, subst)
    assert old_len == len(replaced)

    open(file, 'w').write(replaced)


def build_subst_table(file):
    f = open(file, 'r')
    table = {}

    for line in f:
        old, new = line.rstrip('\r\n').split(': ')
        table[re.compile(r'\b%s\b' % old)] = new

    return table


def main():
    project_path = sys.argv[1]
    replacements_file = sys.argv[2]

    substitutions = build_subst_table(replacements_file)

    for file in all_files(project_path):
        validate_file(file, substitutions)


if __name__ == '__main__':
    sys.exit(main())
