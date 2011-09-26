#!/usr/bin/env python
"""Fix gcov output."""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>'

import os.path
import sys

def main():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'py_lib'))
    sys.path.insert(0, path)
    import seqan.fixgcov
    return seqan.fixgcov.main()

if __name__ == '__main__':
    sys.exit(main())
