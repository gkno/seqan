#!/bin/sh
rm html/*
./main.py ../projects/library/seqan -d concepts -d pages $@
exit $?
