#!/bin/sh
rm dddoc_cache.bin
rm -rf html/*
../util/bin/dddoc.py -d concepts -d pages ../core ../extras $@
exit $?
