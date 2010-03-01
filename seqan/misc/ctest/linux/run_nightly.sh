#!/bin/bash
cd $HOME/Nightly
ctest -S Linux_Nightly.cmake -VV -d 2>&1 | tee CTEST_nightly.log
