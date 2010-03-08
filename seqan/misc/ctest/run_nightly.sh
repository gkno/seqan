#!/bin/bash
cd $HOME/Nightly
ctest -S Seqan_Nightly.cmake -VV -d 2>&1 > CTEST_nightly.log
