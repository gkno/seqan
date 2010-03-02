#!/bin/bash
cd $HOME/Nightly
/opt/local/bin/ctest -S Mac_Nightly.cmake -VV -d 2>&1 > CTEST_nightly.log

