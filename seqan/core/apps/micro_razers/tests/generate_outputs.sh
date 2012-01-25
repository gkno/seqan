#!/bin/sh
#
# Output generation script for micro_razers

MICRO_RAZERS=../../../../build/Release/core/apps/micro_razers/micro_razers

# ============================================================
# First Section
# ============================================================



${MICRO_RAZERS} adeno-genome.fa adeno-reads36_1.fa -o se-adeno-reads36_1_default.out > se-adeno-reads36_1_default.stdout

for sl in 14 15 16 17 18 19 20 ; do
	${MICRO_RAZERS} -sL $sl adeno-genome.fa adeno-reads36_1.fa -o se-adeno-reads36_1_sl${sl}.out > se-adeno-reads36_1_sl${sl}.stdout
	${MICRO_RAZERS} -sL $sl -sE adeno-genome.fa adeno-reads36_1.fa -o se-adeno-reads36_1_sl${sl}_se.out > se-adeno-reads36_1_sl${sl}_se.stdout
done


${MICRO_RAZERS} -sL 18 -m 20 -pa adeno-genome.fa adeno-reads36_1.fa -o se-adeno-reads36_1_sl18_m20_pa.out > se-adeno-reads36_1_sl18_m20_pa.stdout
