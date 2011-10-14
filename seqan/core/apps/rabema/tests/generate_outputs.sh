#!/bin/sh
#
# Expected output generation for rabema.

# We use the current trunk version of 2011-10-13 (r10463) for building the
# reference.
RAZERS=../../../../build/Release/core/apps/razers2/razers2
RABEMA=../../../../build/Release/core/apps/rabema/rabema

# ============================================================
# Map reads for gold standard and with too low error rate.
# ============================================================

${RAZERS} -m 10000 -vv -of 4 -ds -i 92     -o gold-adeno-hamming-08.sam adeno-genome.fa reads.fasta
${RAZERS} -m 10000 -vv -of 4 -ds -i 92 -id -o gold-adeno-edit-08.sam    adeno-genome.fa reads.fasta

${RAZERS} -vv -of 4 -ds -i 92     -o razers2-adeno-hamming-08.sam adeno-genome.fa reads.fasta
${RAZERS} -vv -of 4 -ds -i 92 -id -o razers2-adeno-edit-08.sam    adeno-genome.fa reads.fasta

${RAZERS} -vv -of 4 -ds -i 96     -o razers2-adeno-hamming-04.sam adeno-genome.fa reads.fasta
${RAZERS} -vv -of 4 -ds -i 96 -id -o razers2-adeno-edit-04.sam    adeno-genome.fa reads.fasta

# ============================================================
# Build Gold Standard
# ============================================================

${RABEMA} build_standard -d hamming -e 8 -o gold-adeno-hamming-08.wit adeno-genome.fa gold-adeno-hamming-08.sam > gold-adeno-hamming-08.stdout
${RABEMA} build_standard -d edit    -e 8 -o gold-adeno-edit-08.wit    adeno-genome.fa gold-adeno-edit-08.sam > gold-adeno-edit-08.stdout

# ============================================================
# Compare Against Gold Standard
# ============================================================

${RABEMA} compare -d hamming -e 8 adeno-genome.fa razers2-adeno-hamming-08.sam gold-adeno-hamming-08.wit > razers2-adeno-hamming-08.stdout
${RABEMA} compare -d hamming -e 8 adeno-genome.fa razers2-adeno-hamming-04.sam gold-adeno-hamming-08.wit > razers2-adeno-hamming-04.stdout
${RABEMA} compare -d edit    -e 8 adeno-genome.fa razers2-adeno-edit-08.sam    gold-adeno-edit-08.wit > razers2-adeno-edit-08.stdout
${RABEMA} compare -d edit    -e 8 adeno-genome.fa razers2-adeno-edit-04.sam    gold-adeno-edit-08.wit > razers2-adeno-edit-04.stdout
