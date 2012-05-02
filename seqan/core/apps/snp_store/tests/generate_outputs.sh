#!/bin/sh
#
# Output generation script for snp_store

SNP_STORE=../../../../build/Debug/core/apps/snp_store/snp_store

# ============================================================
# First Section
# ============================================================

genome=human-chr22-inf2.fa 
readsGff=human-reads2.gff
readsSam=human-reads2.sam

#echo "${SNP_STORE} $genome $readsGff -o snps_default.out -id indels_default.out > snp_store_default.stdout"
${SNP_STORE} $genome $readsGff -o snps_default.out -id indels_default.out > snp_store_default.stdout

# command options often in use: sam input and do realignment
#echo "${SNP_STORE} $genome $readsSam -if 1 -re -o snps_realign.out -id indels_realign.out > snp_store_realign.stdout"
${SNP_STORE} $genome $readsSam -if 1 -re -o snps_realign.out -id indels_realign.out > snp_store_realign.stdout

# orientation aware and pile up correction, threshold method, hide qualites, indel thershold 1
#echo "${SNP_STORE} $genome $readsSam -if 1 -it 1 -re -oa -mp 1 -m 1 -hq -o snps_realign_m0mp1oa.out -id indels_realign_m0mp1oa.out > snp_store_realign_m0mp1oa.stdout"
${SNP_STORE} $genome $readsSam -if 1 -it 1 -re -oa -mp 1 -m 1 -hq -o snps_realign_m0mp1oa.out -id indels_realign_m0mp1oa.out > snp_store_realign_m0mp1oa.stdout


