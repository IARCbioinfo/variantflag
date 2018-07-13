#!/bin/sh

data_PATH=~/data/DATA-1/MMB_MEF_WGS_signatureRatio/
grep "#CHROM" ${data_PATH}raw_VCF_mutect/MNU-1_calls.vcf | cut -f 1-5 > ${data_PATH}header.txt
