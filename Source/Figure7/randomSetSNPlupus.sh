#!/bin/bash
# plink --bfile EUR_phase3_autosomes --recode --out EUR_phase3_autosomes_snp
# cut -f 2 EUR_phase3_autosomes_snp.map > snps.map
# grep -v -w -f SNPpleio.txt snps.map > snpsNOpleio.map
DIRIN=../../Data/Fig7b
DIROUT=../../Data/Fig7b/randomSetLupus
for ((i=1;i<=1020;i++));
do
	date
	FILEOUT=$DIROUT/randomSet$i.txt
	echo $i
	shuf -n 917 $DIRIN/snpsNOpleio.map > $FILEOUT
done
