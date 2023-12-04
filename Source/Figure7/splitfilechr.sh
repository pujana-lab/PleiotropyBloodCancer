#!/bin/bash

for i in {1..22}
do plink2 --bfile EUR_phase3_autosomes --chr $i --make-bed --out EUR_phase3_chr$i
done

