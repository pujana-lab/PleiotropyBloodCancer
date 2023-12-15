# Decompress pgen.zst to pgen 
../plink2dir/plink2 --zst-decompress ../1000G/all_phase3.pgen.zst > ../1000G/all_phase3.pgen

# "vzs" modifier to directly operate with pvar.zst
# "--chr 1-22" excludes all variants not on the listed chromosomes
# "--output-chr 26" uses numeric chromosome codes
# "--max-alleles 2": PLINK 1 binary does not allow multi-allelic variants
# "--rm-dup" removes duplicate-ID variants
# "--set-missing-var-id" replaces missing IDs with a pattern
mkdir autosomes
../plink2dir/plink2 --pfile ../1000G/all_phase3 vzs \
       --chr 1-22,X,Y,MT \
       --output-chr 26 \
       --max-alleles 2 \
       --rm-dup exclude-mismatch \
       --set-missing-var-ids '@_#_$1_$2' \
       --make-pgen \
       --new-id-max-allele-len 662 \
       --out autosomes/all_phase3_autosomes

# Prepare sub-population filter file
mkdir EUR
awk 'NR == 1 || $5 == "EUR" {print $1}' ../1000G/all_phase3.psam > EUR/EUR_1kg_samples.txt

# Generate sub-population fileset
../plink2dir/plink2 --pfile autosomes/all_phase3_autosomes \
       --keep EUR/EUR_1kg_samples.txt \
       --make-pgen \
       --out EUR/EUR_phase3_autosomes

# pgen to bed (MAF>0.05)
../plink2dir/plink2 --pfile EUR/EUR_phase3_autosomes \
       --maf 0.05 \
       --snps-only \
       --make-bed \
       --out EUR/EURmaf0_05_phase3_autosomes
       
# Split bed/bim/fam by chromosome (MAF>0.05)
mkdir EURmaf0_05
for i in {1..24}
do ../plink2dir/plink2 --bfile EUR/EURmaf0_05_phase3_autosomes --chr $i --make-bed --out EURmaf0_05/EURmaf0_05_phase3_chr$i
done
../plink2dir/plink2 --bfile EUR/EURmaf0_05_phase3_autosomes --chr 26 --make-bed --out EURmaf0_05/EURmaf0_05_phase3_chr26

# pgen to bed (MAF>0.01)
../plink2dir/plink2 --pfile EUR/EUR_phase3_autosomes \
       --maf 0.01 \
       --snps-only \
       --make-bed \
       --out EUR/EURmaf0_01_phase3_autosomes
       
# Split bed/bim/fam by chromosome (MAF>0.01)
mkdir EURmaf0_01
for i in {1..24}
do ../plink2dir/plink2 --bfile EUR/EURmaf0_01_phase3_autosomes --chr $i --make-bed --out EURmaf0_01/EURmaf0_01_phase3_chr$i
done
../plink2dir/plink2 --bfile EUR/EURmaf0_01_phase3_autosomes --chr 26 --make-bed --out EURmaf0_01/EURmaf0_01_phase3_chr26
