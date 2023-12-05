# Decompress pgen.zst to pgen 
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

# "vzs" modifier to directly operate with pvar.zst
# "--chr 1-22" excludes all variants not on the listed chromosomes
# "--output-chr 26" uses numeric chromosome codes
# "--max-alleles 2": PLINK 1 binary does not allow multi-allelic variants
# "--rm-dup" removes duplicate-ID variants
# "--set-missing-var-id" replaces missing IDs with a pattern
plink2 --pfile all_phase3 vzs \
       --chr 1-22,X \
       --output-chr 26 \
       --max-alleles 2 \
       --rm-dup exclude-mismatch \
       --set-missing-var-ids '@_#_$1_$2' \
       --make-pgen \
       --out autosomes/all_phase3_autosomes

# Prepare sub-population filter file
awk 'NR == 1 || $5 == "EUR" {print $1}' all_phase3.psam > EUR/EUR_1kg_samples.txt

# Generate sub-population fileset
plink2 --pfile autosomes/all_phase3_autosomes \
       --keep EUR/EUR_1kg_samples.txt \
       --make-pgen \
       --out EUR/EUR_phase3_autosomes

# pgen to bed
plink2 --pfile EUR/EUR_phase3_autosomes \
       --maf 0.01 \
       --snps-only \
       --make-bed \
       --out EUR/EUR_phase3_autosomes
       
# Split bed/bim/fam by chromosome
for i in {1..23}
do plink2 --bfile EUR/EUR_phase3_autosomes --chr $i --make-bed --out EUR/EUR_phase3_chr$i
done

