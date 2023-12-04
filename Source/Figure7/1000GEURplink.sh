# https://www.cog-genomics.org/plink/2.0/
# https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
# https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html
# Decompress pgen.zst to pgen 
# (dataExt/1000G)
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

# "vzs" modifier to directly operate with pvar.zst
# "--chr 1-22" excludes all variants not on the listed chromosomes
# "--output-chr 26" uses numeric chromosome codes
# "--max-alleles 2": PLINK 1 binary does not allow multi-allelic variants
# "--rm-dup" removes duplicate-ID variants
# "--set-missing-var-id" replaces missing IDs with a pattern
# (dataExt/1000Gplink/autosomes)
plink2 --pfile ../1000G/all_phase3 vzs \
       --chr 1-22 \
       --output-chr 26 \
       --max-alleles 2 \
       --rm-dup exclude-mismatch \
       --set-missing-var-ids '@_#_$1_$2' \
       --make-pgen \
       --out autosomes/all_phase3_autosomes

# Prepare sub-population filter file
# (dataExt/1000Gplink/EUR)
awk 'NR == 1 || $5 == "EUR" {print $1}' all_phase3.psam > EUR_1kg_samples.txt

# Generate sub-population fileset
# (dataExt/1000Gplink/EUR)
plink2 --pfile autosomes/all_phase3_autosomes \
       --keep EUR/EUR_1kg_samples.txt \
       --make-pgen \
       --out EUR/EUR_phase3_autosomes

plink2 --pfile EUR_phase3_autosomes --ld rs2493214 rs35731977

# pgen to bed
# "--maf 0.005" remove most monomorphic SNPs 
# (still may have some when all samples are heterozyguous -> maf=0.5)
# (dataExt/1000Gplink/EURmaf0_05)
plink2 --pfile EUR/EUR_phase3_autosomes \
       --maf 0.05 \
       --snps-only \
       --make-bed \
       --out EURmaf0_05/EURmaf0_05_phase3_autosomes
       
# Split bed/bim/fam by chromosome
for i in {1..22}
do plink2 --bfile EURmaf0_05/EURmaf0_05_phase3_autosomes --chr $i --make-bed --out EURmaf0_05/EURmaf0_05_phase3_chr$i
done

plink --bfile EUR_phase3_chr1 --ld rs2493214 rs35731977
