## .................................................................................
## Purpose: Prepare PfCross for hmmIBD
##
## Notes:
## .................................................................................
library(tidyverse)
library(vcfR)
library(vcfRmanip)

#............................................................
# read in and process cross
#...........................................................
# bcftools view SNP_INDEL_WG.combined.filtered.crosses.vcf.gz --include 'FILTER="PASS" && N_ALT=1 && TYPE="snp" && AC>0' --output-type z --output-file SNP_INDEL_WG.combined.filtered.crosses.seg.biallelicsnps.vcf.gz
pfcross <- vcfR::read.vcfR("vcfdata/MalariaGen_Crosses/3d7_hb3.combined.final.seg.biallelicsnps.vcf.gz")

pfcrosschrompos <- tibble::tibble(CHROM = vcfR::getCHROM(pfcross),
                                  POS = vcfR::getPOS(pfcross))
pfcrossgt <- vcfRmanip::gtmat012(pfcross)
pfcrossgt[pfcrossgt == 1] <- -1
pfcrossgt[pfcrossgt == 2] <- 1
pfcrossgt[is.na(pfcrossgt)] <- -1

pfout <- cbind.data.frame(pfcrosschrompos, pfcrossgt)
pfout$CHROM <- as.integer(as.factor(pfout$CHROM)) # coerce on purpose
pfout$POS <- as.integer(pfout$POS)

#......................
# write out
#......................
readr::write_tsv(pfout, file = "analyses/pfcross_contrast/hmmibd_inputs/PfCross.tab.txt")

#......................
# run hmmIBD
#......................
outdir <- "analyses/pfcross_analysis/hmmibd_outputs/"
dir.create(outdir, recursive = TRUE)
setwd(outdir)
system("../../run_hmmIBD/copy-hmmIBD-master/hmmIBD -i ../hmmibd_inputs/PfCross.tab.txt -o PfCross_3d7hb3_hmmibd_out")
# clean up
setwd(here::here())
