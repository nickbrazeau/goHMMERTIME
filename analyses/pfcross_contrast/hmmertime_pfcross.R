## .................................................................................
## Purpose: Run PfCross w/ HMMERTIME
##
## Notes:
## .................................................................................
library(vcfR)
library(HMMERTIME)

#............................................................
# read in and process cross
#...........................................................
# bcftools view 3d7_hb3.combined.final.vcf.gz --include 'FILTER="PASS" && N_ALT=1 && TYPE="snp" && AC>0' --output-type z --output-file S3d7_hb3.combined.final.seg.biallelicsnps.vcf.gz
pfcross <- vcfR::read.vcfR("vcfdata/MalariaGen_Crosses/3d7_hb3.combined.final.seg.biallelicsnps.vcf.gz")
pfcross_subset <- vcfRmanip::select_samples(pfcross, smplvctr = c("3D7/PG0051-C/ERR019061", "HB3/PG0052-C/ERR019054",
                                                                  "C04/PG0061-C/ERR019059", "C02/PG0053-C/ERR019067"))
# hmmertime run
ret <- HMMERTIME::runMCMC(vcfRobj = pfcross_subset, # vcfR object we simulated
                          vcfploid = 1, # ploidy of VCF
                          m_max = 10, # max COI to consider
                          rho = 7.4e-7, # recombination rate
                          k_max = 10, # max switch rate to consider
                          e1 = 0.05, # error for going from homozygous to heterozygous
                          e2 = 0.05, # error for going from heterozygous to homozygous
                          burnin = 1e4,
                          samples = 1e4,
                          reportIteration = 1e3,
                          verbose = TRUE,
                          parallelize = TRUE)
#......................
# save out
#......................
dir.create("analyses/pfcross_analysis/HMMERTIME_outputs/", recursive = TRUE)
saveRDS(ret, "PfCross_3d7hb3_HMMERTIME_out.RDS")
