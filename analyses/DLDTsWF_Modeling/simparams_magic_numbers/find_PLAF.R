## .................................................................................
## Purpose: Get Shape parameters for beta PLAF
## Notes:
## .................................................................................

vcf <- vcfR::read.vcfR("inst/lowmod_neut_Pf7k_nuc_pass_biallelic_segsites_finalsmpls_diversesites.snpeff.vcf.gz")
PLAF <- rowSums(vcfRmanip::gtmat012(vcf)/2, na.rm = T)/(ncol(vcf@gt) - 1)
hist(PLAF)

betaparams <- fitdistrplus::fitdist(PLAF, "beta")
hist(rbeta(1e3, 0.98, 2.29))
saveRDS(object = betaparams, "analyses/DLDTsWF_Modeling/simparams_magic_numbers/optim_betashapes.RDS")
