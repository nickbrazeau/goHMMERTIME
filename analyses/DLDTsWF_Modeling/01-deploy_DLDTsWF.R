## .................................................................................
## Purpose: Run a set of genetic simulations under the Discrete-Loci Discrete-Time
##          structured Wright Fisher Model for Plasmodium falciparum genetics (PMC7192906).
##          with
##
## Date: 18 November, 2020
##
## Notes: Transmission/COI occurs through the mosquito,
##        which is the unit for the binomial sampling -- see the polySimIBD vignettes
##        for more details
## .................................................................................
set.seed(48)
#remotes::install_github("nickbrazeau/polySimIBD")
library(polySimIBD)
library(furrr)
library(tidyverse)

#............................................................
# Magic Numbers
#...........................................................
# Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
# Aimee gets this number by taking the inverse of Mile's estimate of the CO recombination rate of 13.5 kb/cM
rho <- 7.4e-7
# going to assume we can only detect things 10 generations ago
tlim <- 10

# approximate average of Pf3d7 Chromosome Lengths
pflen <- 1.664e6

# assuming single chromosome for ease
# assuming 5e2 loci
pos <- sort(sample(1.664e6, 5e2))
CHROMPOS <- tibble::tibble(CHROM = "CHROM1",
                           POS = pos)

# Deme Size: we will consider hosts (N) on a logarithmic-base 10 scale
N <- round(10^seq(1, 3, l = 11))
# Mean COIs for lambda from Verity, Aydemir, Brazeau et al. 2020 (PMC7192906)
coilamdas <- readRDS("analyses/DLDTsWF_Modeling/simparams_magic_numbers/optim_lambda.RDS")
# various levels of M -- "migration" which drives mosquito-deme exchange
M <- c(0, 0.25, 0.5, 1)

#............................................................
# Parameter Map
#...........................................................
# expand out combinations
paramsdf <- expand.grid(N, coilamdas, M) %>%
  tibble::as_tibble() %>%
  magrittr::set_colnames(c("N", "mean_coi", "m"))

# when M = 0, there is no relatedness because there is no genetic exchange
# when Mean COI: fixed at 1, there is going to be either an "all" or nothing phenomenon
# due to the binomial sampling in WF
# so split out for and run more interations of the "interesting" paramsets
m1 <- unique(paramsdf$mean_coi)[1]
baseparamsdf <- paramsdf %>%
  dplyr::filter(m == 0 | mean_coi == m1) %>%
  dplyr::mutate(lvl = "base")
interestparamsdf <- paramsdf %>%
  dplyr::filter(m != 0) %>%
  dplyr::filter(mean_coi != m1) %>%
  dplyr::mutate(lvl = "interest")
# 10 more iteration of these more interesting ones
baseparamsdf <- lapply(1:10, function(x){baseparamsdf}) %>%
  dplyr::bind_rows()

# 100 more iteration of these more interesting ones
interestparamsdf <- lapply(1:100, function(x){interestparamsdf}) %>%
  dplyr::bind_rows()

# combine
paramsdf <- dplyr::bind_rows(baseparamsdf, interestparamsdf) %>%
  dplyr::arrange(N, mean_coi, m)

# add details for run sim
paramsdf <- paramsdf %>%
  dplyr::mutate(pos = list(pos),
                rho = rho,
                tlim = tlim,
                hosts = list(1:2),
                simnum = paste0("sim_", 1:dplyr::n()))

#............................................................
# save parameter map out
#...........................................................
dir.create("results/cluster_results/DLDTsWFsims/", recursive = TRUE)
saveRDS(paramsdf,
        "results/cluster_results/DLDTsWFsims/param_map.RDS")


run_dldtswf <- function(N, mean_coi, m, lvl, pos, rho, tlim, hosts, simnum) {

  # make interesting
  make_btwness <- TRUE
  while (make_btwness) {

    #......................
    # run structured WF
    #......................
    swfsim <- polySimIBD::sim_swf(pos =       pos,
                                  N =         N,
                                  m =         m,
                                  rho =       rho,
                                  mean_coi =  mean_coi,
                                  tlim =      tlim)
    # extract ARG and down sample arg
    ARG <- polySimIBD::get_arg(swfsim, host_index = hosts)
    this_coi <- swfsim$coi[hosts]
    # extract Haplotype Matrix from arg
    hapmat <- polySimIBD::get_haplotype_matrix(ARG)

    #......................
    # get realized results
    #......................
    pairwiseIBD <- polySimIBD::get_realized_pairwise_ibd(swf = swfsim, host_index = hosts)
    wthnIBD1 <- polySimIBD::get_within_host_IBD(swf = swfsim, host_index = hosts[[1]])
    wthnIBD2 <- polySimIBD::get_within_host_IBD(swf = swfsim, host_index = hosts[[2]])
    effCOI1 <- polySimIBD::get_realized_coi(swf = swfsim, host_index = hosts[[1]])
    effCOI2 <- polySimIBD::get_realized_coi(swf = swfsim, host_index = hosts[[2]])
    # out
    realized <- list(
      pairwiseIBD = pairwiseIBD,
      withinIBD_host1 = wthnIBD1,
      withinIBD_host2 = wthnIBD2,
      effCOI_host1 = effCOI1,
      effCOI_host2 = effCOI2)

    #............................................................
    # catch interesting
    #...........................................................
    if (lvl == "interest" & any(pairwiseIBD %in% c(0,1)) ) {
      # if interesting, re-run until there is some between IBD
      make_btwness <- TRUE
    } else {
      make_btwness <- FALSE # break loop as goal is accomplished or is known uninteresting
    }
  } # end while loop catch

  #......................
  # simulate reads
  #......................
  WSAF.list <- polySimIBD::sim_biallelic(COIs = this_coi,
                                         haplotypematrix = hapmat,
                                         shape1 = 0.9836757, # magic number from spRF
                                         shape2 = 2.1892612, # magic number from spRF
                                         coverage = 100,
                                         alpha = 1,
                                         overdispersion = 1e-10,
                                         epsilon = 0.025)

  #......................
  # convert to vcfRobj
  #......................
  meta <- paste("##fileformat=VCFv4.3", "##Simulated with polySimIBD", "##ploidy=2", collapse = "\n")

  # get Fix:  must be in this order and only these
  fix <- data.frame(CHROM = "CHROM1",
                    POS = pos,
                    ID = NA,
                    REF = "A", # arbitrary choice
                    ALT = "T", # arbitrary choice
                    QUAL = NA,
                    FILTER = "PASS",
                    INFO = NA)
  # get GT
  # here we will need to make a choice on the "threshold" for "calling" a homozyg-ref, heterozyg, or homozyg-alt allele
  gtmatsim <- matrix("0/1", nrow = nrow(WSAF.list$NRWSAcounts), ncol = ncol(WSAF.list$NRWSAcounts))
  gtmatsim[WSAF.list$NRWSAF > 0.95] <- "1/1"
  gtmatsim[WSAF.list$NRWSAF < 0.05] <- "0/0"
  # now that we have genotype calls, we can combine these into depths
  gt <- matrix(NA, nrow = nrow(gtmatsim), ncol = ncol(gtmatsim))
  # quick for loop to protect against vector setting
  for (i in 1:nrow(gtmatsim)) {
    for (j in 1:ncol(gtmatsim)) {
      AD <- paste(WSAF.list$WS.coverage[i,j] - WSAF.list$NRWSAcounts[i,j], WSAF.list$NRWSAcounts[i,j], sep = ",")
      gt[i, j] <- paste(gtmatsim[i,j], AD, WSAF.list$WS.coverage[i,j], sep = ":")
    }
  }
  gt <- cbind(FORMAT = "GT:AD:DP", gt)

  # write out new vcfRobj
  require(vcfR) # need this for class
  vcfRobj <- new("vcfR", meta = meta, fix = as.matrix(fix), gt = gt)

  #......................
  # send out
  #......................
  out <- list(realized = realized,
              ARG = ARG,
              swfsim = swfsim,
              hosts = hosts,
              this_coi = this_coi,
              WSAF.list = WSAF.list,
              vcfRobj = vcfRobj)

  #......................
  # save out on my local slurm machine -- future users will need to change this
  #......................
  # out
  outpath = paste0("results/cluster_results/DLDTsWFsims/",
                   simnum, ".RDS")
  saveRDS(out, file = outpath)
  return(0)
}


#............................................................
# run
#...........................................................
furrr::future_pmap(paramsdf, run_dldtswf)


