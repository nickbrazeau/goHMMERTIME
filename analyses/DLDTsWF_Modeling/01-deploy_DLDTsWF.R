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

#............................................................
# Running the Simulation
#...........................................................
get_truth_from_pairwise_arg <- function(arg, this_coi){
  # arg doesn't store host information so need to carry COI information
  if (!length(this_coi) == 2) {
    stop("must be a pairwise comparison")
  }
  # get connections
  conn <- purrr::map(arg, "c")
  # get timing of connections
  tm <- purrr::map(arg, "t")
  # find pairwise
  get_pairwise_ibd <- function(conni, tmi, this_coi) {
    smpl1con <- conni[1:this_coi[1]]
    smpl2con <- conni[(this_coi[1]+1):(cumsum(this_coi)[2])]

    #......................
    # get IBD
    #......................
    # connections between 1 and 2
    pwconn <- which(smpl2con %in% 0:(this_coi[1]-1) )
    locimatches <- rep(1, length(pwconn))
    # note we are 0 based in connections
    # note bvtrees always point left
    # catch if there are multiple matches within sample 2 to the pairwise
    # this is a coalescent true that looks like below if host COI is 2,2
    # c: -1 -1 1 2
    # t: -1 -1 5 1
    if (length(pwconn) != 0) {
      for (i in 1:length(pwconn)) {
        haplotypeindex <- this_coi[1] + pwconn[i] - 1 # -1 for 0-based
        internalconn <- which(smpl2con %in% haplotypeindex )
        if (length(internalconn) != 0) {
          for (i in 1:length(internalconn)) {
            internalhaplotypeplace <- this_coi[1] + internalconn[i] # here 1-based in R
            if (tmi[internalhaplotypeplace] < tmi[this_coi[1] + internalconn[i]]) { # here 1-based in R
              locimatches[i] <- locimatches[i] + 1
            }
          }
        }
      }
    }
    #......................
    # calculating within sample IBD as the number of strains within a host that are identical
    # this means there are COI - 1 strains that can be identical
    # this arises in two scenarios in bvtrees (always point left):
    #       (1) samples coalesce within the same host at some T
    #       (2) samples coalesce to the same strain in a different host at different time Ts
    #......................
    # within sample1 is easy since we start at 0
    withinIBD_smpl1 <- sum(smpl1con %in% (1:this_coi[1]-1)) + sum(duplicated(smpl1con[smpl1con != -1]))
    # within sample2 adjust slightly for "offset"
    withinIBD_smpl2 <- sum(smpl2con %in% this_coi[1]:(length(conni)-1)) + sum(duplicated(smpl2con[smpl2con != -1]))

    # return
    out <- list(pairwiseIBD = sum(locimatches),
                withinIBD_smpl1 = withinIBD_smpl1,
                withinIBD_smpl2 = withinIBD_smpl2,
                eff_pairwiseIBD = eff_pairwiseIBD)
    return(out)
  }
  # calculate
  numerator <- purrr::map2(.x = conn, .y = tm,
                           .f = get_pairwise_ibd, this_coi = this_coi)

  #......................
  # IBD between
  # If between IBD exceed the minimum COI, then
  # cap at minimum COI (i.e. the between amount of IBD cannot be greater than the
  # number of strains that are within the smallest host COI)
  #......................
  pairwiseIBDvec <- purrr::map_dbl(numerator, "pairwiseIBD")
  pairwiseIBDvec <- ifelse(pairwiseIBDvec > min(this_coi), min(this_coi), pairwiseIBDvec)
  pairwiseIBD <- sum(pairwiseIBDvec)/(min(this_coi) * length(conn)) # min combn * num Loci

  #......................
  # within
  #......................
  # -1 here for the SELF comparison
  withinIBD_smpl1 <- sum(purrr::map_dbl(numerator, "withinIBD_smpl1")) / ((this_coi[1]-1) * length(conn))
  withinIBD_smpl2 <- sum(purrr::map_dbl(numerator, "withinIBD_smpl2")) / ((this_coi[2]-1) * length(conn))

  # catch when MOI = 1 and no within possible, so set to 0
  withinIBD_smpl1 <- ifelse(is.nan(withinIBD_smpl1), 0, withinIBD_smpl1)
  withinIBD_smpl2 <- ifelse(is.nan(withinIBD_smpl2), 0, withinIBD_smpl2)

  #......................
  # effective between
  #......................
  get_effective_coi <- function(arg, hostcoi_index) {
    # get connections
    conn <- purrr::map(arg, "c")
    # get connections for this specific host
    conn <- lapply(conn, function(x)return(x[hostcoi_index]))
    conn <- unique(conn)
    connmat <- matrix(NA, ncol = length(hostcoi_index), nrow = length(conn))
    for (i in 1:nrow(connmat)) {
      connmat[i,] <- conn[[i]]
    }
    # look to see if all loci are coalesced w/in and only w/in for each strain
    clonecount <- apply(connmat, 2, function(x) {all(x %in% (hostcoi_index-1))})
    return(length(hostcoi_index) - sum(clonecount))
  }
  # run
  effcoi1 <- get_effective_coi(arg = arg, hostcoi_index = 1:this_coi[1])
  effcoi2 <- get_effective_coi(arg = arg, hostcoi_index = (this_coi[1]+1):sum(this_coi))
  eff_coi <- c(effcoi1, effcoi2)

  # return
  ret <- list(pairwiseIBD = pairwiseIBD,
              eff_coi = eff_coi,
              withinIBD_smpl1 = withinIBD_smpl1,
              withinIBD_smpl2 = withinIBD_smpl2)
  return(ret)
}


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
    # get True IBD
    #......................
    trueIBD <- get_truth_from_pairwise_arg(arg = ARG, this_coi = this_coi)

    #............................................................
    # catch interesting
    #...........................................................
    if (lvl == "interest" & any(trueIBD$pairwiseIBD %in% c(0,1)) ) {
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
  out <- list(trueIBD = trueIBD,
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


