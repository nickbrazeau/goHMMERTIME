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
#remotes::install_github("nickbrazeau/HMMERTIME")
library(HMMERTIME)
library(tidyverse)
library(drake)

#............................................................
# make sim map
#...........................................................
simruns <- list.files("results/cluster_results/DLDTsWFsims/", pattern = "sim_[0-9]+.RDS",
                      full.names = TRUE)
simmap <- tibble::tibble(nm = sub(".RDS", "", basename(simruns)),
                         path = simruns)

#............................................................
# HMMERTIME run function
#...........................................................
dir.create("/pine/scr/n/f/nfb/Projects/goHMMERTIME/run_HMMERTIME/",
           recursive = TRUE)
run_HMMERTIME <- function(path) {
  simret <- readRDS(path)

  # hmmertime run
  ret <- HMMERTIME::runMCMC(vcfRobj = simret$vcfRobj, # vcfR object we simulated
                            vcfploid = 2, # ploidy of VCF
                            PLAF = 1-simret$WSAF.list$rbetaPLAF,
                            m_max = 15, # max COI to consider
                            rho = 7.4e-7, # recombination rate
                            k_max = 25, # max switch rate to consider, NB it affect poisson of layering on ibd lvls
                            e1 = 0.05, # error for going from homozygous to heterozygous
                            e2 = 0.05, # error for going from heterozygous to homozygous
                            burnin = 1e4,
                            samples = 1e4,
                            reportIteration = 1e3,
                            verbose = TRUE,
                            parallelize = TRUE)

  #......................
  # save out on my local slurm machine -- future users will need to change this
  #......................
  # out
  outpath = paste0("/pine/scr/n/f/nfb/Projects/goHMMERTIME/run_HMMERTIME/",
                   sub(".RDS", "", basename(path)), ".hmmertime_ret.RDS")
  saveRDS(ret, file = outpath)
  return(0)
}

#............................................................
# Make Drake Plan
#...........................................................
plan <- drake::drake_plan(
  fits = target(
    run_HMMERTIME(path),
    transform = map(
      .data = !!simmap
    )
  )
)


#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "drake_clst/slurm_clustermq_LL.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = nrow(simmap),
     log_make = "hmmertime_deploy_drake.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")

