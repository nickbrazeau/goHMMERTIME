## .................................................................................
## Purpose: Bring together results from HMMERTIME run
##
##
## Notes:
## .................................................................................

library(tidyverse)
#............................................................
# read in DLDTsWF Simulations and trueIBD
#...........................................................
simtrueIBD <- readRDS("results/tidied_results/DLDTsWF_simulation_runs_trueIBD.RDS")

#............................................................
# read in results
#...........................................................
ret.paths <- tibble::tibble(simnum = sub(".hmmertime_ret.RDS", "", list.files("results/cluster_results/run_HMMERTIME",
                                                                pattern = "sim_[0-9]+.hmmertime_ret.RDS")),
                            paths = list.files("results/cluster_results/run_HMMERTIME",
                                               pattern = "sim_[0-9]+.hmmertime_ret.RDS", full.names = T))

# get hmmertime results
read_hmmertime_quants <- function(path) {
  ret <- readRDS(path)
  return(ret$mcmcout[[1]]$summary$quantiles)
}

ret.paths$hmmertimequants <- purrr::map(ret.paths$paths, read_hmmertime_quants)

#......................
# save out
#......................
dplyr::left_join(simtrueIBD, ret.paths, by = "simnum") %>%
  dplyr::select(-c("paths")) %>%
  saveRDS(., "results/tidied_results/hmmertime_results_from_dtdlswfsims.RDS")
