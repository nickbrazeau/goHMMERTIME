## .................................................................................
## Purpose: Tidy up and par down the results from the dldtswf
##
## Notes:
## .................................................................................
library(tidyverse)
#............................................................
# read in DLDTsWF simulations
#...........................................................
# run maps from dldtswf sims
runmap <- readRDS("results/cluster_results/DLDTsWFsims/param_map.RDS")
# results from dldtswf
run_rets <- tibble::tibble(simnum = sub(".RDS", "", list.files("results/cluster_results/DLDTsWFsims/",
                                                               pattern = "sim_[0-9]+.RDS")),
                           paths = list.files("results/cluster_results/DLDTsWFsims/",
                                              pattern = "sim_[0-9]+.RDS", full.names = T))
read_in_dldtswf_trueIBD <- function(path) {
  ret <- readRDS(path)$realized
  return(ret)
}
# get trueIBD
run_rets <- run_rets %>%
  dplyr::mutate(realized = purrr::map(paths, read_in_dldtswf_trueIBD))

#......................
# bring together and tidy up
# dldtswf simulations
#......................
out <- dplyr::left_join(runmap, run_rets, by = "simnum") %>%
  dplyr::select(c("N", "mean_coi", "m", "lvl", "simnum", "realized"))

# save out
saveRDS(out, "results/tidied_results/DLDTsWF_simulation_runs_trueIBD.RDS")
