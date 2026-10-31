# Reproduce all results and figures. Run from the repository root:
#   Rscript run_all.R
source("R/setup.R")
source("R/model.R")
source("R/data.R")
source("R/calibration.R")
source("R/simulate.R")
source("R/figures.R")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
message(">>> Done.")
