# Runs the full analysis. Open the .Rproj file in RStudio, then run this script.
# Required packages: install.packages(c("dlnm", "ggplot2"))

source("R/01_descriptive.R")
source("R/02_main_model.R")
source("R/03_sensitivity.R")
source("R/04_dlnm.R")

writeLines(capture.output(sessionInfo()), "results/session_info.txt")
