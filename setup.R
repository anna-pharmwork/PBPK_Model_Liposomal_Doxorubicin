# Packages and plotting theme ---------------------------------------------------
required_packages <- c("deSolve", "ggplot2", "dplyr", "scales")
missing_packages  <- required_packages[!(required_packages %in% rownames(installed.packages()))]
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

theme_thesis <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 9.5, hjust = 0.5, color = "gray30"),
      axis.title       = element_text(face = "bold", size = 11),
      legend.title     = element_text(face = "bold", size = 10),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
}
