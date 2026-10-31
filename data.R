# Digitized calibration data ----------------------------------------------------
# Source: Luo R. et al., Int J Pharm 519 (2017) 1-10, Fig. 1A (mouse plasma,
# 5 mg/kg IV, total doxorubicin). Digitized with WebPlotDigitizer.
# NOTE: points are digitization nodes along the plotted curves, not independent
# observations (the paper reports 13 sampling times, n = 3 per time).

load_digitized <- function(path) {
  d <- read.csv(path)
  stopifnot(all(c("time_h", "conc_ng_ml") %in% names(d)))
  d %>%
    transmute(time = as.numeric(time_h),
              C_obs_ng_ml = as.numeric(conc_ng_ml),
              C_obs = C_obs_ng_ml / 1000) %>%     # ng/mL -> ug/mL
    arrange(time)
}

exp_data_plasma_NP   <- load_digitized("data/luo2017_fig1A_doxil_plasma_digitized.csv")     # Doxil
exp_data_plasma_free <- load_digitized("data/luo2017_fig1A_free_dox_plasma_digitized.csv")  # free DOX

stopifnot(
  nrow(exp_data_plasma_NP) == 69L, nrow(exp_data_plasma_free) == 75L,
  !anyNA(exp_data_plasma_NP$time), !anyNA(exp_data_plasma_free$time),
  all(exp_data_plasma_NP$C_obs > 0), all(exp_data_plasma_free$C_obs > 0)
)
