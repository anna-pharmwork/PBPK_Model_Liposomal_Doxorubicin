# Simulations with calibrated parameters ---------------------------------------
# Dense grid early on (0-2 h) to capture peaks.
times_sim <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2, seq(2.5, 72, by = 0.5))

run_ode <- function(y0, pars, atol, rtol) {
  as.data.frame(ode(y = y0, times = times_sim, func = pbpk_system, parms = pars,
                    method = "lsoda", atol = atol, rtol = rtol, maxsteps = 1000000))
}

sim_base <- run_ode(y0_pld, calibrated_pars, atol = 1e-4, rtol = 1e-4)
sim_base$C_NP_pl <- sim_base$A_NP_pl / calibrated_pars["V_pl"]
sim_base$C_heart <- sim_base$A_heart / calibrated_pars["V_heart"]

sim_free <- run_ode(y0_free, calibrated_pars, atol = 1e-10, rtol = 1e-8)
sim_free$C_DOX_pl <- sim_free$A_DOX_pl / calibrated_pars["V_pl"]
sim_free$C_heart  <- sim_free$A_heart  / calibrated_pars["V_heart"]

# One-at-a-time sensitivity to k_release (Fig 3.3: time courses) ---------------
k_rel_values <- c(0.001, 0.003, 0.008, 0.020, 0.050)
df_sens_plot <- bind_rows(lapply(k_rel_values, function(kr) {
  p <- calibrated_pars
  p["k_release"] <- kr; p["k_release_tum"] <- kr
  o <- run_ode(y0_pld, p, atol = 1e-4, rtol = 1e-4)
  data.frame(time = o$time, C_heart = o$A_heart / calibrated_pars["V_heart"],
             k_rel = paste0("k_rel = ", kr, " h\u207b\u00b9"))
}))

# Cmax(heart) vs k_release (Fig 3.4) ------------------------------------------
k_range <- seq(0.005, 0.50, by = 0.02)
cmax_results <- data.frame(k_release = k_range, Cmax_heart = NA_real_)
for (i in seq_along(k_range)) {
  p <- calibrated_pars
  p["k_release"] <- k_range[i]; p["k_release_tum"] <- k_range[i]
  o <- run_ode(y0_pld, p, atol = 1e-4, rtol = 1e-4)
  cmax_results$Cmax_heart[i] <- max(o$A_heart / calibrated_pars["V_heart"])
}
write.csv(cmax_results, "results/cmax_heart_vs_k_release.csv", row.names = FALSE)

# Peak heart concentration, free DOX bolus vs liposomal dose
heart_cmax <- data.frame(
  formulation = c("free_DOX_bolus", "liposomal_calibrated"),
  Cmax_heart_ug_per_mL = c(max(sim_free$C_heart), max(sim_base$C_heart))
)
write.csv(heart_cmax, "results/heart_cmax_free_vs_liposomal.csv", row.names = FALSE)
