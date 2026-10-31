# Calibration ---------------------------------------------------------------------
# Objective: mean squared error in log10 space for two profiles:
#   (1) total DOX in plasma after Doxil          (liposomal + free, model arm "pld")
#   (2) total DOX in plasma after free DOX bolus (model arm "free")
# Optimised (log-scale): Vmax, Km, k_pt, k_tp, k_tp_DOX, k_deg, k_release
# (k_release_tum is tied to k_release).

cost_function_balanced <- function(opt_pars, return_parts = FALSE) {

  fixed_only <- params_fixed[!names(params_fixed) %in% kinetic_names]

  pars <- c(
    setNames(as.numeric(fixed_only), names(fixed_only)),
    Vmax          = unname(exp(opt_pars[["log_Vmax"]])),
    Km            = unname(exp(opt_pars[["log_Km"]])),
    k_pt          = unname(exp(opt_pars[["log_kpt"]])),
    k_tp          = unname(exp(opt_pars[["log_ktp"]])),
    k_tp_DOX      = unname(exp(opt_pars[["log_ktp_dox"]])),
    k_deg         = unname(exp(opt_pars[["log_kdeg"]])),
    k_release     = unname(exp(opt_pars[["log_krel"]])),
    k_release_tum = unname(exp(opt_pars[["log_krel"]]))
  )

  if (any(!is.finite(pars)) ||
      any(pars[c("V_pl", "V_heart", "V_tum", "Kp_heart", "Km")] <= 0) ||
      pars["CL_tot"] < 0 || pars["Q_heart"] < 0) {
    return(1e20)
  }

  times <- sort(unique(c(0, exp_data_plasma_NP$time, exp_data_plasma_free$time)))

  result <- tryCatch({
    list(
      pld  = as.data.frame(deSolve::ode(y = y0_pld,  times = times, func = pbpk_system, parms = pars)),
      free = as.data.frame(deSolve::ode(y = y0_free, times = times, func = pbpk_system, parms = pars))
    )
  }, error = function(e) NULL)
  if (is.null(result)) return(1e20)

  total_pl <- function(o) (pmax(o$A_NP_pl, 0) + pmax(o$A_DOX_pl, 0)) / pars["V_pl"]

  pred_NP   <- approx(result$pld$time,  total_pl(result$pld),  xout = exp_data_plasma_NP$time)$y
  pred_free <- approx(result$free$time, total_pl(result$free), xout = exp_data_plasma_free$time)$y

  if (any(!is.finite(pred_NP)) || any(!is.finite(pred_free))) return(1e20)

  # Floor keeps log10() finite; it does not fix the model, it only allows a
  # finite error to be computed.
  epsilon <- 1e-8
  err_NP   <- mean((log10(pmax(pred_NP,   epsilon)) - log10(exp_data_plasma_NP$C_obs))^2)
  err_free <- mean((log10(pmax(pred_free, epsilon)) - log10(exp_data_plasma_free$C_obs))^2)

  if (return_parts) {
    return(list(err_NP = err_NP, err_free = err_free,
                pred_NP = pred_NP, pred_free = pred_free))
  }

  ssr <- err_NP + err_free   # equal weight of both profiles
  if (!is.finite(ssr)) return(1e20)
  ssr
}

initial_calib <- c(
  log_Vmax = log(2.0), log_Km = log(20.0), log_kpt = log(0.05), log_ktp = log(0.01),
  log_ktp_dox = log(0.04), log_kdeg = log(0.03), log_krel = log(0.003)
)

message(">>> Calibration (Nelder-Mead, single start) ...")
fit <- optim(initial_calib, cost_function_balanced, method = "Nelder-Mead",
             control = list(maxit = 5000, reltol = 1e-8))

# Calibrated parameter set on the natural scale
fixed_calibrated <- params_fixed[!names(params_fixed) %in% kinetic_names]
calibrated_pars <- c(
  fixed_calibrated,
  Vmax          = unname(exp(fit$par[["log_Vmax"]])),
  Km            = unname(exp(fit$par[["log_Km"]])),
  k_pt          = unname(exp(fit$par[["log_kpt"]])),
  k_tp          = unname(exp(fit$par[["log_ktp"]])),
  k_tp_DOX      = unname(exp(fit$par[["log_ktp_dox"]])),
  k_deg         = unname(exp(fit$par[["log_kdeg"]])),
  k_release     = unname(exp(fit$par[["log_krel"]])),
  k_release_tum = unname(exp(fit$par[["log_krel"]]))
)
stopifnot(!anyDuplicated(names(calibrated_pars)))

# Fit diagnostics ------------------------------------------------------------------
parts <- cost_function_balanced(fit$par, return_parts = TRUE)
fe    <- parts$pred_NP / exp_data_plasma_NP$C_obs      # fold error, Doxil arm

fit_metrics <- data.frame(
  metric = c("optim_convergence_code", "objective_total", "log10_MSE_doxil_arm",
             "log10_MSE_free_arm", "doxil_median_fold_error",
             "doxil_fraction_FE_0.5_to_2", "doxil_fraction_FE_0.3_to_3",
             "doxil_RMSE_ug_per_mL", "doxil_MAPE_percent"),
  value  = c(fit$convergence, fit$value, parts$err_NP, parts$err_free, median(fe),
             mean(fe > 0.5 & fe < 2), mean(fe > 0.3 & fe < 3),
             sqrt(mean((parts$pred_NP - exp_data_plasma_NP$C_obs)^2)),
             mean(abs(parts$pred_NP - exp_data_plasma_NP$C_obs) / exp_data_plasma_NP$C_obs) * 100)
)

dir.create("results", showWarnings = FALSE)
write.csv(data.frame(parameter = names(calibrated_pars), value = as.numeric(calibrated_pars)),
          "results/calibrated_parameters.csv", row.names = FALSE)
write.csv(fit_metrics, "results/fit_metrics.csv", row.names = FALSE)
print(fit_metrics, digits = 4)
