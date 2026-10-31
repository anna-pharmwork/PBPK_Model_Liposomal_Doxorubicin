# Model definition: fixed parameters and the 6-state ODE system -----------------
# Units: amounts in ug, volumes in mL, concentrations in ug/mL, time in h.
# Species: 25 g mouse, 5 mg/kg IV dose.

params_fixed <- c(
  V_pl          = 0.85,    # mL     plasma volume
  V_heart       = 0.125,   # mL     myocardium volume
  V_tum         = 0.100,   # mL     tumour volume (not used in the equations)
  Q_heart       = 30.45,   # mL/h   plasma flow through the heart
  Kp_heart      = 3.5,     # -      heart/plasma partition coefficient (adjusted)
  CL_tot        = 35.0,    # mL/h   systemic clearance of free DOX (adjusted)
  Dose          = 125.0,   # ug     5 mg/kg x 0.025 kg

  # Starting values of the kinetic constants (overwritten by calibration,
  # except where noted in R/calibration.R)
  Vmax          = 5.0,     # ug/h   max MPS uptake rate
  Km            = 15.0,    # ug/mL  Michaelis constant of MPS uptake
  k_pt          = 0.02,    # 1/h    plasma -> tumour transfer of liposomes (EPR)
  k_tp          = 0.005,   # 1/h    tumour -> plasma transfer of liposomes
  k_tp_DOX      = 0.01,    # 1/h    free DOX tumour -> plasma
  k_deg         = 0.02,    # 1/h    liposome degradation in MPS (releases DOX)
  k_release     = 0.015,   # 1/h    DOX release from liposomes in plasma
  k_release_tum = 0.015    # 1/h    DOX release from liposomes in tumour
)

kinetic_names <- c("Vmax", "Km", "k_pt", "k_tp", "k_tp_DOX",
                   "k_deg", "k_release", "k_release_tum")

# State vector (all amounts in ug):
#   A_NP_pl   liposomal DOX in plasma     A_DOX_pl  free DOX in plasma
#   A_MPS     liposomal DOX in MPS        A_heart   free DOX in heart
#   A_tum_NP  liposomal DOX in tumour     A_tum     free DOX in tumour
pbpk_system <- function(t, state, parameters) {

  model_values <- c(as.list(state), as.list(parameters))

  required_parameters <- c(
    "V_pl", "V_heart", "V_tum", "Q_heart", "Kp_heart", "CL_tot",
    "Vmax", "Km", "k_pt", "k_tp", "k_tp_DOX", "k_deg",
    "k_release", "k_release_tum"
  )
  missing_parameters <- setdiff(required_parameters, names(model_values))
  if (length(missing_parameters) > 0) {
    stop("Missing parameters: ", paste(missing_parameters, collapse = ", "))
  }

  with(model_values, {

    # Guard against tiny negative values from solver round-off
    A_NP_pl_safe  <- max(0, A_NP_pl)
    A_MPS_safe    <- max(0, A_MPS)
    A_tum_NP_safe <- max(0, A_tum_NP)
    A_DOX_pl_safe <- max(0, A_DOX_pl)
    A_heart_safe  <- max(0, A_heart)
    A_tum_safe    <- max(0, A_tum)

    # Concentrations, ug/mL
    C_NP_pl   <- A_NP_pl_safe / V_pl
    C_DOX_pl  <- A_DOX_pl_safe / V_pl
    C_heart   <- A_heart_safe / V_heart
    C_v_heart <- C_heart / Kp_heart

    # Saturable liposome uptake by the MPS, ug/h
    v_mps <- Vmax * C_NP_pl / (Km + C_NP_pl)

    # Perfusion-limited exchange of free DOX plasma <-> heart, ug/h
    J_heart <- Q_heart * (C_DOX_pl - C_v_heart)

    dA_NP_pl  <- - k_release * A_NP_pl_safe - v_mps -
                   k_pt * A_NP_pl_safe + k_tp * A_tum_NP_safe
    dA_MPS    <- v_mps - k_deg * A_MPS_safe
    dA_tum_NP <- k_pt * A_NP_pl_safe - k_tp * A_tum_NP_safe -
                 k_release_tum * A_tum_NP_safe
    dA_DOX_pl <- k_release * A_NP_pl_safe + k_deg * A_MPS_safe +
                 k_tp_DOX * A_tum_safe -
                 (CL_tot / V_pl) * A_DOX_pl_safe - J_heart
    dA_heart  <- J_heart
    dA_tum    <- k_release_tum * A_tum_NP_safe - k_tp_DOX * A_tum_safe

    list(c(dA_NP_pl, dA_MPS, dA_tum_NP, dA_DOX_pl, dA_heart, dA_tum))
  })
}

# Initial states: liposomal (Doxil-like) dose vs free-DOX bolus
y0_pld <- c(A_NP_pl = as.numeric(params_fixed["Dose"]),
            A_MPS = 0, A_tum_NP = 0, A_DOX_pl = 0, A_heart = 0, A_tum = 0)
y0_free <- c(A_NP_pl = 0, A_MPS = 0, A_tum_NP = 0,
             A_DOX_pl = as.numeric(params_fixed["Dose"]), A_heart = 0, A_tum = 0)
