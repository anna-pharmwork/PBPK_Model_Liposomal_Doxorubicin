# Figures ---------------------------------------------------------------------------
ug_ml <- "\u00b5g/mL"
dir.create("figures", showWarnings = FALSE)

# Fig 3.1: calibration, liposomal (encapsulated) DOX in plasma ----------------------
plot_fig1 <- ggplot() +
  geom_line(data = sim_base, aes(x = time, y = C_NP_pl, color = "PBPK model"), linewidth = 1.1) +
  geom_point(data = exp_data_plasma_NP, aes(x = time, y = C_obs, shape = "Experiment (Doxil)"),
             size = 2.5, color = "navy") +
  scale_y_log10() +
  scale_color_manual(name = "", values = c("PBPK model" = "navy")) +
  scale_shape_manual(name = "", values = c("Experiment (Doxil)" = 16)) +
  labs(title = "A. Encapsulated DOX/Doxil (plasma)",
       x = "Time (h)", y = paste0("Concentration (", ug_ml, ", log)")) +
  theme_thesis()

# Fig 3.2: heart kinetics, free DOX bolus vs liposomal dose -------------------------
df_compare_heart <- rbind(
  data.frame(time = sim_free$time, Concentration = sim_free$C_heart, Formulation = "Free DOX (IV bolus)"),
  data.frame(time = sim_base$time, Concentration = sim_base$C_heart, Formulation = "PLD (Doxil / Caelyx)")
) %>% filter(Concentration > 0)   # log axis: drop zero/negative round-off values

plot_fig2 <- ggplot(df_compare_heart, aes(x = time, y = Concentration, color = Formulation)) +
  geom_line(linewidth = 1.1) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c("Free DOX (IV bolus)" = "#D95F02", "PLD (Doxil / Caelyx)" = "#1B9E77")) +
  labs(title = "Fig. 3.2. Comparison of kinetics in the myocardium (semi-logarithmic scale)",
       subtitle = "Reduction of peak cardiac concentration (Cmax) for PLD",
       x = "Time after administration (h)",
       y = paste0("DOX concentration in the heart (", ug_ml, ", log)"),
       color = "Dosage form:") +
  theme_thesis()

# Fig 3.3: sensitivity to k_release, heart time courses -----------------------------
plot_fig3 <- ggplot(df_sens_plot, aes(x = time, y = C_heart, color = k_rel)) +
  geom_line(linewidth = 1.1) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Fig. 3.3. Influence of the release rate (k_release) on accumulation in the myocardium",
       subtitle = "Influence of the release rate of liposomal DOX on its concentration in the myocardium",
       x = "Time (h)", y = paste0("DOX concentration in the heart (", ug_ml, ")"),
       color = "Release rate constant k_release:") +
  theme_thesis()

# Fig 3.4: Cmax(heart) vs k_release -------------------------------------------------
plot_fig4 <- ggplot(cmax_results, aes(x = k_release, y = Cmax_heart)) +
  geom_line(color = "#2C7BB6", linewidth = 1.3) +
  geom_point(color = "#2C7BB6", size = 2.5) +
  labs(title = "Fig. 3.4. Dependence of the peak concentration in the myocardium (Cmax) on k_release",
       subtitle = "Dependence of peak myocardial DOX concentration on the release rate",
       x = "Release rate constant, k_release (h\u207b\u00b9)",
       y = paste0("Peak concentration in the heart, Cmax (", ug_ml, ")")) +
  theme_thesis()

message(">>> Saving figures ...")
ggsave("figures/Figure_3_1_Calibration_Validation.png", plot_fig1, width = 9.5, height = 4.8, dpi = 300)
ggsave("figures/Figure_3_2_Free_vs_PLD_Heart.png",      plot_fig2, width = 7.5, height = 4.8, dpi = 300)
ggsave("figures/Figure_3_3_Sensitivity_Dynamic.png",    plot_fig3, width = 7.5, height = 4.8, dpi = 300)
ggsave("figures/Figure_3_4_Cmax_vs_krelease.png",       plot_fig4, width = 7.5, height = 4.8, dpi = 300)
