# ========================================================================#
#  - Create figure that shows a picture of curves varying
#  only in amplitude and phase.
#  - Simulate the data using the registr package for simulated curves
#=========================================================================#

# 1) Import Packages ------------------------------------------------------
library(tidyverse)  # CRAN v1.3.0
library(registr)    # CRAN v1.0.0
library(tikzDevice) # CRAN v0.12.3.1


# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-02", "figures")

source(file.path(functions_path, "theme_gunning.R"))
theme_gunning()
theme_update(strip.text = element_text(size = 11),
             axis.text = element_text(size = 10),
             axis.title = element_text(size = 11))


# 3) Simulation  ----------------------------------------------
set.seed(1996)
# Generate functions differing only in amplitude and phase
Amplitude <- simulate_unregistered_curves(I = 20, lambda = 8, D = 100, phase_variation = F)
Phase <- simulate_unregistered_curves(I = 20, lambda = 0, D = 100, phase_variation = T)
# checks:
head(Amplitude)
head(Phase)
stopifnot(dim(Amplitude) == dim(Phase))


# 4) Data Wrangling -------------------------------------------------------

# Select only variables needed
Amplitude <- Amplitude %>% select(id, index, amp_var = latent_mean)
Phase <- Phase %>% select(id, index, phase_var = latent_mean)

# Join them up
amp_phase_data <- inner_join(Amplitude, Phase, by = c("id", "index")) %>%
  rename("Amplitude Variation" = 3, "Phase Variation" = 4) %>%
  gather(-id, -index, key = "var_type", value = "x_t")

# Calculate a mean to overlay on the data
amp_phase_mean <- amp_phase_data %>%
  group_by(index, var_type) %>%
  summarise(mean_fun = mean(x_t)) %>%
  ungroup()


# 5) Plot -----------------------------------------------------------------
# Plot data
amp_phase_plot <- ggplot(amp_phase_data) +
  aes(x = index, y = x_t, group = id, col = factor(id)) +
  facet_wrap( ~ var_type) +
  geom_line(alpha = 0.75) +
  geom_line(inherit.aes = F, data = amp_phase_mean,
            aes(x = index, y =mean_fun), size = 1.5) +
  scale_x_continuous(expand = c(0, 0.05)) +
  labs(y = "$x(t)$", x = "$t$") +
  theme(legend.position = "none")

# View Plot
amp_phase_plot

# Save Plots
# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "amp-phase-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 0.5 *  doc_width_inches,
     standAlone = TRUE)
amp_phase_plot
dev.off()

tinytex::lualatex(file.path(plots_path, "amp-phase-plot.tex"))
