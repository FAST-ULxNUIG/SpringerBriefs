library(data.table)
library(fda)
library(tidyverse)
library(tikzDevice)
source(here::here("functions", "functional_boxplot_gg.R"))
source(here::here("functions", "theme_gunning.R"))

# Some settings for the Figure: -------------------------------------------
theme_gunning()
theme_update(strip.text = element_text(size = 10),
             axis.text = element_text(size = 9),
             axis.title = element_text(size = 10),
             plot.title = element_text(size = 11))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
plots_path <- here::here("chapter-06", "figures")

# Read in the data:

# -------------------------------------------------------------------------
interpolated_data <- readRDS(here::here("chapter-06", "data", "interpolated-data.rds"))
GRF_dataset_PRO_meta <- interpolated_data$GRF_dataset_PRO_meta
bspl_35 <- interpolated_data$bspl_35

# -------------------------------------------------------------------------
GRF_dataset_PRO_meta_V <- GRF_dataset_PRO_meta[component=="vertical"]
GRF_dataset_PRO_meta_V_fd <- fd(coef = t(as.matrix(GRF_dataset_PRO_meta_V[, paste0("bspl4.", 1:35)])),
                                           basisobj = bspl_35)

GRF_dataset_PRO_meta_ML <- GRF_dataset_PRO_meta[component=="medio_lateral"]
GRF_dataset_PRO_meta_ML_fd <- fd(coef = t(as.matrix(GRF_dataset_PRO_meta_ML[, paste0("bspl4.", 1:35)])),
                                           basisobj = bspl_35)

GRF_dataset_PRO_meta_AP <- GRF_dataset_PRO_meta[component=="anterior_posterior"]
GRF_dataset_PRO_meta_AP_fd <- fd(coef = t(as.matrix(GRF_dataset_PRO_meta_AP[, paste0("bspl4.", 1:35)])),
                                           basisobj = bspl_35)

par(mfrow = c(1, 3))
p1 <- fda::boxplot.fd(GRF_dataset_PRO_meta_V_fd,
                      xlab = "Normalised Time",
                      ylab = "Force (Proportion of Body Weight)",
                      main = "Vertical")
p2 <- fda::boxplot.fd(GRF_dataset_PRO_meta_ML_fd,
                      xlab = "Normalised Time",
                      ylab = "Force (Proportion of Body Weight)",
                      main = "Medio-lateral")
p3 <- fda::boxplot.fd(GRF_dataset_PRO_meta_AP_fd,
                      xlab = "Normalised Time",
                      ylab = "Force (Proportion of Body Weight)",
                      main = "Anterior-Posterior")

V_fboxplot <-create_functional_boxplot_gg(time_grid = 0:100, fd_obj = GRF_dataset_PRO_meta_V_fd)$ggplot + 
  labs(title = "Vertical")
AP_fboxplot <- create_functional_boxplot_gg(time_grid = 0:100, fd_obj = GRF_dataset_PRO_meta_AP_fd)$ggplot +
  labs(title = "Anterior-Posterior")
ML_fboxplot <- create_functional_boxplot_gg(time_grid = 0:100, fd_obj = GRF_dataset_PRO_meta_ML_fd)$ggplot +
  labs(title = "Medio-Lateral")




(combined_plot <- list(V_fboxplot, AP_fboxplot, ML_fboxplot) %>%
lapply(function(x) {
         x + labs(x = "Normalised Time",
                  y = "Force (BW)")
       }) %>%
  ggpubr::ggarrange(plotlist = ., ncol = 3 , nrow = 1))


tikz(file.path(plots_path, "functional_boxplots.tex"),
     width = 1.4 * doc_width_inches, 
     height = 0.475 *  doc_width_inches, standAlone = TRUE)
combined_plot
dev.off()
tinytex::lualatex(file.path(plots_path, "functional_boxplots.tex"))



