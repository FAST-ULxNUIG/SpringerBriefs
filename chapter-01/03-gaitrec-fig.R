# Load Packages: ----------------------------------------------------------
library(tidyverse)  # CRAN v1.3.1
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
source(here::here("functions", "theme_gunning.R"))


# Some settings for the Figure: -------------------------------------------
theme_gunning()
theme_update(strip.text = element_text(size = 10),
             axis.text = element_text(size = 9),
             axis.title = element_text(size = 10),
             plot.title = element_text(size = 11))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
plots_path <- here::here("chapter-01", "figures")



# Read in and Reshape Data: -----------------------------------------------
GRF_dataset_PRO_meta <- readRDS(file = here::here("chapter-06",
                                                  "data",
                                                  "GRF_dataset_PRO_meta.rds"))
# do reshaping with data.table because it's a little easier
GRF_dataset_PRO_meta <- data.table(GRF_dataset_PRO_meta)

GRF_dataset_PRO_meta <- melt.data.table(GRF_dataset_PRO_meta,
                                        measure.vars = paste0("time_", 0:100),
                                        variable.name = "time",
                                        value.name = "force")

GRF_dataset_PRO_meta[, time := as.numeric(stringr::str_remove(time, "time_"))]

GRF_dataset_PRO_meta[, obs_id := paste(SUBJECT_ID, SESSION_ID, TRIAL_ID, sep = "_")]


# Take Sample of Observations ---------------------------------------------
set.seed(1)
random_sample_ids <- sample(unique(GRF_dataset_PRO_meta$obs_id), size = 100)

GRF_dataset_PRO_meta_sample <- GRF_dataset_PRO_meta[obs_id %chin% random_sample_ids]

GRF_dataset_PRO_meta_sample[, side := fifelse(side=="left", "Left Side", "Right Side")]
GRF_dataset_PRO_meta_sample[, component := factor(component,
                                                  levels = c("vertical", "anterior_posterior", "medio_lateral"),
                                                  labels = c("Vertical", "Anterior-Posterior", "Medio-Lateral"))]


# Make Figure and Save: ---------------------------------------------------
(p <- ggplot(data = GRF_dataset_PRO_meta_sample[side == "Left Side"]) +
   aes(x = time, y = force,
       group = obs_id) +
   facet_wrap(~ component, scales = "free_y", nrow = 1) +
   geom_line(alpha = 0.45, colour = "darkgrey") +
   geom_line(data = . %>% filter(SUBJECT_ID == 1003 & TRIAL_ID == 3), colour = "red4", lwd = 0.8) +
   geom_line(data = . %>% filter(SUBJECT_ID == 99 &  TRIAL_ID == 1), colour = "blue4", lwd = 0.8) +
   labs(x = "Normalised Time ($\\%$ of Stance Phase)",
        y = "Force (BW)") +
   theme(legend.position = "none",
         plot.subtitle = element_text(size = 10, hjust = 0.5)))

tikz(file.path(plots_path, "gaitrec-intro-plot.tex"),
     width = 1 * doc_width_inches, 
     height = (1/3) *  doc_width_inches, standAlone = TRUE)
p
dev.off()

tinytex::lualatex(file.path(plots_path, "gaitrec-intro-plot.tex"))
