library(data.table) # CRAN v1.14.2
library(fda)        # CRAN v5.5.1
library(tidyverse)  # CRAN v1.3.1
library(tikzDevice) # CRAN v0.12.3.1
library(refund)     # CRAN v0.1-26
library(RColorBrewer) # CRAN v1.1-2
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



# Read in Data: -----------------------------------------------------------
data_path <- here::here("chapter-06", "data", "interpolated-data.rds")
interpolated_data <- readRDS(data_path)
GRF_dataset_PRO_meta <- interpolated_data$GRF_dataset_PRO_meta
bspl_35 <- interpolated_data$bspl_35
GRF_dataset_PRO_meta[, uniqueN(SESSION_ID), by = SUBJECT_ID]

# remove discrete values from dataset now we just working with basis coefficients.
GRF_dataset_PRO_meta[, paste0("time_",0:100) := NULL]
GRF_dataset_PRO_averages <- GRF_dataset_PRO_meta[,
                                                 as.list(apply(.SD, 2, mean)), # average basis coefficients of all trials
                                                 by = .(SUBJECT_ID, SESSION_ID, side, component, CLASS_LABEL, CLASS_LABEL_DETAILED, SEX, AGE, HEIGHT, 
                                                        BODY_WEIGHT, BODY_MASS, SHOE_SIZE, AFFECTED_SIDE, SHOD_CONDITION, # defines averaging
                                                        ORTHOPEDIC_INSOLE, SPEED),
                                                 .SDcols = paste0("bspl4.",1:35)] # says which columns to average




# Different Filtering: ----------------------------------------------------
data <- GRF_dataset_PRO_averages[component == "anterior_posterior"]
data <- data[(CLASS_LABEL == "HC" & side == "right") | (AFFECTED_SIDE == 0 & side == "left") | (AFFECTED_SIDE == 1 & side == "right") | (AFFECTED_SIDE == 2 & side == "right")]
(N <- nrow(data))

table(data$CLASS_LABEL)

# create fd object defined by coefficients and basis object
fdobj <- fd(coef = t(as.matrix(data[, paste0("bspl4.",1:35)])),
                                                basisobj = bspl_35)
# and do FPCA:
fpca <- pca.fd(fdobj = fdobj, nharm = 35)
cumsum(fpca$varprop)



# Figure for Anterior-Posterior: ----------------------------------------------------
fpca_eval <- eval.fd(evalarg = 0:100, fpca$harmonics)
fpca_dt <- data.table(time = 0:100, fpca_eval)
fpca_dt_long <- melt.data.table(fpca_dt, id.vars = "time", measure.vars = paste0("PC",1:35))
mean_dt <- data.table(time = 0:100, mean = c(eval.fd(0:100, fpca$meanfd)))
fpca_dt_long <- merge.data.table(x = fpca_dt_long, 
                                 y = mean_dt, 
                                 by = "time",
                                 all.x = TRUE)
var_explained_dt <- data.table(variable = paste0("PC", 1:35), 
                                                  constant = (2 * sqrt(fpca$values))[1:35], # to add to fpcs
                                                  var_explained = paste0("$", round(fpca$varprop * 100, 1), "\\%$"))

fpca_dt_long <- merge.data.table(x = fpca_dt_long, 
                                 var_explained_dt,
                                 by = "variable", 
                                 all.x = TRUE)



fpca_dt_long[, facet_label := paste0("F", variable, " (", var_explained, " of Variance)")]
anterior_posterior_plot <- ggplot(fpca_dt_long[variable %in% paste0("PC", 1:3)]) +
  aes(x = time) +
  facet_wrap(~ facet_label, ncol = 3) +
  geom_line(aes(y = mean)) +
  geom_point(aes(y = mean + constant * value), pch = "+", size = 1.5) +
  geom_point(aes(y = mean - constant * value), pch = "-", size = 3) +
  labs(x = "Normalised Time", y = "Force (BW)", title = "Anterior-Posterior")
anterior_posterior_plot


tikz(file.path(plots_path, "fpcs-case-study-02.tex"),
     width = 1 * doc_width_inches,
     height = 0.35 *  doc_width_inches,
     standAlone = TRUE)
anterior_posterior_plot
dev.off()
tinytex::lualatex(file.path(plots_path, "fpcs-case-study-02.tex"))



# -------------------------------------------------------------------------
boxplots_dt <- data.table(class_label = data$CLASS_LABEL,
                          sex = data$SEX,
                          fpca$scores[, 1:3])
names(boxplots_dt)[-c(1:2)] <- paste0("FPC", 1:3)
boxplots_dt_lng <- melt.data.table(boxplots_dt, id.vars = c("class_label", "sex"))
boxplots_dt_lng[, 
                class_label_fill := factor(class_label, # re-label facets for strip texts.
                       levels = c("HC", "A", "K", "H", "C"),
                       labels = c("Healthy Control", "Ankle", "Knee", "Hip", "Calcaneous"))]
boxplots_dt_lng[, class_label := factor(class_label, levels = c("HC", "A", "K", "H", "C"))]

(boxplot_scores <- ggplot(boxplots_dt_lng) +
  aes(x = class_label, y = value, fill = class_label_fill) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_grid(~ variable) +
  geom_boxplot() +
  labs(y = "FPC Score", x = "Impairement", fill = "Impairement:") +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold")))

tikz(file.path(plots_path, "fpc-scores-case-study-02.tex"),
     width = 1 * doc_width_inches,
     height = 0.45 *  doc_width_inches,
     standAlone = TRUE)
boxplot_scores
dev.off()
tinytex::lualatex(file.path(plots_path, "fpc-scores-case-study-02.tex"))


saveRDS(list(data = data,
             fd_obj = fdobj),
        file = here::here("chapter-06", "data",
        "function-on-scalar-data.rds"))

