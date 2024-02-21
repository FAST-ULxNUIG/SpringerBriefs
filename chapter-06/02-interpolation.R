# Load Packages: ----------------------------------------------------------
library(tidyverse)  # CRAN v1.3.1
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(fda)        # CRAN v5.5.1
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


# Read in and Reshape Data: -----------------------------------------------
GRF_dataset_PRO_meta <- readRDS(file = here::here("chapter-06",
                                                  "data",
                                                  "GRF_dataset_PRO_meta.rds"))
nrow(GRF_dataset_PRO_meta)


k_seq <- seq(15, 80, by = 5)
SSE_vec <- vector("numeric", length = length(k_seq))
y <- t(GRF_dataset_PRO_meta[, paste0("time_", 0:100)])
SSE_mat <- matrix(NA, nrow = ncol(y), ncol = length(k_seq))
for(k in seq_along(k_seq)) {
  # loop through different values of k
  print(paste("iteration", k, "of", length(k_seq)))
  bspl_k <- create.bspline.basis(rangeval = c(0, 100), nbasis = k_seq[k], norder = 4)
  GRF_dataset_PRO_meta_smooth_basis <- smooth.basis(argvals = 0:100,
                                                    y = y,
                                                    fdParobj = bspl_k)
  yhat <- predict(GRF_dataset_PRO_meta_smooth_basis)
  SSE_mat[,k] <- apply((y - yhat)^2, 2, sum) # calculate (discrete approx. to) integrated squared error
}

apply(SSE_mat, 2, mean)
plot(k_seq, apply(SSE_mat, 2, mean))
names(SSE_mat) <- k_seq
boxplot(sqrt(SSE_mat),
        xlab = "Number of Basis Functions",
        ylab = "Root SSE",
        main = "Approximation Error from Interpolation")



# Publication Plot: -------------------------------------------------------
SSE_dt <- data.table(SSE_mat)
colnames(SSE_dt) <- paste(k_seq)
SSE_dt$obs <- factor(seq_len(nrow(SSE_dt)))
SSE_dt_lng <- melt.data.table(data = SSE_dt,
                              id.vars = "obs",
                              measure.vars = paste(k_seq),
                              variable.name = "k", 
                              value.name = "SSE")

(choose_k_plot <- ggplot(data=SSE_dt_lng) +
  aes(x = k, y = sqrt(SSE), colour = (k==35)) +
  geom_boxplot(outlier.size = 0.4) +
  scale_color_manual(values = c("black", "darkgreen")) +
  labs(x = "Number of Basis Functions",
       y = "Root SSE",
       title = "Approximation Error from Interpolation") +
  theme(legend.position = "none",
        plot.margin = margin(r = 40, l = 40)) +
  geom_vline(xintercept = 35, lty = 2, colour = "darkgreen"))

tikz(file.path(plots_path, "interpolation-plot-trial-K.tex"),
     width = 0.75 * doc_width_inches, 
     height = 0.5 *  doc_width_inches, standAlone = TRUE)
choose_k_plot
dev.off()

tinytex::lualatex(file.path(plots_path, "interpolation-plot-trial-K.tex"))

# from this, choose basis

# Do basis representation: ------------------------------------------------
bspl_35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 35, norder = 4)
GRF_dataset_PRO_meta_smooth_basis <- smooth.basis(argvals = 0:100,
                                                  y = y,
                                                  fdParobj = bspl_35)
yhat <- predict(GRF_dataset_PRO_meta_smooth_basis)


GRF_dataset_PRO_meta <- data.table(GRF_dataset_PRO_meta)
GRF_dataset_PRO_meta[, obs_id := paste(SUBJECT_ID, SESSION_ID, TRIAL_ID, sep = "_")]

t_finer <- seq(0, 100, by = 0.025)
y_hat_finer <- eval.fd(evalarg = t_finer, GRF_dataset_PRO_meta_smooth_basis$fd)
set.seed(1)
random_sample_obs_id <- sample(GRF_dataset_PRO_meta$obs_id, 3)
GRF_dataset_PRO_meta_random_sample <- GRF_dataset_PRO_meta[obs_id %in% random_sample_obs_id]
yhat_random_sample <- t(y_hat_finer[, GRF_dataset_PRO_meta[,which(obs_id %in% random_sample_obs_id)]])
colnames(yhat_random_sample) <- paste0("yhat_time_", t_finer)
GRF_dataset_PRO_meta_random_sample <- cbind(GRF_dataset_PRO_meta_random_sample, yhat_random_sample)
GRF_dataset_PRO_meta_random_sample_lng <- melt.data.table(GRF_dataset_PRO_meta_random_sample,
                                        measure.vars = c(paste0("time_", 0:100), paste0("yhat_time_", 0:100)),
                                        variable.name = "time",
                                        value.name = "force")
GRF_dataset_PRO_meta_random_sample_lng[, object := fifelse(stringr::str_detect(time, "yhat"),
                                                           "yhat", "y")]

GRF_dataset_PRO_meta_random_sample_lng[, time := stringr::str_remove(time, "(|yhat_)time_")]
GRF_dataset_PRO_meta_random_sample_lng[, time := as.numeric(time)]

GRF_dataset_PRO_meta_random_sample_lng[, component := factor(component,
                                                  levels = c("vertical", "anterior_posterior", "medio_lateral"),
                                                  labels = c("Vertical", "Anterior-Posterior", "Medio-Lateral"))]
fitted_vs_observed_plot <- ggplot(data = GRF_dataset_PRO_meta_random_sample_lng[side=="left"]) +
  aes(x = time, y = force, group = obs_id, colour = obs_id) +
  facet_wrap(~ component, scales = "free_y", ncol=3) +
  geom_point(data =  . %>% filter(object == "y"), size = 1, alpha = 0.3) +
  geom_line(data =  . %>% filter(object == "yhat"), size = 0.65) +
  theme(legend.position = "none") +
  labs(x = "Normalised Time",
       y = "Force (BW)")
  

tikz(file.path(plots_path, "interpolation-plot-fitted.tex"),
     width = 1 * doc_width_inches, 
     height = 0.35 *  doc_width_inches, standAlone = TRUE)
fitted_vs_observed_plot
dev.off()
tinytex::lualatex(file.path(plots_path,  "interpolation-plot-fitted.tex"))



# Save for Later: ---------------------------------------------------------
GRF_dataset_PRO_meta_basis_coefs <- t(GRF_dataset_PRO_meta_smooth_basis$fd$coefs)
colnames(GRF_dataset_PRO_meta_basis_coefs) <- paste0("bspline_coef", 1:35)
GRF_dataset_PRO_meta <- cbind(GRF_dataset_PRO_meta, t(GRF_dataset_PRO_meta_smooth_basis$fd$coefs))
saveRDS(object = list(GRF_dataset_PRO_meta = GRF_dataset_PRO_meta,
                      bspl_35 = bspl_35),
        here::here("chapter-06", "data", "interpolated-data.rds"))



