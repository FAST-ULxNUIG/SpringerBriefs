# Packages: ---------------------------------------------------------------
library(readr) # CRAN v1.3.1
library(fda) # CRAN v5.1.9
library(tidyverse) # CRAN v1.3.0

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


# Some Helper Functions: --------------------------------------------------
source(here::here("functions/center_fd_around_new_mean.R"))
source(here::here("functions/project_data_onto_fpcs.R"))

# Data Import: ------------------------------------------------------------

# From:
# https://springernature.figshare.com/collections/GaitRec_A_large-scale_ground_reaction_force_dataset_of_healthy_and_impaired_gait/4788012/1

# Citation:
############
# Horsak, Brian; Slijepcevic, Djordje; Raberger, Anna-Maria; Schwab, Caterine; Worisch, Marianne; Zeppelzauer, Matthias (2020).
# GaitRec: A large-scale ground reaction force dataset of healthy and impaired gait. figshare. Collection. 
# https://doi.org/10.6084/m9.figshare.c.4788012.v1
############

# Download the following Ground Reaction Forces Data:
# Force (F)
# vertical, anterior-posterior and medio-lateral force components (V, AP, ML)
# Legs: Both Left (left) and Right (right)
# Processing: Smoothed and Normalized (PRO)
# This is a big data set, it might take a while:

# Vertical Force:
GRF_F_V_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063191")
GRF_F_V_PRO_right <- read_csv("https://springernature.figshare.com/ndownloader/files/22063119")

# Medio-lateral Force:
GRF_F_ML_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063113")
GRF_F_ML_PRO_right <- read_csv("https://springernature.figshare.com/ndownloader/files/22063086")

# Anterior-Posterior Force:
GRF_F_AP_PRO_left <- read_csv("https://springernature.figshare.com/ndownloader/files/22063185")
GRF_F_AP_PRO_right <- read_csv("https://springernature.figshare.com/ndownloader/files/22063101")

# And associated MetaDeta, e.g., session and subject information:
GRF_metadata <- read_csv("https://ndownloader.figshare.com/files/22062960")
#.. you might have to wait.. it Will download!




# Data Wrangling: ---------------------------------------------------------
# Do this manually, possibly could do it not manually.
GRF_F_V_PRO_left$side <- "left"
GRF_F_V_PRO_right$side <- "right"
GRF_F_ML_PRO_left$side <- "left"
GRF_F_ML_PRO_right$side <- "right"
GRF_F_AP_PRO_left$side <- "left"
GRF_F_AP_PRO_right$side <- "right"

GRF_F_V_PRO_left$component <- "vertical"
GRF_F_V_PRO_right$component <- "vertical"
GRF_F_ML_PRO_left$component <- "medio_lateral"
GRF_F_ML_PRO_right$component <- "medio_lateral"
GRF_F_AP_PRO_left$component <- "anterior_posterior"
GRF_F_AP_PRO_right$component <- "anterior_posterior"

# renaming the grf columns:
GRF_F_V_PRO_left <- rename_at(.tbl = GRF_F_V_PRO_left,
                              .vars = vars(paste0("F_V_PRO_", 1:101)), 
                              .funs = ~ paste0("time_", 0:100))

GRF_F_V_PRO_right <- rename_at(.tbl = GRF_F_V_PRO_right,
                               .vars = vars(paste0("F_V_PRO_", 1:101)), 
                               .funs = ~ paste0("time_", 0:100))

GRF_F_ML_PRO_left <- rename_at(.tbl = GRF_F_ML_PRO_left,
                               .vars = vars(paste0("F_ML_PRO_", 1:101)), 
                               .funs = ~ paste0("time_", 0:100))

GRF_F_ML_PRO_right <- rename_at(.tbl = GRF_F_ML_PRO_right,
                                .vars = vars(paste0("F_ML_PRO_", 1:101)), 
                                .funs = ~ paste0("time_", 0:100))

GRF_F_AP_PRO_left <- rename_at(.tbl = GRF_F_AP_PRO_left,
                               .vars = vars(paste0("F_AP_PRO_", 1:101)), 
                               .funs = ~ paste0("time_", 0:100))

GRF_F_AP_PRO_right <- rename_at(.tbl = GRF_F_AP_PRO_right,
                                .vars = vars(paste0("F_AP_PRO_", 1:101)), 
                                .funs = ~ paste0("time_", 0:100))


stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_V_PRO_right)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_ML_PRO_left)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_ML_PRO_right)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_AP_PRO_left)))
stopifnot(all.equal(names(GRF_F_V_PRO_left), names(GRF_F_AP_PRO_right)))

# create full dataset for analysis:

GRF_dataset_PRO <- rbind(
  GRF_F_V_PRO_left,
  GRF_F_V_PRO_right,
  GRF_F_AP_PRO_left,
  GRF_F_AP_PRO_right,
  GRF_F_ML_PRO_left,
  GRF_F_ML_PRO_right
)


inner_join(x = GRF_dataset_PRO, y = GRF_metadata)

GRF_dataset_PRO_meta <- merge.data.frame(
  x = GRF_dataset_PRO, 
  y = GRF_metadata,
  by = c("SUBJECT_ID", "SESSION_ID"))

stopifnot(nrow(GRF_dataset_PRO) == nrow(GRF_dataset_PRO_meta))

GRF_dataset_PRO_meta <- relocate(GRF_dataset_PRO_meta,
                                 time_0:time_100,
                                 .after = last_col())

# For now work only with balanced, controlled subset of the data:
GRF_dataset_PRO_meta_test <- GRF_dataset_PRO_meta %>%
  filter(SHOD_CONDITION == 1 & SPEED == 2 & TEST == 1 & CLASS_LABEL == "HC" & side == "right")

y <- t(GRF_dataset_PRO_meta_test[, paste0("time_", 0:100)])
bspl_35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 35, norder = 4)
GRF_dataset_PRO_meta_test_smooth_basis <- smooth.basis(argvals = 0:100,
                                                  y = y,
                                                  fdParobj = bspl_35)

GRF_dataset_PRO_meta_test_basis_coefs <- t(GRF_dataset_PRO_meta_test_smooth_basis$fd$coefs)
GRF_dataset_PRO_meta_test <- data.table(cbind(GRF_dataset_PRO_meta_test, t(GRF_dataset_PRO_meta_test_smooth_basis$fd$coefs)))


GRF_dataset_PRO_meta_test[, paste0("time_",0:100) := NULL]
GRF_dataset_PRO_averages_test <- GRF_dataset_PRO_meta_test[,
                                                 as.list(apply(.SD, 2, mean)), # average basis coefficients of all trials
                                                 by = .(SUBJECT_ID,component),
                                                 .SDcols = paste0("bspl4.",1:35)] # says which columns to average

GRF_dataset_PRO_averages_vertical_test <- GRF_dataset_PRO_averages_test[component == "vertical"]
# create fd object defined by coefficients and basis object
GRF_fdobj_PRO_averages_vertical_test <- fd(coef = t(as.matrix(GRF_dataset_PRO_averages_vertical_test[, paste0("bspl4.",1:35)])),
                                      basisobj = bspl_35)

GRF_dataset_PRO_averages_anterior_posterior_test <- GRF_dataset_PRO_averages_test[component == "anterior_posterior"]
# create fd object defined by coefficients and basis object
GRF_fdobj_PRO_averages_anterior_posterior_test <- fd(coef = t(as.matrix(GRF_dataset_PRO_averages_anterior_posterior_test[, paste0("bspl4.",1:35)])),
                                                basisobj = bspl_35)

anterior_posterior_fd_eval_test <- eval.fd(evalarg = 0:100, fdobj = GRF_fdobj_PRO_averages_anterior_posterior_test)
max_anterior_posterior_test <- apply(anterior_posterior_fd_eval_test, 2, max)



# Extract Model Results: --------------------------------------------------

model_results <- readRDS(file = here::here("chapter-06", "data", "model_results.rds"))
GRF_fpca_averages_vertical <- model_results$fpca_vertical
fpcr_best <- model_results$fpcr
fRegress_best_fit <- model_results$fregress  
pfr <- model_results$pfr  



# FPCR --------------------------------------------------------------------


GRF_fdobj_PRO_averages_vertical_test_cent <- center_fd_around_new_mean(fdobj = GRF_fdobj_PRO_averages_vertical_test, GRF_fpca_averages_vertical$meanfd)
GRF_fpc_scores_test <- project_data_onto_fpcs(fdobj = GRF_fdobj_PRO_averages_vertical_test_cent,
                                              pca.fd_obj = GRF_fpca_averages_vertical)

test_df_fpcr <- data.frame(max_anterior_posterior = max_anterior_posterior_test,
                           GRF_fpc_scores_test)
names(test_df_fpcr)[-1] <- paste0("vertical_fpca_scores", 1:35)

plot(predict(object = fpcr_best, newdata = test_df_fpcr),
     max_anterior_posterior_test)

fpcr_test_yhat <- predict(object = fpcr_best, newdata = test_df_fpcr)
plot(fpcr_test_yhat, fregress_test_yhat)





# fRegress() --------------------------------------------------------------
# Need to create constant fd object for scalar predictor to work in fRegress.CV(), weird, I know.
constant_basis <- create.constant.basis(c(0,100))
constant_fd <- fd(coef = matrix(1, nrow = 1, ncol = ncol(GRF_fdobj_PRO_averages_vertical_test$coefs)), basisobj = constant_basis)
# List containing predictors (jntercept and functional covariate)
xfd_list_test <- list(constant_fd, # for intercept
                 GRF_fdobj_PRO_averages_vertical_test) 
fregress_test_yhat <- predict.fRegress(object = fRegress_best_fit,
                 newdata = xfd_list_test
                 )

cor(fregress_test_yhat, max_anterior_posterior_test)

boxplot(abs(fpcr_test_yhat - max_anterior_posterior_test),
        abs(fregress_test_yhat - max_anterior_posterior_test))


# pfr ---------------------------------------------------------------------

vertical_fd_eval_test <- t(eval.fd(evalarg = 0:100, fdobj = GRF_fdobj_PRO_averages_vertical_test))
pfr <- pfr(max_anterior_posterior ~ lf(X = vertical_fd_eval, bs = "bs", k = 35, argvals = 0:100))
pfr_test_yhat <- predict(object = pfr, newdata = list(vertical_fd_eval = vertical_fd_eval_test))


abs_errors_dt <- data.table(
  ind = seq_len(length(max_anterior_posterior_test)),
  fpcr = abs(fpcr_test_yhat - max_anterior_posterior_test),
  fRegress =  abs(c(fregress_test_yhat) - max_anterior_posterior_test),
  pfr = abs(pfr_test_yhat - max_anterior_posterior_test))

abs_errors_dt_lng <- melt.data.table(data = abs_errors_dt, 
                                     id.vars = "ind")

abs_errors_dt_lng[, method := fcase(
  variable == "pfr", "\\texttt{pfr()}",
  variable == "fRegress", "\\texttt{fRegress()}",
  variable == "fpcr", "FPCR"
)]

error_plot <- ggplot(abs_errors_dt_lng) +
  aes(x = method, y = value, group = method, fill = method) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Method",
       y = "Absolute Error", 
       title = "Out-of-Sample Prediction Errors") +
  theme(legend.position = "none")

tikz(file = file.path(plots_path, "sofr-plot-error.tex"),
     width = 0.6 * doc_width_inches, 
     height = 0.5 * doc_width_inches, 
     standAlone = TRUE)
error_plot
dev.off() 
tinytex::lualatex(file.path(plots_path, "sofr-plot-error.tex"))
