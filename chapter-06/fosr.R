library(data.table) # CRAN v1.14.2
library(fda)        # CRAN v5.5.1
library(tidyverse)  # CRAN v1.3.1
library(tikzDevice) # CRAN v0.12.3.1
library(refund)     # CRAN v0.1-26

source(here::here("functions", "theme_gunning.R"))
# Some settings for the Figure: -------------------------------------------
theme_gunning()
theme_update(strip.text = element_text(size = 10),
             axis.text = element_text(size = 9),
             axis.title = element_text(size = 10),
             plot.title = element_text(size = 11))

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

data_stored <- readRDS(file = "chapter-06",
                       "data",
                       "function-on-scalar-data.rds")

dt <- data_stored$data
fdobj <- data_stored$fd_obj
N <- nrow(dt) # sample size

# Pointwise Minimisation: -------------------------------------------------
dt[, CLASS_LABEL := factor(CLASS_LABEL,
                           levels = c("HC", "A", "K", "H", "C"),
                           labels = c("Healthy Control", "Ankle", "Knee", "Hip", "Calcaneous"))]
fdobj_eval <- eval.fd(0:100, fdobj)



beta_mat <- matrix(NA, nrow = 101, ncol = 5)
for(tind in seq_along(0:100)) {
  print(paste0("t = ", c(0:100)[tind]))
  df_lm <- data.frame(force_t = fdobj_eval[tind,], class_label = dt$CLASS_LABEL)
  contrasts(df_lm$class_label) <- contr.sum(5)
  test_lm <- lm(force_t ~ class_label, data = df_lm)
  beta_mat[tind, ] <- coef(test_lm)
}

# and last column is obtained by sum of other predictors (constraint)
beta_mat <- cbind(beta_mat, - apply(beta_mat[,-1], 1, sum))
colnames(beta_mat) <- c("Intercept", rownames(contrasts(df_lm$class_label)))

# can also do pointwise regressions with refund; this smoooths not interpolates though.
fosr2s <- fosr2s(Y = t(fdobj_eval),
                 X = model.matrix(~ class_label, data = df_lm),
                 argvals = 0:100,
                 nbasis = 35)

fosr2s$est.func <- cbind(fosr2s$est.func, - apply(fosr2s$est.func[,-1], 1, sum))
colnames(fosr2s$est.func) <- colnames(beta_mat) 
colnames(fosr2s$se.func) <- colnames(fosr2s$est.func)[1:5]

par(mfrow = c(2, 3))
for(j in 1:5) {
  plot(fosr2s$est.func[,j] - beta_mat[, j], type = "l") # plot difference
} # identical


# Fosr --------------------------------------------------------------------
modmat <- cbind(Intercept = 1, model.matrix(~ factor(dt$CLASS_LABEL) - 1))
(constraints = matrix(c(0, 1, 1, 1, 1, 1), 1)) # all betas sum to zero
colnames(modmat)[-1] <- stringr::str_remove(colnames(modmat)[-1], pattern =  "factor\\(dt\\$CLASS_LABEL\\)")
fosr <- fosr(fdobj = fdobj,
             X = modmat, # design matrix
             con = constraints, # sum-to-zero constraints
             method = "GLS", # (Reiss et al., 2010)  
             argvals = 0:100)

colnames(fosr$se.func) <- colnames(fosr$est.func) <- colnames(modmat)


# pffr() -------------------------------------------------------------------
dt_copy <- copy(dt)
contrasts(dt_copy$CLASS_LABEL) <- contr.sum(n = 5)
X_constraint <- model.matrix(~ CLASS_LABEL, data = dt_copy)
colnames(X_constraint) <- c("Intercept", rownames(contrasts(dt_copy$CLASS_LABEL))[1:4])
pffr_df <- as.data.frame(X_constraint)
names(pffr_df)[1:5]
# just in case we could fit smooth error:
pffr_df$id <- factor(seq_len(N))
pffr_df$Y <- t(fdobj_eval)


# pffr_fit <- pffr(formula = Y ~ class_label1 + class_label2 + class_label3 + class_label4 + s(id, bs = "re"),
#      data = pffr_df,
#      yind = 0:100,
#      algorithm = "bam",
#      bs.yindex = list(bs="ps", k=10, m=c(2, 1)),
#      bs.int  = list(bs="ps", k=35, m=c(2, 1)))

names(pffr_df)[2] <- "Healthy_Control"

pffr_fit <- pffr(formula = Y ~ Healthy_Control + Ankle + Knee + Hip,
                 data = pffr_df,
                 yind = 0:100,
                 algorithm = "bam",
                 bs.yindex = list(bs="ps", k=35, m=c(2, 1)),
                 bs.int  = list(bs="ps", k=35, m=c(2, 1)))


# fRegress: ---------------------------------------------------------------

# Three (perhaps strange) notes on this implementation: ---------------------

# 1)
# We could use the contr.sum() version of the design matrix to construct a 
# suitable design matrix (same as last fit) to enforce sum-to-zero constraint.
# however, we use the "hack" provided in Ramsay, Hooker and Graves p. 148

# 2)
# For fitting the fosr model, xfdlist is typically a list with each element
# containing the vector that would be a column of the standard design matrix.
# However, to manually use predict.fRegress(), for some reason there we need
# to create each scalar covariate as a constant "function" using a constant basis.

# 3) 
# We manually do 10-fold cross-validation to choose smoothing parameter.
# The reason being is that fRegress.CV only does leave-one-out which takes
# too long to run on this dataset


impairement_classes <- unique(dt$CLASS_LABEL)
constant_basis <- create.constant.basis(rangeval = c(0, 100))

# manually setting up xfdlist with "hack" to enforce sum to zero:
p <- length(impairement_classes) + 1
xfd_list <- vector(mode = "list", length = p)
names(xfd_list) <- c("Intercept", as.character(impairement_classes))
xfd_list_constant_coefs <- c(rep(1, times = N), 0)
xfd_list[[1]] <- fd(coef = matrix(xfd_list_constant_coefs, nrow = 1, ncol = N + 1), basisobj = constant_basis)
for (j in 2:p) {
  xj <- (dt$CLASS_LABEL == impairement_classes[j-1])
  xfd_list[[j]] <- fd(coef = matrix(data  = c(xj, 1), nrow = 1, ncol = N + 1), basisobj = constant_basis)
}
coefs_tmp <- fdobj$coefs
coefs_augmented <- cbind(coefs_tmp, matrix(data = 0, nrow = fdobj$basis$nbasis, 1))
Y_fd_augmented <- fd(coef = coefs_augmented, 
                     basisobj =  fdobj$basis)


# (RHG p. 148)
# Set up basis for predictor and regression functions:
# use same basis as data
betas_basis <- fdobj$basis

# Let's do CV to choos smoothing parameter:
log_lambda_seq <- seq(from = -10, to = 10, by = 1)
SSE_cv_vec <- vector(mode = "numeric", length = length(log_lambda_seq))

# Do shuffling of dataset prior to cross-validation
ind_obs <- seq_len(N)
set.seed(1)
ind_obs_shuffled <- sample(ind_obs)
Y_fd_augmented_shuffled <- Y_fd_augmented[c(ind_obs_shuffled, N + 1),]
x_fd_list_shuffled <- lapply(xfd_list, function(x) {x[c(ind_obs_shuffled, N + 1)]})
# Create folds:
nfolds <- 10
cv_obs <- cut(seq_len(N), labels = FALSE, breaks = nfolds)
folds <- unique(cv_obs)
stopifnot(seq_len(nfolds) == folds)

SSE_mat <- matrix(NA, nrow = nfolds, ncol = length(log_lambda_seq))

for(k in seq_along(log_lambda_seq)) {
  
  # Create parameters using same lambda value:
  cat("-------------------------------------------\n")
  print(paste0("lambda = ", 10^log_lambda_seq[k]))
  cat("-------------------------------------------\n")
  beta_fdPar_k <- fdPar(fdobj = betas_basis, Lfdobj = 2, lambda = 10^log_lambda_seq[k])
  beta_list_k <- replicate(n = p, expr = beta_fdPar_k, simplify = FALSE)
  
  print("Doing Cross-Validation:")
  for(fold_ind in folds) {
    
    print(paste("Fold", fold_ind, "of", nfolds))
    
    # Do test:train split.
    train_inds <- which(cv_obs != fold_ind)
    test_inds <- which(cv_obs == fold_ind)
    
    # Do training:
    Y_fd_augmented_train <- Y_fd_augmented_shuffled[c(train_inds, N+1),]
    x_fd_list_train <- lapply(x_fd_list_shuffled, function(x) {x[c(train_inds, N+1)]})
    fRegress_fit <- fRegress(Y_fd_augmented_train,
                             xfdlist = x_fd_list_train,
                             betalist = beta_list_k)
    # Do testing
    Y_fd_test <- Y_fd_augmented_shuffled[test_inds,]
    x_fd_list_test <- lapply(x_fd_list_shuffled, function(x) {x[c(test_inds)]})
    Y_fd_test_hat <- predict.fRegress(object = fRegress_fit, newdata = x_fd_list_test)
    err_fd_test <- Y_fd_test - Y_fd_test_hat
    
    # Store Cross-Validated Errors
    SSE_mat[fold_ind, k] <- sum(inprod(err_fd_test^2))
  }
}

SSE_vec <- apply(SSE_mat, 2, sum) # sum SSE over folds
# plot log(lambda) vs. SSE
par(mfrow = c(1, 2))
plot(log_lambda_seq, log(SSE_vec), type = "b", pch = 20)
# zoom in:
plot(log_lambda_seq[1:13], log(SSE_vec)[1:13], type = "b", pch = 20)

# take log(lambda) that minimises SSE CV"
(log_lambda_best <- log_lambda_seq[which.min(SSE_vec)])

# Do a model fit with best lambda:
beta_fdPar_best <- fdPar(fdobj = betas_basis, Lfdobj = 2, lambda = 10^log_lambda_best)
beta_list_best <- replicate(n = p, expr = beta_fdPar_best, simplify = FALSE)
names(beta_list_best) <- names(xfd_list)
fRegress_best_fit <- fRegress(Y_fd_augmented,
                              xfdlist = xfd_list,
                              betalist = beta_list_best)


# And construct confidence intervals
y2cMap <- smooth.basis(argvals = 0:100,
                       y = eval.fd(0:100, fdobj[1]),
                       fdParobj = fdPar(fdobj$basis, Lfdobj = 2, lambda = 10^-10))$y2cMap # didn't have smooth.basis stored.

# Calculate residual error covariance for getting CIs
err_fd <- fRegress_best_fit$yhatfdobj[1:N, ] - fRegress_best_fit$yfdobj[1:N, ]
err_fd_cov_bifd <- var.fd(fdobj1 = err_fd, fdobj2 = err_fd)
err_fd_cov_eval <- eval.bifd(0:100, 0:100, err_fd_cov_bifd)

par(mfrow = c(1, 2))
plot(fosr$resid, type = "l")
plot(err_fd)


# Get Standard Errors for fRegress()
# can skip this for now: takes a while!
fRegress_best_fit_stderrList <- fRegress.stderr(y = fRegress_best_fit,
                                                y2cMap = y2cMap,
                                                SigmaE = err_fd_cov_eval)



# This approach DOES NOT WORK for getting CIs: ----------------------------
# do it by constructing design matrix:
p2 <- length(impairement_classes)
xfd_list_2 <- vector(mode = "list", length = p2)
names(xfd_list_2) <- colnames(X_constraint)

for (j in 1:p2) {
  xfd_list_2[[j]] <- fd(coef = matrix(data  = X_constraint[, j], nrow = 1, ncol = N), basisobj = constant_basis)
}

beta_fdPar_best <- fdPar(fdobj = betas_basis, Lfdobj = 2, lambda = 10^log_lambda_best)
beta_list_best_2 <- replicate(n = p2, expr = beta_fdPar_best, simplify = FALSE)
names(beta_list_best_2) <- names(xfd_list_2)


fRegress_best_fit_2 <- fRegress(fdobj,
                                xfdlist = xfd_list_2,
                                betalist = beta_list_best_2)


for(j in c("Intercept", "Healthy Control", "Ankle", "Knee", "Hip")) {
  plot(fRegress_best_fit_2$betaestlist[[j]]$fd, col = 2)
  lines(fRegress_best_fit$betaestlist[[j]]$fd, col = 2)
} # close but still different

# Calculate residual error covariance for getting CIs
err_fd_2 <- fRegress_best_fit_2$yhatfdobj - fRegress_best_fit_2$yfdobj
err_fd_cov_bifd_2 <- var.fd(fdobj1 = err_fd_2, fdobj2 = err_fd_2)
err_fd_cov_eval_2 <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, err_fd_cov_bifd_2)
# Get Standard Errors for fRegress()
# can skip this for now: takes a while!
fRegress_best_fit_stderrList_2 <- fRegress.stderr(y = fRegress_best_fit_2,
                                                  y2cMap = y2cMap,
                                                  SigmaE = err_fd_cov_eval_2)
names(fRegress_best_fit_stderrList_2$betastderrlist) <- names(fRegress_best_fit_2$betaestlist)



# Comparison: -------------------------------------------------------------


pffr_coefs <- sapply(coef(pffr_fit)[["smterms"]], function(x) x[["value"]])
colnames(pffr_coefs) <- stringr::str_remove(colnames(pffr_coefs), "\\(yindex\\)")
colnames(pffr_coefs)[colnames(pffr_coefs) == "Healthy_Control"] <- "Healthy Control"
pffr_coefs <- cbind(pffr_coefs, Calcaneous = - apply(pffr_coefs[,-1], 1, sum))
pffr_yind <- coef(pffr_fit)[["smterms"]][[1]][["coef"]][["yindex.vec"]]


par(mfrow = c(2, 3))
for(j in c("Intercept", "Healthy Control", "Ankle", "Knee", "Hip", "Calcaneous")) {
  plot(0:100, beta_mat[,j], type = "l", xlab = expression(t), ylab = expression(beta(t)))
  title(j)
  lines(fRegress_best_fit$betaestlist[[j]]$fd, col = 2)
  lines(fRegress_best_fit_2$betaestlist[[j]]$fd, col = 3)
  lines(0:100, fosr$est.func[, j], col = 4)
  lines(pffr_yind, pffr_coefs[,j], col = 5)
}
# all reasonably similar!

# Let's extract point wise confidence intervals
j <- "Hip"

pffr_se_j <- coef(pffr_fit)[["smterms"]][[paste0(j, "(yindex)")]][["se"]]
plot_hip_dt <- data.table(t = c(0:100, 0:100, 0:100, pffr_yind),
           model = c(rep("P-GLS (\\texttt{fosr()})", 101),
                     rep("Two-Step (\\texttt{fosr2s()})", 101),
                     rep("P-OLS (\\texttt{fRegress()})", 101),
                     rep("FAMM (\\texttt{pffr()})", length(pffr_yind))
                     ),
           point_est = c(fosr$est.func[, j],
                         fosr2s$est.func[, j],
                         eval.fd(0:100, fRegress_best_fit_2$betaestlist[[j]]$fd),
                         pffr_coefs[,j]),
           lower = c(fosr$est.func[, j] - 2 * fosr$se.func[, j],
                     fosr2s$est.func[, j] - 2 * fosr2s$se.func[, j],
                     eval.fd(0:100, fRegress_best_fit_2$betaestlist[[j]]$fd) - 2 * eval.fd(0:100, fRegress_best_fit_stderrList_2$betastderrlist[[j]]),
                     pffr_coefs[, j] - 2 * pffr_se_j),
           upper = c(fosr$est.func[, j] + 2 * fosr$se.func[, j],
                     fosr2s$est.func[, j] + 2 * fosr2s$se.func[, j],
                     eval.fd(0:100, fRegress_best_fit_2$betaestlist[[j]]$fd) + 2 * eval.fd(0:100, fRegress_best_fit_stderrList_2$betastderrlist[[j]]),
                     pffr_coefs[, j] + 2 * pffr_se_j))


compare_plot <- ggplot(data = plot_hip_dt) +
  aes(x = t, y = point_est) +
  facet_wrap(~ model, nrow = 2) +
  geom_line() +
  geom_line(aes(y= lower), linetype = 3, size = 0.65) +
  geom_line(aes(y= upper), linetype = 3, size = 0.65) +
  labs(y = "$\\beta_4 (t)$ Impaired: Hip",
       x = "Normalised Time ($\\%$ of Stance)")
compare_plot

tikz(here::here("chapter-06", "figures", "fosr-coefs-different-methods.tex"),
     width = 1 * doc_width_inches, standAlone = TRUE,
     height = 0.9 *  doc_width_inches)
compare_plot
dev.off()


tinytex::lualatex(here::here("chapter-06", "figures", "fosr-coefs-different-methods.tex"))

# -------------------------------------------------------------------------




# Do bootstrap: -----------------------------------------------------------

# choose to use P-GLS for this:
# Fosr --------------------------------------------------------------------
inds_for_bootstrap <- seq_len(nrow(modmat))
B <- 1000 # number of bootstrap samples.
bootstrap_estimates <- array(data = NA, dim = c(101, 6, B))
set.seed(96)
for(b in seq_len(B)) {
  print(paste0("Bootstrap Iteration ", b, " of ", B))
  # select bootstrap sample of observations:
  bootstrap_inds <- sample(inds_for_bootstrap, replace = TRUE)
  # fit model on these observations:
  fosr_b <- fosr(fdobj = fdobj[bootstrap_inds, ], # functional response sub setted for bootstrap obs
                 X = modmat[bootstrap_inds, ], # design matrix sub setted for bootstrap obs
                 con = constraints, # sum-to-zero constraints
                 method = "GLS", # (Reiss et al., 2010)  
                 argvals = 0:100)
  bootstrap_estimates[,,b] <- fosr_b$est.func
}


se_mat <- apply(bootstrap_estimates, MARGIN = c(1, 2), FUN = sd)

plot(fosr$est.func)

par(mfrow = c(2, 3))
for(j in 1:6) {
  ylim <- range(
    fosr$est.func[,j] + 2 * se_mat,
    fosr$est.func[,j] - 2 * se_mat,
    fosr$est.func[,j] + 2 * fosr$se.func,
    fosr$est.func[,j] + 2 * fosr$se.func)
  plot(0:100, fosr$est.func[, j], type = "l", ylim = ylim)
  lines(0:100, fosr$est.func[,j] + 2 * se_mat[,j], col = 2, lty = 2)
  lines(0:100, fosr$est.func[,j] - 2 * se_mat[,j], col = 2, lty = 2)
  lines(0:100, fosr$est.func[,j] + 2 * fosr$se.func[,j], col = 3, lty = 2)
  lines(0:100, fosr$est.func[,j] - 2 * fosr$se.func[,j], col = 3, lty = 2)
}
# bootstrap and analytic intervals v similar!

# Do the Crainiceanu et al. (2012) /  Ruppert et al.  (2003) Bands: -------

q_vector <- vector(mode = "numeric", length = 6)
n_samples <- 100000 # no. of samples from mvt normal.
for(j in seq_len(6)) {
  print(j)
  Beta_hat <- fosr$est.func[,j]
  se_boot_hat <- se_mat[,j]
  Cov_Beta_s_Beta_t <- cov(t(bootstrap_estimates[,j,]))
  Beta_samples <- t(mvrnorm(n = n_samples, mu = Beta_hat, Sigma = Cov_Beta_s_Beta_t))
  Beta_samples_centered <- sweep(Beta_samples, MARGIN = c(1), STATS = Beta_hat, FUN = "-")
  t_stat_samples <- sweep(Beta_samples_centered, MARGIN = c(1), STATS = se_boot_hat, FUN = "/")
  t_stat_max_abs_samples <- apply(t_stat_samples, 2, function(x) max(abs(x)))
  q_0.95 <- quantile(t_stat_max_abs_samples, probs = 0.95)
  q_vector[j] <- q_0.95
}

lower_pw <- upper_pw <- lower_sim <- upper_sim <- matrix(NA, nrow = 101, ncol = 6)
for(j in seq_len(j)) {
  print(j)
  Beta_hat <- fosr$est.func[,j]
  se_boot_hat <- se_mat[,j]
  lower_pw[, j] <- Beta_hat - 2 * se_boot_hat
  upper_pw[, j] <- Beta_hat + 2 * se_boot_hat
  lower_sim[, j] <- Beta_hat - q_vector[j] * se_boot_hat
  upper_sim[, j] <- Beta_hat + q_vector[j] * se_boot_hat
}




plot_dt <- data.table(t = rep(0:100),
                      parameter = rep(c("Intercept", "Healthy Control", "Ankle", "Knee", "Hip", "Calcaneous"), each = 101),
                      beta_hat = c(fosr$est.func),
                      lower_pw = c(lower_pw),
                      upper_pw = c(upper_pw),
                      lower_sim = c(lower_sim),
                      upper_sim = c(upper_sim))


plot_dt[, parameter := factor(parameter,
                              levels =  c("Intercept", "Healthy Control", "Ankle", "Knee", "Hip", "Calcaneous"),
                              labels = c("$\\beta_0 (t)$ Overall Mean",
                                         "$\\beta_1 (t)$ Healthy Control",
                                         "$\\beta_2 (t)$ Impaired: Ankle",
                                         "$\\beta_3 (t)$ Impaired: Knee",
                                         "$\\beta_4 (t)$ Impaired: Hip",
                                         "$\\beta_5 (t)$ Impaired: Calcaneus"))]

coef_plot <- ggplot(data = plot_dt) +
  facet_wrap(~ parameter, scales = "free_y") +
  aes(x = t, colour = parameter, fill = parameter) +
  geom_hline(yintercept = 0, colour = 'grey') +
  geom_line(aes(y = beta_hat)) +
  geom_point(data = . %>% filter(parameter != "Intercept"), # TRICK :-)
             inherit.aes = F,
             aes(x = 50, y = -0.06), color = NA) +
  geom_point(data = . %>% filter(parameter != "Intercept"),
             inherit.aes = F, aes(x = 50, y = 0.06), color = NA) +
  geom_line(aes(y = lower_pw), linetype = "dotted") +
  geom_line(aes(y = upper_pw), linetype = "dotted") +
  geom_ribbon(mapping = aes(ymin = lower_sim, ymax = upper_sim), alpha = 0.25, col = NA) +
  labs(x = "Normalised Time ($\\%$ of Stance)",
       y = "$\\widehat{\\beta} (t)$") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", scales::hue_pal()(5))) +
  scale_fill_manual(values = c("black", scales::hue_pal()(5)))

coef_plot

tikz(here::here("chapter-06", "figures", "fosr-coefs-pw-sim.tex"),
     width = 1 * doc_width_inches, standAlone = TRUE,
     height = (2.1/3) *  doc_width_inches)
coef_plot
dev.off()


tinytex::lualatex(here::here("chapter-06", "figures", "fosr-coefs-pw-sim.tex"))


# -------------------------------------------------------------------------

set.seed(1)
Fperm_fd <- Fperm.fd(yfdPar = fdobj,
                     xfdlist = xfd_list_2,
                     betalist = beta_list_best_2, 
                     nperm = 400, 
                     argvals = 0:100)
fosr_perm <- fosr.perm(fdobj = fdobj,
                       X = modmat, # design matrix
                       con = constraints, # sum-to-zero constraints
                       method = "GLS", # (Reiss et al., 2010)  
                       argvals = 0:100, 
                       nperm = 400)

tikz(here::here("chapter-06", "figures", "fosr-perm-tests.tex"),
     width = 1 * doc_width_inches, standAlone = TRUE,
     height = 0.5 *  doc_width_inches)
par(mfrow = c(1, 2))
matplot(Fperm_fd$Fnullvals,
        col = "grey",
        type = "l",
        lty = 1,
        ylab = "F statistics",
        xlab = "Normalised Time ($\\%$ of Stance)",
        ylim = c(-0.0005, 0.45),
        xlim = c(0, 100))
title("\\texttt{Fperm.fd()}")
lines(Fperm_fd$Fvals, type = "l", col = "blue")
abline(h = Fperm_fd$qval, col = "red", lty = 2)
plot(fosr_perm, 
     xlabel = "Normalised Time ($\\%$ of Stance)",
     ylim = c(0, 70))
title("\\texttt{fosr.perm()}")
dev.off()

tinytex::lualatex(here::here("chapter-06", "figures", "fosr-perm-tests.tex"))
