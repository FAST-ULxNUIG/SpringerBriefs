# =======================================================================#
# Make a Figure that demonstrates landmark registration being applied
# in a tep by step fashion

# Use the childrens' gait data set and add some phase variation so
# the effect of registration is clear
# ========================================================================#

# 1) Load Packages ------------------------------------------------------
library(fda)        # CRAN v5.1.9
library(tidyverse)  # CRAN v1.3.0
library(ggpubr)     # CRAN v0.4.0
library(registr)    # CRAN v1.0.0
library(Polychrome) # CRAN v1.2.6
library(tikzDevice) # CRAN v0.12.3.1 

# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
plots_path <- here::here("chapter-02", "figures")
results_path <- here::here("chapter-02", "results")
source(file.path(functions_path, "theme_gunning.R"))
source(file.path(functions_path, "landmarkreg-fixed.R"))
theme_gunning()
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12),
             axis.ticks = element_blank(),
             plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
             plot.subtitle = element_text(size = 12, hjust = 0.5))


# 3) Get gait data -----------------------------------------------------------

data("gait")
gait_time <- seq(0, 1, length.out = 20)
smth_gait <- smooth.basis(argvals = gait_time,
                          y = gait, 
                          fdParobj = fdPar(create.bspline.basis(rangeval = range(gait_time), nbasis = 20, norder = 4),
                                           lambda = 10^-5))
# we could possibly use order = 5 and choose a penalty that penalises 3rd derivative
# so that the 1st derivative we use is mooth. But for the purpose of the demonstration this suffices

eval_100 <- seq(0, 1, length.out = 100)
eval_gait <- eval.fd(evalarg = eval_100, fdobj = smth_gait$fd)


# 4) Take a random sample of 12 curves ---------------------------------------
set.seed(1995)
ss <- 12
sample_12 <- sample(1:39, replace = F, size = ss)
matplot(eval_gait[, sample_12, 2], type  = "l")


# 5) Simulate phase variation ------------------------------------------------
set.seed(1996)
phase_var <- runif(n = ss, min = 0.6, max = 1.4)
eval_100_mat <- matrix(eval_100, nrow = nrow(eval_gait), ncol = ncol(eval_gait))
phase_evals <-   matrix(data = NA, nrow = nrow(eval_100_mat), ncol = ncol(eval_100_mat))
for(i in 1:ss){
  phase_evals[, i] <- eval_100_mat[, i] ^ phase_var[i]
}


sim_phase.fd <- smooth.basis(argvals = phase_evals, y = eval_gait[,sample_12,],
                             fdParobj = fdPar(create.bspline.basis(rangeval = range(gait_time), nbasis = 20, norder = 4),
                                              lambda = 10^-10))

phase_varying_eval <- eval.fd(sim_phase.fd$fd, evalarg = eval_100)
phase_varying_deriv_eval <- eval.fd(sim_phase.fd$fd, evalarg = eval_100, Lfdobj = 1)
matplot(phase_varying_deriv_eval[,, 2], type = "l")
matplot(phase_varying_eval[,, 2], type = "l")


# 6) Extract the landmarks ---------------------------------------------------
# we're going to point and click on the zero crossing of each derivative:
# DO NOT RUN THIS, I have already point and clicked. We can just load the landamrks
# Don'run unless you want to pick out the zero crossings again:

sim_phase.fd$fd$fdnames <- smth_gait$fd$fdnames
# veloctiyfun <- deriv.fd(sim_phase.fd$fd[, 2],1)
# landmarks <- matrix(0, nrow = ss, ncol=1)
# par(mfrow=c(1,1), ask=T)
# icase=1
# for(icase in 1:ss){
#   velveci <- predict(veloctiyfun[icase], eval_100)
#   plot(eval_100, velveci, "l", main=paste("case", icase), xlab="time", ylab="angle")
#   abline(h = 0)
#   my_loc = locator(n=1)
#   landmarks[icase]=my_loc$x[1]
# }
# 
# # can also get their corresponding values while we're at it (in both the function and the derivative)
# deriv_vals <- sapply(X = 1:ss, function(y){
#   predict(veloctiyfun[y], landmarks[y])
# })
# 
# fun_vals <- sapply(X = 1:ss, function(y){
#   predict(sim_phase.fd$fd[, 2][y], landmarks[y])
# })
# 
# landmarks_df <- data.frame(subject = paste0("sub", 1:ss), time = landmarks, deriv_vals, fun_vals)
# 
# saveRDS(object = landmarks_df, file = file.path(results_path, "landmarks_df.rds")) # so we don't have to re click each time
# saveRDS(object = landmarks, file = file.path(results_path, "landmarks.rds"))


# can skip to here:
landmarks_df <- readRDS(file = file.path(results_path, "landmarks_df.rds"))
landmarks <- readRDS(file = file.path(results_path, "landmarks.rds"))

# 7) Perform landmark registration -------------------------------------------
my_landmarks_mean = mean(landmarks)

#### set up basis for the warping functions
wbasisLM <- create.bspline.basis(rangeval = range(eval_100), nbasis=4, norder = 3,
                                 breaks = c(0, my_landmarks_mean, 1))

WfdLM = fd(matrix(0,4,1), wbasisLM) # are these initial guesses
WfdParLM <- fdPar(WfdLM, 1, 1e-10)

regListLM <- landmarkreg(fdobj = sim_phase.fd$fd[, 2],
                         ximarks = landmarks,
                         x0marks = my_landmarks_mean,
                         WfdPar = WfdParLM,
                         monwrd = TRUE)

registered_eval <- eval.fd(evalarg = eval_100, fdobj = regListLM$regfd)
warping_eval <- eval.fd(evalarg = eval_100, fdobj = regListLM$warpfd)
fun_eval <- eval.fd(evalarg = eval_100, fdobj = sim_phase.fd$fd[, 2])
deriv_eval <-   eval.fd(evalarg = eval_100, fdobj = veloctiyfun)

colnames(registered_eval) <- colnames(warping_eval) <- colnames(fun_eval) <- colnames(deriv_eval) <- paste0("sub", 1:ss)

# 8) Create plots ------------------------------------------------------------
# we make 4 different plots and combine them:
par(ask = F)
set.seed(2021)
my_pallette_12 <- createPalette(N = 12,
                                seedcolors = c("#ff0000", "#00ff00")) %>% as.vector()


r1 <- data.frame(time = eval_100, fun_eval) %>%
  gather(-time, key = "subject", value = "x_t") %>%
  ggplot() +
  aes(x = time, y = x_t, group = subject, color = subject) +
  geom_line(alpha = 0.8) +
  geom_point(data = landmarks_df, inherit.aes = F,
             mapping = aes(x = time, y = fun_vals, color = subject),
             size = 2.25) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1),labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_color_manual(values = my_pallette_12) +
  theme(legend.position = "none") +
  labs(title = "Unregistered Data $x(t)$",
       x = "Normalised Time",
       y = "Angle ($^{\\circ}$)") +
  annotate(geom = "text", x = 0.25, y = 70, size = (5/ 14) * 9,
           label = "Peak knee flexion angle \n varies from curve to curve") +
  geom_curve(mapping = aes(x = 0.5, xend = 0.6, y= 70, yend = 75),
             curvature = - 0.4, color = "black",
             arrow = arrow(length = unit(0.1, "inches")))



r2 <- data.frame(time = eval_100, deriv_eval) %>%
  gather(-time, key = "subject", value = "x_t") %>%
  ggplot() +
  aes(x = time, y = x_t, group = subject, color = subject) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.7) +
  geom_point(data = landmarks_df, inherit.aes = F,
             mapping = aes(x = time, y = deriv_vals, color = subject),
             size = 2.25, alpha = 0.7) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1),labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_color_manual(values = my_pallette_12) +
  theme(legend.position = "none") +
  labs(title = "Unregistered First Derivative $x'(t)$",
       x = "Normalised Time",
       y = "$x'(t)$") +
  annotate(geom = "text", x = 0.26, y = 400, size = (5/ 14) * 9,
           label = "The peak timings are picked \n out manually from the zero \n crossings of the derivative") +
  geom_curve(mapping = aes(x = 0.3, xend = 0.625, y = 250, yend = 0),
             curvature = 0.4, color = "black",
             arrow = arrow(length = unit(0.15, "inches"))) +
  annotate(geom = "text", x = 0.4, y = -250, size = (5/ 14) * 9,
           label = "A target landmark time is \n chosen to register the curves to") +
  geom_curve(mapping = aes(x = 0.63, xend = mean(landmarks), y = - 250, yend = -50),
             curvature = 0.2, color = "black",
             arrow = arrow(length = unit(0.15, "inches"))) +
  geom_vline(xintercept = mean(landmarks)) +
  geom_point(mapping = aes(x= mean(landmarks), y = 0), inherit.aes = F, pch = 13, size = 7)



r3 <- data.frame(time = eval_100, registered_eval) %>%
  gather(-time, key = "subject", value = "x_t") %>%
  ggplot() +
  aes(x = time, y = x_t, group = subject, color = subject) +
  geom_vline(xintercept = mean(landmarks)) +
  geom_line(alpha = 0.8) +
  geom_point(data = landmarks_df, inherit.aes = F,
             mapping = aes(x = mean(landmarks), y = fun_vals, color = subject),
             size = 2.25) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1),labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_color_manual(values = my_pallette_12) +
  theme(legend.position = "none") +
  labs(title = "Registered Data $x \\circ h (t)$",
       x = "Normalised Time",
       y = "Angle ($^{\\circ}$)") +
  annotate(geom = "text", x = 0.34, y = 60, size = (5/ 14) * 9,
           label = "The curves are now aligned \n at peak knee flexion angle") +
  geom_curve(mapping = aes(x = 0.5, xend = 0.64, y = 65, yend = 75),
             curvature = -0.4, color = "black",
             arrow = arrow(length = unit(0.1, "inches")))


r4 <- data.frame(time = eval_100, warping_eval) %>%
  gather(-time, key = "subject", value = "x_t") %>%
  ggplot() +
  aes(x = time, y = x_t, group = subject, color = subject) +
  geom_vline(xintercept = mean(landmarks)) +
  geom_line(alpha = 0.8) +
  geom_point(data = landmarks_df, inherit.aes = F,
             mapping = aes(y = time, x = rep(mean(landmarks), length(landmarks)), color = subject),
             size = 2.25) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_color_manual(values = my_pallette_12) +
  theme(legend.position = "none") +
  labs(title = "Warping functions $h(t)$",
       x = "Registered Time",
       y = "Observed Time") +
  annotate(geom = "text", x = 0.325, y = 0.85, size = (5/ 14) * 9,
           label = "Construct warping functions \n that map the target time \n to the peak timings")+
  geom_curve(mapping = aes(x = 0.55, xend = 0.7, y = 0.875, yend = 0.825),
             curvature = -0.2, color = "black",
             arrow = arrow(length = unit(0.1, "inches")))

comb_r <- ggarrange(r1, r2, r4, r3, ncol = 2, nrow = 2, align = "hv")
comb_r

# Save Plots 
# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "landmark-reg-plot.tex"),
     width = 1.3 * doc_width_inches, 
     height = 1.25 *  doc_width_inches,
     standAlone = TRUE)
comb_r
dev.off()

tinytex::lualatex(file.path(plots_path, "landmark-reg-plot.tex"))

