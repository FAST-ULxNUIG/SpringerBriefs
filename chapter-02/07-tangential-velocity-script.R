# 1) Import Packages ------------------------------------------------------
library(fda)        # CRAN v5.1.9
library(tidyverse)  # CRAN v1.3.0
library(ggpubr)     # CRAN v0.4.0
library(registr)    # CRAN v1.0.0
library(Polychrome) # CRAN v1.2.6
library(ggtext)     # [github::wilkelab/ggtext] v0.1.0
library(rmatio)     # CRAN v0.18.0
library(tikzDevice)

# 2) Settings -------------------------------------------------------------
functions_path <- here::here("functions")
data_path <- here::here("chapter-02", "juggling-data")
plots_path <- here::here("chapter-02", "figures")
source(file.path(functions_path, "landmarkreg-fixed.R"))

theme_set(new = theme_bw())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12),
             axis.ticks = element_blank(),
             plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
             plot.subtitle = element_text(size = 12, hjust = 0.5))



# 3) Import Data -------------------------------------------------------------
# Juggling data should be downloaded and unzipped from the following link:
# https://archive.mbi.ohio-state.edu/vod/mbi-media/e-162/1400789954_juggling.zip
# the stored in a folder called Juggling within the working directory

## 3.1) From Text Files -------------------------------------------------------
x <- read_table2(file.path(data_path,"Xcoord.txt"), col_names = FALSE) # x coordinates
y <- read_table2(file.path(data_path,"Ycoord.txt"), col_names = FALSE) # y coordinates
z <- read_table2(file.path(data_path,"Zcoord.txt"), col_names = FALSE) # z co-ordinates
ni <- as.vector(read_csv(file.path(data_path,"ni.txt"), col_names = FALSE))$X1 # no. of cycles in trial
seqn <- as.vector(t(read_table2(file.path(data_path,"seqn.txt"), col_names = FALSE)[1, ])) # no. of observations

## 3.2) From Matlab files -----------------------------------------------------
landmarkregCell <- R.matlab::readMat(file.path(data_path,"landmarkregCell.mat"))
recordfd0Cell <- R.matlab::readMat(file.path(data_path,"recordfdCell0.mat"))
Afeaturetimes <- rmatio::read.mat(file.path(data_path,"Afeaturetimes.mat"))
Vfeaturetimes <- rmatio::read.mat(file.path(data_path,"Vfeaturetimes.mat"))
# Fix clumsiness of reading Matlab files into R

Afeaturetimes <- Afeaturetimes$Afeaturetimes
Vfeaturetimes <- Vfeaturetimes$Vfeaturetimes
#  The tangential velocity features for each cycle as follows:
#  1. a minor local maximum in tangential velocity positioned between 
#     the minimum and the primary maximum (launch)
#  2. a major local maximum in tangential velocity (drop), and
#  3. a primary local minimum in tangential velocity (handoff)

dimnames(Vfeaturetimes)[[1]] <- c("launch", "drop", "handoff")

# From Matlab File:
#  The tangential acceleration features for each cycle are as follows:
#  1. a major acceleration peak 
#  2. a minor peak
#  3. a minimum
#  4. a small maximum, but small only because the ball is also being
#     accelerated
#  5. a minimum
dimnames(Afeaturetimes)[[1]] <- c("major_peak", "minor_peak", "minimum", "small_maximum", "minimum")


# 4) Copy of data processing from Matlab Scripts -------------------------------------------------
ni_min1 <- ni -1
timei <- 5 * ni_min1
(timepercycle = timei/seqn)  #  mean duration per cycle
(meantimepercycle = mean(timepercycle)) # average over records




# 5) Perform landmark registration on combined data ----------------------------------------
# on COMBINED data (background, grey in plot)
# do this for each separately for each trial (1:10)
# through a loop
# using the landmark timings for the launch velocity that are given
landmark_registered_TV <- vector(mode = "list", length = 10)
unregistered_TV <- vector(mode = "list", length = 10)
launch_peak_1 <- vector(length = 10)
for(i in 1:10){
  xyz_mat <- as.matrix(cbind(x[1:ni[i], i], y[1:ni[i], i], z[1:ni[i], i]), ncol = 3)
  time_1 <- 5 * 0: (ni[i]-1)
  time_1 <- time_1/1000
  
  fd_rec1 <- smooth.basis(argvals = time_1,
                          y = xyz_mat, returnMatrix = F,
                          fdParobj = create.bspline.basis(rangeval = range(time_1), nbasis = seqn[i] * 10))
  
  vel_eval <- eval.fd(evalarg = time_1, fdobj = fd_rec1$fd, Lfdobj = 1)
  
  TV_rec1 <- sqrt(apply(X = vel_eval^2, MARGIN = 1, FUN = sum))
  
  wbasis <- create.bspline.basis(rangeval = range(time_1), nbasis = 113, norder = 6)
  Wfd0   <- fd(matrix(data = 0, wbasis$nbasis, 1), wbasis)
  WfdPar <- fdPar(fdobj = Wfd0, Lfdobj = 4, 10^-8)
  
  
  launch_times <- Vfeaturetimes[1,,i]
  
  while(tail(launch_times, 1) == 0){
    launch_times <- head(launch_times, -1)
  }
  
  equally_spaced <- seq(launch_times[1], by = meantimepercycle/1000, length.out = length(launch_times))
  
  unregistered_TV[[i]] <- data.frame(time = time_1 - launch_times[1], unreg = TV_rec1) %>% filter(time >= 0)
  
  plot(time_1, TV_rec1, type = "l", lty = 1, lwd = 2)
  abline(v = launch_times)
  
  
  wbasis <- create.bspline.basis(rangeval = range(time_1), nbasis = (length(launch_times)+1) * 9, norder = 6)
  Wfd0   <- fd(matrix(data = 0, wbasis$nbasis, 1), wbasis)
  WfdPar <- fdPar(fdobj = Wfd0, Lfdobj = 4, 10^-8)
  
  fd_rec1$fd$basis$rangeval
  
  landmark_reg <- landmarkreg(fdobj = fd_rec1$fd,
                              ximarks = matrix(rep(launch_times, each = 3), ncol = length(launch_times), nrow = 3),
                              WfdPar =  WfdPar,
                              x0marks = equally_spaced)
  
  vel_reg_eval <- eval.fd(evalarg = time_1, fdobj = landmark_reg$regfd, Lfdobj = 1)
  TV_rec1_reg <- sqrt(apply(X = vel_reg_eval^2, MARGIN = 1, FUN = sum))
  landmark_registered_TV[[i]] <- data.frame(time = time_1 - launch_times[1], reg = TV_rec1_reg) %>% filter(time >= 0)
}




# 6) Look at trial 1 individually -----------------------------------------------------
i = 1
xyz_mat <- as.matrix(cbind(x[1:ni[i], i], y[1:ni[i], i], z[1:ni[i], i]), ncol = 3)
time_1 <- 5 * 0: (ni[i]-1)
time_1 <- time_1/1000

fd_rec1 <- smooth.basis(argvals = time_1,
                        y = xyz_mat, returnMatrix = F,
                        fdParobj = create.bspline.basis(rangeval = range(time_1), nbasis = seqn[i] * 10))

vel_eval <- eval.fd(evalarg = time_1, fdobj = fd_rec1$fd, Lfdobj = 1)

TV_rec1 <- sqrt(apply(X = vel_eval^2, MARGIN = 1, FUN = sum))

wbasis <- create.bspline.basis(rangeval = range(time_1), nbasis = 113, norder = 6)
Wfd0   <- fd(matrix(data = 0, wbasis$nbasis, 1), wbasis)
WfdPar <- fdPar(fdobj = Wfd0, Lfdobj = 4, 10^-8)

launch_times <- Vfeaturetimes[1,,i]

while(tail(launch_times, 1) == 0){
  launch_times <- head(launch_times, -1)
}

equally_spaced <- seq(launch_times[1], by = meantimepercycle/1000, length.out = length(launch_times))


wbasis <- create.bspline.basis(rangeval = range(time_1), nbasis = (length(launch_times)+1) * 9, norder = 6)
Wfd0   <- fd(matrix(data = 0, wbasis$nbasis, 1), wbasis)
WfdPar <- fdPar(fdobj = Wfd0, Lfdobj = 4, 10^-8)

fd_rec1$fd$basis$rangeval

landmark_reg <- landmarkreg(fdobj = fd_rec1$fd,
                            ximarks = matrix(rep(launch_times, each = 3), ncol = length(launch_times), nrow = 3),
                            WfdPar =  WfdPar,
                            x0marks = equally_spaced)

vel_reg_eval <- eval.fd(evalarg = time_1, fdobj = landmark_reg$regfd, Lfdobj = 1)
TV_rec1_reg <- sqrt(apply(X = vel_reg_eval^2, MARGIN = 1, FUN = sum))


reg_fd_eval <- eval.fd(evalarg = time_1, fdobj = landmark_reg$regfd, Lfdobj = 0)
unreg_fd_eval <- eval.fd(evalarg = time_1, fdobj =  fd_rec1$fd, Lfdobj = 0)

TV_launc_unreg <- sqrt(apply(eval.fd(launch_times, fdobj = fd_rec1$fd, Lfdobj = 1)^2,
                             MARGIN = 1, FUN = sum))
TV_launc_reg <- sqrt(apply(eval.fd(equally_spaced, fdobj = landmark_reg$regfd, Lfdobj = 1)^2,
                           MARGIN = 1, FUN = sum))
launch_points_unreg <- data.frame(time = launch_times - launch_times[1], tan_vel = TV_launc_unreg)
launch_points_reg <-  data.frame(time = equally_spaced - launch_times[1], tan_vel =  TV_launc_reg)
launch_points_unreg$registration = "unreg"
launch_points_reg$registration = "reg"
launch_points <- bind_rows(launch_points_reg, launch_points_unreg) %>%
  mutate(registration = factor(x = registration, levels = c("unreg", "reg"), labels = c("Unregistered", "Landmark Registered")))
launch_points$trial <- 1



# Plot of landmark registered overlaid on unregistered -------------------------------------------------------------------
# the 1st trial will be huighlighted
unreg_TV <- data.frame(data.table::rbindlist(l = unregistered_TV, idcol = "trial"))
landreg_TV <-  data.frame(data.table::rbindlist(l = landmark_registered_TV, idcol ="trial"))

flabels <- data.frame(x = rep(8.89, 2), y = c(2.4,2.4),
                      label = c(" <span style='color:red'>**Launch Peak**</span> <br> not aligned <br> to target ",
                                " <span style='color:red'>**Launch Peak**</span> <br> aligned <br> to target "),
                      registration =  factor(c("Unregistered", "Landmark Registered"), levels = c("Unregistered", "Landmark Registered")))

tv_plot <- inner_join(unreg_TV, landreg_TV, by = c("trial", "time")) %>%
  gather(-trial, -time, key = "registration", value = "tan_vel") %>%
  mutate(registration = factor(x = registration, levels = c("unreg", "reg"), labels = c("Unregistered", "Landmark Registered"))) %>%
  ggplot() +
  geom_vline(xintercept = equally_spaced - launch_times[1], col = "black", lty =2, alpha = 0.7) +
  facet_wrap(~registration,
             nrow = 2) +
  aes(x = time, y = tan_vel, group = trial) +
  geom_line(color = "grey") +
  geom_line(data = . %>% filter(trial == 1), col = "black", lwd = 0.8) +
  geom_point(data = launch_points, pch = 16, size = 2, col = "red") +
  scale_x_continuous(expand = c(0.01, 0.01))+
  labs(y = "$\\sqrt{x'(t)^2 + y'(t)^2 + z'(t)^2}$",
       x = "Time $t$ in seconds")+
  geom_richtext(mapping = aes(x =x, y =y, label = label),
                data = flabels, inherit.aes = F,
                size = (5/ 14) * 8, label.padding = unit(0.1, units ="line"))+
  geom_curve(mapping = aes(x = 8.3, xend = 7.95, y = 2.5, yend = 2.1),
             curvature = 0.35, color = "black",
             arrow = arrow(length = unit(0.05, "inches"), type = "closed"), size = 0.5) +
  theme(strip.text = element_text(face = "bold"))


tv_plot

tikz(file.path(plots_path, "tv_plot.tex"),
     width = 7.87402, 
     height = 7.87402/2,
     standAlone = TRUE)
tv_plot
dev.off()

tinytex::lualatex(file.path(plots_path, "tv_plot.tex"))

# ggsave(plot = tv_plot, filename = "Figurestv_plot.png",
#        device = "png", units = "cm", width = 20, height = 10)


# Others ------------------------------------------------------------------

# some other options for plots that were not considered:
colnames(reg_fd_eval) <- colnames(unreg_fd_eval) <- c("x", "y", "z")


reg_lng <- data.frame(time = time_1, reg_fd_eval) %>%
  gather(-time, key = "coord", value = "reg")

unreg_lng <- data.frame(time = time_1, unreg_fd_eval) %>%
  gather(-time, key = "coord", value = "unreg")

inner_join(reg_lng, unreg_lng, by = c("time", "coord")) %>%
  gather(-time, -coord, key = "registered", value = "value") %>%
  ggplot() +
  facet_wrap( ~ coord, nrow = 3, scales = "free_y") +
  aes(x = time, y = value, color = registered, group = registered) +
  geom_line(lwd = 0.7) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("red4", "cornflowerblue")) +
  theme(legend.position = "none")



data.frame(time_1, reg = TV_rec1_reg, unreg = TV_rec1) %>%
  gather(-time_1, key = "registered", value = "tv") %>%
  ggplot() +
  geom_vline(xintercept = equally_spaced, col = "darkgrey", lty =2) +
  aes(x = time_1, group = registered, color = registered, y = tv)+
  geom_line()  +
  scale_x_continuous(expand = c(0, 0))



