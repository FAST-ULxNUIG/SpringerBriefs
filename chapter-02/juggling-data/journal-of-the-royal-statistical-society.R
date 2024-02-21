# 1) Import Packages ------------------------------------------------------
library(fda)
library(tidyverse)
library(ggpubr)
library(registr)
library(Polychrome)

# 2) Settings -------------------------------------------------------------
theme_set(new = theme_bw())
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 12),
             strip.text = element_text(size = 12),
             axis.ticks = element_blank(),
             plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
             plot.subtitle = element_text(size = 12, hjust = 0.5))

# 3) Read Data From Matlab files -----------------------------------------------------
landmarkregCell <- R.matlab::readMat("juggling/landmarkregCell.mat")
recordfd0Cell <- R.matlab::readMat("juggling/recordfdCell0.mat")
Afeaturetimes <- rmatio::read.mat(filename = "juggling/Afeaturetimes.mat")
Vfeaturetimes <- rmatio::read.mat(filename = "juggling/Vfeaturetimes.mat")
ni <- as.vector(read_csv("juggling/ni.txt", col_names = FALSE)) # no. of cycles in trial
seqn <- as.vector(t(read_table2("juggling/seqn.txt", col_names = FALSE)[1, ])) # no. of observations
n_trials <- 10
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

# Copy from Matlab Script -------------------------------------------------
ni_min1 <- ni -1
timei <- 5 * ni_min1
(timepercycle = timei/seqn)  #  mean duration per cycle
(meantimepercycle = mean(timepercycle$X1)) # average over records


# let's loop through and landmark register for each trial:
par(mfrow = c( 6, 1))

TV_list <- vector(mode = "list", length = 10)
for(i in 1:10){
  #======================================================#
  # Extract Smoothed Data from recordfd0Cell$recordfdCell0
  
  # recordfd0Cell$recordfdCell0 elements 21:30 contain 
  # fd objects storing the xyz coordinates:
  # it should be an S3 class but some problems with basis
  # being recognised because of how it's set up in R
  # we can set it back up as an FD object in R:
  record_xyz.fd <- recordfd0Cell$recordfdCell0[[(20 + i)]][[1]] # extract
  coef_r <- record_xyz.fd$coef # basis coefficients
  basis_r <- record_xyz.fd$basisobj # basisobj, this will need reconstructing from what we know
  basis_r_reconstruct <- create.bspline.basis(rangeval = basis_r$rangeval, # rangeval given in object
                                               nbasis = basis_r$nbasis, # nbasis given
                                               norder = 6, # this is a guess since 1) it satisfies breaks and nbasis and 2) it facilitates acceleration estimates.
                                               breaks = c(basis_r$rangeval[1], basis_r$params, basis_r$rangeval[2])) # internal knots given
  record_xyz.fd_recon <- fd(coef = coef_r, basisobj = basis_r_reconstruct) # create fd object.
  
  #======================================================#
  # Landmark register each sycle to launch velocity
  
  # 1:10 is the landmark registered data for trials 1:10.
  landmark_r <- landmarkregCell$landmarkregCell[[ i ]][[1]]
  landmark_r.fd <- fd(coef = landmark_r$coef,
                      basisobj = create.bspline.basis(rangeval = landmark_r$basisobj$rangeval,
                                                      nbasis = landmark_r$basisobj$nbasis,
                                                      norder = 6,
                                                      breaks = c(landmark_r$basisobj$rangeval[1],
                                                                 landmark_r$basisobj$params,
                                                                 landmark_r$basisobj$rangeval[2])))
  


  # We could just use these for the plots
  # but we'll actually perform the landmark registration ourselves:
  
  # use the same warping functions as they used
  # these are stored in 11:20
  warp_coef <- landmarkregCell$landmarkregCell[[10 + i]][[1]]$coef
  warp_basis <- landmarkregCell$landmarkregCell[[10 + i]][[1]]$basisobj
  stopifnot(warp_basis$rangeval == record_xyz.fd_recon$basis$rangeval)
  wbasis <- create.bspline.basis(rangeval = warp_basis$rangeval,
                                 nbasis = warp_basis$nbasis, norder = 6)

  stopifnot(all(warp_basis$params == warp_basis$params))
  
  Wfd0   <- fd(matrix(data = 0, warp_basis$nbasis, 1), wbasis)
  WfdPar <- fdPar(fdobj = Wfd0, Lfdobj = 4, 10^-8)
  Vfeaturetimes[,, i]
  Vfeaturetimes[,, i]
  unequally_spaced_points <- matrix(Vfeaturetimes[1,2:n_cycle, i], nrow = 1, ncol = (n_cycle-1))
  equally_spaced_points <- seq(0.0471 + meantimepercycle/1000, length.out = (n_cycle -1), by =  meantimepercycle/1000)
  landmark_reg <- landmarkreg(fdobj = record_xyz.fd_recon,
                              ximarks = unequally_spaced_points,
                              WfdPar =  WfdPar,
                              x0marks = equally_spaced_points)
  plot(landmark_r.fd, main = paste0("-------------------------", "trial", i))
  lines(landmark_reg$regfd[3], col = "red")
  
  # Tangential velocity is the sqrt(x'(t)^2 + y'(t) ^2 + x'(t) ^ 2)
  # first evaluate velocity:
  velocity_eval_seq <- seq(landmark_r.fd$basis$rangeval[1],
                           landmark_r.fd$basis$rangeval[2],
                           length.out = 2000)
  velocity_eval_reg <- eval.fd(evalarg = velocity_eval_seq,
                           fdobj = landmark_reg$regfd,
                           Lfdobj = 1)
  velocity_eval_unreg <-  eval.fd(evalarg = velocity_eval_seq,
                              fdobj = record_xyz.fd_recon,
                              Lfdobj = 1)
  velocity_launch_reg <- eval.fd(evalarg = c(0, equally_spaced_points),
                                       fdobj = landmark_reg$regfd,
                                       Lfdobj = 1)
  
  velocity_launch_unreg <- eval.fd(evalarg = unequally_spaced_points[1, ],
                                 fdobj = record_xyz.fd_recon,
                                 Lfdobj = 1)
  
  velocity_eval_ramsay <- eval.fd(evalarg = velocity_eval_seq,
                               fdobj = landmark_r.fd,
                               Lfdobj = 1)
  
  
  
  # Calculate sqrt(x'(t)^2 + y'(t) ^2 + x'(t) ^ 2) by summing over columns:
  TV_reg <- sqrt(apply(X = velocity_eval_reg^2, MARGIN = 1, FUN = sum))
  TV_ramsay <- sqrt(apply(X = velocity_eval_ramsay^2, MARGIN = 1, FUN = sum))
  
  plot(TV_reg, type = "l")
  lines(TV_ramsay)
  
  TV_unreg <- sqrt(apply(X = velocity_eval_unreg^2, MARGIN = 1, FUN = sum))

  TV_launch_reg <- sqrt(apply(X = velocity_launch_reg^2, MARGIN = 1, FUN = sum))
  TV_launch_unreg <- sqrt(apply(X = velocity_launch_unreg^2, MARGIN = 1, FUN = sum))
  
  
  
 
  
  TV_rec <- data.frame(t = velocity_eval_seq, unregistered = TV_unreg, reg = TV_reg)
  TV_list[i] <- TV_rec
  
  }

