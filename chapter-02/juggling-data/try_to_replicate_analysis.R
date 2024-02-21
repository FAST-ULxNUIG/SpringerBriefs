# 1) Import Packages ------------------------------------------------------
library(fda)
library(tidyverse)


# 2) Import Data -------------------------------------------------------------


## 2.1) From Text Files -------------------------------------------------------
x <- read_table2("juggling/Xcoord.txt", col_names = FALSE) # x coordinates
y <- read_table2("juggling/Ycoord.txt", col_names = FALSE) # y coordinates
z <- read_table2("juggling/Zcoord.txt", col_names = FALSE) # z co-ordinates
ni <- as.vector(read_csv("juggling/ni.txt", col_names = FALSE)) # no. of cycles in trial
seqn <- as.vector(t(read_table2("juggling/seqn.txt", col_names = FALSE)[1, ])) # no. of observations


## 2.2) From Matlab files -----------------------------------------------------
landmarkregCell <- R.matlab::readMat("juggling/landmarkregCell.mat")
recordfd0Cell <- R.matlab::readMat("juggling/recordfdCell0.mat")
Afeaturetimes <- rmatio::read.mat(filename = "juggling/Afeaturetimes.mat")
Vfeaturetimes <- rmatio::read.mat(filename = "juggling/Vfeaturetimes.mat")
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



# 3) Explore Data --------------------------------------------


## 3.1) Low Dimensional Representation Data ------------------------------------------

# From Matlab Script:
#  Cell array ContuousRegCell stores smooths of the the data by 
#  unpenalized least squares using a minimal basis.
#  The part of each record before the first velocity feature time for
#  the first cycle is removed.

#  The basis for this smooth is constructed by using nine knots
#  per cycle, and the smoothing is by unpenalized least squares.
#  The derivative estimates are unstable at the end points, but this
#  basis is much more suitable for registering the data
#  Because the records vary in length, this basis must be set up
#  separately for each record

#  We will extract some data for the first record: 

# let's subset the feature vectors to contain only the first record first:
Afeatures_1 <- Afeaturetimes[,,1]
Vfeatures_1 <- Vfeaturetimes[,,1]
# recordfd0Cell$recordfdCell0 elemennts 1:10 contain 
# the time argumens for trials 1:10
record_1_time <- recordfd0Cell$recordfdCell0[[1]][[1]] %>% as.vector
recordfd0Cell$recordfdCell0[[8]][[1]][1]


plot(record_1_time)
# recordfd0Cell$recordfdCell0 elemennts 11:20 contain 
# evaluations of x y z coordinates for trials 1:10
record_1_xyz <- recordfd0Cell$recordfdCell0[[11]][[1]]
par(mfrow = c(3, 1))
for(p in 1:3){
  plot(record_1_time, record_1_xyz[, p], type = "l", lty = 1)
  abline(v = Vfeatures_1[1, 1:12])
}

# recordfd0Cell$recordfdCell0 elements 21:30 contain 
# fd objects storing the xyz coordinates:
# it should be an S3 class but some problems with basis
# being recognised because of how it's set up in R
# we can set it back up as an FD object in R:
record_1_xyz.fd <- recordfd0Cell$recordfdCell0[[21]][[1]] # extract
coef_r1 <- record_1_xyz.fd$coef # basis coefficients
basis_r1 <- record_1_xyz.fd$basisobj # basisobj, this will need reconstructing from what we know
create.bspline.basis(rangeval = c(0, 100), norder = 4, breaks = seq(0, 100, by = 1))
basis_r1_reconstruct <- create.bspline.basis(rangeval = basis_r1$rangeval, # rangeval given in object
                                             nbasis = basis_r1$nbasis, # nbasis given
                                             norder = 6, # this is a guess since 1) it satisfies breaks and nbasis and 2) it facilitates acceleration estimates.
                                             breaks = c(basis_r1$rangeval[1], basis_r1$params, basis_r1$rangeval[2])) # internal knots given
record_1_xyz.fd_recon <- fd(coef = coef_r1, basisobj = basis_r1_reconstruct) # create fd object.
par(mfrow = c(3, 1))
plot(record_1_xyz.fd_recon[1]) # from FD object
lines(record_1_time, record_1_xyz[, 3], type = "l", lty = 1, col = "red") # from evaluations to check if it worked:


# recordfd0Cell$recordfdCell0 elements 31:40 contain 
# fd objects storing the tangential velocity:
record_1_TV.fd <- recordfd0Cell$recordfdCell0[[31]][[1]] # extract
coef_r1_TV <- record_1_TV.fd$coef
basis_r1_TV <- record_1_TV.fd$basisobj
# I don't think we can reconstruct this basis....



## 3.2) Landmark Registration Data-------------------------------------------------

# From Matlab files:

#  Landmark registration of each curve using first velocity peak.
#  The first velocity peak is at or just before the throw, and used
#  to define the beginning and end of each cycle.  This registration
#  approximately normalizes the cycle durations.

#  load cellarray containing results

landmarkregCell <- landmarkregCell$landmarkregCell

landmarkregCell[[21]]

# 1:10 is the landmark registered data for trials 1:10.
landmark_r1 <- landmarkregCell[[1]][[1]]
landmark_r1.fd <- fd(coef = landmark_r1$coef,
                     basisobj = create.bspline.basis(rangeval = landmark_r1$basisobj$rangeval,
                                   nbasis = landmark_r1$basisobj$nbasis,
                                   norder = 6,
                                   breaks = c(landmark_r1$basisobj$rangeval[1],
                                              landmark_r1$basisobj$params,
                                              landmark_r1$basisobj$rangeval[2])
                                   ))

plot(landmark_r1.fd)

landmark_r1_eval <- eval.fd(evalarg = record_1_time, fdobj = landmark_r1.fd)

### Recreate Plot (a) -------------------------------------------------------


# recreate Figure 1 (a) of Ramsay, Gribble and Kurtek (2014).
# We'll add lines for the observed and target landmarks.
par(mfrow = c(3, 1))
for(dim in 1:3){
  plot(record_1_time, record_1_xyz[, dim], type = "l")
  lines(record_1_time, landmark_r1_eval[,, dim], type = "l", col = "red")
  abline(v = Afeatures_1[1, 1:12], col = "black", lty = 2)
  abline(v = seq(Afeatures_1[1, 1], length.out = 12, by =  713/1000), col = "blue", lty = 2)
}


tangential_registered_record1 <- eval.fd(record_1_time, fdobj = landmark_r1.fd)
TV_rec1 <- sqrt(apply(X = tangential_registered_record1^2, MARGIN = 1, FUN = sum))
tangential_unregistered_record1 <- eval.fd(record_1_time, fdobj = record_1_xyz.fd_recon)
TV_rec1 <- sqrt(apply(X = tangential_registered_record1^2, MARGIN = 1, FUN = sum))
TV_rec1_unreg <- sqrt(apply(X = tangential_unregistered_record1^2, MARGIN = 1, FUN = sum))
plot(record_1_time, TV_rec1, type = "l")
abline(v = seq(Vfeatures_1[1, 1], length.out = 12, by =  713/1000) + , col = "blue", lty = 2)
lines(record_1_time, TV_rec1_unreg)
abline(v = Vfeatures_1[1, 1:12], col = "red", lty = 2)

plot(record_1_time, TV_rec1_unreg, type = "l")
abline(v = Vfeatures_1[1, 1:12] + record_1_time[1], col = "red", lty = 2)

### Recreate Plot (b) -------------------------------------------------------
# Let's look at the warping functions:
# for trials 1:10 they're stored in 11:20
warp_coef <- landmarkregCell[[11]][[1]]$coef
warp_basis <- landmarkregCell[[11]][[1]]$basisobj

warp_fd <- fd(create.bspline.basis(rangeval = warp_basis$rangeval, nbasis = warp_basis$nbasis, norder = 6,  breaks = c(warp_basis$rangeval[1],
                                                                                            warp_basis$params, 
                                                                                            warp_basis$rangeval[2])), coef = warp_coef)
plot(record_1_time, eval.fd(evalarg = record_1_time, fdobj = warp_fd) - record_1_time, type = "l", ylim = c(-0.03, 0.13))
#abline(v = as.vector(seq(Vfeatures_1[1, 1], length.out = 12, by =  713/1000)))
#abline(v = Vfeatures_1[1, 1:12], col = "red")

# let's try the landmark Registration ourselves
# need the full trials:






wbasis <- create.bspline.basis(rangeval = record_1_xyz.fd_recon$basis$rangeval,
                               nbasis = 113, norder = 6)
Wfd0   <- fd(matrix(data = 0, warp_basis$nbasis, 1), wbasis)
WfdPar <- fdPar(fdobj = Wfd0, Lfdobj = 4, 10^-8)



landmark_reg <- landmarkreg(fdobj = record_1_xyz.fd_recon,
                            ximarks = matrix(c(Vfeatures_1[1, 2:12]) - Vfeatures_1[1, 1], nrow = 1, ncol = 11),
                            WfdPar =  WfdPar,
                            x0marks = seq(meantimepercycle/1000, length.out = 11, by =  meantimepercycle/1000))

plot(record_1_time, eval.fd(evalarg = record_1_time, fdobj = warp_fd) - record_1_time,type = "l", xlab = "Seconds", ylab = "Deformation Function = h(t) - t", col = "red", ylim = c(-0.03, 0.13))
lines(record_1_time, eval.fd(evalarg = record_1_time, fdobj = landmark_reg$warpfd) - record_1_time)
abline( h = 0, lty = 2)
# Almost the same, we'll take it.


# Continuous Registration -------------------------------------------------

# Register the Tangential Velocity to ewually spaced points:
nbasis0 = 3 #  number of basis functions per cycle for warping
norder  = 6
nknots0 = 9
period  = meantimepercycle


# Tangential velocity is the sqrt(x'(t)^2 + y'(t) ^2 + x'(t) ^ 2)
# first evaluate velocity:
velocity_eval_seq <- seq(landmark_r1.fd$basis$rangeval[1],
                         landmark_r1.fd$basis$rangeval[2],
                         length.out = 2000)
velocity_eval <- eval.fd(evalarg = velocity_eval_seq,
                         fdobj = landmark_reg$regfd,
                         Lfdobj = 1)
# Calculate sqrt(x'(t)^2 + y'(t) ^2 + x'(t) ^ 2) by summing over columns:
tangential_velocity <- sqrt(apply(X = velocity_eval^2, MARGIN = 1, FUN = sum))
plot(x = velocity_eval_seq, y = tangential_velocity, type = "l")
abline(v = seq(0, length.out = 13, by =  meantimepercycle)/1000, col = "red")




Ti = period*seqn[1] # end of interval
rngi = c(0,Ti)   #  range for registered functions
rngi/1000
startime = 0
tvec <- velocity_eval_seq
indi <- which(tvec <= max(rngi)/1000)
tvec <- tvec[indi]
ntvec <- length(tvec)
tvec <- seq(0, Ti[1], length.out = ntvec)

sclfac   = Ti/(tvec[ntvec]-tvec[1])
tvec     = tvec*sclfac;


# Non-peridoic Smooth -----------------------------------------------------
Nnbasis   = seqn[1] * nknots0 + 5
Nbasisobj = create.bspline.basis(rangeval = rngi, nbasis = Nnbasis, norder = norder)
NfdParobj = fdPar(Nbasisobj)
data     = tangential_velocity[indi]
Nfdobj   = smooth.basis(argvals = tvec,y = data, fdParobj = NfdParobj)
plot(Nfdobj)
abline(v = seq(0, length.out = 13, by =  meantimepercycle), col = "red")



# Periodic Smooth ---------------------------------------------------------
Pnbasis = 10

Pbasisobj = create.fourier.basis(rangeval = c(0, Ti), nbasis = Pnbasis, period = period)
PfdParobj = fdPar(Pbasisobj)
# periodic smooth of the data
Pfdobj   = smooth.basis(argvals = tvec,y = data, fdParobj = PfdParobj);
lines(Pfdobj)

#  set up a basis for the function W defining warping function h
wbasisi = create.bspline.basis(rangeval = rngi, norder = 4, nbasis = seqn[1] * nbasis0 + 3)

Wfd0   <- fd(matrix(data = 0, wbasisi$nbasis, 1), wbasisi)
WfdPari = fdPar(fdobj = Wfd0,Lfdobj = 2, lambda = 10^-8)



Rfd_p <- register.fd(y0fd = Pfdobj$fd, yfd = Nfdobj$fd, WfdParobj = WfdPari)

plot(Rfd_p$regfd)
lines(Pfdobj$fd, col)

plot(tvec, eval.fd(fdobj = Rfd_p$warpfd, tvec) - tvec, type = "l")
abline( h = 0)

# ok...... that didn't work.