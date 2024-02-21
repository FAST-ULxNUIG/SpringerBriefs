# -------------------------------------------------------------------------
center_fd_around_new_mean <- function(fdobj, mean.fd.obj)
{
  #  center functional data around a different mean
  #  useful for doing training and testing with functional data
  # i.e., for centering a test set of functional data around
  # the training set's mean.
  
  if (!(is.fd(fdobj) || is.fdPar(fdobj))) 
    stop("First argument is neither an fd or an fdPar object.")
  if (is.fdPar(fdobj)) fdobvj = fdobj$fd
  
  if (!(is.fd(mean.fd.obj) || is.fdPar(mean.fd.obj))) 
    stop("Second argument is neither an fd or an fdPar object.")
  if (is.fdPar(mean.fd.obj)) fdobvj = mean.fd.obj$fd
  
  coef     <- as.array(fdobj$coefs)
  coefd    <- dim(coef)
  ndim     <- length(coefd)
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  coefmean <- mean.fd.obj$coefs
  
  if(!(length(dim(coefmean)) == length(dim(coef)))) {
    stop("Dimensions of Mean Function and data don't match.")
  }
  if(length(dim(coefmean)) > 2 || length(dim(coef)) > 2) {
    if(! all((dim(coefmean)[c(1, 3)] == dim(coef)[c(1, 3)]))) {
      stop("Dimensions of Mean Function and data don't match.")
    }
  }
  
  if(!(dim(coefmean)[2] == 1)) {
    stop("Mean fd object contains replicate observations.")
  }
  
  if (ndim == 2) {
    coef     <- sweep(coef,1,coefmean[,1])
  } else {
    nvar <- coefd[3]
    for (j in 1:nvar) {
      coef[,,j] <- sweep(coef[,,j, drop = FALSE],1, coefmean[,1,j]) # drop = FALSE added to stop bug when only one observation.
    }
  }
  fdnames      <- fdobj$fdnames
  fdnames[[3]] <- paste("Centered", fdnames[[3]])
  centerfdobj  <- fd(coef, basisobj, fdnames)
  return(centerfdobj)	
}

