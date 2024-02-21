project_data_onto_fpcs <- function(fdobj, pca.fd_obj) {
  # adapted code from `fda` package for projecting
  # individual functions onto FPCs.
  # useful for when pca is calculated on training set but
  # we want fpc scores for the test set.
  
  
  stopifnot(class(pca.fd_obj) == "pca.fd")
  # NOTE: ASSUMES DATA REPRESENTED ON SAME BASIS AS FPCS
  # FOR MULTIVARIATE FUNCTIONAL DATA.
  stopifnot(fdobj$basis == pca.fd_obj$harmonics$basis)
  # this could be improved in time.
  
  
  harmfd <- pca.fd_obj$harmonics
  stopifnot(is.fd(harmfd))
  
  stopifnot(is.fd(fdobj))
  
  fd_coefs <- fdobj$coefs
  
  ndim  <- length(dim(fd_coefs))
  nrep <- dim(fd_coefs)[2]
  nharm <- dim(harmfd$coefs)[2]
  
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  type     <- basisobj$type
  
  if(ndim == 2) {
    nvar <- 1
  } else if(ndim == 3) {
    nvar <- dim(fd_coefs)[3]
  } else stop("Dimension of coefficient matrix wrong!")
  
  
  if (nvar == 1) {
    harmscr  <- inprod(fdobj, harmfd)
  } else {
    harmscr  <- array(0, c(nrep, nharm, nvar))
    coefarray <- fdobj$coefs
    harmcoefarray <- harmfd$coefs
    for (j in 1:nvar) {
      fdobjj  <- fd(as.matrix(coefarray[,,j]), basisobj)
      harmfdj <- fd(as.matrix(harmcoefarray[,,j]), basisobj)
      harmscr[,,j] <- inprod(fdobjj, harmfdj)
    }
  }
  harmscr
}



