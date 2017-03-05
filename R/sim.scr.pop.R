#                           --- Simulate SCR data ---

library(secr)
library(fields)

create.thinned.capthist <- function (capthist) {
  # attributes for thinned capture history
  attr <- attributes(capthist)
  nzeros <- summary(capthist)$zeros
  attr$dim[1] <- attr$dim[1] - nzeros 
  
  # remove individuals not detected and adjust attributes accordingly
  which.zero <- apply(capthist, 1, function(x) all(x==0))
  det.index <- names(which.zero[which.zero==F])
  attr$dimnames[[1]] <- det.index
  capthist <- capthist[as.numeric(det.index),,]
  attributes(capthist) <- attr
  
  return(capthist)
}

#---------------

create.zero.capthist <- function (capthist) {
  # attributes for zero capture history
  attr.zero <- attributes(capthist)
  ndets <- dim(capthist)[3]
  attr.zero$dim[1] <- 2
  attr.zero$dim[2] <- 1
  attr.zero$dimnames[[1]] <- c("1", "2")
  attr.zero$dimnames[[2]] <- "1"
  attr.zero$covariates <- NULL
  
  # create zero capthist where zero.capthist is a 2 by 1 by ndets matrix of zeros
  zero.capthist <- array(0, dim=c(2, 1, ndets))
  attributes(zero.capthist) <- attr.zero
  
  return(zero.capthist)
}

#---------------

sim.open.scr.pop <- function(lambda0, sigma, occasions, phi, N, D,
                             seed=12345) {
  
  # create the detectors
  detype <- "count"
  spacing <- sigma
  dets <- make.grid(nx=7, ny=7, spacing=spacing, detector=detype)
  
  # create the mesh
  mesh <- make.mask(dets, buffer=4*sigma, nx=64, ny=64, type="trapbuffer")
  
  # simulate a population from this intensity surface
  pop <- sim.popn(D=D, core=mesh, buffer=0, Nbuffer=N, buffertype = "convex", 
                  seed=seed)
  
  # multi-occasion capture history
  capthist <- sim.capthist(dets, popn=pop, detectpar=list(lambda0=lambda0, 
                                                          sigma=sigma), 
                           detectfn="HHN", noccasions=occasions, nsessions=1, 
                           seed=seed)
  
  # randomly generate survival times with known survival probability
  n <- dim(capthist)[1]
  set.seed(seed)
  survival.occ <- rgeom(n, 1-phi)
  
 
  # thin capture history 
  for (i in 1:n) {
    if (survival.occ[i] < occasions) {
      capthist[i, (survival.occ[i]+1):occasions, ] <- 0
    }
  }
  
  # remove capture history of individuals not detected after thinning
  thin.capthist <- create.thinned.capthist(capthist)
  
  # capture history for an individual not detected 
  zero.capthist <- create.zero.capthist(capthist)
  
  # plot simulated population and capture history
  plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
  points(pop$x,pop$y,pch=19,cex=0.25)
  plot(dets,add=TRUE)
  plot(thin.capthist, border=sigma, tracks=TRUE, gridlines=FALSE, rad=3, add=TRUE)
  
  
  return(list("capthist"=thin.capthist, "zero.capthist"=zero.capthist, "mesh"=mesh))
}

#---------------

# simulate an example population
lambda0 <- 0.1 ; sigma <- 500 ; occasions <- 5 ; phi <- 0.9 ; D <- 1; N <- 500
sim <- sim.open.scr.pop(lambda0, sigma, occasions, phi, N, D)