library(secr)
library(fields)
setwd("C:/Users/Matthew/Desktop/MMath project/R/R")

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

jag.sim.open.pop <- function(lambda0, sigma, occasions, phi, N, D, mesh, dets,
                             seed=12345) {
  
  # simulate a population from this intensity surface
  pop <- sim.popn(D=D, core=mesh, Nbuffer=N, buffer=0, poly=mesh, 
                  seed=seed)
  
  # trap usage
  det.usage <- attributes(dets)$usage[,1]
  usage.mat <- as.matrix(attributes(dets)$usage)
  for (i in 2:occasions) {
    usage.mat <- cbind(usage.mat, round(det.usage*phi**(i-1)))
  }
  dimnames(usage.mat)[[2]] <- as.character(1:occasions)
  attributes(dets)$usage <- usage.mat
  
  # multi-occasion capture history
  capthist <- sim.capthist(dets, popn=pop, detectpar=list(lambda0=lambda0, sigma=sigma), 
                           detectfn="HHN", noccasions=occasions, nsessions=1, 
                           seed=seed)
  
  # plot simulated population and capture history
  plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
  points(pop$x,pop$y,pch=19,cex=0.4)
  plot(dets, add=TRUE)
  plot(capthist, border=sigma, tracks=TRUE, gridlines=FALSE, rad=3, add=TRUE)
  
  # randomly generate survival times with known survival probability
  n <- dim(capthist)[1]
  set.seed(seed)
  survival.occ <- rgeom(n, 1-phi)
  
  # thin capture history 
  for (j in 1:n) {
    if (survival.occ[j]+1 < occasions) {
      capthist[j, (survival.occ[j]+2):occasions, ] <- 0
    }
  }
  
  # remove capture history of individuals not detected after thinning
  thin.capthist <- create.thinned.capthist(capthist)
  
  # capture history for an individual not detected 
  zero.capthist <- create.zero.capthist(capthist)
  
  # plot simulated population and capture history
#  plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
#  points(pop$x,pop$y,pch=19,cex=0.4)
#  plot(dets, add=TRUE)
#  plot(thin.capthist, border=sigma, tracks=TRUE, gridlines=FALSE, rad=3, add=TRUE)
  
  
  return(list("capthist"=thin.capthist, "zero.capthist"=zero.capthist))
}

#---------------

load(file="jagCH1.rda")
load(file="jagmask.rda") # loads mask object jagmask
dets <- traps(jagCH1)

# simulate jaguar capture history
lambda0 <- 0.024 ; sigma <- 3250 ; occasions <- 5 ; phi <- 0.1; 
D <- 1; N <- 75
sim <- jag.sim.open.pop(lambda0, sigma, occasions, phi, N, D, jagmask, dets,
                        seed=12345)

dist <- distances(dets, jagmask) # calculate trap distances 
pars <- c(lambda0=log(lambda0), sigma=log(sigma), phi=qlogis(phi)) # starting values
