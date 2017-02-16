#                           --- Simulate SCR data ---

library(secr)
library(fields)

sim.open.scr.pop <- function(lambda0, sigma, occasions, phi, D = 0.1, 
                             seed=12345) {
  
  # create the detectors
  detype <- "count"
  spacing <- sigma
  dets <- make.grid(nx=7, ny=7, spacing=spacing, detector=detype)
  
  # create the mesh
  mesh <- make.mask(dets, buffer=4*sigma, nx=64, ny=64, type="trapbuffer")
  
  # simulate a population from this intensity surface
  pop <- sim.popn(D=D, core=mesh, buffer=0, seed=seed)
  
  # multi-occasion capture history
  capthist <- sim.capthist(dets, popn=pop, detectpar=list(lambda0=lambda0, sigma=sigma), 
                           detectfn="HHN", noccasions=occasions, nsessions=1, 
                           seed=seed)
  
  # plot simulated population and capture history
  plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
  points(pop$x,pop$y,pch=19,cex=0.25)
  plot(dets,add=TRUE)

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
  return(list("capthist"=capthist, "mesh"=mesh))
}

# simulate an example population
lambda0 <- 0.1 ; sigma <- 800 ; occasions <- 20 ; phi <- 0.9 ; D <- 1
sim <- sim.open.scr.pop(lambda0, sigma, occasions, phi)
capthist.sim <- sim$capthist
mesh <- sim$mesh

plot(capthist.sim, border=sigma, tracks=TRUE, gridlines=FALSE, rad=3, add=TRUE)
