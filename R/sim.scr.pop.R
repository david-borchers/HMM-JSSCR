#                           --- Simulate SCR data ---

library(secr)
library(fields)

sim.scr.pop <- function(lambda0, sigma, occasions, seed=123) {
  
  # create the detectors
  detype <- "count"
  spacing <- sigma
  dets <- make.grid(nx=7, ny=7, spacing=spacing, detector=detype)
  
  # create the mesh
  mesh <- make.mask(dets, buffer=4*sigma, nx=64, ny=64, type="trapbuffer")
  
  # set density via abundance on mesh
  Nmesh <- 100
  M <- dim(mesh)[1]
  D <- Nmesh / (M * attributes(mesh)$area)
  
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
  plot(capthist, border=sigma, tracks=TRUE, gridlines=FALSE,rad=3,add=TRUE)
  
  return(list(mesh=mesh, capthist=capthist))
  
}

# simulate an example population
lambda0 <- 2 ; sigma <- 20 ; occasions <- 5
scr.pop <- sim.scr.pop(lambda0=lambda0, sigma=sigma, occasions=occasions)
pop.capthist <- scr.pop$capthist
mesh <- scr.pop$mesh

# randomly generate survival times with known survival probability
n <- dim(pop.capthist)[1]
phi <- 0.5
set.seed(123)
survival.occ <- rgeom(n, phi)

# thin capture history 
for (i in 1:n) {
  if (survival.occ[i] < occasions) {
    pop.capthist[i, (survival.occ[i]+1):occasions, ] <- 0
  }
}

