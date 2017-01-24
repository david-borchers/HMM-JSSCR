#                         --- HMM Likelihood Functions ---

# calculate distance between two pairs of coordinates
distances <- function (X, Y) {
  ## X and Y are 2-column matrices of coordinates
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}

# ---------------------
# calculates P_it(s) for SCR likelihood calculation.
calc.p.Pit.s <- function(capthist, mesh, pars, dist) {
  
  # calculate distances from each detector to each mesh point
  if(is.null(dist)) dist <- distances(traps(capthist), mesh)
  
  # Proximity detector with count data
  if (detector(traps(capthist))=="count") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    g0 <- exp(pars["g0"])
    sigma <- exp(pars["sigma"])
    er <- g0 * exp(-dist**2 / 2 / sigma**2) # Poisson encounter rate at all M distances from each detector
    
    # calculate log(P_it(s)): n x nt x M  matrix
    log.Pit.s <- count.log.Pit.s(capthist, er) # log(P_{it}(s))
    
  }
  
  return(exp(log.Pit.s))
  
}

# ---------------------
# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of count capture history (ith individual's capture history)
# of length K, containing capture frequencies at each detector
# er is K by M matrix containing encounter rates for each trap at each mesh point
count.log.Pi.si = function(wi,er) {
  one=rep(1,length(wi))
  log.Pi.si = wi %*% log(er)  - one %*% er - sum(lfactorial(wi)) # log(\prod_k Poisson(n_{ks}))
  return(log.Pi.si)
}

# ---------------------
# Calculates log(P_i(s)) for one occasion for all individuals for all s in mesh
# ch1occ is n by K matrix of counts
# returns an n by M matrix
count.log.Pi.s = function(capthist1,er) {
  log.Pi.s=apply(capthist1, 1, count.log.Pi.si,er=er)
  return(t(log.Pi.s))
}

# ---------------------
# Calculates log(P_it(s)) for all occasions for all individuals for all s in mesh
# capthist is n by nt by K matrix of counts
# returns an n by M matrix
count.log.Pit.s <- function(capthist, er) {
  n=dim(capthist)[1]
  nt=dim(capthist)[2]
  M=dim(er)[2]
  log.Pit.s=array(dim=c(n,nt,M))
  for(i in 1:nt) { # calculate log(P_i(s)) for each occasion:
    log.Pit.s[,i,]=count.log.Pi.s(capthist[,i,],er)
  }
  return(log.Pit.s)
}

# ---------------------
calc.P.occ.mat <- function(capthist, Pit.s, anim, occ, nmesh) {
  P <- matrix(nrow=nmesh, ncol=2)
  P[,1] <- Pit.s[anim, occ,]
  if (all(capthist[anim, occ,] == 0)) {
    P[,2] <- 1
  } else {
    P[,2] <- 0
  }
  return(P)
}

# ---------------------
calc.P.mat <- function(capthist, Pit.s, anim, occ, nmesh) {
  P <- matrix(nrow=nmesh, ncol=2)
  P[,1] <- Pit.s[anim, occ,]
  if (all(capthist[anim,,] == 0)) {
    P[,2] <- 1
  } else {
    P[,2] <- 0
  }
  return(P)
}

# ---------------------
hmm.negloglike <- function (capthist, mesh, pars, dist) {
  n <- dim(capthist)[1]
  nt <- dim(capthist)[2]
  M <- dim(mesh)[1]
  a <- rep(attributes(mesh)$area, M) # area of each mesh cell
  D <- rep(exp(pars["D"]), M) # Density at each mesh point
  phi <- plogis(pars["phi"])
  
  Pi <- matrix(c(1,0), nrow=M, ncol=2, byrow = T)
  tpm <- matrix(c(phi, 0, 1-phi, 1) ,nrow=2, ncol=2)
  lscale <- 0
  
  # Calculate log(P_i(s)) at each mesh point
  Pit.s <- calc.p.Pit.s(capthist, mesh, pars, dist = NULL)
  
  for (i in 1:n) { # loop over individuals 
    for (j in 1:nt) { # loop over occasions
      P <- calc.P.mat(capthist, Pit.s, i, j, M)
      alpha <- (Pi * P) %*% tpm
      lscale <- lscale + log(sum(alpha))
      alpha <- alpha/sum(alpha)
    }
  }
  return(loglik = -lscale)
}

# ---------------------
# Example
D <- scr.pop$D
pars <- c(log(g0), log(sigma), log(D), qlogis(phi)) # starting values
names(pars) <- c("g0","sigma","D", "phi")
dist <- distances(traps(pop.capthist),mesh) # calculate trap distances 

est <- optim(pars, hmm.negloglike, capthist=pop.capthist, mesh=mesh, dist=dist, control=list(trace=5))
plogis(est$par["phi"])
# true phi = 0.45, estimate phi = 0.788
# scaling not carried out correctly?
# calc.P.occ.mat or calc.P.mat?
