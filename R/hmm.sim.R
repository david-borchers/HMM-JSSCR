#                         --- HMM SIMULATION ---

# load parallel library
library(parallel)

# Calculate the number of cores
no.cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no.cores)
clusterEvalQ(cl, c(library(secr), library(fields)))
clusterExport(cl, c('sim.open.scr.pop',
                    'lambda0', 'sigma', 'occasions', 'phi'))


# create tables
n <- 100
seeds <- sample(1:1000000, n, replace=T)
store.capthists <- list(length = n)
store.capthists <- parLapply(cl, seeds, function(x) 
  sim.open.scr.pop(lambda0, sigma, occasions, phi, seed=x))

stopCluster(cl)

#---

# Initiate cluster
cl <- makeCluster(no.cores)
clusterEvalQ(cl, c(library(secr), library(fields)))
clusterExport(cl, c('calc.P.mat', 'calc.p.Pit.s', 'count.log.Pi.s', 
                    'count.log.Pi.si', 'count.log.Pit.s', 'distances', 
                    'hmm.negloglike', 'dist', 'pars', 'mesh'))

results <- list(length = n) 
results <- parLapply(cl, store.capthists, function(x)
  optim(pars, hmm.negloglike, capthist=x$capthist, mesh=x$mesh, dist=dist, 
        hessian = T, control=list(trace=1)))


stopCluster(cl)

#---