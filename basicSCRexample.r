library(secr)

# set detector type
# detype="multi"
detype="count"

# create the detectors
# --------------------
spacing=20
dets = make.grid(nx=7, ny=7, spacing=spacing,detector=detype)
plot(dets,border=0)
# set normal detection function parameters
g0=2.5; sigma=spacing
# create the mesh
mesh = make.mask(dets, buffer=4*sigma, nx=64, ny=64, type="trapbuffer")

# set density via abundance on mesh
# ---------------------------------
Nmesh = 100
M=dim(mesh)[1]
D = Nmesh / (M * attributes(mesh)$area); D

# simulate a population from this intensity surface
# -------------------------------------------------
pop=sim.popn(D=D, core=mesh, buffer=0, seed=12345)
# plot mesh with individuals' locations and detectors overlaid
plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
dim(pop)[1] # check simulated population size

# Generate capture histories
# --------------------------
# single-occasion capture history
capthist1=sim.capthist(dets,popn=pop, detectpar=list(g0=g0,sigma=sigma), detectfn="HHN", noccasions=1, nsessions=1,seed=12345)
# multi-occasion capture history
nt=5
capthist=sim.capthist(dets,popn=pop, detectpar=list(g0=g0,sigma=sigma), detectfn="HHN", noccasions=nt, nsessions=1,seed=12345)
summary(capthist)
n=dim(capthist)[1]
plot(mesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets,add=TRUE)
plot(capthist, border=sigma, tracks=TRUE, gridlines=FALSE,rad=3,add=TRUE)

# Find the MLEs
# =============
pars=c(log(g0), log(sigma), log(D))
names(pars) = c("g0","sigma","D") # need this ` cause scr.negloglik uses these names
dist=distances(traps(capthist),mesh) # calculated distances only once, then pass
# multiple-occasions capture history:
est=optim(pars,scr.negloglik, capthist=capthist, mesh=mesh, dist=dist, control=list(trace=5))
screst=secr.fit(capthist,mask=mesh,detectfn="HHN",model=list(lambda0~1)) # do secr fit for comparison
screst;exp(est$par);exp(pars["D"]) # compare our estimates to those from secr (and true density)

# single-occasion capture history:
est1=optim(pars,scr.negloglik, capthist=capthist1, mesh=mesh, dist=dist, control=list(trace=5))
screst1=secr.fit(capthist1,mask=mesh,detectfn="HHN",model=list(lambda0~1)) # do secr fit for comparison
screst1;exp(est1$par);exp(pars["D"]) # compare our estimates to those from secr (and true density)

# Plot some stuff to see if it makes sense:
plot.p..(capthist,mesh,est$par,dist)
plot(dets,add=TRUE,detpar=list(col="white")) 

i=1
# execute the line below repeatedly to plot estimated location of each detection
plot.Pi(i,capthist,mesh,est$par,dist,contour=FALSE);i=i+1


