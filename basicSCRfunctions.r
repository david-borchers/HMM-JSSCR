#' @title SCR negative log-likelihood
#'
#' @description  Calculates the log-likelihood for an SCR model. Currently only multi-catch traps and 
#' count detectors implemented, and no detection function or density model covariates implemented.
#' 
#' @return Value of negative log-likelihood
#' 
#' @param pars Vector of paramters, named "g0", "sigma" and "D".
#' @param capthist Capture history object, or library \code{secr} class 'capthist.
#' @param mesh Integration mesh, of library \code{secr} class 'mask': an Mx2 matrix.
#' @param dist KxM matrix of distances between the K detectors in traps(capthist) and the M mesh points.
#' 
#' @export
scr.negloglik=function(pars, capthist, mesh, dist=NULL) {

  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.Pis(capthist, mesh, pars, dist=dist)

# E[n(s)] = \int a(s) * D(s) * p..(s) ds:
  En.s = a * D * scr$p..s # E[n(s)] = a(s) * D(s) * p..(s)
  En = sum(En.s)
  
  # calculate a(s) * D(s) * P_i(s) at all mesh points
  aDPi.s=t(a * D * t(scr$Pi.s)) # transposing to multiply by row, not colum
  Pch.i=apply(aDPi.s,1,sum) # Prob of observed capture histories (vector of length n)
  
  # log-likelihood
  n = length(Pch.i)
  l = sum(log(Pch.i)) - En - lfactorial(n)
  
  return(loglik = -l)
}


#' @title Calculates p..(s) and P_i(s) for SCR likelihood calculation.
#'
#' @description  Calculates p..(s) (a vector of length n) and P_i(s) an nxM matrix of probabilities for 
#' SCR likelihood calculation.
#'
#' @param pars Vector of paramters, named "g0", "sigma" and "D".
#' @param capthist Capture history object, or library \code{secr} class 'capthist.
#' @param mesh Integration mesh, of library \code{secr} class 'mask': an Mx2 matrix.
#' @param dist KxM matrix of distances between the K detectors in traps(capthist) and the M mesh points.
#' 
#' @return A list with these elements: 
#' \itemize{
#'  \item{p..s}{Probability of detection by any detector, for activity centres at each mesh point 
#'  (vecotr of length M)}
#'  \item{Pi.s}{For each of n detected individuals, the rows of this matrix contain the probability of 
#'  obtaining the observed capture history, for all M possible activity centre locations (mesh points) 
#'  (nxM matrix).}
#'  }
#' 
#' @export
calc.p.Pis=function(capthist, mesh, pars, dist) {
  
  # calculate distances from each detector to each mesh point
  if(is.null(dist)) dist=distances(traps(capthist),mesh)
  
  # Proximity detector with count data
  if (detector(traps(capthist))=="count") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    g0=exp(pars["g0"])
    sigma=exp(pars["sigma"])
    er <- g0 * exp(-dist^2 / 2 / sigma^2) # Poisson encounter rate at all M distances from each detector
    
    # probability of being caught at least once if at mesh vertex s of M
    p.ts <- 1 - exp(- apply(er,2,sum)) # over a single occasion
    p..s <- 1 - exp(- apply(er * dim(capthist)[2],2,sum)) # over all occasions
    
    # calculate log(P_i(s)): nxM matrix
    #comb.capthist=apply(capthist,c(1,3),sum)
    #log.Pi.s = counts1.log.Pi.s(comb.capthist,er)
    log.Pit.s=counts.log.Pi.s(capthist,er) # log(P_{it}(s))
    log.Pi.s = apply(log.Pit.s,c(1,3),sum) # P_i(s) = \sum_t log(P_{it}(s)) over occasions
  }
  # Multi-catch trap
  if (detector(traps(capthist))=="multi") {
    # Probabilities for each detector and each mesh point (assumed same at all times)
    g0=plogis(pars["g0"])
    sigma=exp(pars["sigma"])
    gtk <- g0 * exp(-dist^2 / 2 / sigma^2)
    log.gtk <- log(gtk) # for use below
    log.gtk1 <- log(1-gtk) # for use below
    
    # probability of being caught at least once if at mesh vertex s of M
    p..s <- 1 - apply(1-gtk, 2, prod) ^ dim(capthist)[2] 
    
    # calculate log(P_i(s)): nxM matrix
    log.Pi.s=multi.log.Pi.s(capthist,log.gtk,log.gtk1)
  }
  
  # put all this crap in a list to pass
  return(list(p..s=p..s, Pi.s=exp(log.Pi.s)))
}



# Calucates distances between two sets of coordinates
# X is a by 2 matrix, Y is b by 2 matrix, output is a by b matrix
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

# MULTI-CATCH detectors (detector(traps(capthist))=="multi")
# ---------------------
# Calculates log(P_i(s)) for a single individual for all s in mesh
# wi is ith row of multi-catch trap capture history (ith individual's capture history)
# of length nt, and with detector number in non-zero elements
multi.log.Pi.si = function(wi,log.gtk,log.gtk1) {
  delta=rep(0,dim(log.gtk)[1])
  delta[wi[wi>0]]=1
  log.Pi.si <- delta %*% log.gtk  + (1-delta) %*% log.gtk1 # log(\prod_k g_itk(s))
  return(log.Pi.si)
}

# Calculates log(P_i(s)) for all individuals for all s in mesh
# capthist is n by nt multi-catch trap capture history matrix
# log.gtk and log.gtk1 are K by M matrices
# returns an n by M matrix
multi.log.Pi.s = function(capthist1,log.gtk,log.gtk1) {
  log.Pi.s=apply(capthist1, 1, multi.log.Pi.si,log.gtk=log.gtk,log.gtk1=log.gtk1)
  return(t(log.Pi.s))
}


# COUNT detectors (detector(traps(capthist))=="count")
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

# Calculates log(P_i(s)) for one occasion for all individuals for all s in mesh
# ch1occ is n by K matrix of counts
# returns an n by M matrix
count.log.Pi.s = function(capthist1,er) {
  log.Pi.s=apply(capthist1, 1, count.log.Pi.si,er=er)
  return(t(log.Pi.s))
}

# Calculates log(P_i(s)) for all occasions for all individuals for all s in mesh
# capthist is n by nt by K matrix of counts
# returns an n by M matrix
counts.log.Pi.s = function(capthist,er) {
  n=dim(capthist)[1]
  nt=dim(capthist)[2]
  M=dim(er)[2]
  log.Pit.s=array(dim=c(n,nt,M))
  for(i in 1:nt) { # calculate log(P_i(s)) for each occasion:
    log.Pit.s[,i,]=count.log.Pi.s(capthist[,i,],er)
  }
  return(log.Pit.s)
}



# ================== Some plotting functions =====================

plot.p..=function(capthist, mesh, pars, dist=NULL, contour=TRUE, key=TRUE, ...) {
  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.logPis(capthist, mesh, pars, dist=dist)
  
  covariates(mesh)$p..s = p..s
  plotcovariate(mesh,covariate="p..s",main="p..s",contour=contour, key=key, ...)
}




plot.Pi=function(i,capthist, mesh, pars, dist=NULL, contour=TRUE, key=TRUE, ...) {
  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  dets=traps(capthist) # detectors
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.logPis(capthist, mesh, pars, dist=dist)
  
  covariates(mesh)$Pi.s = scr$Pi.s[i,]/sum(scr$Pi.s[i,])
  plotcovariate(mesh,covariate="Pi.s",main=expression(f(s[i]*"|"*capthist[i])),contour=contour, key=key, ...)
  plot(dets,add=TRUE,detpar=list(pch=19,col="black"))
  if(dim(ch)[2]==1) freq = apply(apply(ch,c(2,3),sum),2,sum)
  else freq=apply(capthist[i,,],2,sum)
  detected=which(freq>0)
  text(dets$x[detected],dets$y[detected],labels=freq[detected],cex=0.75,col="white")
}




#' @title Plots image (and optionally contours) of mask covariate value
#' @param mask is an object of class `mask'
#' @param covariate is a character variable with the name of one of the covariates in mask (the one to plot)
#' @param contour is a logical, TRUE if want contour plots on image
#' @param ... other arguments to be passed to \code{prep4image}
plotcovariate=function(mask,covariate,contour=TRUE,key=TRUE, ...) {
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  dat=data.frame(x=mask$x,y=mask$y,z=covariates(mask)[[cnum]])
  prep4image(dat,contour=contour,key=key,...)
}

#' @title Prepares data frame for plotting with image/contour/persp.
#'   
#' @description From an input data frame with columns x, y and z, this function 
#'   creates a list with elements x, y and z in a format suitable for passing to
#'   functions \code{\link{image}}, \code{\link{contour}} or 
#'   \code{\link{persp}}. The coordinates in \code{data} are assumed to come
#'   from a 2D grid of points.
#'   
#' @param data a data frame with columns x and y being Cartesian coordinates, 
#'   and z being the values of some variable at each coordinate.
#' @param plot if \code{TRUE} then an image plot will be drawn using 
#'   \code{\link{image.plot}}
#' @param contour if \code{TRUE} then contours will be added (only used when 
#'   \code{plot=TRUE})
#' @param key logical for whether or not to include key when \code{plot = TRUE} (\code{\link{image.plot}} is used when \code{key = TRUE}, \code{\link{image}} is used when \code{key = FALSE})
#' @param ... other arguments to pass to \code{\link{image}} or \code{\link{image.plot}} (only used 
#'   when \code{plot=TRUE})
#'   
#' @details Sorts z on values of x first, then y, then creates a matrix of 
#'   z-values from this. Returns a list with elements x (unique values of x, in 
#'   increasing order), y (unique values of y, in increasing order) and z 
#'   (matrix of z-values in appropriate order for image/contour/persp). 
#'   
#'   If the original z is a factor variabele, the z returned is a matrix of integers 
#'   between 1 and length(levels(z)) and the output list has an attributes called 
#'   ``facnames'' that is a character vector containing the levels as factor 
#'   variables, with z=1 corresponding to the first name, z=2 to the second, etc.
#' @export
#' @importFrom fields image.plot
prep4image = function(data, plot = TRUE, contour = TRUE, key = TRUE, ...){
  
  # convert factor data$z to integer:
  zfactor=FALSE
  if(is.factor(data$z)) {
    zfactor=TRUE
    fac=data$z
    facnames=levels(fac)
    nlevels=length(facnames)
    data$z=rep(0,length(fac))
    got=rep(FALSE,nlevels)
    for(i in 1:nlevels){
      j=which(fac==facnames[i])
      if(length(j)>0) got[i]=TRUE
      data$z[j]=i
    }
    facnames=facnames[got] # remove factor names not in mask
  }
  data = as.matrix(data)
  
  x = sort(unique(data[,"x"]))
  y = sort(unique(data[,"y"]))
  
  z = matrix(NA, nr = length(x), nc = length(y))
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      m = which(data[,"x"] == x[i] & data[,"y"] == y[j]) ; m
      z[i,j] = if(length(m) == 0) NA else data[,"z"][m]
    }
  }
  
  if(plot){
    if(key){
      image.plot(x, y, z, ...)
    }else{
      image(x, y, z, ...)
    }
    if(contour) contour(x, y, z, add = TRUE)
  }
  
  outlist=list(x = x, y = y, z = z)
  if(zfactor) attributes(outlist)$facnames=facnames
  
  invisible(outlist)
  
}


#' @title add new covariate consisting of distances to points that have covariate==distance.to
#' @param mask is object of class "mask".
#' @param covariate is name of one of the elements of covariates(mask)
#' @param distance.to is the target value of covariate; distances to mask points with this value are calculated.
#' @param dname is the name you want to give to the new covariate
add.dist.to=function(mask,covariate,distance.to,dname=NULL,overwrite=FALSE){
  covs=names(covariates(mask))
  if(is.null(covs)) stop("No covariates in mask. Can't add distance to anything")
  ncov=length(covs)
  if(is.null(dname)) dname=paste("dist2",covariate,"=",distance.to,sep="")
  if(is.element(covariate,covs)) {
    if(dname %in% names(covariates(mask))) {
      if(overwrite) {
        warning(paste("Covariate called",dname,"already existed. It has been overwritten."))
        covi=which(names(covariates(mask))==dname)
      } else stop(paste("Covariate called",covariate,"already exists. Use overwrite=TRUE if you want to overwrite it."))
    } else {
      covi=ncov+1
    }
    cov=covariates(mask)[which(covs==covariate)]
    if(is.null(cov)) stop(paste("No covariate called",covariate,"in mask."))
    targets=which(cov==distance.to)
    distances=rdist(mask,mask[targets,])
    d2=apply(distances,1,min)
    covariates(mask)[[covi]]=d2
    names(covariates(mask))[[covi]]=dname
    return(mask)
  } else stop(paste("No covariate called",covariate,"in mask."))
}

