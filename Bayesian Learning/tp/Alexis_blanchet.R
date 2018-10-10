#install.packages("MASS")
#install.packages("abind")
#install.packages("mnormt")
#install.packages("LaplacesDemon")
#install.packages("coda")
library(MASS)
library(abind)
library(mnormt)
library(LaplacesDemon)
library(coda)
source('GM_functions.R')

########################################
########################################
d <- 2 ## dimension
K <- 3 ## number of mixture components
N <- 500 ## sample size
p <- c(3/10, 2/10, 5/10) ## weight parameter
Mu <- rbind(c(-1, -1), c(0, 1), c(1, -1)) ## centers
set.seed(1)
NN <- rmultinom(n = 1, size = N, prob = p)
## NN: number of points in each cluster
Sigma <- array(dim = c(d, d, K)) ## the covariance matrices
X <- matrix(ncol = d, nrow = 0) ## the dataset
for(j in 1:K){
  Sigma[, , j] <- rwishart(nu = 5, S = 0.05*diag(d))
}
for(j in 1:K){
  X <- rbind(X, mvrnorm(n=NN[j], mu = Mu[j, ], Sigma=Sigma[, , j]))
}
labs <- rep(0, N) ##vector of labels
count=1
for(j in 1:K){
  labs[count:(count+NN[j]-1)] <- j
  count <- count + NN[j]
}
#' Plot the labelled data
plot(X[,1], X[,2], col=labs)


Kfit <- 9## try with Kfit= 2,3,4,5 ...
outputem <- emalgo(x=X,k=Kfit, tol=1e-6)
## inspect the objective function (stopping criterion)
length(outputem$objective)
## Plot the (labelled) data
plot(X[,1], X[,2], col = labs, pch = 19)
## Add the starting points (from kmeans) to the plot
Init <- initem(X,Kfit)
points(Init$Mu[,1],Init$Mu[,2], col="orange",
       pch=18,cex=10*Init$p)
## add the centers from EM
points(outputem$last$Mu[,1],outputem$last$Mu[,2], col="blue",
       pch=18,cex=10*outputem$last$p)
## add the true centers
points(Mu[,1],Mu[,2], col="black",pch=8,cex=10*p)
for(j in 1:Kfit){
  ellips <- draw_sd(outputem$last$Mu[j,], outputem$last$Sigma[,,j])
  lines(ellips[1,], ellips[2,], col='blue')
}
for(j in 1:K){
  ellips <- draw_sd(Mu[j,], Sigma[,,j])
  lines(ellips[1,], ellips[2,], col='black')
}

plot(outputem$objective)

#########################'
#' 3. VB
#########################' 
#' Bayesian model: 
#' p ~ dirichlet(alpha);  alpha = (alpha0, ... , alpha0)
#' [ xi | p ] ~ Multinomial(p)
#' [ mu_j | Lambda_j ] ~ Normal(m0, beta0 Lambda_j^(-1))
#' Lambda_j ~ Wishart(W0, nu0)
#' [ X| xi=j, mu, Lambda ] ~ Normal (mu_j, Lambda_j^(-1))

#' hyper-parameters : to be varied
Kfit <- 5 ## try with Kfit= 2,3,4,5 ... 
alpha0 <- 0.1
m0 <- rep(0,2)
beta0 <- 0.1
W0 <- 1*diag(2)
nu0 <- 10
#' Run VB 
#'
seed <- 10
set.seed(seed)
outputvb <- vbalgo(x=X,k=Kfit, alpha0 = alpha0, W0inv = solve(W0),
                   nu0 = nu0, m0 = m0, beta0=beta0, tol=1e-6)

#' plot the lowerbound over iterations 
plot(outputvb$lowerbound)
#' show a summary of VB's output
T <- ncol(outputvb$alphamat)
outputvb$alphamat[,T]
outputvb$Marray[,,T]
dim(outputvb$Marray)
dim(outputvb$Winvarray)
print(T)
#' Visual summary of VB's output :
#' posterior expectancy of each parameter
p_vb <- outputvb$alphamat[,T]/sum(outputvb$alphamat[,T])## complete the code
print(p_vb)
## (variational posterior expectancy of mixture weights)
Mu_vb <-outputvb$Marray[,,T]## complete the code
## (variational posterior expectancy of mixture centers)
Sigma_vb <- array(dim=c(d,d,Kfit))
for(j in 1:Kfit)
{
  Sigma_vb[,,j] <- outputvb$Winvarray[,,j,T]/outputvb$Numat[j,T]## complete the code
  ## (variational posterior expectancy of mixture covariances)
}

## show the data, true centers and initial positions from K-means
graphics.off()
plot(X[,1], X[,2], col=labs)
points(Mu[,1],Mu[,2], col="black",pch=8,cex=10*p) 
set.seed(seed)
Init <-  initem(X,Kfit)
points(Init$Mu[,1],Init$Mu[,2], col="orange",pch=18,cex = 10*Init$p)
## Add a  summary of the VB solution
nonneg <- which(p_vb>0.001)
print(nonneg)
for(j in nonneg){
  points(Mu_vb[j,1], Mu_vb[j,2], col="blue",
         pch=18,cex= 10 * p_vb[j])
  ellips <- draw_sd(mu = Mu_vb[j,], 
                    sigma = Sigma_vb[,,j])
  lines(ellips[1,], ellips[2,], col='blue')
}
####################################################'
####' Metropolis-Hastings
####################################################'
#' Basic testing for the MH sampler
Kmc <- Kfit ## try with different values
init <- initem(x=X, k=Kmc)
hpar <- list( alpha0=rep(alpha0, Kmc),
              m0 = rep(0, d), beta0 = beta0, 
              W0 = W0, nu0 = nu0)

ppar <- list(var_Mu = 0.001,
             nu_Sigma = 500,
             alpha_p = 500) 
nsample = 5000

set.seed(1)
pct <- proc.time()
outputmh <- MHsample(x=X, k=Kmc, nsample= nsample,
                     init=init, hpar=hpar, ppar=ppar)
newpct <- proc.time()
elapsed <- newpct - pct
elapsed
outputmh$naccept ## should not be ridiculously low. 
length(outputmh$lpostdens)

outputmh_1 <- mcmc(cdfTrace(x=c(-1,1), outputmh))
outputmh_2 <- mcmc(cdfTrace(x=c(-1,-1), outputmh))
outputmh_3 <- mcmc(cdfTrace(x=c(-1,0), outputmh))
heidel.diag(outputmh_3)
plot(outputmh_1)

gelman.diag(mcmc.list(outputmh_1,outputmh_2,outputmh_3))
gelman.plot(mcmc.list(outputmh_1,outputmh_2,outputmh_3))




#' Predictive density
xx <- seq(-2,2,length.out=20)
yy <- xx
dtrue <- outer(X= xx, Y=yy,
               FUN = function(x,y){
                 wrapper(x=x, y=y,
                         FUN=function(u,v){
                           exp(gmllk(x = c(u,v), Mu = Mu,
                                     Sigma = Sigma, p = p))
                         })
               })

dpredmh <-  outer(X= xx, Y=yy,FUN = function(x,y){wrapper(x = x, y = y,FUN =function(u,v){MHpredictive(c(u,v),outputmh)})})

breaks <- c(seq(0.01,0.09, length.out=5),seq(0.1,0.3,length.out=5))
nbreaks <- length(breaks)
contour(xx,yy, z = dtrue, nlevels=nbreaks, levels = breaks)
contour(xx,yy, z = dpredmh,  nlevels=nbreaks, levels = breaks,
                            add=TRUE, col='red')




                    
 
I <- dim(outputem$Muarray)[3]
I
dim(outputem$pmat)

Pexcess <- rep(0,10)
Pexcess_em <- Pexcess; Pexcess_vb <- Pexcess; Pexcess_mh <- Pexcess
thres_vect <- seq(-3, 3, length.out=30)
for(i in seq_along(thres_vect))
  {
  
  threshold <- rep(thres_vect[i], 2)
  Pexcess[i] <- 1 - gmcdf(x = threshold, Mu = Mu, Sigma=Sigma, p=p)
  Pexcess_em[i] <- 1-gmcdf(threshold,outputem$Muarray[,,I],outputem$Sigmaarray[,,,I],outputem$pmat)## complete the code:
    ##maximum likelihood estimator using EM output
    
  Pexcess_vb[i] <- 1-vbPredictiveCdf(threshold,outputvb$alphamat[,T],outputvb$Betamat[,T],outputvb$Marray[,,T],outputvb$Winvarray[,,,T],outputvb$Numat[,T])## complete the code:
    ## posterior predictive estimator using VB output:
    ## use vbPredictiveCdf
    
  Pexcess_mh[i] <- MHpredictiveCdf(threshold,outputmh)## complete the code:
    ## posterior predictive estimator using MH output:
    ## use MHpredictiveCdf.
    
  }
ylim <- range(Pexcess, Pexcess_em,Pexcess_vb)
plot(thres_vect,Pexcess, ylim = ylim)
lines(thres_vect, Pexcess_vb, col='red')
lines(thres_vect, Pexcess_em, col='blue')
lines(thres_vect, Pexcess_mh, col='green')



Pexcess <- rep(0,10)
Pexcess_em <- Pexcess; Pexcess_vb <- Pexcess; Pexcess_mh <- Pexcess
thres_vect <- seq(1, 5, length.out=30)
for(i in seq_along(thres_vect))
{
  
  threshold <- rep(thres_vect[i], 2)
  Pexcess[i] <- 1 - gmcdf(x = threshold, Mu = Mu, Sigma=Sigma, p=p)
  Pexcess_em[i] <- 1-gmcdf(threshold,outputem$Muarray[,,I],outputem$Sigmaarray[,,,I],outputem$pmat)## complete the code:
  ##maximum likelihood estimator using EM output
  
  Pexcess_vb[i] <- 1-vbPredictiveCdf(threshold,outputvb$alphamat[,T],outputvb$Betamat[,T],outputvb$Marray[,,T],outputvb$Winvarray[,,,T],outputvb$Numat[,T])## complete the code:
  ## posterior predictive estimator using VB output:
  ## use vbPredictiveCdf
  
  Pexcess_mh[i] <- MHpredictiveCdf(threshold,outputmh)## complete the code:
  ## posterior predictive estimator using MH output:
  ## use MHpredictiveCdf.
  
}
ylim <- range(Pexcess, Pexcess_em,Pexcess_vb)
plot(thres_vect,Pexcess, ylim = ylim)
lines(thres_vect, Pexcess_vb, col='red')
lines(thres_vect, Pexcess_em, col='blue')
lines(thres_vect, Pexcess_mh, col='green')
