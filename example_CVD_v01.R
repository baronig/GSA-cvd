# ==============================================
# Example script for implementing the global sensitivity analysis based on
# the combined variance and distribution based strategy
# 
# 12 February 2020
# 
# Gabriele Baroni
# g.baroni@unibo.it
# University of Bologna (Italy)

# Till Francke
# francke@uni-potsdam.de
# university of Potsdam (Germany)
# 
# ==============================================
# The script is divided in 5 sections
# 1. define the model, parameters and factors distribution
# 2. sampling
# 3. run the model
# 4. look at the input-output space
# 5. estimate and visualize the indices
# ==============================================

# install.packages("randtoolbox")
# install.packages("ks")

rm(list = ls())

library(ks)
source("Sobol_sampling.R")
source("CVD_function.R")

# ==============================================
# ==============================================
# 1. chose functions and settings
# Ishigami-Homma (3 factors)

ai  = c(2,1)         # parameters
n.k = 3              # number of factors
lb  = c(-pi,-pi,-pi) # lower boundary
ub  = c(pi,pi,pi)    # upper boundary

e.sam = 14      # exponental of the number of samples e.g., = 14 
N     = 2^e.sam # number of samples N = as power of 2 as in Sobol seq
n.m   = 10      # number of slides m

# ==============================================
# Initialize vector and matrix

vals = c() 
for (k in 1:n.k){
  # k=1
  temp <- list(list(var=paste("x",k,sep=""), dist="unif", params = list(min = lb[k], max = ub[k])))
  vals = c(vals, temp) 
}

n.seed   = 1548674 
n.scramb = 1 # type of scrambling

# ==============================================
# ==============================================
# 2. Sobol sampling

X.sam <- makeMCSample(N,vals,n.scramb,n.seed)

# ==============================================
# ==============================================
# 3 run the model

Y.sam = sin(X.sam[, 1]) + ai[1] * sin(X.sam[, 2])^2 + ai[2] * X.sam[, 3]^4 * sin(X.sam[, 1])

# ==============================================
# ==============================================
# 4 look at input-output space
layout(matrix(c(1:3), nr=1, byrow=F))
par(oma=c(4,3,1,1))
par(mar=c(0.5,2.5,1,0))
par(mgp=c(0.5,0.8,0))

col.sel = c("#386cb0","#f0027f","#7fc97f")

for (k in 1:n.k){
  if (k==1){
    plot(X.sam[,k],Y.sam,ylab="",xlab="", xlim=range(X.sam[,k]),ylim=range(Y.sam),col=col.sel[1],pch=1)
    mtext("Y",side=2,line = 2)
    mtext(expression("x"[1]),side =1,line = 2.5)
    abline(h=0,lty=1,col=1)
  } else {
    plot(X.sam[,k],Y.sam,yaxt="n",ylab="",xlab="",xlim=range(X.sam[,k]),ylim=range(Y.sam),col=col.sel[1],pch=1)
    abline(h=0,lty=1,col=1)
    mtext(bquote(x[.(k)]),side =1,line = 2.5)
  }
}

# ==============================================
# ==============================================
# 5: sensitivity analysis

index = CVD(X.sam,Y.sam,n.m)

# plot indices
par(mfrow=c(1,1))
plot(index[,1],index[,2],type="p",xlim=c(0,0.8),ylim=c(0,0.8),pch=19,col=col.sel,ylab="",cex = 1.5)
abline(v=0.1,col="gray",lty=2)
abline(h=0.1,col="gray",lty=2)
mtext(expression("Main effect S "[i]),1,line=2.5)
mtext(expression("Interaction I "[i]),2,line=2.5)
legend(0.6,0.8,expression(x[1],x[2],x[3]),col=col.sel,pch=19,cex = 1.5)
