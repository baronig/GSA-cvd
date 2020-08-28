# ==============================================
# Example script for running the global sensitivity analysis
# based on the combined variance and distribution based strategy
# (CVD): Baroni and Francke, under review
# 
# This script also compares a standard approach that
# combines the main effect and the distribution-based measure
# PAWN
# 
# 16 July 2020
# 
# Gabriele Baroni
# g.baroni@unibo.it
# University of Bologna (Italy)
# 
# Till Francke
# francke@uni-potsdam.de
# University of Potsdam (Germany)
# 
# ==============================================
# The script is divided in 5 sections
# 1. define the model, parameters and factors distribution
# 2. sampling
# 3. run the model
# 4. look at the input-output space
# 5. estimate and visualize the indices

# 0. initialisation ==============================================
rm(list = ls())

# install packages, if missing, load them
install_load = function(pckg_name)
{
  if (!require(pckg_name, character.only = TRUE))
  {
    install.packages(pckg_name)
    library(pckg_name, character.only = TRUE)
  }
} 

install_load(pckg_name = "randtoolbox")
install_load(pckg_name = "ks")
install_load(pckg_name = "plyr")
install_load(pckg_name = "lattice")
install_load(pckg_name = "Matrix")

source("Sobol_sampling.R")
source("CVD_function.R")
source("SCA_function.R")

# 1. chose functions and settings ==============================================
# Ishigami-Homma function (3 factors)

ai  = c(2, 1)          # fixed parameters
n_k = 3                # number of factors
lb  = c(-pi, -pi, -pi) # lower boundary
ub  = c(pi, pi, pi)    # upper boundary

e_sam = 14      # exponent of the number of samples e.g., = 14 
N     = 2^e_sam # number of samples N = as power of 2 as in Sobol seq
n_m   = 10      # number of slides m

Ishigami_Homma_3 = function(x)
{  
  sin(x[,1]) + ai[1] * sin(x[,2])^2 + ai[2] * x[,3]^4 * sin(x[,1])
}

model_f = Ishigami_Homma_3 #choose function
  
# Initialize vector and matrix

vals = c() #specifications for Sobol sampling
for (k in 1:n_k){
  # k=1
  temp = list(list(var=paste("x", k, sep=""), dist="unif", params = list(min = lb[k], max = ub[k])))
  vals = c(vals, temp) 
}

n_seed   = 1548674 
n_scramb = 1 # type of scrambling

# 2. Sobol sampling ==============================================

X_sam = makeMCSample(N, vals, n_scramb, n_seed)


# 3 run the model ==============================================

Y_sam = model_f(X_sam)

# 4 visualize input-output space ==============================================
layout(matrix(c(1:3), nr=1, byrow=F))
par(oma=c(4, 3, 1, 1))
par(mar=c(0.5, 2.5, 1, 0))
par(mgp=c(0.5, 0.8, 0))

col_sel = c("#386cb0", "#f0027f", "#7fc97f")

for (k in 1:n_k){
  if(k==1) yaxt=NULL else yaxt="n"
  plot(X_sam[, k], Y_sam, yaxt=yaxt, ylab="", xlab="", xlim=range(X_sam[, k]), ylim=range(Y_sam), col=col_sel[1], pch=1)
    if (k==1) mtext(expression(italic("Y")), side=2, line = 2)
  mtext(bquote(italic("x")[.(k)]), side =1, line = 2.5)
  abline(h=0, lty=1, col=1)
}

# 5: sensitivity analysis ==============================================

CVD_index = CVD(X_sam, Y_sam, n_m)
SCA_index = SCA(X_sam, Y_sam, n_m)

# plot indices
par(mfrow=c(1, 2))
par(oma=c(4, 3, 1, 1))
par(mar=c(0.5, 2, 1, 2))
par(mgp=c(0.5, 0.8, 0))

plot(SCA_index[, 1], SCA_index[, 2], type="p", xlim=c(0, 0.8), ylim=c(0, 0.8), pch=19, col=col_sel, ylab="", cex = 1.5)
abline(v=0.1, col="gray", lty=2)
abline(h=0.1, col="gray", lty=2)
mtext(expression("Main effect "*italic("S"[i])), 1, line=2.5)
mtext(expression(italic("PAWN"[i])), 2, line=2)
legend(0.5, 0.8, expression(italic("x"[1]), italic("x"[2]), italic("x"[3])), col=col_sel, pch=19, cex = 1)
mtext("Standard methods", 3, line=0.5)

plot(CVD_index[, 1], CVD_index[, 2], type="p", xlim=c(0, 0.8), ylim=c(0, 0.8), pch=19, col=col_sel, ylab="", cex = 1.5)
abline(v=0.1, col="gray", lty=2)
abline(h=0.1, col="gray", lty=2)
mtext(expression("Main effect "*italic("S"[i])), 1, line=2.5)
mtext(expression("Interaction "*italic("I"[i])), 2, line=2)
legend(0.5, 0.8, expression(italic("x"[1]), italic("x"[2]), italic("x"[3])), col=col_sel, pch=19, cex = 1)
mtext("CVD strategy", 3, line=0.5)
