# ===================================================
# ===================================================
# Global sensitivity analysis based on a standard 
# combined use of variance and distribution-based methods.
# Specifically, the main effect is estimated based on the DLR
# method presented by Kucherenko and Song (2017) and the 
# distribution based index based on PAWN index by 
# Pianosi and Wagener (2018) 
# The function estimate main and interaction index
# based on a generic sample N and interval m.
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
# References
# Pianosi, Francesca, and Thorsten Wagener. "Distribution-Based Sensitivity Analysis from a Generic Input-Output Sample." Environmental Modelling & Software 108 (October 1, 2018): 197-207. https://doi.org/10.1016/j.envsoft.2018.07.019.
# Kucherenko, S., and S. Song. "Different Numerical Estimators for Main Effect Global Sensitivity Indices." Reliability Engineering & System Safety 165 (September 2017): 222-38. https://doi.org/10.1016/j.ress.2017.04.003.

# ===================================================
#  X = input
#  Y = output
#  m = number of slides
# ===================================================

SCA = function(X, Y, m) {

  # initialize  
  Si      = rep(NA, n_k)                   # main effect Si
  PAWNi   = rep(NA, n_k)                   # PAWN index 
  temp_PAWN = array(NA,c(n_m,n_k))         # PAWN
  MCY     = array(NA, c(n_m, n_k))         # initialize mean of the conditional Y
  y_F     = array(NA, c(N, n_m, n_k))      # conditional samples
  
  # calculations
  VUY = var(Y) # variance of the unconditional Y
  
  # ============================================
  # to estimate main effect
  for (k in 1:n_k){
    # k=1
    qua_x = quantile(X[, k], probs = seq(0, 1, (1/n_m)), na.rm=TRUE)
    
    for (m in 1:n_m) {
      # m=1
      y_F_temp = Y[X[, k] >= qua_x[m] & X[, k] < qua_x[m+1]]
      
      # calculate mean for Si
      MCY[m, k]  = mean(y_F_temp, na.rm=T)
      
      # calculate KD for PAWN
      temp_ks = ks.test(y_F_temp,Y)
      temp_PAWN[m,k] = temp_ks$statistic  
      
    } #end loop over n_m

    Si[k]  = var(MCY[,k])/VUY
    PAWNi[k] = median(temp_PAWN[,k])
      
  } # loop over k
  
  index = cbind(Si, PAWNi)
  return(index)

} # end function

