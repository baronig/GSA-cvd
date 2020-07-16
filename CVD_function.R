# ===================================================
# ===================================================
# Global sensitivity analysis based on combined
# variance and distribution-based approach (CVD).
# Baroni and Francke, under review
# 
# The function estimate main and interaction index
# based on a generic sample N and interval m.

# 25 February 2020
# 
# Gabriele Baroni
# g.baroni@unibo.it
# University of Bologna (Italy)
# 
# Till Francke
# francke@uni-potsdam.de
# University of Potsdam (Germany)

# ===================================================
#  X = input
#  Y = output
#  m = number of slides
# ===================================================

CVD = function(X, Y, m) {

  # initialize  
  Si      = rep(NA, n_k)                   # Si based on CVD
  Ii      = rep(NA, n_k)                   # interaction terms 
  MCY     = array(NA, c(n_m, n_k))         # initialize mean of the conditional Y
  n_comb  = factorial(n_m)/(factorial(2)
                      *factorial(n_m-2))   # number of combinations
  Ii_temp = array(NA, c(n_comb, n_k))      # interaction indices
  y_F     = array(NA, c(N, n_m, n_k))      # conditional samples

  # calculations
  VUY = var(Y) # variance of the unconditional Y
  
  # ============================================
  # to estimate main effect
  for (k in 1:n_k){
    # k=1
    qua_x = quantile(X[, k], probs = seq(0, 1, (1/n_m)), na.rm=TRUE)
    
    for (m in 1:n_m) {
  
      y_F_temp = Y[X[, k] >= qua_x[m] & X[, k] < qua_x[m+1]]
      MCY[m, k]  = mean(y_F_temp, na.rm=T)                        # calculate conditional mean
      y_F[(1:length(y_F_temp)), m, k] = y_F_temp - MCY[m, k]       # store the centered conditional distributions
      
    } #end loop over n.m

    s      = smooth.spline(X[, k], Y, nknots = n_m)
    Si[k]  = var(s$y, na.rm=T)/VUY
      
  } # loop over k
  
  # ============================================
  # Interaction index based on conditional distribution
  for (k in 1:n_k){

    j = 1
    for (m1 in 1:(n_m-1)) {
      # m1=1
      y_F1 = y_F[, m1, k]
      y_F1 = y_F1[!is.na(y_F1)]
      
      for (m2 in (m1+1):n_m) {
        
        y_F2 = y_F[, m2, k]
        y_F2 = y_F2[!is.na(y_F2)]
        
        ub = max(y_F1, y_F2)
        lb = min(y_F1, y_F2)
        df_F1 = kde(x = y_F1, xmin = lb, xmax = ub)
        df_F2 = kde(x = y_F2, xmin = lb, xmax = ub)
        Ii_temp[j, k] = sum(abs(df_F1$estimate/sum(df_F1$estimate)-df_F2$estimate/sum(df_F2$estimate)))
  
        j = j +1
      }
    }
    
    Ii[k] = 0.5*mean(Ii_temp[, k], na.rm = T)      
 
  } # loop over k  
  
  index = cbind(Si, Ii)
  return(index)

} # end function

