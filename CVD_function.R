# ===================================================
# ===================================================
# Global sensitivity analysis based on combined
# variance and distribution-based approach (CVD).
# The function estimate main and interaction index
# based on a generic sample N and interval m.
#  
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
#  y.sam = output
#  X.sam = input
#  m = number of slides
# ===================================================

CVD = function(X, Y, m) {

  # initialize  
  Si      = rep(NA, n.k)                   # Si based on CVD
  Ii      = rep(NA, n.k)                   # interaction terms 
  MCY     = array(NA, c(n.m, n.k))         # initialize mean of the conditional Y
  n.comb  = factorial(n.m)/(factorial(2)
                      *factorial(n.m-2))   # number of combinations
  Ii.temp = array(NA, c(n.comb, n.k))      # interaction indices
  y.F     = array(NA, c(N, n.m, n.k))      # conditional samples

  # calculations
  VUY = var(Y.sam) # variance of the unconditional Y
  
  # ============================================
  # to estimate main effect
  for (k in 1:n.k){
    # k=1
    qua.x <- quantile(X[, k], probs=seq(0, 1, (1/n.m)), na.rm=TRUE)
    
    for (m in 1:n.m) {
  
      y.F.temp = Y[X.sam[, k] >= qua.x[m] & X.sam[, k] < qua.x[m+1]]
      MCY[m, k]  = mean(y.F.temp, na.rm=T)                        # calculate conditional mean
      y.F[(1:length(y.F.temp)), m, k] = y.F.temp - MCY[m, k]       # store the centered conditional distributions
      
    } #end loop over n.m

    s      = smooth.spline(X[, k], Y, nknots = n.m)
    Si[k]  = var(s$y, na.rm=T)/VUY
      
  } # loop over k
  
  # ============================================
  # Interaction index based on conditional distribution
  for (k in 1:n.k){

    j = 1
    for (m1 in 1:(n.m-1)) {
      # m1=1
      y.F1 = y.F[, m1, k]
      y.F1 = y.F1[!is.na(y.F1)]
      
      for (m2 in (m1+1):n.m) {
        
        y.F2 = y.F[, m2, k]
        y.F2 = y.F2[!is.na(y.F2)]
        
        ub = max(y.F1, y.F2)
        lb = min(y.F1, y.F2)
        df.F1 = kde(x=y.F1, xmin=lb, xmax=ub)
        df.F2 = kde(x=y.F2, xmin=lb, xmax=ub)
        Ii.temp[j, k] = sum(abs(df.F1$estimate/sum(df.F1$estimate)-df.F2$estimate/sum(df.F2$estimate)))
  
        j = j +1
      }
    }
    
    Ii[k] = 0.5*mean(Ii.temp[, k], na.rm = T)      
 
  } # loop over k  
  
  index = cbind(Si, Ii)
  return(index)

} # end function

  
