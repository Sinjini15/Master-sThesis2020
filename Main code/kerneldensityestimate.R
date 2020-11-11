kerneldensityestimate <- function(data,h, H) {
  library(kdensity)
  library(ks)
  
  data <- as.matrix(data)
  kerd <- kde(data, H=H, gridsize = 300)
  xygrid <- kerd$eval.points
  p_hat <- kerd$estimate
  x_points <- xygrid[[1]]
  y_points <- xygrid[[2]]
  
  

  
  #==============Hypothesis testing section===========# 
  #==============Generate the z_(i) values such that==# 
  #z_(i) = p_hat_h(X_i) where X_i belongs to [X_1,X_2,...,X_n]
  # Here each X_i = (t_i,R_i) and n = total number of data points
  
  kerd_hypothesis <- kde(data,H=H, eval.points= data)
  
  # Order the z_i values in increasing order
  est_z <- kerd_hypothesis$estimate
  
  h1 <- kerd_hypothesis$H
  z_ord <- sort(est_z, decreasing = FALSE)
  
  # Compute k as 
  k <- floor((n+1)*alpha)
  
  # Here alpha = 0.05 since we want to predict within 95# 
  confidencek <-  floor((n+1)*alpha)
  
  # Fetch the z_(k) value from the sorted z_k set
  if (k != 0){    
    z_k <- z_ord[k]} else {    
      z_k <- z_ord[1]}
  
  # Compute t_k = z_(k) - (K(0)/(n*h^d))
  
  org <- as.matrix(cbind(min(x_points),min(y_points)))
  
  K_0 <- kde(data, H=H, eval.points=org)
  c_k <- z_k - ((K_0$estimate)/(n*H[1]))

  
outlist <- list (p_hat,c_k, xygrid,kerd)
return(outlist)
}