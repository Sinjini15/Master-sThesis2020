kerneldensity <- function( data, h, eval.points, type) {
  
  xseq <- eval.points[,1]
  yseq <- eval.points[,2]
  tp <- data[,1]
  rp <- data[,2]
  n <- length(tp)
  neval <- length(xseq)
  
  #==== Choosing the kernel type =====#
  
  if(type == "gaussian") {
    
    W1 <- matrix(rep(xseq, rep(n, neval)), ncol = n, byrow = TRUE)
    Y1 <- matrix(rep(tp, neval), ncol = n, byrow = TRUE)
    W1 <- W1- Y1
    Wx <- exp(-0.5 * (W1/h)^2)/(h*sqrt(2*pi))
    W2 <- matrix(rep(yseq, rep(n, neval)), ncol = n, byrow = TRUE)
    Y2 <- matrix(rep(rp, neval), ncol = n, byrow = TRUE)
    W2 <- W2-Y2
    Wy <- exp(-0.5 * (W2/h)^2)/(h*sqrt(2*pi))
    est2d <- (Wx %*% t(Wy))/n
  }
  
  else if (type == "epanechnikov"){
    W1 <- t(matrix(rep(tp, neval), ncol = n, byrow = TRUE))
    Y1 <- matrix(rep(xseq, n), ncol = neval, byrow = TRUE)
    W1 <- W1- Y1
    Wx <- 0.75*(1-(W1/h)^2)/h
    Wx <- colSums(Wx)/(n*h)
    #reject non negative values
    ind <- which(Wx<0, arr.ind = T)
    for (i in 1:length(ind)) {
      val <- ind[i]
      Wx[val] <- 0
    }
    W2 <- t(matrix(rep(rp,neval), ncol = n, byrow = TRUE))
    Y2 <- matrix(rep(yseq,n), ncol = neval, byrow = TRUE)
    W2 <- W2-Y2
    Wy <- 0.75*(1-(W2/h)^2)/h
    Wy <- colSums(Wy)/(n*h)
    #reject non negative values
    ind <- which(Wy<0, arr.ind = T)
    for (i in 1:length(ind)) {
      val <- ind[i]
      Wy[val] <- 0
    }
    est2d <- Wx%*% t(Wy)
  } else if (type == "rectangular"){
    
    W1 <- matrix(rep(tp,neval), ncol = n, byrow = TRUE)
    Y1 <- t(matrix(rep(xseq, n), ncol = neval, byrow = TRUE))
    W1 <- (W1- Y1)/h
    indx1 <- which(W1 < -1, arr.ind = T)
    if (length(indx1) == 0) {
      W1 <- W1
    } else {
      for (i in 1:length(indx1[,1])) {
        a <- indx1[i,1]
        b <- indx1[i,2]
        W1[a,b] <- 0
      }
      
    }
    indx2 <- which(W1 > 1, arr.ind = T)
    if (length(indx2) == 0){
      W1 <- W1
    } else {
      for (i in 1:length(indx2[,1])) {
        a <- indx2[i,1]
        b <- indx2[i,2]
        W1[a,b] <- 0
      }
    }
    #replace all these indices with 0.5
    
    indx3 <- which(W1!=0, arr.ind = T)
    if(length(indx3)==0){
      for (i in 1:length(W1[,1])) {
        W1[i,1] <- 0.5
        W1[i,2] <- 0.5
      }
    } else {
      for (i in 1:length(indx3[,1])) {
        a <- indx3[i,1]
        b <- indx3[i,2]
        W1[a,b] <- 0.5
      }
    }
    
    Wx <- rowSums(W1)
    Wx <- Wx/(n*h)
    
    
    W2 <- matrix(rep(rp, neval), ncol = n, byrow = TRUE)
    Y2 <- t(matrix(rep(yseq, n), ncol = neval, byrow = TRUE))
    W2 <- (W2-Y2)/h
    indy1 <- which(W2 < -1, arr.ind = T)
    if (length(indy1)==0) {
      W2 <- W2
    } else {
      for (j in 1:length(indy1[,1])) {
        a <- indy1[j,1]
        b <- indy1[j,2]
        W2[a,b] <- 0
      }
      
    }
    indy2 <- which(W2 > 1, arr.ind = T)
    if (length(indy2)==0){
      W2 <- W2
    } else {
      for (j in 1:length(indy2[,1])) {
        a <- indy2[j,1]
        b <- indy2[j,2]
        W2[a,b] <- 0
      }
    }
    
    #replace all these indices with 0.5
    
    indy3 <- which(W2!=0, arr.ind = T)
    if (length(indy3)==0){
      for (j in 1:length(W2[,1])) {
        W2[j,1] <-0.5
        W2[j,2] <- 0.5
      }
    } else {
      for (j in 1:length(indy3[,1])) {
        a <- indy3[j,1]
        b <- indy3[j,2]
        W2[a,b] <- 0.5
      }
      
    }
    
    Wy <- rowSums(W2)
    Wy<- Wy/(n*h)
    
    est2d <- (Wx %*% t(Wy)) }
  
  else{
    
    W1 <- matrix(rep(tp, neval), ncol = n, byrow = TRUE)
    Y1 <- t(matrix(rep(xseq, n), ncol = neval, byrow = TRUE))
    W1 <- (W1- Y1)/h
    Wx <- (pi/4) * cos((pi/2) * W1)
    
    x  <- which(Wx < 0, arr.ind = T)
    
    if (length(x)!=0) {
      for (i in 1:length(x[,2])) {
        x_val <- x[i,1]
        y_val <- x[i,2]
        Wx[x_val, y_val] <- 0
      }
    } else {
      Wx <- Wx
    }
    
    
    Wx <- rowSums(Wx)/(n*h)
    
    
    W2 <- matrix(rep(rp, neval), ncol = n, byrow = TRUE)
    Y2 <- t(matrix(rep(yseq, n), ncol = neval, byrow = TRUE))
    W2 <- (W2-Y2)/h
    Wy <- (pi/4) * cos(pi/2 * W2)
    y  <- which(Wy < 0, arr.ind = T)
    
    if (length(y)!=0) {
      for (i in 1:length(y[,2])) {
        x_val <- y[i,1]
        y_val <- y[i,2]
        Wy[x_val, y_val] <- 0
      }
    } else {
      Wy <- Wy
    }
    Wy <- rowSums(Wy)/(n*h)
    est2d <- (Wx %*% t(Wy))/(h*n)

  }
 
  outlist <- list(h,xseq,yseq,est2d)
  return(outlist)
}

