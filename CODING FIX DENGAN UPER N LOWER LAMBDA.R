library(splines)
MPL<-function(x,eps=1e-009)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1]) %*% t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)] %*% diag(1/diago) %*% t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}

#x<-c(8827,8854,12427,10916,14909,13185,13709,15997,12999,24221,18199,17765,21357,19201,10860,9210,8244,10588,8227,8177,9428,9620,10791,9950,10776,10166,11294,12193,11927,11757,9044,9289,12058,11229,17639,12087,16239,15926,12500,10578,10967,10904,11905,10277,9776,9282,10671,11108,10011,13250,12522,11841,9780,11798,13052,10610,10067,10937,10948,11609,10913,10698,12448,9773,11999,9972,10707,8994,10020,10514,12816,15463,16351,16047,13158,13455,10511,16002,9874,16438,19319,9184,10199,10042,11162,11001,11565,10326,9466,9840,12320,10851,10263,11254,10726,14808,13051,11579,12349,11848,12031,11563,10323,10703,11648,13384,8971,8944,8967,9388,12762,14058,16897,12571,13803,14054,16503,18345,13094)
#y<-c(65.84,64.71,72.97,67.75,78.9,73.95,72.98,81.95,72.79,85.21,83.45,82.11,82.51,80.81,71.2,67.64,65.94,73.16,67.41,66.84,71.45,70.16,70.06,68.56,72.69,68.55,69.87,71.56,71.74,75.22,69.04,69.03,77.17,75.4,82.5,75.89,82.46,81.86,78.77,73.83,72.55,70.99,73.17,69.54,68.61,70.79,73.6,68.89,70.85,74.97,76.95,77.94,71.04,76.58,74.65,70.97,69.95,71,73.14,75.89,73.15,73.36,74.67,70.77,73.19,69.45,70.81,67.19,69.53,67.03,80.39,83.08,84.35,84.08,75.46,80.69,70.96,84.31,87.69,69.37,71.87,71,74.06,71.86,73.46,71.38,66.95,67.97,71.94,67.31,68.25,66.96,69.68,81.02,74.89,74.05,72.93,72.39,74.85,71.75,70.12,69.67,74.02,77.16,65.05,63.39,66.99,67.87,79.59,79.93,82.71,74.56,76.54,79.32,82.01,82.74,77.22)
x<-c(0.55,-0.05,0.15,0.10,0.08,0.46,0.56,-0.26,-0.11,0.13,0.46,0.57,0.42,-0.08,0.26,0.46,0.42,0.25,0.05,0.07,-0.07,0.18,0.31,0.46,0.27,0.40,0.07,-0.24,0.22,0.08,-0.08,-0.04,0.03,0.08,0.13,0.48,0.54,0.14,0.08,0.01,0.07,0.05,0.11,0.05,-0.17,0.24,0.45,0.71,0.59,0.05,0.77,1.14,0.75,0.52,0.47,-0.12,1.05,0.11,0.32,0.65)
y<-c(-0.05,0.15,0.10,0.08,0.46,0.56,-0.26,-0.11,0.13,0.46,0.57,0.42,-0.08,0.26,0.46,0.42,0.25,0.05,0.07,-0.07,0.18,0.31,0.46,0.27,0.40,0.07,-0.24,0.22,0.08,-0.08,-0.04,0.03,0.08,0.13,0.48,0.54,0.14,0.08,0.01,0.07,0.05,0.11,0.05,-0.17,0.24,0.45,0.71,0.59,0.05,0.77,1.14,0.75,0.52,0.47,-0.12,1.05,0.11,0.32,0.65,0.17)


gcv1 <- function(x, y, m, k, lambda_lower = 0.001, lambda_upper = 0.99) 
  {
  n <- length(y)
  gcv <- 10^10
  aic <- 10^10
  t1 <- seq(min(x), max(x), length.out = 150)
  p <- rep(0, (n - 2))
  for (z in 1:(n - 2)) 
    {
    p[z] <- t1[z + 1]
    }
  comb1 <- combn(p, 1, FUN = NULL)
  c1 <- t(comb1)
  for (m in 2:4) 
    {
    k1 <- c1[, 1]
    K1 <- length(k1)
    for (j1 in 1:K1) 
      {
      bs1 <- bs(x, df = NULL, knot = k1[j1], degree = m - 1, intercept = TRUE, Boundary.knots = range(x))
      B <- cbind(bs1)
      D <- diag(ncol(bs1))
      BtB <- t(B) %*% B
      
      #Mencari lambda optimal
      best_lambda <- NULL
      best_gcv <- Inf
      for (lambda in seq(lambda_lower, lambda_upper, length.out = 150)) 
        {
        beta <- solve(BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
        L <- B %*% solve(BtB + lambda * t(D) %*% D) %*% t(B)
        yhat <- B %*% beta
        MSE <- t(y - yhat) %*% (y - yhat) / n
        I <- diag(n)
        k <- sum(diag(I - L))
        GCV <- MSE / (k/n)^2
        
        if (GCV < best_gcv) 
          {
          best_gcv <- GCV
          best_lambda <- lambda
          knot1 <- k1[j1]
          deter <- det(BtB)
          orde1 <- m
        }
      }
      cat("orde = ", orde1, "knot 1 = ", knot1, "Optimal lambda = ", best_lambda, "GCV = ", best_gcv, "determinant = ", deter, "\n")
    }
  }
}

gcv2<-function(x, y, m, k, lambda_lower = 0.001, lambda_upper = 0.99)
{
  n=length(y)
  gcv=10^10
  aic=10^10
  t1=seq(min(x),max(x),length.out=500)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,2,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    K1=length(k1)
    for (j1 in 1:K1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      D = diag(ncol(bs1))
      wtw = t(w) %*% w
      #Mencari lambda optimal
      best_lambda <- NULL
      best_gcv <- Inf
      for (lambda in seq(lambda_lower, lambda_upper, length.out = 150)) 
      {
        beta <- solve(BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
        L <- B %*% solve(BtB + lambda * t(D) %*% D) %*% t(B)
        yhat <- B %*% beta
        MSE <- t(y - yhat) %*% (y - yhat) / n
        I <- diag(n)
        k <- sum(diag(I - L))
        GCV <- MSE / (k/n)^2
        
        if (GCV < best_gcv) 
        {
          best_gcv <- GCV
          best_lambda <- lambda
          knot1 <- k1[j1]
          knot2 <- k2[j2]
          deter <- det(BtB)
          orde1 <- m
        }
      }
      cat("orde = ", orde1, "knot 1 = ", knot1, "Optimal lambda = ", best_lambda, "GCV = ", best_gcv, "determinant = ", deter, "\n")
    }
  }
}
  
gcv3<-function(x, y, m, k, lambda_lower = 0.001, lambda_upper = 0.99)
{
  n=length(y)
  gcv=10^10
  aic=10^10
  t1=seq(min(x),max(x),length.out=150)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,3,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    k3=c1[,3]
    K1=length(k1)
    for (j1 in 1:K1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1],k3[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1) #membuat matriks B
      D = diag(ncol(bs1))
      wtw = t(w) %*% w
      #Mencari lambda optimal
      best_lambda <- NULL
      best_gcv <- Inf
      for (lambda in seq(lambda_lower, lambda_upper, length.out = 150)) 
      {
        beta <- solve(BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
        L <- B %*% solve(BtB + lambda * t(D) %*% D) %*% t(B)
        yhat <- B %*% beta
        MSE <- t(y - yhat) %*% (y - yhat) / n
        I <- diag(n)
        k <- sum(diag(I - L))
        GCV <- MSE / (k/n)^2
        
        if (GCV < best_gcv) 
        {
          best_gcv <- GCV
          best_lambda <- lambda
          knot1 <- k1[j1]
          knot2 <- k2[j2]
          knot3 <- k3[j3]
          deter <- det(BtB)
          orde1 <- m
        }
      }
      cat("orde = ", orde1, "knot 1 = ", knot1, "Optimal lambda = ", best_lambda, "GCV = ", best_gcv, "determinant = ", deter, "\n")
    }
  }
}
gcv4<-function(x, y, m, k, lambda_lower = 0.001, lambda_upper = 0.99)
{
  n=length(y)
  gcv=10^10
  aic = 10^10
  t1=seq(min(x),max(x),length.out=150)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,4,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    k3=c1[,3]
    k4=c1[,4]
    K1=length(k1)
    for (j1 in 1:K1) 
      {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1],k3[j1],k4[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      D = diag(ncol(bs1))
      wtw = t(w) %*% w
      #Mencari lambda optimal
      best_lambda <- NULL
      best_gcv <- Inf
      for (lambda in seq(lambda_lower, lambda_upper, length.out = 150)) 
      {
        beta <- solve(BtB + lambda * t(D) %*% D) %*% (t(B) %*% y)
        L <- B %*% solve(BtB + lambda * t(D) %*% D) %*% t(B)
        yhat <- B %*% beta
        MSE <- t(y - yhat) %*% (y - yhat) / n
        I <- diag(n)
        k <- sum(diag(I - L))
        GCV <- MSE / (k/n)^2
        
        if (GCV < best_gcv) 
        {
          best_gcv <- GCV
          best_lambda <- lambda
          knot1 <- k1[j1]
          knot2 <- k2[j2]
          knot3 <- k3[j3]
          knot4 <- k4[j4]
          deter <- det(BtB)
          orde1 <- m
        }
      }
      cat("orde = ", orde1, "knot 1 = ", knot1, "Optimal lambda = ", best_lambda, "GCV = ", best_gcv, "determinant = ", deter, "\n")
    }
  }
}



pspline<-function(x,y,m,k,lambda_lower = 0.001, lambda_upper = 0.99)
  {
  n<-length(y)
  print(n)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  bs1 = bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  w = cbind(bs1)
  D = diag(ncol(bs1))
  wtw = t(w) %*% w
  z = MPL(wtw)
  beta= MPL(wtw + lamda*t(D) %*% D) %*% (t(w) %*% y)
  cat("Nilai parameter adalah", "\n", beta, "\n")
  L = w %*% MPL(wtw + lamda*t(D) %*% D) %*% t(w)
  yhat = (w %*% beta)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  I <- matrix(0, ncol=n, nrow = n)
  for (j in 1:n) 
    I[j,j] = 1
  l = sum(diag(I-L))
  GCV=(MSE/(l/n)^2)
  AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(3+m))))
  cat("NIlai GCV adalah ","\n",GCV,"\n")
  cat("NIlai AIC adalah ","\n",AIc,"\n")
}



