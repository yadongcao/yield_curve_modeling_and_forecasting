
# https://www.r-bloggers.com/2021/05/arbitrage-free-nelson-siegel-model-with-r-code/
setwd("C:/D/AI_Workspace/yield_curve_modeling_and_forecasting/AFNS")
getwd()

#install.packages('expm')

library(readxl)
library(expm)
library(numDeriv)

# DNS factor loading matrix
NS.B<-function(lambda, tau) {
  col1 <- rep.int(1,length(tau))
  col2 <- (1-exp(-lambda*tau))/(lambda*tau)
  col3 <- col2 - exp(-lambda*tau) 
  return(cbind(col1,col2,col3))
}

# yield adjustment term in AFNS
AFNS.C<-function(sigma, lambda, tau) {
  s <- sigma; la <- lambda; t <- tau
  
  s11<-s[1,1]; s12<-s[1,2]; s13<-s[1,3]
  s21<-s[2,1]; s22<-s[2,2]; s23<-s[2,3]
  s31<-s[3,1]; s32<-s[3,2]; s33<-s[3,3]
  
  A<-s11^2+s12^2+s13^2; D<-s11*s21+s12*s22+s13*s23
  B<-s21^2+s22^2+s23^2; E<-s11*s31+s12*s32+s13*s33 
  C<-s31^2+s32^2+s33^2; F<-s21*s31+s22*s32+s23*s33
  
  r1<--A*t^2/6; 
  r2<--B*(1/(2*la^2)-(1-exp(-la*t))/(la^3*t)+
            (1-exp(-2*la*t))/(4*la^3*t))
  r3<--C*(1/(2*la^2)+exp(-la*t)/(la^2)-t*exp(-2*la*t)/(4*la)-
            3*exp(-2*la*t)/(4*la^2)-2*(1-exp(-la*t))/(la^3*t)+
            5*(1-exp(-2*la*t))/(8*la^3*t))
  r4<--D*(t/(2*la)+exp(-la*t)/(la^2)-(1-exp(-la*t))/(la^3*t))
  r5<--E*(3*exp(-la*t)/(la^2)+t/(2*la)+t*exp(-la*t)/(la)-
            3*(1-exp(-la*t))/(la^3*t))
  r6<--F*(1/(la^2)+exp(-la*t)/(la^2)-exp(-2*la*t)/(2*la^2)-
            3*(1-exp(-la*t))/(la^3*t)+3*(1-exp(-2*la*t))/(4*la^3*t))
  return(r1+r2+r3+r4+r5+r6)
}

# parameter restrictions
trans<-function(b) {
  bb <- b
  bb[1]  <- 1/(1+exp(b[1]))  # kappa11
  bb[13] <- b[13]^2          # lambda
  bb[14:npara] <- b[14:npara]^2          # measurement error
  return(bb)
}

# log likelihood function
loglike<-function(para_un,m.spot) {
  # parameter restrictions
  para <- trans(para_un)
  
  # restricted parameters
  kappa  <- rbind(c(para[1],0,0),
                  c(0,para[2],0),
                  c(0,0,para[3]))
  sigma  <- rbind(c(para[4],0,0),
                  c(para[5],para[6],0),
                  c(para[7],para[8],para[9]))
  theta  <- para[10:12]
  lambda <- para[13]
  H      <- diag(para[14:npara])
  
  B  <- NS.B(lambda,v.mat); tB <- t(B) # factor loading matrix
  C  <- AFNS.C(sigma,lambda,v.mat)     # yield adjustment
  
  # Conditional and Unconditional covariance matrix : Q, Q0
  m    <- eigen(kappa) 
  eval <- m$values 
  evec <- m$vectors; ievec<-solve(evec)
  Smat <- ievec%*%sigma%*%t(sigma)%*%t(ievec)
  Vdt  <- Vinf <- matrix(0,nf,nf)
  
  for(i in 1:nf) { for(j in 1:nf) {
    Vdt[i,j] <-Smat[i,j]*(1-exp(-(eval[i]+eval[j])*dt))/
      (eval[i]+eval[j]) # conditional
    Vinf[i,j]<-Smat[i,j]/(eval[i]+eval[j]) # unconditional
  }}
  
  # Analytical Covariance matrix
  # Q : conditional, Q0 : unconditional
  Q  <- evec%*%Vdt%*%t(evec)
  Q0 <- evec%*%Vinf%*%t(evec)
  
  # initialzation of vector and matrix
  prevX <- theta; prevV <- Q0
  Phi1  <- expm(-kappa*dt)
  Phi0  <- (diag(nf)-Phi1)%*%theta
  loglike <- 0 # log likelihood function
  
  for(i in 1:nobs) {
    Xhat <- Phi0+Phi1%*%prevX        # predicted state
    Vhat <- Phi1%*%prevV%*%t(Phi1)+Q # predicted cov
    
    y        <- m.spot[i,] # the observed yield
    yimplied <- B%*%Xhat+C # the model-implied yields
    er       <- y-yimplied # prediction error
    
    # updating 
    ev <- B%*%Vhat%*%tB+H; iev<-solve(ev)
    KG <- Vhat%*%tB%*%iev # Kalman Gain
    
    prevX <- Xhat+KG%*%er       # E[X|y_t]   updated state 
    prevV <- Vhat-KG%*%B%*%Vhat # Cov[X|y_t] updated cov
    
    # log likelihood function
    loglike <- loglike - 0.5*(nmat)*log(2*pi)-
      0.5*log(det(ev))-0.5*t(er)%*%iev%*%er
    gm.factor[i,] <<- prevX
  }
  
  return(-loglike)
}