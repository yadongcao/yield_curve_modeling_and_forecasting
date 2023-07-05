

library(tidyverse)
# 读取数据
# https://www.r-bloggers.com/2021/05/arbitrage-free-nelson-siegel-model-with-r-code/
filename <- "../中国利率市场/data/ytm_data.csv"
raw_df <- read.csv(file = filename, encoding = "UTF-8")

data <- as_tibble(raw_df)
# 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, and 20 years.



colnames(data)

# M03
0.25 * 12

all_data <- data 


colnames(all_data) <- c("Date", "0.5", "1", "2","3", "4","5", "6","7", "8","9","10", "15","20", "30")



colnames(all_data)

library(lubridate)
all_data[['Date']] <- as.POSIXct(all_data[['Date']],
                                   format = "%Y-%m-%d")
# 选择
need_data <- all_data[all_data$Date>=ymd(20120101) 
                           & all_data$Date<=ymd(20221231),]

need_data <- drop_na(need_data)
need_data <- need_data %>% 
  arrange(Date)


  
dt <- 1/360 #1/52 # weekly
nf <- 3

df.spot <- need_data

# divide date and data
v.ymd <- df.spot[,1]
v.mat <- as.numeric(colnames(df.spot)[-1])
m.spot <- as.matrix(df.spot[,-1])

nmat <- length(v.mat)
nobs <- nrow(m.spot)

# factor estimates
gm.factor <- matrix(0, nobs, nf)

#—————————————————–
# initial guess for unconstrained parameters(para_un)
#—————————————————–
init_para_un <- c(
  1.226637,  0.840692,  0.603496, # kappa
  0.006327,  
  -0.005464,  0.003441,  
  -0.000707, -0.003399,  0.010891, # sigma
  0.032577, -0.012536, -0.002748, # theta 
  0.5    ,                       # lambda
  rep(0.000524,nmat)              # measurement error
)

npara <- length(init_para_un) # # of observations

m<- optim(init_para_un,loglike,
         control = list(maxit=100, trace=2),
         method=c("Nelder-Mead"),m.spot=m.spot)
prev_likev <- m$value
v.likev    <- m$value


i <- 1
repeat {
  print(paste(i,"-th iteration"))
  m<- optim(m$par,loglike,
           control = list(maxit=100, trace=2),
           method=c("Nelder-Mead"),m.spot=m.spot)
  
  v.likev <- cbind(v.likev, m$value)
  print(paste(m$value," <- likelihood value"))
  
  print(c("result",i,  m$value, prev_likev))    
  if (abs(m$value - prev_likev) < 0.000000001) {
      break
  }
  prev_likev <- m$value
  i <- i + 1
  
  
  print(v.likev)
}

name_theta <- c("kappa_1", "kappa_2", "kappa_3",
  "sigma_11", "sigma_21", "sigma_22", 
  "sigma_31", "sigma_32", "sigma_33",
  "theta_1",  "theta_2", "theta_3", 
  "lambda", 
  "epsilon_1", "epsilon_2", "epsilon_3", 
  "epsilon_4", "epsilon_5", "epsilon_6", 
  "epsilon_7", "epsilon_8", "epsilon_9", 
  "epsilon_10", "epsilon_11", "epsilon_12", "epsilon_13", "epsilon_14")


# draw 3 factor estimates
x11(width=6, height=5)
matplot(gm.factor,type="l", ylab="L,S,C", lty = 1, 
        main = "AFNS 3 Factor Estimates (L,S,C)", lwd=2)




# Delta method for statistical inference
grad    <- jacobian(trans, m$par)
hess    <- hessian(func=loglike, x=m$par,m.spot=m.spot)
vcv_con <- grad%*%solve(hess)%*%t(grad)

# parameter | std.err | t-value | p-value
theta   <- trans(m$par)
stderr  <- sqrt(diag(vcv_con))
tstat   <- theta/stderr
pvalue  <- 2*pt(-abs(tstat),df=nobs-npara)
df.est  <- cbind(theta, round(stderr,4),
                 round(tstat,4), round(pvalue,4))

rownames(df.est) <- name_theta # parameter name
colnames(df.est) <-c("parameter","std.err","t-stat","p-value")
print(df.est)

para_un <- c(2.946924769, 0.845110501, 0.634485561,
             -0.087337617, 0.165153975, 0.406481743,
             0.119850488, -0.169886793, 0.285420706,
             3.209325303, 0.007235799, 0.030360527,
             0.380626177,
             0.193047489, 0.146635489, 0.088179026, 
             -0.030440797, 0.015689684, 0.021855681,
             0.025738274, 0.056673128, 0.025700232,
             0.052354582, 0.115847195, 0.097461404,
             0.097718369, 0.094590434)

# 实际值跟因子计算值
build_curve<-function(para_un, factor, mat_vec) {
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
  
  B  <- NS.B(lambda,mat_vec); tB <- t(B) # factor loading matrix
  C  <- AFNS.C(sigma,lambda,mat_vec)     # yield adjustment
  
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
  
  prevX <- factor
  prevV <- Q0
  Phi1  <- expm(-kappa*dt)
  Phi0  <- (diag(nf)-Phi1)%*%theta
  loglike <- 0 # log likelihood function
  
  
  Xhat <- Phi0+Phi1%*%prevX        # predicted state
  
  
  Vhat <- Phi1%*%prevV%*%t(Phi1)+Q # predicted cov
  
  # 
  #y        <- m.spot[1,] # the observed yield
  yimplied <- B%*%Xhat+C # the model-implied yields
  #er       <- y-yimplied # prediction error
  
  return(yimplied)
  
  
}

# para_un -> m$par
# factor -> gm.factor
# AFNS参数
para_un <- m$par
# DNS因子
factor<-t(t(gm.factor[1,]))
# 实际值
y        <- m.spot[1,]

mat_vec <- c(0.5, 1, 2,  3,  4,  5,  6, 7, 8, 9, 10, 15, 20, 30)
test <- build_curve(para_un, factor, mat_vec)

# 保存数据
# 真实的利率数据，加上日期
# 因子数据，加上日期
# 模型计算出来的利率数据，加上日期


col_name_vec <- c('6M', '1Y', '2Y', '3Y', '4Y', '5Y', '6Y', '7Y', '8Y', '9Y',
                  '10Y', '15Y', '20Y', '30Y')
colnames(m.spot) <- col_name_vec

temp_real_yield_df <- as_tibble(m.spot)

real_yield_df <- temp_real_yield_df %>%
  add_column(need_data[,1],
             .before = "6M") 


# 因子数据
factor_df <- as_tibble(gm.factor, digits=6)%>%
                add_column(need_data[,1],
                           .before = "V1") 
colnames(factor_df) <- c("Date", "Level", "Slope", "Curvature")

# 模型计算
rows <- dim(factor_df)[1]
temp_yimplied_m <- matrix(0, 0, 14)

for(i in seq_along(1:rows)) {
  temp_factor_df <- factor_df[i,]
  factor <- c(temp_factor_df[[2]], temp_factor_df[[3]], temp_factor_df[[4]])
  #para_un <- m$par
  yimplied <- build_curve(para_un, factor, mat_vec)
  
  temp_yimplied_m <- rbind(temp_yimplied_m, t(yimplied))

}

temp_yimplied_df <- as_tibble(temp_yimplied_m)
colnames(temp_yimplied_df) <- col_name_vec

yimplied_df <- temp_yimplied_df %>%
  add_column(need_data[,1],
             .before = "6M") 


write.table(real_yield_df, file = "CN_Real_Yield.csv", sep=",", row.names = FALSE)

write.table(factor_df, file = "CN_AFNS_Factor.csv", sep=",", row.names = FALSE)

write.table(yimplied_df, file = "CN_AFNS_Model_Yield.csv", sep=",", row.names = FALSE)

# CN_AFNS_Pred_Stats_Yield.csv
# 1.读取预测因子

pred_factor <- read.csv(file = "../中国利率市场/data/ML_Pred_AFNS_Factor_1Y.csv", encoding = "UTF-8")
pred_factor <- as_tibble(pred_factor)



para_un <- c(2.90615233, 0.82607139, 0.63300710,
             -0.09104841, 0.11148613, 0.45139625, 
             0.19204882, -0.14628164, -0.18943611,
             2.64910629, -0.02652632, -0.11065036,  
             0.40565026,
             0.16145567, 0.10306592, 0.09454710,
             0.03178410, 0.01677096, 0.02772500,  
             0.02621404, 0.07172742, -0.02737894,
             0.08785550, 0.12496281, 0.10120062, 0.09316415, 0.08639887)

mat_vec <- c(0.5, 1, 2,  3,  4,  5,  6, 7, 8, 9, 10, 15, 20, 30)


rows <- dim(pred_factor)[1]
temp_pred_yimplied_m <- matrix(0, 0, 14)

for(i in seq_along(1:rows)) {
  temp_factor_df <- pred_factor[i,]
  factor <- c(temp_factor_df[[2]], temp_factor_df[[3]], temp_factor_df[[4]])
  para_un <- para_un
  yimplied <- build_curve(para_un, factor, mat_vec)
  
  temp_pred_yimplied_m <- rbind(temp_pred_yimplied_m, t(yimplied))
  
}

temp_pred_yimplied_df <- as_tibble(temp_pred_yimplied_m)
colnames(temp_pred_yimplied_df) <- col_name_vec

pred_yimplied_df <- temp_pred_yimplied_df %>%
  add_column(pred_factor[,1],
             .before = "6M") 

colnames(pred_yimplied_df)
names(pred_yimplied_df)[1] <- 'Date'

write.table(pred_yimplied_df, file = "CN_AFNS_Pred_ML_Yield_1Y.csv", sep=",", row.names = FALSE)