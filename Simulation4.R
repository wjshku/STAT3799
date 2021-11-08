
sigmax <- 5
p <- function(x) 1/sqrt(2*pi*sigmax^2)*exp(-0.5*(x-0)^2/sigmax^2)
p_1 <- Deriv(p)
p_2 <- Deriv(p_1)
q <- function(x) 1/sqrt(2*pi*sigmax^2)*exp(-0.5*(x-0)^2/sigmax^2)
q_1 <- Deriv(q)
q_2 <- Deriv(q_1)
mx_1 <- Deriv(mx)
mx_2 <- Deriv(mx_1)

sigmak <- integrate(function(x) k(x)^2*x^2, -Inf, Inf)
rk <- integrate(function(x) k(x)^2, -Inf, Inf)
muk <- integrate(function(x) k(x)*x^2, -Inf, Inf)
#rp <- integrate(function(x) p_2(x)^2, -Inf, Inf)
#rq <- integrate(function(x) q_2(x)^2, -Inf, Inf)

sigmax <- 5
theta22 <- integrate(function(x) mx_2(x)^2*p(x), -Inf, Inf)
sigmaep <- 1


simu <- function(multi, n, m, sigmax = 5, sigmaep = 1){
  
  X <- rnorm(n*multi, sd = sigmax)
  Y <- mx(X) + rnorm(n, sd = sigmaep)
  
  return(data.frame(X, Y))
}

nw_simu <- function(data, rounds, testd, n, h, k){
  
  esti_nw <- testd$Y
  esti <- data.frame()
  X <- data$X
  Y <- data$Y
  real <- testd$Y
  test <- testd$X
  
  for(j in 1:length(test)){
    #set.seed(i) # random, may exist better solutions
    esti_nw[j] <- mhat(n,h,test[j],X,Y,k)
  }
  
  MSE_NW <- mean((real - esti_nw)^2)
  #browser()
  
  esti <- rbind(esti,data.frame(X = test,esti_nw))
  
  return(list(MSE = MSE_NW, esti = esti))
  
}

#Self-Supervised
ss_simu <- function(data, data_un, rounds, testd, n, h, m, g, k){
  esti <- data.frame()
  esti_ss <- testd$Y
  
  X <- data$X
  Y <- data$Y
  Xu <- data_un$X
  
  test <- testd$X
  real <- testd$Y
  
  for(j in 1:length(test)){
    #set.seed(i) # random, may exist better solutions
    esti_ss[j] <- rhat(n,h,test[j],X,Y,k,m,g,Xu)
  }
  
  MSE_SS <- mean((real - esti_ss)^2)
  #browser()
  
  esti <- rbind(esti,data.frame(X = test,esti_ss))
  
  
  return(list(MSE = MSE_SS, esti = esti))
  
}

hy_simu <- function(rounds, testd, lambda, esti_nw, esti_ss){
  
  esti <- data_frame()
  
  real <- testd$Y
  
  esti_hy <- lambda * esti_nw + (1-lambda) * esti_ss
  
  esti <- rbind(esti,data.frame(X = testd$X ,esti_ss))
  
  MSE_HY <- mean((real - esti_hy)^2)
  
  return(list(MSE = MSE_HY, esti = esti))
}

#Grid search for bandwidth

para <- paraset3(64,10/19,sigmaep,rk$value, muk$value,theta22$value)

data <- simu(5, para$n, para$m, sigmax = sigmax, sigmaep = sigmaep)

rounds <- 50

ratioh <- seq(0.5,4,0.5)
ratiog <- seq(3,10,1)

mMSE_NW <- 0
mMSE_HY <- 0

set.seed(20)
record_nw <- 1 : rounds#
record_hy <- 1 : rounds#

for(exp in 1:rounds){
  sample <- sample.int(n = nrow(data), size = floor(.75*nrow(data)), replace = FALSE)
  train_label <- data[sample[1:para$n],]
  train_un <- data[sample[(1+para$n):(para$m+para$n)],]
  #test <- data[-sample,]
  test <- data.frame(X = 2, Y = mx(2))
  result_nw <- data.frame()
  for(i in 1:length(ratioh)){
    h <- ratioh[i]*para$h
    MSE_NW <- nw_simu(train_label, rounds, test, para$n,h, k)$MSE
    
    result_nw <- rbind(result_nw, data.frame(mMSE_NW = MSE_NW, h))
  }
  record_nw[exp] <- result_nw$mMSE_NW#
  mMSE_NW <- mMSE_NW + result_nw$mMSE_NW
  
  result_hy <- data.frame()
  for(i in 1:length(ratioh)){
    for(j in 1:length(ratiog)){
      h <- ratioh[i]*para$h
      g <- ratiog[j]*para$g
      para$lambda <- 1 + h^2/g^2
      nw_res <- nw_simu(train_label, rounds, test, para$n,h, k)
      ss_res <- ss_simu(train_label, train_un, rounds, test, para$n,h, para$m,g, k)
      
      hy_res <- hy_simu(rounds, test, para$lambda, nw_res$esti$esti_nw, ss_res$esti$esti_ss)
      
      result_hy <- rbind(result_hy, data.frame(mMSE_HY = hy_res$MSE, h, g))
    }
  }
  record_hy[exp] <- result_hy$mMSE_HY#
  mMSE_HY <- mMSE_HY + result_hy$mMSE_HY
}

hist(record_nw)
mean(record_nw)

hist(record_hy)
mean(record_hy)

result_nw$mMSE_NW <- mMSE_NW/rounds
result_hy$mMSE_HY <- mMSE_HY/rounds

min(result_nw$mMSE_NW,na.rm = TRUE)
min(result_hy$mMSE_HY,na.rm = TRUE)

match(min(result_nw$mMSE_NW),result_nw$mMSE_NW)
match(min(result_hy$mMSE_HY,na.rm = TRUE),result_hy$mMSE_HY)

plot(mMSE_NW~h,data = result_nw,col = "blue")

plot(mMSE_HY~h,data = filter(result_hy, g == result_hy$g[36]),main = "g = 0.87")
plot(mMSE_HY~g,type = "b",data = filter(result_hy, h == result_hy$h[36]),main = "h = 0.64")


plot(mMSE_HY~h,data = filter(result_hy, g == result_hy$g[36]), type="b", pch=19, col="red", xlab="h", ylab="MSE",main = "g = 2.27")
# Add a line
lines(mMSE_NW~h,data = result_nw, pch=18, col="blue", type="b", lty=2)
# Add a legend
legend("topright", legend=c("MSE_HY", "MSE_NW"),col=c("red", "blue"), lty=1:2, cex=0.6)
```
