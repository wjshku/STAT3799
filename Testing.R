set.seed(30)
rounds <- 500
nw_res <- data.frame()
for(i in 1:rounds){
  X <- rnorm(32, 0, 5)
  Y <- mx(X) + rnorm(32, 0, 1)
  train <- data.frame(X,Y)
  nw_res <- rbind(nw_res, nw_simu(train, rounds, test, 32 , ratioh[3], k)$esti)
}
var(nw_res$esti_nw)
hist(nw_res$esti_nw)
quantile(nw_res$esti_nw,c(.025, .975))



x <- test$X
n <- 500
h <- n^(-1/5)
X <- rnorm(n, 0, 5)
Y <- mx(X) + rnorm(n, 0, 1)
mhat(n,h,x,X,Y,k)

mhat_1(n,h,x,X,Y,k,k_1)

hist(coverage$lower,main = "h = 0.1, g = 1, n = 1024")
var(coverage$lower, na.rm = T)

n <- para$n
h <- result_hy$h[HY_opt]
m <- para$m
g <- result_hy$g[HY_opt]
X <- train_label$X
Y <- train_label$Y
mx_deriv1 <- mhat_1(n,h,test$X,X,Y,k,k_1)
px <- phat(n,h,test$X,X,k)
qx <- phat(m,g,test$X,train_un$X,k)
esti_var(n,h,m,g,mx_deriv1,sigmaep,px,px,rk$value,sigmak$value,E_H_ratio = 1)


