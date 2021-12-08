for(i in 1:10){
  train_label <- simu_norm(para$n, 1, meanx = 0, sigmax = sigmax, sigmaep = sigmaep)
  h_opt <- nw_cv(hlist = seq(0.1,2,0.1),train_label)
  #mx_deriv1 <- mhat_1(para$n,h_opt,test$X,train_label$X,train_label$Y,k,k_1)
  
  width <- 0.005
  
  mx_1_esti <- (mhat(para$n,h_opt,test$X*(1 + width),train_label$X,train_label$Y,k) - mhat(para$n,h_opt,test$X*(1 - width),train_label$X,train_label$Y,k))/(test$X*2*width)
  print(str_c("mx_deriv1",mx_deriv1))
  print(str_c("difference",mx_1_esti))
  print(str_c("h_opt: ",h_opt))
}

p_1_esti <- function(para, train_label, test_x, width){
  esti <- (phat(para$n, para$h, test_x*(1 + width),train_label$X, k) - phat(para$n, para$h, test_x*(1 - width),train_label$X, k))/(test_x*2*width)
  return(esti)
}

p_2_esti <- function(para, train_label, test_x, width){
  esti <- (p_1_esti(para, train_label, test_x * (1 + width), width) - p_1_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

p_3_esti <- function(para, train_label, test_x, width){
  esti <- (p_2_esti(para, train_label, test_x * (1 + width), width) - p_2_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

record <- 1:50
for(i in 1:50){
  train_label <- simu_norm(para$n, 1, meanx = 0, sigmax = sigmax, sigmaep = sigmaep)
  record[i] <- (mx_3_esti(para, train_label, test_x, width))
}
hist(record)

mx_1_esti <- function(para, train_label, test_x, width){
  esti <- (mhat(para$n, para$h, test_x*(1 + width),train_label$X, train_label$Y, k) - mhat(para$n, para$h, test_x*(1 - width),train_label$X, train_label$Y, k))/(test_x*2*width)
  return(esti)
}

mx_2_esti <- function(para, train_label, test_x, width){
  esti <- (mx_1_esti(para, train_label, test_x * (1 + width), width) - mx_1_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

mx_3_esti <- function(para, train_label, test_x, width){
  esti <- (mx_2_esti(para, train_label, test_x * (1 + width), width) - mx_2_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

p_1(test_x)

pesti <- phat(para$n, para$h, test$X ,train_label$X, k)
p1esti <- p_1_esti(para, train_label, test_x, width)

p2esti <- p_2_esti(para, train_label, test_x, 0.01)

p3esti <- p_3_esti(para, train_label, test_x, width)

mxesti <- mhat(para$n, para$h, test$X ,train_label$X, train_label$Y, k)
mx1esti <- mx_1_esti(para, train_label, test_x, width)
mx2esti <- mx_2_esti(para, train_label, test_x, width)
mx3esti <- mx_3_esti(para, train_label, test_x, width)


C1 <- (mx1esti*p1esti/p(x) + mx_2(x)*p_2(x)/p(x) + mx_3(x)*p_1(x)/p(x) - p_3(x)*p_1(x)/p(x)^2*mx_2(x))/4
C2 <- rk$value*sigmaep^2/p(x)
C3 <- sigmak$value * mx_1(x)^2/p(x)

mse <- function(hg){
  value <- C1 * hg[1]^4 * hg[2]^4 + C2/(para$n * hg[1]) + C3 * hg[1]^4 /(para$m * hg[2]^2)
  return(value)
}

optim(c(0.1,0.5),mse)


p_2 <- Deriv(p_1)
p_3 <- Deriv(p_2)
mx_2 <- Deriv(mx_1)
mx_3 <- Deriv(mx_2)
p_2 <- Deriv(p_1)
p_3 <- Deriv(p_2)
mx_2 <- Deriv(mx_1)
mx_3 <- Deriv(mx_2)

p_1_esti <- function(para, train_label, test_x, width){
  esti <- (phat(para$n, para$h, test_x*(1 + width),train_label$X, k) - phat(para$n, para$h, test_x*(1 - width),train_label$X, k))/(test_x*2*width)
  return(esti)
}

p_2_esti <- function(para, train_label, test_x, width){
  esti <- (p_1_esti(para, train_label, test_x * (1 + width), width) - p_1_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

p_3_esti <- function(para, train_label, test_x, width){
  esti <- (p_2_esti(para, train_label, test_x * (1 + width), width) - p_2_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

mx_1_esti <- function(para, train_label, test_x, width){
  esti <- (mhat(para$n, para$h, test_x*(1 + width),train_label$X, train_label$Y, k) - mhat(para$n, para$h, test_x*(1 - width),train_label$X, train_label$Y, k))/(test_x*2*width)
  return(esti)
}

mx_2_esti <- function(para, train_label, test_x, width){
  esti <- (mx_1_esti(para, train_label, test_x * (1 + width), width) - mx_1_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

mx_3_esti <- function(para, train_label, test_x, width){
  esti <- (mx_2_esti(para, train_label, test_x * (1 + width), width) - mx_2_esti(para, train_label, test_x *(1 - width), width))/(test_x*2*width)
  return(esti)
}

mse <- function(hg){
  value <- C1 * hg[1]^4 * hg[2]^4 + C2/(para$n * hg[1]) + C3 * hg[1]^4 /(para$m * hg[2]^2)
  return(value)
}

hg_optsearch <- function(para, train_label, width){
  pesti <- phat(para$n, para$h, test$X ,train_label$X, k)
  p1esti <- p_1_esti(para, train_label, test_x, width)
  p2esti <- p_2_esti(para, train_label, test_x, width)
  p3esti <- p_3_esti(para, train_label, test_x, width)
  
  mxesti <- mhat(para$n, para$h, test$X ,train_label$X, train_label$Y, k)
  mx1esti <- mx_1_esti(para, train_label, test_x, width)
  mx2esti <- mx_2_esti(para, train_label, test_x, width)
  mx3esti <- mx_3_esti(para, train_label, test_x, width)
  
  
  C1 <- (mx1esti*p1esti/pesti + mx2esti*p2esti/pesti + mx3esti*p1esti/pesti - p3esti*p1esti/pesti^2*mx2esti)/4
  C2 <- rk$value*sigmaep^2/pesti
  C3 <- sigmak$value * mx1esti^2/pesti
  
  mse <- function(hg){
    value <- C1 * hg[1]^4 * hg[2]^4 + C2/(para$n * hg[1]) + C3 * hg[1]^4 /(para$m * hg[2]^2)
    return(value)
  }
  
  hgopt <- optim(c(0.1,0.5),mse)$par
  return(hgopt)
}



