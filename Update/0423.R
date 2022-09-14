data_prep("DOW30", 1)

# 기본 input 정의 # window : 250
window <- 250
start <- 10
end <- start+window
base <- data1[start:end,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용


daily_return <- matrix(0,window, dim(base)[2])
for ( j in 1:ncol(base)) {
  for ( i in 2:nrow(base)) {
    daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
  }
} # i 인덱스가 1~T까지 돌아감


R <- t(daily_return) # dimension : p*T

# r-bar
rb <- apply(R,mean,MARGIN = 1) 


lamb_seq(30)

#####################

ww <-  matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
rp <- 0.005 ######
t0 <- proc.time()[3]
for (nth in 1:ncol(lambda)){
  # algorithm
  # 초기값 지정
  alpha <- rep(0, p) # theta는 shortfall에서 나온 새로운 모수
  beta <- theta <- 0
  gamma <- z <- rep(0,window)
  w <- v <- rep(0,p)
  e <- rep(1,p)
  II <- diag(p)
  q <- 0.1 # quantile값 #######
  KK <- floor(q*window) # big K
  kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
  eta <- 4  ####### #0.8 #0.7 #0.5 #0.03 #70 #100 #5 #30 #0.01 #2 # 10 #0.02 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-5) # 이거에 따라서도 차이가 심함 check 필요!!
  lambda_1 <- lambda[,nth]
  
  while (TRUE) {
    
    w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
    v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
    z <- proxZ(R,w,gamma,eta,kvec)
    #z <- fastproxZ(R,w,gamma,eta,kvec)
    alpha <- alpha + c(eta)*(w-v)
    beta <- beta + eta*(t(e)%*%w - 1)
    theta <- theta + eta*(t(w)%*%rb-rp)
    gamma <- gamma + c(eta)*(z - t(R)%*%w) # ; gamma2 <- gamma2 + c(eta)*(z2 - t(R)%*%w)
    
    
    G3 <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 # primal-dual gap으로 계산해야함!
   
    if (G3 < tau) break
    
    j <- j+1
    
    if (j>=5000) break # 반복수도 생각해야함..
  }
  ww[,nth] <- w
  print(nth)
  if(nth == ncol(lambda)) { real_viz(ww, eta, rp, tau); t0 <- proc.time()[3]-t0 }
}