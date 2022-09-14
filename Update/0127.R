# 기본 작업
setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴
library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid

###############################################
# simulation data generating
set.seed(123)
t <- 50 ; r<-3 ; k <- 12 

# risk factors
risk_factors <- rmvnorm(50, rep(0,r), diag(r))

# vectors of error term
epsilon<- rmvnorm(50, rep(0,k), 0.05*diag(k)) ## 논문에서는 diag(r)로 나와있음,,! 해결 필요

# loading matrix
col1 <- c(0.77,0.64,0)
col2 <- c(0.9,0,-0.42)
col3 <- c(0, 0.31,0.64)

BB <- matrix(0,3,12)
for ( i in (1:4)) {
  BB[,i] = t(t(col1))
  BB[,i+4] = t(t(col2))
  BB[,i+8] = t(t(col3))
}
BB

# return matrix
RR <- risk_factors%*%BB + epsilon

# mu vector
mu <- apply(RR, 2, mean)
mu

# covariance matrix
sigma <- t(BB)%*%BB + 0.05*diag(k)
sigma_hat <- cov(RR)

#############################################################
# lambda grid 생성
k <- length(mu)
q <- 0.01
qi <- rep(0,k)
for (i in 1:k) qi[i] <- i*(q/(2*k))
aa <- lseq(from = 10^-5, to = 10^2, length=30)/qnorm(1-qi[1]) # lseq 사용해야함
lambda <- matrix(0,k,30)
for ( j in 1:length(aa)) {
  for ( i in 1:k) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}

###########################################
# weight 값 도출 / 함수로 만들려면 필요한 input : mu, sigma, theta, eta, lambda

ww <- vv <- matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
jj <- rep(0, ncol(lambda))
for ( z in 1:100 ){
  # algorithm
  # 초기값 지정
  alpha <- rep(0, k)
  beta <- 0
  w <- rep(0,k)
  v <- rep(1/12,k)
  e <- rep(1,k)
  II <- diag(k)
  theta <- 40 # 상대위험회피정도 : theta > 0 (여기서는 임의의 값으로 두었음,,)
  eta <- 10 # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-6)
  lambda_1 <- lambda[,z]
  
  while (TRUE) {
    # update
    w <- solve(theta*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
    v <- proxSortedL1(w+(1/eta)*alpha, lambda_1)
    alpha <- alpha + eta*(w-v)
    beta <- beta + eta*(t(e)%*%w - 1)
    
    # dual gap 계산 및 체크
    G <- abs(-t(alpha + as.vector(beta%*%e))%*%w + beta + sum(lambda_1*sort(abs(v), decreasing = TRUE)))
    if (G < tau) break
    
    # 다음 단계 넘어가기 전에 전 단계 값을 저장
    preA <- alpha
    preB <- beta
    preW <- w
    preV <- v
    
    j <- j+1
    if (j>=10000) break
  }
  ww[,z] <- w
  vv[,z] <- v
  jj[z] <- j
}

# 그래프
plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=1, main=paste("data = SP500, theta = ", theta, ", eta = ", eta))
for(i in 2:4) lines(ww[i,], type='l', col=1, lwd=2)
for(i in 5:8) lines(ww[i,], type='l', col=2, lwd=2)
for(i in 9:12) lines(ww[i,], type='l', col=3, lwd=2)


# long-only 영역 표시하기
c <- rep(0,ncol(ww))
for ( i in 1:ncol(ww)) {
  if (sum(ww[,i]>=0) == nrow(ww)) c[i] <- ncol(ww)-i
}
for ( i in 2:ncol(ww)-2) {
  if (c[i+1]==0) c[i] <- 0
}
abline(v=ncol(ww)-max(c))


plot(1:40, ww[1,1:50], ylim=c(-0.3,0.3), type='l', col=1, lwd=1, main="data = DOW30, theta =15, eta=1")
for(i in 2:30) lines(ww[i,1:50], type='l', col=i)
# THETA : 1번째 iteration의 각 weight values' range를 결정
# 상대위험회피계수 높을 수록 위험을 기피한다는 뜻이므로, 너무 큰 값의 매수나 매도가 발생하지 않음..?
# eta : 수렴하는 속도를 조절 : eta값이 클수록 더 큰 람다가 주어져야 그룹화 되는 양상ㅇ
# 오늘 목표 : 
# 1) 실제 데이터에는 eta = 1.5로 고정하고, theta value를 계산할 수 있는 수식이 있을지 더 develop
# 2) 실제로 섹터가 나뉘어있던 데이터를 기반으로 해보기,,(SGLasso 논문 참고)
