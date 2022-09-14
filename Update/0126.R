# 0126
# Daily Return matrix 생성
library(readxl)
Dow30 <- read_excel("Dow30.xlsx",sheet = 1, col_names = TRUE)
Dow30 <- as.matrix(Dow30[,seq(from=1,to=ncol(Dow30),2)])

# 기본 input 정의 # window : 250
window <- 250
base <- Dow30[1:250,1:31]

daily_return <- matrix(0,dim(base)[1], dim(base)[2])
for ( j in 1:ncol(base)) {
  for ( i in 2:nrow(base)) {
    daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i,j]
  }
}

mu <- as.vector(mean_cov_comp(daily_return[,1:30])$Mean)
sigma <- mean_cov_comp(daily_return[,1:30])$Covariance_matrix
la <- mean_cov_comp(daily_return[,1:30])$la


# 위험회피계수 계산(RRA : relative risk aversion) :   
w <- rep(1/30,30) # 초기값 기준임..
rf <- mean(daily_return[,31])# 무위험 자산 수익률 -> 여기서는 index의 수익률을 사용해보자..!
risk <-  t(w)%*%sigma%*%w
return <- sum(w*mu)
theta <- ( return - rf ) / risk # 위험회피계수
theta <-as.numeric(abs(theta)) # 논문에는 >0으로 되어있는데 음수가 나와서 여기선 임의로 절대값 씌워서 사용


##################################
library(emdbook)

# k <- ncol(daily_return)-1
k <- length(mu)
q <- 0.01 # 논문에서 주어진 값
qi <- rep(0,k)
for (i in 1:k) qi[i] <- i*(q/(2*k))
aa <- lseq(10^-5, 10^2, 100)/qnorm(1-qi[1])
# aa <- exp(seq(from = log(10^-5), to = log(10^2), length=100))/qnorm(1-qi[1]) # lambda sequence
lambda <- matrix(0,k,100)
for ( j in 1:length(aa)) {
  for ( i in 1:k) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}

# algorithm
ww <- vv <- matrix(0, k,ncol(lambda))
jj <- rep(0,ncol(lambda))

for( i in 1:ncol(lambda)) {
  # 초기값 지정
  alpha <- rep(0, k)
  beta <- 0
  w <- rep(0,k)
  v <- rep(1/k,k)
  e <- rep(1,k)
  II <- diag(k)
  theta <- 3 #theta # 상대위험회피정도 : theta > 0 (미국기준 range가 보통 -0.5에서 3인듯!)
  eta <- 0.1 # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-6)
  lambda_1 <- lambda[,i]
  
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
  ww[,i] <- w
  vv[,i] <- v
  jj[i] <- j
}

#################
plot(1:100, ww[1,], ylim=c(min(ww), max(ww)), type='l', main='eta = 0.1 / theta = 3 / j >= 10000')
for ( i in 2:30) lines(ww[i,], type='l', col=i+1)
for ( i in 2:4) lines(ww[i,], type='l', col=1)
for ( i in 5:8) lines(ww[i,], type='l', col=2)
for ( i in 9:12) lines(ww[i,], type='l', col=3)

