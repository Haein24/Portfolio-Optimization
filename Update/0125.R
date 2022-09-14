# 0125
# Dow30 데이터 불러오기
library(readxl)

setwd("Desktop/portfolio_selection")
Dow30 <- read_excel("Dow30.xlsx",sheet = 1, col_names = TRUE)
Dow30 <- Dow30[,seq(from=1,to=ncol(Dow30),2)] # 멀쩡한 데이터 만들기,,

# 기본 input 정의 # window : 250
data <- Dow30[1:100,1:30]
data <- as.matrix(data) # dataframe을 넣기 전에 matrix로 변형해서 넣어야함
mu <- mean_cov_comp(data, 0.3)$Mean
#sigma <- cov(data)
sigma <- mean_cov_comp(data, 0.3)$Covariance_matrix


############ Algorithm 적용 ################
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴

# lambda grid 생성
k <- ncol(data)
q <- 0.01
qi <- rep(0,k)
for (i in 1:k) qi[i] <- i*(q/(2*k))
aa <- seq(from = 10^-5, to = 10^2, length=100)/qnorm(1-qi[1])
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
  theta <- 10 # 상대위험회피정도 : theta > 0 (여기서는 임의의 값으로 두었음,,)
  eta <- 1 # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-6)
  lambda_1 <- lambda[,i]
  
  while (TRUE) {
    # update
    w <- solve(theta*sigma + eta*(II+e%*%t(e)))%*%(t(mu)-alpha-as.vector(beta%*%e)+eta*(v+e))
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
