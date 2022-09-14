# lambda grid 생성
# 한 '행'이 하나의 람다 시퀀스를 지칭

k <- 12
jth <- seq(from = 10^-5, to = 10^2, length=100)
jlist <- 1:100
ilist <- 1:12
lambda_mat <- matrix(0, 100, 12)
for (j in jlist) {
  alpha <- jth[j]/qnorm(1-0.01/(2*k))
  for (i in ilist) {
    lambda_mat[j,i] <- alpha*(qnorm(1-(0.01*i)/(2*k)))
  }
}

# w 초기값 설정 ..?
w0 <- runif(11,-0.2, max = 0.2)
w0[12] <- 1-sum(w0)
w0


repp <- matrix(0, nrow = 101, ncol = 12)
jj <- rep(0,100)
repp[1,]<-rep(0,k)

for ( z in 2:nrow(lambda_mat)+1 ) {
  w <- v <- rep(0,k)
  alpha <- rep(0,k)
  beta <- 0
  j <- 0 # iteration counter
  e <- rep(1,length(mu))
  lambda <- lambda_mat[z-1,]
  eta <- 0.005
  theta <- 1 #0.03 # prefixed value # 찾아야하는 값..!!
  tau <- 1e-6 # 정의하기 : 현재 임의로 지정한 값임
  #gg <- 1 # 초기값을 크게 줘서 while문이 돌아가도록 설정..! 했는데 이값때문에 문제가 되는 듯..?! 대체 왜..?
  rm(gg)
  
  # updating 부분
  while ( TRUE) { # stopping rule 관련
    # 현재 값을 저장
    preW <- w
    preV <- v
    preA <- alpha
    preB <- beta
    
    # update
    w <- solve(theta*sigma + eta*(diag(k) + e%*%t(e)))%*%(mu-preA-as.vector(preB%*%e)+eta*(preV+e))
    v <- proxSortedL1(w+(1/eta)*preA, lambda)
    alpha <- preA + eta*(w-v)
    beta <- preB + eta*(t(e)%*%w -1)
    
    # gg값 업데이트
    # gg <- abs(t(preA)%*%(preW-preV)+preB*(t(e)%*%preW -1) + (eta/2)*(t(preW-preV)%*%(preW-preV) +(t(e)%*%preW -1)^2))
     gg <- abs(t(alpha)%*%(w-v)+beta*(t(e)%*%w -1) + (eta/2)*(t(w-v)%*%(w-v) +(t(e)%*%w -1)^2))
    # gg <- g(w,v, alpha, beta) # 논문에서 나온 값 이용한 경우
    
    # iteration +1
    j <- j+1
    
    if (j>=2000) break
    if (gg<=tau) break
  }
  repp[z,] <- w
  jj[z] <- j
}



#########################################################
mu <- mu
sigma <- sigma
k <- length(mu)

repp <- matrix(0, nrow = 101, ncol = 12)
jj <- rep(0,100)
repp[1,]<-rep(0,k)

for ( z in 2:nrow(lambda_mat)+1 ) {
  w <- v <- rep(1/12,k)
  alpha <- rep(2,k)
  beta <- 2
  j <- 0 # iteration counter
  e <- rep(1,length(mu))
  lambda <- lambda_mat[z-1,]
  eta <- 0.005
  theta <- 1 #0.03 # prefixed value # 찾아야하는 값..!!
  tau <- 1e-6 # 정의하기 : 현재 임의로 지정한 값임
  #gg <- 1 # 초기값을 크게 줘서 while문이 돌아가도록 설정..! 했는데 이값때문에 문제가 되는 듯..?! 대체 왜..?
  #rm(G)
  
  while ( TRUE ) {
    preW <- w
    preV <- v
    preA <- alpha
    preB <- beta
    
    w <- solve(sigma + eta*(diag(12)+e%*%(t(e))))%*%(mu-preA-as.vector(preB%*%e) + eta*(preV+e))
    #v <- sortedL1Prox(w + (1/eta)*preA, lambda/eta)
    v <- proxSortedL1(w + (1/eta)*preA, lambda/eta)
    alpha <- preA + eta*(w-v)
    beta <- preB + eta*(t(e)%*%w -1)
    
   # G <- t(alpha)%*%(w-v)+beta*(t(e)%*%w -1) + (eta/2)*(t(w-v)%*%(w-v) +(t(e)%*%w -1)^2) # 방향은 얘가 맞는 것 같은데 아닌 것 같음..
   G <-  -t(alpha + as.vector(beta%*%e))%*%w + beta + sum(lambda*sort(abs(v), decreasing=TRUE)) # 논문상
    
    j <- j+1
    if (G<=tau) break
    if (j>=2000) break
  }
  repp[z,] <- preW
  jj[z-1] <- j
}


#####################################
# 21.01.24
dyn.load('Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
source('Desktop/portfolio_selection/SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴

set.seed(123)
t <- 50 ; r<-3 ; k <- 12 

library(mvtnorm)
# risk factors
risk_factors <- rmvnorm(50, rep(0,r), diag(r))

# vectors of error term
epsilon<- rmvnorm(50, rep(0,k), 0.05*diag(k)) ## 논문에서는 diag(r)로 나와있음,,

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


# lambda grid 생성
q <- 0.01
qi <- rep(0,12)
for (i in 1:12) qi[i] <- i*(q/(2*k))
aa <- seq(from = 10^-5, to = 10^2, length=100)/qnorm(1-qi[1])
lambda <- matrix(0,12,100)
for ( j in 1:length(aa)) {
  for ( i in 1:12) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}



# algorithm
# 초기값 지정
alpha <- rep(0, k)
beta <- 0
w <- rep(0,k)
v <- rep(1/12,k)
e <- rep(1,k)
II <- diag(k)
theta <- 10 # 상대위험회피정도 : theta > 0 (여기서는 임의의 값으로 두었음,,)
eta <- 1 # eta > 0 (여기서는 임의의 값으로 두었음,,)
j <- 0
tau <- 10^(-6)
lambda_1 <- lambda[,1]

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

