setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴
library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid

# 앞으로 데이터 셋 불러올때 일반화 해서 불러오기 위한 작업
data_set <- c("DOW30", "DAX30", "SP100", "FTSE100", "FTSE250", "SP500")

# Daily Return matrix 생성
library(readxl)
#setwd('Desktop/portfolio_selection/data')
title <- paste("data/",data_set[3],".xlsx", sep="")
print(title)
data1 <- read_excel(title ,sheet = 1, col_names = TRUE)
data1 <- as.matrix(data1[,seq(from=1,to=ncol(data1),2)])

# 결측치 제거
h <- c()
for ( i in 1:ncol(data1)) { if (sum(is.na(data1[,i])>0)) h <- append(h, i) }
data1 <- data1[,-h] # 결측치가 있는 열의 경우 삭제
dim(data1)

# 기본 input 정의 # window : 250
window <- 250
base <- data1[500:750,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용

daily_return <- matrix(0,dim(base)[1], dim(base)[2])
for ( j in 1:ncol(base)) {
  for ( i in 2:nrow(base)) {
    daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i,j]
  }
}

mu <- as.vector(mean_cov_comp(daily_return[,1:ncol(base)])$Mean)
sigma <- mean_cov_comp(daily_return[,1:ncol(base)])$Covariance_matrix
la <- mean_cov_comp(daily_return[,1:ncol(base)])$la

########## lambda sequence 생성
k <- length(mu)
ll <- 100 # 생성할 람다 시퀀스 개수
q <- 0.01
qi <- rep(0,k)
for (i in 1:k) qi[i] <- i*(q/(2*k))
aa <- lseq(from = 10^-5, to = 10^2, length=ll)/qnorm(1-qi[1]) # lseq 사용해야함
lambda <- matrix(0,k,ll)
for ( j in 1:length(aa)) {
  for ( i in 1:k) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i])
  }
}

#####################
t0 <- proc.time()[3]
ww <- vv <- matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
gg <- jj <- rep(0, ncol(lambda))
for ( z in 1:ncol(lambda) ){
  # algorithm
  # 초기값 지정
  alpha <- rep(0, k)
  beta <- 0
  w <- rep(0,k)
  v <- rep(1/k,k)
  e <- rep(1,k)
  II <- diag(k)
  phi <- 40 # 상대위험회피정도 : phi > 0 (여기서는 임의의 값으로 두었음,,)
  eta <-  15 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-4)
  lambda_1 <- lambda[,z]
  
  while (TRUE) {
    # update
    w <- solve(phi*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
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
  gg[z] <- G
  print(z)
  # if (z == ncol(lambda)) t0 <- proc.time()[3] - t0
}


plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=2, main="data = simulation, phi = 40, eta = 15", xlab=xlb, ylab=ylb) # round(eta,3)))
for(i in 2:4) lines(ww[i,], type='l', col=1, lwd=2)
for(i in 5:8) lines(ww[i,], type='l', col=2, lwd=2)
for(i in 9:12) lines(ww[i,], type='l', col=3, lwd=2)


plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=1, main=paste("data = ", title, ", phi = ", phi, ", eta = ",eta)) # round(eta,3)))
for(i in 2:ncol(base)) lines(ww[i,], type='l', col=1)

# long-only 영역 표시하기
c <- rep(0,ncol(ww))
for ( i in 1:ncol(ww)) {
  if (sum(round(ww[,i],5)>=0) == nrow(ww)) c[i] <- ncol(ww)-i
}
for ( i in 2:ncol(ww)-2) {
  if (c[i+1]==0) c[i] <- 0
}
abline(v=ncol(ww)-max(c), col=8)


## 매우작은 값은 0으로 만들어서 sparse함이 잘 드러나게 만듦
for ( i in 1:81) ww[i,abs(ww[i,])<10^-6] <- 0 

# distance
ii <- 8 # 몇번째 람다 시퀀스 이용할 것인지 지정
abline(v=ii, col=6) # 사용할 람다 시퀀스에 표시

distance <- dist(ww[,ii])
distanced <- hclust(distance, method = "single")
plot(distanced, hang=-1, main=paste(ii,"th lambda",sep=""))

groups <- cutree(distanced, k=24) # 이것도 기준이 있긴해야할텐데 맞지,,?
groups

# write.csv(ww, file='weights/FTSE250.csv', row.names = FALSE)



