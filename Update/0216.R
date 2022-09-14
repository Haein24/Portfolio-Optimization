dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴

library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid
library(readxl)

simul <- function(number) {
  title <- paste("data/",data_set[number],".xlsx", sep="")
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
  start <- 100
  end <- start+window
  base <- data1[start:end,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용
  
  daily_return <- matrix(0,window, dim(base)[2])
  for ( j in 1:ncol(base)) {
    for ( i in 2:nrow(base)) {
      daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
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
  aa <- lseq(from = 10^-5, to = 10, length=ll)/qnorm(1-qi[1]) # lseq 사용해야함 / real data에 대해서는 범위가 다름
  lambda <- matrix(0,k,ll)
  for ( j in 1:length(aa)) {
    for ( i in 1:k) {
      lambda[i,j] <- aa[j]*qnorm(1-qi[i])
    }
  }
  
  #####################
  ww <<- vv <- matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
  gg <- jj <- rep(0, ncol(lambda))
  for ( z in 1:ncol(lambda) ){
    # algorithm
    # 초기값 지정
    alpha <- rep(0, k)
    beta <- 0
    w <- rep(0,k)
    v <- rep(0,k)
    e <- rep(1,k)
    II <- diag(k)
    phi <- 4 #15 #40 #k # 상대위험회피정도 : phi > 0 (여기서는 임의의 값으로 두었음,,)
    eta <- 0.02  #15 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
    j <- 0
    tau <- 10^(-5)
    lambda_1 <- lambda[,z]
    
    while (TRUE) {
      # update
      w <- solve(phi*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
      v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
      alpha <- alpha + eta*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      
      # dual gap 계산 및 체크
      G <- abs(-t(alpha + as.vector(beta%*%e))%*%w + beta + sum(lambda_1*sort(abs(v), decreasing = TRUE)))
      #G <- t(w-v)%*%(w-v) + (t(e)%*%w-1)^2 
      if (G < tau) break
      
      # 다음 단계 넘어가기 전에 전 단계 값을 저장
      #preA <- alpha
      #preB <- beta
      #preW <- w
      #preV <- v
      
      j <- j+1
      if (j>=10000) break
    }
    ww[,z] <<- w
    print(z)
  }
  #write.csv(ww, file=paste("result(3,10)/result_",data_set[number],'.csv',sep=""), row.names=FALSE)
  plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=1, main=paste("data = ", title, ", phi = ", phi, ", eta = ",eta))
  for(i in 2:nrow(ww)) lines(ww[i,], type='l') #, col=ddd[i])
}


simul(1)

for ( number in 2:5 ) {#length(data_set) ) {
  simul(number)
}
ww6 <- ww


### 시각화
# write.csv(ww, file='weights/FTSE250.csv', row.names = FALSE)
ww <- as.matrix(read.csv("result(k,15)/result_SP500.csv"))
plot(ww[1,], ylim=c(min(ww),max(ww)), type='l', lwd=1)#, main=paste("data = ", Dow, ", phi = ", phi, ", eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) lines(ww[i,], type='l') #, col=ddd[i])

ncol(ww)

table(round(a,5))

# 그룹화
length(table(rr))

jj <- round(ww[,8],3)
length(table(jj))
       
c <- 1
dd <- seq(0, length(jj))
for (i in 1:length(unique(jj))) {
  for ( j in 1:length(jj)) {
    if (jj[j]==unique(jj)[i]) dd[j] <- c
  }
  c <- c+1
}
