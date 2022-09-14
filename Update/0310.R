setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # 여기서 ㅊcproxSortedL1을 이용해서 함수를 생성하고, 읽어옴
library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid
library(readxl)
library(quadprog)
# 앞으로 데이터 셋 불러올때 일반화 해서 불러오기 위한 작업
data_set <- c("DOW30", "DAX30", "SP100", "FTSE100", "FTSE250", "SP500")

number <- 1
title <- paste("data/",data_set[number],".xlsx", sep="")
print(title)
data1 <- read_excel(title ,sheet = 1, col_names = TRUE)
data1 <- as.matrix(data1[,seq(from=1,to=ncol(data1),2)])

# 결측치 제거
h <- c()
for ( i in 1:ncol(data1)) { if (sum(is.na(data1[,i])>0)) h <- append(h, i) }
data1 <- data1[,-h] # 결측치가 있는 열의 경우 삭제
dim(data1)

############## 위는 더 이상 사용하지 않음


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
rb <- apply(R,mean,MARGIN = 1) # 자산별 250일 평균 계산
# rb <- (R%*%rep(1,window))/window # 같은값 계산하는 코드

# simulation data용
rb <- mu
R <- t(RR)
rp <- 0.01  #임의로 지정
window <- ncol(R)

# rp : 목표 수익률 -> 포트폴리오 평균 수익률을 써야하나..? -> 사실 이건 hyperparameter가 맞긴함!
#rp_prep <- rep(0,window)
#for( i in 2:nrow(base)) rp_prep[i-1] <- (data1[start:end,ncol(base)][i]-data1[start:end,ncol(base)][i-1])/data1[start:end,30][i-1]
#rp <- mean(rp_prep) # 기간동안 포트폴리오 평균 수익률을 사용
rp <- 0.01 # 임의로 지정하는 경우 

########## lambda sequence 생성

p <- length(rb) # 고려하는 자산의 개수
ll <- 30 # 생성할 람다 시퀀스 개수
q <- 0.01
qi <- rep(0,p)
for (i in 1:p) qi[i] <- i*(q/(2*p))
aa <- lseq(from = 10^-5, to = 10^2, length=ll)/qnorm(1-qi[1]) # 1까지로 설정해서 봄 !!(원래는 10^2까지인것 까먹으면 안됨)
lambda <- matrix(0,p,ll)
for ( j in 1:length(aa)) {
  for ( i in 1:p) {
    lambda[i,j] <- aa[j]*qnorm(1-qi[i]) # lambda[i,j] <- aa[j] : Lasso ver.
  }
}



lamb_seq(30)

#####################
#jj <- rep(0, ncol(lambda))
#GG <- rep(0,5500)
#GG <- c('G1', 'G2', 'G3')
ww <-  matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
rp <- 0.005
t0 <- proc.time()[3]
for (nth in 1:ncol(lambda)){
  # algorithm
  # 초기값 지정
  alpha <- rep(0, p) # theta는 shortfall에서 나온 새로운 모수
  beta <- theta <- 0
  gamma <- z <- rep(0,window)
  w <- v <- rep(0,p) # v <- rep(0,k)
  #v <- rep(1/p,p) #초기값은 다 0으로 놓긴했었음
  e <- rep(1,p)
  II <- diag(p)
  q <- 0.1 # quantile값
  KK <- floor(q*window) # big K
  kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
  eta <- 4 #0.8 #0.7 #0.5 #0.03 #70 #100 #5 #30 #0.01 #2 # 10 #0.02 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
  j <- 0
  tau <- 10^(-5) # 이거에 따라서도 차이가 심함 check 필요!!
  lambda_1 <- lambda[,nth]
  
  while (TRUE) {
    #GG1 <- GG2 <- GG3 <- rep(0,100)
    w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
    v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
    z <- proxZ(R,w,gamma,eta,kvec)
    #z <- fastproxZ(R,w,gamma,eta,kvec)
    alpha <- alpha + c(eta)*(w-v)
    beta <- beta + eta*(t(e)%*%w - 1)
    theta <- theta + eta*(t(w)%*%rb-rp)
    gamma <- gamma + c(eta)*(z - t(R)%*%w) # ; gamma2 <- gamma2 + c(eta)*(z2 - t(R)%*%w)
    
    
    # dual gap 계산 및 체크 ## 아닌것 같음!
    #G1 <- abs(-t(alpha + as.vector(beta%*%e))%*%w + beta + sum(lambda_1*sort(abs(v), decreasing = TRUE)))
    #G2 <- abs(-t(alpha + as.vector(beta%*%e) + as.vector(theta%*%rb))%*%w + beta + sum(lambda_1*sort(abs(v), decreasing = TRUE))) # 계산 필요
    G3 <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 # primal-dual gap으로 계산해야함!
    #
    #if ((G1 < tau)|(G2 < tau)|(G3 < tau)) {
    #  print(GG[which(c(G1,G2,G3)==min(c(G1,G2,G3)))]) 
    #  break
    
    if (G3 < tau) break
    
    j <- j+1
    #GG1[j] <-G1
    #GG2[j] <-G2
    #GG3[j] <-G3
    if (j>=5000) break # 반복수도 생각해야함..
  }
  ww[,nth] <- w
  print(nth)
  if(nth == ncol(lambda)) { real_viz(ww, eta, rp, tau); t0 <- proc.time()[3]-t0 }
}



simul_viz(ww, eta)
real_viz(vv, eta, rp, tau)



# visualization function
real_viz <- function(ww, month, eta, rp, tau){
  plot(ww[1,], ylim=c(min(ww),max(ww)),type='l', lwd=1, main=paste("eta = ",eta," month = ",month," rp = ",rp, "epsil = ", tau)) # round(eta,3)))
  for(i in 2:nrow(ww)) lines(ww[i,], type='l') #, col=ddd[i])
}

simul_viz <- function(ww, eta) {
  plot(ww[1,], ylim=c(min(ww),max(ww)),type='l', col=1, lwd=1, main=paste("data = simulation, eta = ", eta))
  for(i in 2:4) lines(ww[i,], type='l', col=1, lwd=1)
  for(i in 5:8) lines(ww[i,], type='l', col=2, lwd=1)
  for(i in 9:12) lines(ww[i,], type='l', col=3, lwd=1)
} 
 
plot(ww[1,], ylim=c(-1,1), type='l', lwd=1)#, main=paste("data = ", Dow, ", phi = ", phi, ", eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) points(ww[i,],type='l') #, col=ddd[i])



if ( (uu < hh)|(pp<hh|(yy<hh))) print(min(uu,pp,yy))

vv <- ww # 현재는 sp100, eta=0.02



plot(ww[,40], col=2)
points(ww[5:8,41], col=3, lwd=3)
points(ww[9:12, 41], col=1, lwd=2)

abline(v=4.5)
abline(v=8.5)


plot(ww[1,], ylim=c(-1,1),type='l', lwd=1, main=paste("data = ", data_set[number], " eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) lines(ww[i,], type='l') #, col=ddd[i])


# 평가지표 생성
# 다음 한달간의 수익률
Nwindow <- 21
Nstart <- end+1
Nend <- Nstart+Nwindow
base <- data1[Nstart:Nend,1:ncol(data1)-1] # 마지막 인덱스값은 제외하고 사용

N_daily_return <- matrix(0,Nwindow, dim(base)[2])
for ( j in 1:ncol(base)) {
  for ( i in 2:nrow(base)) {
    N_daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
  }
}
NxtR <- N_daily_return%*%ww
dim(NxtR)
max(colMeans(NxtR)) # 각 weight sequence에 대한 다음 30일의 평균수익률

vv <- ww # DAX30 , eta = 0.01 저장 중 (람다 30개)
vv2 <- ww # DAX30, eta = 0.01 람다 50개

plot(ww[1,1:30], ylim=c(-1,1),type='l', lwd=1, main=paste("data = ", data_set[number], " eta = ",eta)) # round(eta,3)))
for(i in 2:nrow(ww)) lines(ww[i,1:30], type='l') #, col=ddd[i])
