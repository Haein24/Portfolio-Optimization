# simuation 데이터용 코드 : 

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
RR <- cbind(RR,rep(0.1, 50))

############  시뮬 데이터 생성 완료 ############
# 시뮬레이션 데이터용 람다 함수
lamb_seq_sm <- function(ll) { # ll : 생성하려는 람다 시퀀스의 개수
  p <<- ncol(data1)-1 # 인덱스 빼고 고려하는 자산의 개수
  q <- 0.01
  qi <- rep(0,p)
  for (i in 1:p) qi[i] <- i*(q/(2*p))
  aa <- lseq(from = 10^-5, to = 10^2, length=ll)/qnorm(1-qi[1])
  lambda <<- matrix(0,p,ll)
  for ( j in 1:length(aa)) {
    for ( i in 1:p) {
      lambda[i,j] <<- aa[j]*qnorm(1-qi[i])
    }
  }
}


daily_return <- data1 <- RR 
lamb_seq_sm(30)
lambda <- cbind(lambda, rep(lambda[1,17],nrow(lambda)))

#################
n_core = detectCores() # 사용가능한 코어가 몇개인지 확인
cl = makeCluster(n_core-3) # 사용가능한 코어중 하나를 빼놓고 할당
cl
registerDoParallel(cl)

##################
# 기본 input 정의
window <- 50
start <- 1
Nwindow <- 21 # monthly update
lenmonth <- 1 # 몇 달치
w_tilde <- rep(0,p)

# mu <- 0.05 # mu = 0.1
########## 최종 결과를 저장할 데이터 프레임 생성
results <- data.frame("month"=0, "lambda"=0,"SLOPE_30"=0, "return_MV"=0,"SLOPE_17"=0, "lasso_17"=0,"index"=0, "EW"=0)
#results_turnover <- data.frame("month"=0, "lambda"=0,"Auto"=0, "SLOPE"=0, "lasso"=0) 
results_30 <- results_MV <- results_w_las <- results_w_17 <- matrix(0,lenmonth,ncol(data1)-1) 

#active_diff <- matrix(0, lenmonth, ncol(daily_return)-1) # setdifference 체크용

#t0 <- proc.time()[3] # 시간 측정
for (month in 1:lenmonth) {
  end <- start+window-1
  
  R <- t(daily_return[start:end,-ncol(daily_return)]) # dimension : p*T / 인덱스 수익률은 빼고 진행
  
  # r-bar
  rb <- apply(R,mean,MARGIN = 1) # 자산별 250일 평균 계산

  
  ww <- foreach(lambda=iter(lambda, by='column'), .combine=cbind,.packages = c("quadprog","CVXR")) %dopar% { #, .export=c("proxSortedL1", "proxZ")) %dopar% { #, .combine=cbind) {
    dyn.load('~/Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
    alpha <- rep(0, p)
    beta <- theta <- 0
    gamma <- z <- rep(0,window)
    w <- v <- rep(0,p)
    e <- rep(1,p)
    II <- diag(p)
    q <- 0.3;     eta <- 0.15;     rp <- 0.005
    KK <- floor(q*window) # big K
    kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
    j <- 0
    tau <- 10^(-5) 
    
    #eq2 <- c(eta/2)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R))
    
    while (TRUE) {
      
      #if (month != 1) {
      #  w <- solve( c(eta)*(c(1+(2*mu/eta))*II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z+c(2*mu/eta)*w_tilde) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      #} else {
      #  w <- solve(c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      #}
      
      w <- solve(c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      v <- proxSortedL1(w+(1/eta)*alpha, lambda/eta)
      z <- proxZ(R,w,gamma,eta,kvec)
      alpha <- alpha + c(eta)*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%rb-rp)
      gamma <- gamma + c(eta)*(z - t(R)%*%w)
      
      
      # dual gap 계산 및 체크 ## 아닌것 같음!
      Gap <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 # primal-dual gap으로 계산해야함!
      
      if (Gap < tau) break
      
      j <- j+1
      if ( j >= 1000 ) break
    }
    return(w) # 결과 내보내기
  }

  # 시각화
  plot(ww[1,1:30], ylim=c(min(ww),max(ww)), type='l', col=1, lwd=2, main="data = simulation")
  for(i in 2:4) lines(ww[i,1:30], type='l', col=1, lwd=2)
  for(i in 5:8) lines(ww[i,1:30], type='l', col=2, lwd=2)
  for(i in 9:12) lines(ww[i,1:30], type='l', col=3, lwd=2)
  

}

stopImplicitCluster()


