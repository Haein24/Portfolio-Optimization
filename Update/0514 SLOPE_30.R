data_prep <- function(data_name, sheet) { # number : 몇번째 데이터셋을 사용할 것인지 지정 
  title <- paste("data/",data_name,".xlsx", sep="")
  print(title)
  data1 <<- read_excel(title ,sheet = sheet, col_names = TRUE) # 두번째 시트사용
  data1 <<- as.matrix(data1[,seq(from=1,to=ncol(data1),2)])
  
  # 결측치 제거
  h <- c()
  for ( i in 1:ncol(data1)) { if (sum(is.na(data1[,i])>0)) h <- append(h, i) }
  data1 <- data1[,-h] # 결측치가 있는 열의 경우 삭제
  data1 <<- data1[order(1:nrow(data1), decreasing=TRUE),] # 역순으로 정렬해야 시간의 순서대로 됨
  print(paste("nrow : ", dim(data1)[1], "ncol : ", dim(data1)[2]))
}

daily_return_prep <- function(data_name) { # number : 몇번째 데이터셋을 사용할 것인지 지정 
  daily_return <<- matrix(0,nrow(data1)-1, dim(data1)[2])
  for ( j in 1:ncol(data1)) {
    for ( i in 2:nrow(data1)) {
      daily_return[i-1,j] <<- (data1[i,j] - data1[i-1,j])/data1[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
    }
  }
}


proxZ <- function(R, w, gamma, eta, kvec) {
  dt <- data.frame(cbind(t(R)%*%w, gamma))
  colnames(dt) <- c('y','gamma')
  idx <- order(dt$y) # original order를 저장
  y <- sort(dt$y)
  
  n_gamma <- dt[idx,]$gamma 
  
  #A <- matrix(0, ncol(R), ncol(R))
  A <- diag(1,ncol(R), ncol(R))
  for ( i in 1:ncol(R)) { A[i, i-1] <- -1 }
  A <- A[2:ncol(R),]
  
  bvec <- rep(0,ncol(R)-1)
  
  dvec <- ( kvec + as.vector(eta)*y - n_gamma)
  Dmat <- diag(eta, ncol(R), ncol(R))
  
  # 계산
  res <- solve.QP(Dmat,dvec, t(A), bvec)$solution
  z[idx] <- res # z를 원래 인덱스로 되돌리기
  return(z)
}
fold_w_proc <- function(R){
  rb <<- apply(R,mean,MARGIN = 1) # 자산별 평균 계산
  window <- ncol(R)
  
  kvec <- c(rep(1/KK, KK), rep(0, window-KK))
  ww<<-  matrix(0, nrow = nrow(lambda), ncol = ncol(lambda)) # reset ww & rww
  for (nth in 1:ncol(lambda)){
    # 초기값으로 reset
    alpha <- rep(0, p) # theta는 shortfall에서 나온 새로운 모수
    beta <- theta <- 0
    gamma <- z <<- rep(0,window) # ""window""가 fold내 데이터수 뺀 값이 되어야함
    w <- v <- rep(0,p)
    j <- 0
    lambda_1 <- lambda[,nth]
    
    while (TRUE) {
      w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      v <- proxSortedL1(w+(1/eta)*alpha, lambda_1/eta)
      z <<- proxZ(R,w,gamma,eta,kvec)
      alpha <- alpha + c(eta)*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%rb-rp)
      gamma <- gamma + c(eta)*(z - t(R)%*%w) # ; gamma2 <- gamma2 + c(eta)*(z2 - t(R)%*%w)
      
      # dual gap 계산 및 체크 ## 아닌것 같음!
      Gap <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2
      
      if (Gap < epsilon) break
      j <- j+1
      if (j>=8000) break
    }
    ww[,nth] <<- w
    # print(nth)
  }
}
lamb_seq <- function(ll) { # ll : 생성하려는 람다 시퀀스의 개수
  p <<- ncol(data1)-1 # 인덱스 빼고 고려하는 자산의 개수
  q <- 0.01
  qi <- rep(0,p)
  for (i in 1:p) qi[i] <- i*(q/(2*p))
  aa <- lseq(from = 10^-5, to = 10, length=ll)/qnorm(1-qi[1])
  lambda <<- matrix(0,p,ll)
  for ( j in 1:length(aa)) {
    for ( i in 1:p) {
      lambda[i,j] <<- aa[j]*qnorm(1-qi[i])
    }
  }
}
setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # cproxSortedL1을 이용해서 함수를 생성
real_viz <- function(ww, month, w_idx){
  plot(ww[1,], ylim=c(min(ww),max(ww)),type='l', lwd=1, main=paste(" month = ",month," selected_idx = ", w_idx,"th")) # round(eta,3)))
  for(i in 2:nrow(ww)) lines(ww[i,], type='l')
} # 시각화용 함수

turnover <- function(results_w) {
  turnover_prep <- rep(0, nrow(results_w)-1)
  for (i in 1:(nrow(results_w)-1)) {
    turnover_prep[i] <- sum(abs(results_w[i+1,]-results_w[i,]))
  }
  print(mean(turnover_prep))
  #turnover_res <<- mean(turnover_prep)
}

library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid
library(readxl)
library(quadprog)
library(doParallel)

######### 병렬작업 준비
n_core = detectCores() # 사용가능한 코어가 몇개인지 확인
cl = makeCluster(n_core-1) # 사용가능한 코어중 하나를 빼놓고 할당(코어 5개 사용)
cl
registerDoParallel(cl)


########## 데이터 준비
data_prep("DOW30", 1)
daily_return_prep(data1)

########## lambda sequence 생성
lamb_seq(30) # 람다 시퀀스 생성 함수 (몇개를 생성할 것인지 지정)

lambda <- cbind(lambda, rep(lambda[1,17],nrow(lambda))) # Lasso를 위해서 31번째에 라쏘를 위한 람다값을 붙여줌

##################
# 기본 input 정의
window <- 250 # 11개월치 사용
start <- 1
Nwindow <- 21 # monthly update
lenmonth <- 12*6 # 몇 년치 데이터를 다룰 것인지

########## 최종 결과를 저장할 데이터 프레임 생성
results_40 <- data.frame("month"=0, "lambda"=0,"less_turnover"=0)
#results_turnover <- data.frame("month"=0, "lambda"=0,"Auto"=0, "SLOPE"=0, "lasso"=0) 
#results_w <- results_w_17 <- results_w_las <- matrix(0,lenmonth,ncol(data1)-1)
results_w_40 <- matrix(0,lenmonth,ncol(data1)-1)

#active_diff <- matrix(0, lenmonth, ncol(daily_return)-1) # setdifference 체크용

t0 <- proc.time()[3] # 시간 측정
for (month in 1:lenmonth) {
  end <- start+window-1
  
  R <- t(daily_return[start:end,-ncol(daily_return)]) # dimension : p*T / 인덱스 수익률은 빼고 진행
  
  # r-bar
  rb <- apply(R,mean,MARGIN = 1) # 자산별 250일 평균 계산
  
  
  ww <- foreach(lambda=iter(lambda, by='column'), .combine=cbind,.packages = "quadprog") %dopar% { #, .export=c("proxSortedL1", "proxZ")) %dopar% { #, .combine=cbind) {
    dyn.load('~/Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
    alpha <- rep(0, p)
    beta <- theta <- 0
    gamma <- z <- rep(0,window)
    w <- v <- rep(0,p)
    e <- rep(1,p)
    II <- diag(p)
    q <- 0.1
    KK <- floor(q*window) # big K
    kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
    eta <- 4 
    j <- 0
    rp <- 0.005
    tau <- 10^(-5) 
    
    while (TRUE) {
      w <- solve( c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
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
      if (j>=5000) break
    }
    return(w) # 결과 내보내기
  }
  
  #w_idx <- which(max(colMeans(t(R)%*%round(ww[,-ncol(ww)],2)))==colMeans(t(R)%*%round(ww[,-ncol(ww)],2))) # 마지막 weights(for lasso) 제외하고 사용
  ## 새로운 기준 : activeset이 30~40%인 ww중에 이전 ww와 turnover가 가장 작은 것으로 선택
  
  w_idx <- which(colSums(round(ww[,-ncol(ww)],2)!=0)/p >= 0.3 & 0.5 >= colSums(round(ww[,-ncol(ww)],2)!=0)/p)
  if (length(w_idx) >= 2 & month >=2) {
    j <- 1
    comp <- rep(0,length(w_idx))
    for(i in w_idx ) {
      comp[j] <- sum(abs(results_w_40[month-1,]-round(ww[,i],2))) # turnover 계산
      j <- j+1
    }
    w_idx <- w_idx[which.min(comp)] # turnover가 가장 작은데 에서 선택
  } else if (length(w_idx) >= 2 & month == 1) { w_idx <- min(w_idx) }
  
  www <- round(ww[,w_idx],2) # 평가에 사용할 weight 저장 (1%단위로 반올림해서 투자)
  
  
  # 평가
  # 결정된 람다를 적용해 다음 21간의 수익률을 계산 및 저장
  ###### 최종적으로 내보낼때는 평가 대상 : 인덱스 수익률 & Equally Weighted 데이터(데이터 프레임에 순서대로 저장)
  N_daily_return <- (data1[(end+Nwindow-1),]-data1[(end),])/data1[(end),]
  # 월수익 계산 및 저장
  nreturn <- N_daily_return[-ncol(data1)]%*%www
  #return_17 <- N_daily_return[-ncol(data1)]%*%round(ww[,17],2)
  #return_lasso <- N_daily_return[-ncol(data1)]%*%round(ww[,31],2)
  #idx <- N_daily_return[ncol(data1)]
  #ew <- N_daily_return[-ncol(data1)]%*%rep(1/p, p) # equally weighted portfolio mean return
  
  
  results_40 <- rbind(results_40, c(month, w_idx, nreturn))
  results_w_40[month,] <- www
  #results_w_17[month,] <- round(ww[,17],2)
  #results_w_las[month,] <- round(ww[,31],2)
  
  start <- start+Nwindow
  
  #######  진행상황 및 grouping 확인용 
  print(month)
  if (month == lenmonth ) { t0 <- proc.time()[3] - t0 }
}

stopImplicitCluster()


turnover(results_w_40)