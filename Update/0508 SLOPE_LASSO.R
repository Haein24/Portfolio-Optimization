###### 0504 : SLOPE & LASSO
# pre-preparation
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
setwd('~/Desktop/portfolio_selection')
dyn.load('SLOPE_code/cproxSortedL1.so')
source('SLOPE_code/proxSortedL1.R') # cproxSortedL1을 이용해서 함수를 생성
real_viz <- function(ww, month, w_idx){
  plot(ww[1,], ylim=c(min(ww),max(ww)),type='l', lwd=1, main=paste(" month = ",month," selected_idx = ", w_idx,"th")) # round(eta,3)))
  for(i in 2:nrow(ww)) lines(ww[i,], type='l')
} # 시각화용 함수

library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid
library(readxl)
library(quadprog)
library(doParallel)


#################
par(mfrow=c(4,1))

######### 병렬작업 준비
n_core = detectCores() # 사용가능한 코어가 몇개인지 확인
cl = makeCluster(n_core-3) # 사용가능한 코어중 하나를 빼놓고 할당(코어 5개 사용)
cl
registerDoParallel(cl)

results_lasso <- data.frame("month"=0, "lambda"=0, "Auto"=0,"SLOPE"=0, "lasso"=0,"index"=0, "EW"=0) # 최종 결과를 저장할 데이터 프레임 생성


########## 데이터 준비
data_prep("DOW30", 1)

########## lambda sequence 생성
lamb_seq(30) # 람다 시퀀스 생성 함수 (몇개를 생성할 것인지 지정)

lambda <- cbind(lambda, rep(lambda[1,17],29)) # Lasso를 위해서 31번째에 라쏘를 위한 람다값을 붙여줌

##################
# 기본 input 정의
window <- 250 # 11개월치 사용
start <- 1
Nwindow <- 21

t0 <- proc.time()[3] # 시간 측정
for (month in 1:12) {
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
  
  
  ww <- foreach(lambda=iter(lambda, by='column'), .combine=cbind,.packages = "quadprog") %dopar% { #, .export=c("proxSortedL1", "proxZ")) %dopar% { #, .combine=cbind) {
    dyn.load('~/Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
    alpha <- rep(0, p)
    beta <- theta <- 0
    gamma <- z <- rep(0,window)
    w <- v <- rep(0,p)
    e <- rep(1,p)
    II <- diag(p)
    q <- 0.1 # quantile값..
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
      if (j>=7000) break
    }
    return(w) # 결과 내보내기
  }
  
  T_daily_return <- daily_return
  
  w_idx <- which(max(colMeans(T_daily_return%*%round(ww[,-ncol(ww)],2)))==colMeans(T_daily_return%*%round(ww[,-ncol(ww)],2))) # 마지막 weights(for lasso) 제외하고 사용
  
  www <- round(ww[,w_idx],2) # 평가에 사용할 weight 저장 (1%단위로 반올림해서 투자)
  
  
  # 평가
  # 결정된 람다를 적용해 다음 21간의 수익률을 계산 및 저장
  Nstart <- end
  Nend <- Nstart+Nwindow
  base <- data1[Nstart:Nend,1:ncol(data1)]
  
  N_daily_return <- matrix(0,Nwindow, dim(base)[2])
  for ( j in 1:ncol(base)) {
    for ( i in 2:nrow(base)) {
      N_daily_return[i-1,j] <- (base[i,j] - base[i-1,j])/base[i-1,j]  # 지금까지 틀렸던거 찾아서 수정..
    }
  }
  
  if (daily_return[250,1]==N_daily_return[1,1]) break
  
  idx_return <- N_daily_return[,ncol(N_daily_return)]
  N_daily_return <- N_daily_return[,-ncol(N_daily_return)]
  
  ###### 최종적으로 내보낼때는 평가 대상 : 인덱스 수익률 & Equally Weighted 데이터(데이터 프레임에 순서대로 저장)
  nreturn <- mean(N_daily_return%*%www) # 평균 수익률
  return_17 <- mean(N_daily_return%*%round(ww[,17],2))
  return_lasso <- mean(N_daily_return%*%round(ww[,31],2))
  idx <- mean(idx_return) # index mean return
  ew <- mean(N_daily_return%*%rep(1/p, p)) # equally weighted portfolio mean return
  
  #######  진행상황 및 grouping 확인용 
  print(paste("month=",month,"& index=",w_idx, "(",sum(www),")"))
  #print(c(nreturn, idx, ew))
  print(c(nreturn, return_17,return_lasso, idx, ew))
  print(table(round(ww[,w_idx],2)))
  print(table(round(ww[,17],2)))
  print(table(round(ww[,31],2)))
  
  #results3 <- rbind(results3, c(month, w_idx, nreturn, idx, ew))
  results_lasso <- rbind(results_lasso, c(month, w_idx, nreturn, return_17, return_lasso, idx, ew))
  start <- start+21
  
  # 시각화
  real_viz(ww[,-ncol(ww)], month, w_idx) # 마찬가지로 
  abline(v=w_idx, col="red", lty=2)
  abline(v=17,col='blue')
  if (month == 12 ) { t0 <- proc.time()[3] - t0 }
}

stopImplicitCluster()

