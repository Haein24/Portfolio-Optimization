library(mvtnorm) # generation simulation data
library(emdbook) # lseq lambda grid
library(readxl)
library(quadprog)
library(doParallel)
library(CVXR)

# 함수 파트
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
      daily_return[i-1,j] <<- (data1[i,j] - data1[i-1,j])/data1[i-1,j]
    }
  }
  #print(paste("nrow : ", dim(daily_return)[1], "ncol : ", dim(daily_return)[2]))
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
  return(mean(turnover_prep))
  #print(mean(turnover_prep))
  #turnover_res <<- mean(turnover_prep)
}

# set_diff : 0에서 값이 생기거나 값이 있다가 0이 된 경우의 개수 평균
setdiff <- function(resultset) {
  set_diff <- 0
  for (i in 2:nrow(resultset)) {
    set_diff <- set_diff + sum((resultset[i,]==0 & resultset[i-1,]!=0) | (resultset[i,]!=0 & resultset[i-1,]==0)) # 계속 누적해서 더해서 계산 
  }
  return(set_diff/nrow(resultset))
}

#par(mfrow=c(4,1))

########## 데이터 준비
# DOW30 DAX30 SP100 FTSE100 FTSE250 SP500
data_prep("DAX30", 1)
daily_return_prep(data1)

########## lambda sequence 생성
lamb_seq(30) # 람다 시퀀스 생성 함수 (몇개를 생성할 것인지 지정)
lambda <- cbind(lambda, rep(lambda[1,17],nrow(lambda))) # Lasso를 위해서 31번째에 라쏘를 위한 람다값을 붙여줌


######### 병렬작업 준비
n_core = detectCores() # 사용가능한 코어가 몇개인지 확인
cl = makeCluster(n_core-3) # 사용가능한 코어중 하나를 빼놓고 할당
cl
registerDoParallel(cl)


##################
# 기본 input 정의
window <- 250 # 11개월치 사용
start <- 250
Nwindow <- 21 # monthly update
lenmonth <- 12*10 # 몇 년치 데이터를 다룰 것인지 : 10년치

# 릿지용
w_tilde <- rep(0,p)
mu <- 0.05 

# 
#q <- 0.1 # 0.15 
#rp <- 0.006 # 0.005
#ridge_YN <- TRUE # FALSE
#
q <- 0.10 ; rp <- 0.005;  ridge_YN <- TRUE

eta <- 6 # 고정
tau <- 10^(-5)
e <- rep(1,p)
II <- diag(p)
KK <- floor(q*window) # big K


########## 최종 결과를 저장할 데이터 프레임 생성
results <- data.frame("month"=0, "lambda"=0,"SLOPE_30"=0, "return_MV"=0,"SLOPE_17"=0, "lasso_17"=0,"index"=0, "EW"=0)
results_30 <- results_MV <- results_w_las <- results_w_17 <- matrix(0,lenmonth,ncol(data1)-1) 


t0 <- proc.time()[3] # 시간 측정
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
    j <- 0
    kvec <- c(rep(1/KK, KK),rep(0, ncol(R)-KK))
    
    
    while (TRUE) {
      
      if (ridge_YN == TRUE) {
        # 릿지 
        w <- solve( c(eta)*(c(1+(2*mu/eta))*II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma -c(2*mu)*w_tilde))
      } else {
        # non-penalty
        w <- solve(c(eta)*(II + e%*%t(e) + rb%*%t(rb) + R%*%t(R) )) %*% ( c(eta)*(v+e+c(rp)*rb +R%*%z) - (rb + alpha + c(beta)*e + c(theta)*rb - R%*%gamma))
      }
      
     
      v <- proxSortedL1(w+(1/eta)*alpha, lambda/eta)
      z <- proxZ(R,w,gamma,eta,kvec)
      alpha <- alpha + c(eta)*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%rb-rp)
      gamma <- gamma + c(eta)*(z - t(R)%*%w)
      
      
      # dual gap 계산 및 체크 ## 아닌것 같음!
      Gap <- t(w-v)%*%(w-v) + t(z-t(R)%*%w)%*%(z-t(R)%*%w) + (t(rb)%*%w-rp)^2 + (t(e)%*%w-1)^2 
      
      if (Gap < tau) break
      
      j <- j+1
      if ( j >= 1000 ) break
    }
    return(w) # 결과 내보내기
  }
  
  
  w_idx <- which(colSums(round(ww[,-ncol(ww)],2)!=0)/p >= 0.3 & 0.4 >= colSums(round(ww[,-ncol(ww)],2)!=0)/p) # 30~40% 정도의 자산을 가지고 운영하겠다
  
  if (length(w_idx) >= 2 & month >=2) {
    j <- 1
    comp <- rep(0,length(w_idx))
    for(i in w_idx ) {
      comp[j] <- sum(abs(results_30[month-1,]-round(ww[,i],2))) # turnover 계산
      j <- j+1
    }
    w_idx <- w_idx[which.min(comp)] # turnover가 가장 작은데 에서 선택
  } else if (length(w_idx) >= 2 & month == 1) { 
    w_idx <- w_idx[which(max(colMeans(t(R)%*%round(ww[,w_idx],2)))==colMeans(t(R)%*%round(ww[,w_idx],2)))] # 가장 return이 높은 지점에서 선택  # sample(w_idx,1) # 임의로 하나 선택
  } else { # 30~40퍼센트 내에 데이터가 없는 경우
    # 조건에 맞는게 없으면 active rate이 중간값인 0.35에서 가장 가까운 지점에서 선택
    w_idx <- which.min(abs(colSums(round(ww[,-ncol(ww)],2)!=0)/p - 0.35))
  }
  
  # 평가
  # 결정된 람다를 적용해 다음 21간의 수익률을 계산 및 저장
  ###### 최종적으로 내보낼때는 평가 대상 : 인덱스 수익률 & Equally Weighted 데이터(데이터 프레임에 순서대로 저장)
  N_daily_return <- (data1[(end+Nwindow-1),]-data1[(end),])/data1[(end),]
  
  # 월수익 계산 및 저장
  return_30 <- N_daily_return[-ncol(data1)]%*%ww[,w_idx] # 투자는 1%대로 반올림해서 사용
  return_17 <- N_daily_return[-ncol(data1)]%*%round(ww[,17],2)
  return_lasso <- N_daily_return[-ncol(data1)]%*%round(ww[,ncol(ww)],2)
  idx <- N_daily_return[ncol(data1)]
  ew <- N_daily_return[-ncol(data1)]%*%rep(1/p, p) # equally weighted portfolio mean return
  
  w_res <- round(ww[,w_idx],2)
  
  ## 
  results_0 <- as.data.frame(table(w_res[w_res != 0]))
  
  for (i in 1:nrow(results_0)) {
    if (results_0[i,2] != 1) {
      results_0$idx[i] <- which.min(((rowSums(R-as.vector(rowMeans(R)))^2)/ncol(R))[which(w_res==results_0[i,1])])
      results_0$new_w[i] <- as.numeric(levels(results_0$Var1)[i])* results_0$Freq[i]
    } else {
      results_0$idx[i] <- which(w_res==results_0[i,1])  
      results_0$new_w[i] <- as.numeric(levels(results_0$Var1)[i])
    }
  }
  
  return_MV <- t(as.vector(N_daily_return[results_0$idx]))%*%results_0$new_w
  ##
  
  results <- rbind(results, c(month, w_idx, return_30, return_MV, return_17, return_lasso, idx, ew)) #, return_17, return_lasso, idx, ew))
  results_30[month,] <- round(ww[,w_idx],2)
  results_MV[month, results_0$idx] <- results_0$new_w #round(ww[results_0$idx, w_idx],2)
  results_w_17[month,] <- round(ww[,17],2)
  results_w_las[month,] <- round(ww[,31],2)
  
  
  
  #######  진행상황 및 grouping 확인용 
  if (month%%12 == 0) {
    print(paste(month*100/lenmonth,"%")) # 진행상황 확인
    
    # 시각화
    #real_viz(ww[,-ncol(ww)], month, w_idx)
    #abline(v=w_idx, col="blue")
    #abline(v=17,col='red', lty=2)
    
  }
  
  ###### 다음으로 넘어가기 위한 작업
  start <- start+Nwindow
  w_tilde <- ww[,w_idx]
  
  if (month == lenmonth ) { t0 <- proc.time()[3] - t0 }
}

stopImplicitCluster()


############### OUTPUT 부분에 사용 ###############
dt_result_D <- data.frame()
port_name <- c('SLOPE_30%','SLOPE_MV' ,'SLOPE_17', 'lasso_17', 'index', 'EW')

for ( i in 3:8) {
  port <- port_name[i-2]
  risk <- round(sd(results[-1,i]),4) # variance risk
  mere <- round(mean(results[-1,i])*12,4)
  if ((i == 3) || (i == 8)) {
    T <- lenmonth # 가지고있는 월별 포트폴리오 수익률의 개수
    KK <- floor(q*T) # big K
    rr <- results[-1,i] # SLOPE_30의 평균 수익률
    shortfall_risk <- mean(rr) - (1/KK)*sum(sort(rr)[1:KK]) # 숏폴 리스크
    shortfall_shr <- round(mean(results[-1,i])/shortfall_risk, 4) # 숏폴 기준 샤프지수
  } else {
    shortfall_risk <- NA
    shortfall_shr <- NA
  }
  shr <- round(mean(results[-1,i])/sd(results[-1,i]),4) # 샤프지수
  dt_result_D <- rbind(dt_result_D, c(port, risk, mere, shr, round(shortfall_risk,4), shortfall_shr))
}

dt_result_D <- cbind(dt_result_D,c(round(turnover(results_30),4),round(turnover(results_MV),4), round(turnover(results_w_17),4), round(turnover(results_w_las),4),NA, NA))
dt_result_D <- cbind(dt_result_D, c(round(setdiff(results_30),4),round(setdiff(results_MV),4), round(setdiff(results_w_17),4), round(setdiff(results_w_las),4), NA, NA))
dt_result_D <- cbind(dt_result_D, c(round(mean(rowSums(results_30 !=0)),4), round(mean(rowSums(results_MV !=0)),4), round(mean(rowSums(results_w_17 !=0)),4), round(mean(rowSums(results_w_las !=0)),4), NA, p)) # round(mean(rowSums(results_MV !=0)),4), 

colnames(dt_result_D) <- c("Portfolio", "Risk", "Mean_Return", "Sharpe_ratio", "shortfall", "shortfall_shr", "Turnover", "Set_Diff", "AP")
View(dt_result_D)