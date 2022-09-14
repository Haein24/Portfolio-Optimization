## mean-variance portfolio ##


########## 데이터 준비
data_prep("FTSE250", 1)
daily_return_prep(data1)

########## lambda sequence 생성
lamb_seq(30) # 람다 시퀀스 생성 함수 (몇개를 생성할 것인지 지정)
# lambda <- cbind(lambda, rep(lambda[1,17],nrow(lambda))) # Lasso를 위해서 31번째에 라쏘를 위한 람다값을 붙여줌


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


########## 최종 결과를 저장할 데이터 프레임 생성
results <- data.frame("month"=0, "lambda"=0,"SLOPE_30"=0, "SLOPE_17"=0, "lasso_17"=0,"index"=0, "EW"=0)
results_30 <- matrix(0,lenmonth,ncol(data1)-1)  # <- results_w_las <- results_w_17 


eta <- 0.01 ; tau <- 10^(-5) 
pi <- 10 # 상대위험회피정도 : pi > 0
rp <- 0.005 # 수정대상1


t0 <- proc.time()[3] # 시간 측정
for (month in 1:lenmonth) {
  end <- start+window-1
  
  R <- t(daily_return[start:end,-ncol(daily_return)]) # dimension : p*T / 인덱스 수익률은 빼고 진행
  
  mu <-  apply(R,mean,MARGIN = 1) # 함수꼴 다시 확인할 필요 있음
  sigma <- cov(t(R)) # p*p
  
  ww <- foreach(lambda=iter(lambda, by='column'), .combine=cbind) %dopar% { #, .export=c("proxSortedL1", "proxZ")) %dopar% { #, .combine=cbind) {
    dyn.load('~/Desktop/portfolio_selection/SLOPE_code/cproxSortedL1.so')
    alpha <- rep(0, p)
    beta <- theta <- 0
    w <- v <- rep(0,p)
    e <- rep(1,p)
    II <- diag(p)
    j <- 0
    
    
    while (TRUE) {
      # update
      #w <- solve(pi*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e)+eta*(v+e))
      w <- solve(pi*sigma + eta*(II+e%*%t(e)))%*%(mu-alpha-as.vector(beta%*%e) + as.vector(theta%*%mu) + eta*(v+e))
      
      v <- proxSortedL1(w+(1/eta)*alpha, lambda/eta)
      alpha <- alpha + eta*(w-v)
      beta <- beta + eta*(t(e)%*%w - 1)
      theta <- theta + eta*(t(w)%*%mu - rp) #추가
      
      # dual gap 계산 및 체크
      Gap <- t(w-v)%*%(w-v) + (t(e)%*%w-1)^2
      if (Gap < tau) break
      
      j <- j+1
      if (j>=1000) break
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
  return_30 <- N_daily_return[-ncol(data1)]%*%round(ww[,w_idx],2) # 투자는 1%대로 반올림해서 사용
  #return_17 <- N_daily_return[-ncol(data1)]%*%round(ww[,17],2)
  #return_lasso <- N_daily_return[-ncol(data1)]%*%round(ww[,ncol(ww)],2)
  idx <- N_daily_return[ncol(data1)]
  ew <- N_daily_return[-ncol(data1)]%*%rep(1/p, p) # equally weighted portfolio mean return
  
  results <- rbind(results, c(month, w_idx, return_30, idx, ew))#, return_17, return_lasso, idx, ew))
  results_30[month,] <- round(ww[,w_idx],2)
  #results_w_17[month,] <- round(ww[,17],2)
  #results_w_las[month,] <- round(ww[,31],2)
  
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


# 최종 output dataframe 만들기
q <- 0.2
KK <- floor(q*window) # big K


dt_result_D <- data.frame()
port_name <- c('SLOPE_30%', 'index', 'EW') # 'SLOPE_17', 'lasso_17', 

for ( i in 3:5) {
  port <- port_name[i-2]
  risk <- round(sd(results[-1,i]),4) # variance risk
  mere <- round(mean(results[-1,i])*12,4)
  if ((i == 3) || (i == 4)|| (i == 5)) {
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

dt_result_D <- cbind(dt_result_D,c(round(turnover(results_30),4), NA, NA))
dt_result_D <- cbind(dt_result_D, c(round(setdiff(results_30),4), NA, NA))
dt_result_D <- cbind(dt_result_D, c(round(mean(rowSums(results_30 !=0)),4), NA, p)) # round(mean(rowSums(results_MV !=0)),4), 

colnames(dt_result_D) <- c("Portfolio", "Risk", "Mean_Return", "Sharpe_ratio", "shortfall", "shortfall_shr", "Turnover", "Set_Diff", "AP")
View(dt_result_D)

###################################
# 최종 output dataframe 만들기
q <- 0.1
KK <- floor(q*window) # big K


dt_result_D01 <- data.frame()
port_name <- c('SLOPE_30%', 'index', 'EW') # 'SLOPE_17', 'lasso_17', 

for ( i in 3:5) {
  port <- port_name[i-2]
  risk <- round(sd(results[-1,i]),4) # variance risk
  mere <- round(mean(results[-1,i])*12,4)
  if ((i == 3) || (i == 4)|| (i == 5)) {
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
  dt_result_D01 <- rbind(dt_result_D01, c(port, risk, mere, shr, round(shortfall_risk,4), shortfall_shr))
}

dt_result_D01 <- cbind(dt_result_D01,c(round(turnover(results_30),4), NA, NA))
dt_result_D01 <- cbind(dt_result_D01, c(round(setdiff(results_30),4), NA, NA))
dt_result_D01 <- cbind(dt_result_D01, c(round(mean(rowSums(results_30 !=0)),4), NA, p)) # round(mean(rowSums(results_MV !=0)),4), 

colnames(dt_result_D01) <- c("Portfolio", "Risk", "Mean_Return", "Sharpe_ratio", "shortfall", "shortfall_shr", "Turnover", "Set_Diff", "AP")
View(dt_result_D01)


