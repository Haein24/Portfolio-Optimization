# grid..

www <- array(dim = c(5,3,12,100))
phi_list <- c(3,12,40)# seq(3,12)
eta_list <-c(0.05, 1, 3, 10, 15) #5 #c(0.05, 0.1, 0.5, 1, 3, 5, 10, 15)

# wwww <-array(dim=c(10,12,30))

for ( x in 1:length(eta_list) ) {
  for ( y in 1:length(phi_list) ) {
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
      phi <- phi_list[y] #k # 상대위험회피정도 : phi > 0 (여기서는 임의의 값으로 두었음,,)
      eta <- eta_list[x] # 15 # round(1/ncol(base),2) # eta > 0 (여기서는 임의의 값으로 두었음,,)
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
      #wwww[y,,z] <- w
      www[x,y, ,z] <- w
      #vv[,z] <- v
      #jj[z] <- j
      #gg[z] <- G
      # print(z)
      #if (z == ncol(lambda)) t0 <- proc.time()[3] - t0
    }
    print(y)
  }
  print(x)
}


#### grid 외부에 저장하기
library(openxlsx)
example <- createWorkbook("example")
for (j in 1:length(phi_list)) {
  for(c in 1:length(eta_list)){
    title <- paste(phi_list[j],"/",eta_list[c], sep="")
    print(title)
    #addWorksheet(example, title)
    #writeDataTable(example, title, data.frame(www[c,j,,]))
  }
}
saveWorkbook(example, file="grid.xlsx")

# 내보냈던 파일 불러오기
library(readxl)
wwww <- array(dim = c(5,3,12,100))
i <- 0
for (x in 1:3) {
  for (y in 1:5) {
    i <- i+1
    wwww[y,x,,] <- as.matrix(read_excel("grid_recovery.xlsx",sheet=paste('복구됨_Sheet',i,sep="")))
  }
}


### 시각화
ylb <- "Weights" ; xlb <- "Grid Lambda"
par(mfrow=c(5,1))
par(mfrow=c(3,5))
www <- wwww
for (j in 1:length(phi_list)) {
  for(c in 1:length(eta_list)){
    title <- paste("phi = ", phi_list[j], ", eta = ", eta_list[c])
    plot(www[c,j,1,], ylim=c(min(www[1,j,,]), max(www[c,j,,])), type='l', col=1, lwd=1, main=title, xlab = xlb, ylab=ylb, cex.main=1, cex.axis=0.7, cex.lab=0.7) # round(eta,3)))
    for(i in 2:4) lines(www[c,j,i,], type='l', col=1, lwd=1)
    for(i in 5:8) lines(www[c,j,i,], type='l', col=2,lwd=1)
    for(i in 9:12) lines(www[c,j,i,], type='l', col=3, lwd=1)
  }
}
