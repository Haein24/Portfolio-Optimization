# z 업데이트를 위한 함수
# proxZ function
library(quadprog)

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

# solve.QP(Dmat,dvec, t(A), bvec)$iteration
