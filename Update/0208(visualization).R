# 그룹화 시각화
setwd('Desktop/portfolio_selection')
data_set <- c("DOW30", "DAX30", "FTSE100", "SP100", "FTSE250", "SP500")
library(dplyr)

# 기본 시각화
eta <- 15 #10
phi <- 40 #3
xlb <- 'nth lambda sequence' ; ylb <- "weight"
for ( a in 1:6) {
  data <- as.matrix(read.csv(paste("result(40,15)/result_",data_set[a],".csv", sep="")))# result(k,15)/result_",data_set[a],".csv", sep="")))
  
  title <- paste("data = ", data_set[a], ", phi = ", phi, ", eta = ", eta)
  
  png(filename=paste("result(40,15)/black/",data_set[a],".png", se=""),width=600,height=395,unit="px") # 395, 525 사용
  plot(data[1,], ylim=c(min(data), max(data)), type='l', col=1, lwd=1, main=title, ylab=ylb, xlab=xlb) # round(eta,3)))
  for(i in 2:nrow(data)) lines(data[i,], type='l')
  dev.off()
}

# 그룹화
# 반올림은 소수점 5번째 자리를 기준으로 사용하고자함
# data_set <- c("DOW30", "DAX30", "FTSE100", "SP100", "FTSE250", "SP500")
ll <- rep(7, 6) # 몇번째 람다를 기준으로 사용할 것인지 지정
ll <- rep(23, 6)
ll <- rep(23,6)

dataset <- data.frame(matrix(ncol=3))
colnames(dataset) <- c("id", "weight", "group")
for ( aa in 3) {
  data <- as.matrix(read.csv(paste("result(3,10)/result_",data_set[aa],".csv", sep="")))
  ww <- round(data[,ll[aa]],5) # 원활한 그룹화를 위해 반올림해 사용 (6~8까지 차이 x)
  # table(ww) # weight기준 몇개의 그룹이 이뤄졌는지
  
  glst <- sort(unique(ww)) # group list : 작->큰 순서로 정렬
  gid <- rep(0, length(ww)) # group id
  for( i in 1:length(glst)) {
    gid[which(ww==glst[i])] <- i # group id
  }
  
  idx <- which(ww!=0)
  title1 <- paste(data_set[aa]," using ",ll[aa],"th lambda / active = ", length(idx), " group = ",length(unique(gid)), sep="")
  # print(title1)
  
  #png(filename=paste('result(k,15)/color/',data_set[aa],".png", se=""),width=600,height=395,unit="px")
  plot(data[idx[1],], ylim=c(min(data),max(data)), type='l', col=gid[idx[1]]+1, lwd=1, main=title1, ylab=ylb, xlab=xlb) # round(eta,3)))
  for(i in 2:length(idx)) lines(data[idx[i],], type='l', col=gid[idx[i]]+1)
  abline(v=ll[aa], type='l', lty='dotdash', col=1)
  #dev.off()
  
  dd <- data.frame(id = rep(data_set[aa],length(ww)),weight = ww, group = gid)
  dataset <- bind_rows(dataset, dd)
  dataset$weight <- as.numeric(dataset$weight)
  dataset$group <- as.numeric(dataset$group)
}
dataset <- dataset[-1,]

View(dataset)

for ( aa in 1:6) {
  pie(table(table(dataset[(dataset$id==data_set[aa])&(dataset$weight!=0),'group'])), main=paste(data_set[aa],"distribution of # of assets in each groups"))
  #barplot(sort(table(dataset[(dataset$id==data_set[aa])&(dataset$weight!=0),'group'])), main=data_set[aa])
}

table(dataset[(dataset$id==data_set[3])&(dataset$weight!=0),'weight'])
      
# 그룹화 후 각 그룹의 weight의 합
for ( i in 3) {
  ff <- data.frame(dataset[dataset$id==data_set[i],] %>%
                        group_by(group) %>%
                        summarise(sum = sum(weight)))
  
  ff$freq <- table(dataset[dataset$id==data_set[i],'group'])

  title <- paste("result(k,15)/ratio/ratio_",data_set[i],'.csv',sep="")
  write.csv(ff[order(-ff$freq),],title,row.names = FALSE)
}



SP100 <- dataset[dataset$id=='SP100', ]
table(SP100$group)
which(SP100$group==12)


###
title <- paste("data/",data_set[4],".xlsx", sep="")
print(title)
data1 <- read_excel(title ,sheet = 1, col_names = TRUE)
data1 <- as.matrix(data1[,seq(from=1,to=ncol(data1),2)])

# 결측치 제거
h <- c()
for ( i in 1:ncol(data1)) { if (sum(is.na(data1[,i])>0)) h <- append(h, i) }
data1 <- data1[,-h]
colnames(data1)[which(SP100$group==12)]



plot(data[1,], ylim=c(-0.5,0.5), type='l', col=1, lwd=1, main=title, ylab=ylb, xlab=xlb) # round(eta,3)))
for(i in 2:nrow(data)) lines(data[i,], type='l')
