data_name <- "DOW30"
sheet <- 2
title <- paste("data/",data_name,".xlsx", sep="")
#print(title)
data1 <<- read_excel(title ,sheet = sheet, col_names = TRUE)

View(data1)
names(data1[,seq(1,62,2)])

dim(daily_return)
years7 <- data1[1709:3019,30]
start <- 1
window <- 21
ro <- 0
n <- 0
while (TRUE) {
  n <- n+1
  ro <- append(ro, ((years7[(start+window-1)] - years7[start])/years7[start]))
  start <- start+window
  if (start >= length(years7)) break 
}
