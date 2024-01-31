
#in R
CG103 <- read.delim("./data/raw/run 91, lane 7_1377, CG103_Nienke and GFP correct reformat.txt")
CG97 <- read.delim("./data/raw/run 73, lane 2_CG97_Nienke correct GFP correction reformat.txt")
CG97 <- CG97[,-2] # select the right column
D <- rbind(CG103,CG97) #merge both datasets
D <- D[which(D$var>=2),] # select clones bigger than 2 cells

#load package
library("plyr")
library("ggplot2")
#make dataset per day
logD56 <- ddply(D, c("mouse", "day"), summarize, var= log2(var))
logD56$var_rounded <- floor(logD56$var)
logD5 <- logD56[which(logD56$day==5),]
logD6 <- logD56[which(logD56$day==6),]
logD7 <- logD56[which(logD56$day==7),]
logD8 <- logD56[which(logD56$day==8),]

empty_table <- data.frame(seq(1,16,1), 2^seq(1, 16, 1), matrix(nrow = 16, ncol = 2))
names_columns <- c("Number_2log", "Number", "Count", "Frequency")
colnames(empty_table) <- names_columns


d5 <- data.frame(table(logD5$var_rounded))
while (nrow(d5) < 16) {
  d5 <- add_row(d5, Var1 = as.factor(nrow(d5)+1), Freq = 0)}
day5 <- empty_table
day5$Count <- d5$Freq
day5$Frequency <- day5$Count/sum(day5$Count)


d6 <- data.frame(table(logD6$var_rounded))
while (nrow(d6) < 16) {
  d6 <- add_row(d6, Var1 = as.factor(nrow(d6)+1), Freq = 0)}
day6 <- empty_table
day6$Count <- d6$Freq
day6$Frequency <- day6$Count/sum(day6$Count)


d7 <- data.frame(table(logD7$var_rounded))
while (nrow(d7) < 16) {
  d7 <- add_row(d7, Var1 = as.factor(nrow(d7)+1), Freq = 0)}
day7 <- empty_table
day7$Count <- d7$Freq
day7$Frequency <- day7$Count/sum(day7$Count)


d8 <- data.frame(table(logD8$var_rounded))
while (nrow(d8) < 16) {
  d8 <- add_row(d8, Var1 = as.factor(nrow(d8)+1), Freq = 0)}
day8 <- empty_table
day8$Count <- d8$Freq
day8$Frequency <- day8$Count/sum(day8$Count)

print("Data distribution day 5 till 8 loaded")

