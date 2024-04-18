
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
logD56$var_rounded <- round(logD56$var)
logD5 <- logD56[which(logD56$day==5),]
logD6 <- logD56[which(logD56$day==6),]
logD7 <- logD56[which(logD56$day==7),]
logD8 <- logD56[which(logD56$day==8),]

empty_table <- data.frame(seq(1,16,1), 2^seq(1, 16, 1))
colnames(empty_table) <- c("Number_2log", "Number")

d5 <- data.frame(table(logD5$var_rounded))
colnames(d5) <- c("Number_2log", "Count")
day5 <- merge(empty_table, d5, by="Number_2log")
day5$Frequency <- day5$Count/sum(day5$Count)

d6 <- data.frame(table(logD6$var_rounded))
colnames(d6) <- c("Number_2log", "Count")
day6 <- merge(empty_table, d6, by="Number_2log")
day6$Frequency <- day6$Count/sum(day6$Count)

d7 <- data.frame(table(logD7$var_rounded))
colnames(d7) <- c("Number_2log", "Count")
day7 <- merge(empty_table, d7, by="Number_2log")
day7$Frequency <- day7$Count/sum(day7$Count)

d8 <- data.frame(table(logD8$var_rounded))
colnames(d8) <- c("Number_2log", "Count")
day8 <- merge(empty_table, d8, by="Number_2log")
day8$Frequency <- day8$Count/sum(day8$Count)

print("Data distribution day 5 till 8 loaded")

