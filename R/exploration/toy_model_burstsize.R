library(magrittr)
library(ggplot2)

families = 500

#-------------------------------------------------------------------------------
#For fixed number of burst divisions (4)

q1_fixed <- sample(1:16, families, replace = T)
p1_fixed <- 16 - q1_fixed

q2_fixed <- c()

for (i in 1:length(q1_fixed)) {
  q2_fixed[i] <- sum(sample(1:16, q1_fixed[i], replace = T))
}
p2_fixed <- (16*q1_fixed) - q2_fixed

q3_fixed <- c()

for (i in 1: length(q2_fixed)) {
  q3_fixed[i] <- sum(sample(1:16, q2_fixed[i], replace = T))
}

p3_fixed <- (16*q2_fixed) - q3_fixed
#-------------------------------------------------------------------------------

number_Q_cells_fixed <- data.frame(q1 = q1_fixed,
                                   p1 = p1_fixed,
                                   q2 = q2_fixed,
                                   p2 = p2_fixed,
                                   q3 = q3_fixed,
                                   p3 = p3_fixed,
                                   mean_q2 = q2_fixed/q1_fixed)
number_Q_cells_fixed %<>% replace(is.na(number_Q_cells_fixed), 0)

cor_prim_sec <- cor.test(number_Q_cells_fixed$p1, number_Q_cells_fixed$p2, method = "spearman")
cor_sec_ter <- cor.test(number_Q_cells_fixed$p2, number_Q_cells_fixed$p3, method = "spearman")

ggplot(number_Q_cells_fixed) + geom_jitter(aes(x = p1, y = p2)) + 
  labs(title = "Primary vs secondary", 
       subtitle = paste("Fixed nr of burst divisions (=4)\n", 
                        "r:", round(cor_prim_sec$estimate, digits = 2)))
ggplot(number_Q_cells_fixed) + geom_jitter(aes(x = p2, y = p3)) + 
  labs(title = "Secondary vs tertiary", 
       subtitle = paste("Fixed nr of burst divisions (=4)\n", 
                        "r:", round(cor_sec_ter$estimate, digits = 2)))

# Here we have a "problem" because the fixed number of burst divisons indeed 
# leads to an inverse relationship between the sizes of families in the primary
# and secondary response. I suspect that adding variation in the number of burst
# divisions will help

#-------------------------------------------------------------------------------
# Model wwith variable number of burst divisions

burst_divs <- sample(c(2,3,4), families, replace = TRUE)
q1 <- c()
p1 <- c()

for (i in 1:families) {
  q1[i] <- sample(1:2^burst_divs[i], 1, replace = T)
  p1[i] <- 2^burst_divs[i] - q1[i]
}

q2 <- rep(0, length(q1))
p2 <- rep(0, length(p1))

for (i in 1:length(q1)) {
  cells <- q1[i]
  while (cells != 0){
    burst_divs_sec <- sample(c(2,3,4), 1, replace = TRUE)
    max_cells <- 2^burst_divs_sec
    nr_Q_formed_per_Q <- sample(1:max_cells, 1, replace = T)
    nr_P_formed_per_Q <- max_cells - nr_Q_formed_per_Q
    q2[i] = q2[i] + nr_Q_formed_per_Q
    p2[i] = p2[i] + nr_P_formed_per_Q
    cells = cells - 1
  }
}

q3 <- rep(0, length(q1))
p3 <- rep(0, length(p1))

for (i in 1:length(q1)) {
  cells <- q2[i]
  while (cells != 0){
    burst_divs_ter <- sample(c(2,3,4), 1, replace = TRUE)
    max_cells <- 2^burst_divs_ter
    nr_Q_formed_per_Q <- sample(1:max_cells, 1, replace = T)
    nr_P_formed_per_Q <- max_cells - nr_Q_formed_per_Q
    q3[i] = q3[i] + nr_Q_formed_per_Q
    p3[i] = p3[i] + nr_P_formed_per_Q
    cells = cells - 1
  }
}

#-------------------------------------------------------------------------------

number_Q_cells <- data.frame(q1 = q1,
                             p1 = p1,
                             q2 = q2,
                             p2 = p2,
                             q3 = q3,
                             p3 = p3,
                             mean_q2 = q2/q1)
number_Q_cells %<>% replace(is.na(number_Q_cells), 0)

cor_prim_sec <- cor.test(number_Q_cells$p1, number_Q_cells$p2, method = "spearman")
cor_sec_ter <- cor.test(number_Q_cells$p2, number_Q_cells$p3, method = "spearman")

ggplot(number_Q_cells) + geom_jitter(aes(x = p1, y = p2)) + 
  labs(title = "Primary vs secondary", 
       subtitle = paste("Variable nr of burst divisions (uniform)\n", 
                        "r:", round(cor_prim_sec$estimate, digits = 2)))
ggplot(number_Q_cells) + geom_jitter(aes(x = p2, y = p3)) + 
  labs(title = "Secondary vs tertiary", 
       subtitle = paste("Variable nr of burst divisions (uniform)\n", 
                        "r:", round(cor_sec_ter$estimate, digits = 2)))

# This does indeed help a bit, but you can clearly see that the first plot of 
# primary vs secondary family size actually consists of multiple negative
# correlations for each of the burst sizes

#-------------------------------------------------------------------------------
# Model wwith variable number of burst divisions with particular chances

burst_divs <- sample(c(2,3,4), families, prob = c(0.2, 0.4, 0.4), replace = TRUE)
q1 <- c()
p1 <- c()

for (i in 1:families) {
  q1[i] <- sample(0:2^burst_divs[i], 1, replace = T)
  p1[i] <- 2^burst_divs[i] - q1[i]
}

q2 <- rep(0, length(q1))
p2 <- rep(0, length(p1))

for (i in 1:length(q1)) {
  cells <- q1[i]
  while (cells != 0){
    burst_divs_sec <- sample(c(2,3,4), 1, prob = c(0.2, 0.4, 0.4), replace = TRUE)
    max_cells <- 2^burst_divs_sec
    nr_Q_formed_per_Q <- sample(0:max_cells, 1, replace = T)
    nr_P_formed_per_Q <- max_cells - nr_Q_formed_per_Q
    q2[i] = q2[i] + nr_Q_formed_per_Q
    p2[i] = p2[i] + nr_P_formed_per_Q
    cells = cells - 1
  }
}

q3 <- rep(0, length(q1))
p3 <- rep(0, length(p1))

for (i in 1:length(q1)) {
  cells <- q2[i]
  while (cells != 0){
    burst_divs_ter <- sample(c(2,3,4), prob = c(0.2, 0.4, 0.4), 1, replace = TRUE)
    max_cells <- 2^burst_divs_ter
    nr_Q_formed_per_Q <- sample(0:max_cells, 1, replace = T)
    nr_P_formed_per_Q <- max_cells - nr_Q_formed_per_Q
    q3[i] = q3[i] + nr_Q_formed_per_Q
    p3[i] = p3[i] + nr_P_formed_per_Q
    cells = cells - 1
  }
}

#-------------------------------------------------------------------------------

number_Q_cells <- data.frame(q1 = q1,
                             p1 = p1,
                             q2 = q2,
                             p2 = p2,
                             q3 = q3,
                             p3 = p3,
                             mean_q2 = q2/q1)
number_Q_cells %<>% replace(is.na(number_Q_cells), 0)

cor_prim_sec <- cor.test(number_Q_cells$p1, number_Q_cells$p2, method = "spearman")
cor_sec_ter <- cor.test(number_Q_cells$p2, number_Q_cells$p3, method = "spearman")

ggplot(number_Q_cells) + geom_jitter(aes(x = p1, y = p2)) + 
  labs(title = "Primary vs secondary", 
       subtitle = paste("Variable nr of burst divisions (with different probabilities)\n", 
                        "r:", round(cor_prim_sec$estimate, digits = 2)))
ggplot(number_Q_cells) + geom_jitter(aes(x = p2, y = p3)) + 
  labs(title = "Secondary vs tertiary", 
       subtitle = paste("Variable nr of burst divisions (with different probabilities)\n", 
                        "r:", round(cor_sec_ter$estimate, digits = 2)))

