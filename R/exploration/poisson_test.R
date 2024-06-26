n <- 50
run_time <- 10
p <- c(b = 2, d = 0.5)

source("C:/Users/madel/OneDrive/Documenten/R/grind.R")

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), { 
    P <- state[1:n]
    P_1 <- c(0,P[1:(n-1)])
    dP <- 2*b*P_1 - (b+d)*P
    return(list(dP))  
  }) 
}
s <- c(1, rep(0,(n-1)))
names(s) <- seq(0,n-1)

output <- run(tmax = run_time)
output <- data.frame(div_index = as.numeric(names(s)), 
                     count = output)

output$corrected_count = (output$count/2^output$div_index)
output$fraction = output$corrected_count/sum(output$corrected_count)

ggplot(output, aes(x = div_index)) + 
  geom_bar(aes(y = fraction), stat = "identity") + 
  geom_point(aes(y = dpois(x = seq(0, 49), lambda = run_time*p[1])))

dpois(x = seq(0, 49), lambda = run_time*p[1])
