#--------------------------------
# Acceptance-Rejection sampling #
#--------------------------------

#----------
# Example 1
#----------

# We use Accept-Reject to estimate the mean and variance of a beta distribution
n = 10000  # number of draws

f = function(x) {dbeta(x,6.3,2.7)}  # target distribution

c = 2.67  # maximum (found using for ex:
optimize(f, interval=c(0, 1), maximum=TRUE)

g = dunif   # proposal distribution

# Accept-Reject algorithm
theta=c()
for (i in 1:n){
  x <- runif(1, 0, 1) 
  y <- runif(1, 0, c) 
  acceptance_rate <- f(x) / (c*g(x)) # acceptance probability
  if (y <= acceptance_rate){
    theta[i] <- x
  }
  else {next}
}

# Plotting the results along with the true distribution
hist(theta, breaks=100, freq=FALSE, main="Random draws from beta(6.3, 2.7)")
lines(x, dbeta(x, 6.3, 2.7))

# and let's see if we get the same results
cat("Estimated mean is equal to", mean(theta, na.rm=TRUE))
quantile(theta, probs=c(0.05,0.95), na.rm=TRUE)
cat("Standard deviation is equal to", sd(theta, na.rm=TRUE))

# check
mean(rbeta(10000, 6.3, 2.7))
sd(rbeta(10000, 6.3, 2.7))
# same results as when using the inbuilt function rbeta()

points = runif(n)
uniforms = runif(n)
accept = uniforms < (f(points)/(c*g(points)))

# plotting
curve(f, lwd=2, main="Acceptance - Rejection method to sample random variates from a Beta(2.7, 6.3)")
points(points, c*uniforms, pch=ifelse(accept,4,4), col=ifelse(accept,"darkred","gray60"), cex=0.5)
legend("bottomleft", c("f","accepted","rejected"), 
       lwd=c(2,NA,NA), col=c("black","darkred","gray60"),
       pch=c(NA,4,4), bg="white")

#----------
# Example 2
#----------

set.seed(2023)
# we want to learn about the mean survival time of butterflies using bayesian analysis
n = 10000  # number of draws
# density of exponential distribution with parameter theta = 1
g <- function(x) return( exp(-x) ) # proposal, easy to sample from

# posterior distribution of the data
# 13 = 6 + 3 + 2 + 2 : number of butterflies alive
# 10 = 4 + 3 + 1 + 0 + 2 : number of butterflies dead
# 1/x is the Jeffrey's prior
# exp(-13*x) * (1-exp(-x))^10 ) is the likelihood of the data
f <- function(x) return( 1/x *exp(-13*x) * (1-exp(-x))^10 ) # target (posterior), hard to sample from

optimize(f, interval=c(0, 100), maximum=TRUE)
c <- 5 * 10^(-7)

# Accept-Reject algorithm
theta=c()
for (i in 1:n){
  x <- runif(1, 0, 1) 
  y <- runif(1, 0, c) 
  acceptance_rate <- f(x) / (c*g(x)) # acceptance probability
  if (y <= acceptance_rate){
    theta[i] <- x
  }
  else {next}
}

# and let's see if we get the same results
cat("Estimated EST is equal to", 1/mean(theta, na.rm=TRUE))
# Estimated EST is equal to 1.946919

points = rexp(n)
uniforms = runif(n, min = 0, max = c*g(points))
accepted = uniforms < f(points)

# plotting
curve(f, lwd=2, main="Acceptance - Rejection method to sample random variates from the target distribution")
points(points, uniforms, pch=ifelse(accepted,4,4), col=ifelse(accepted,"darkred","gray60"), cex=0.5)
legend("bottomleft", c("f","accepted","rejected"), 
       lwd=c(2,NA,NA), col=c("black","darkred","gray60"),
       pch=c(NA,4,4), bg="white")

#----
# end
#----