### Testing new code for cluster processes with G measure for latent
### cluster centers
rm(list=ls())

source("StraussFunctions.R")
source("GmeasureFunctions.R")

true.alpha <- -1
true.beta <- 0.25
true.lambda <- 15
true.theta <- 0.03

true.mu <- true.lambda * abs(true.alpha) * true.beta^{true.alpha}
print(true.mu)
## simulate point process 
pp <- simFiniteGmeasureThomas(true.alpha, true.beta, true.lambda,
                             true.theta)
## plot parents and offsprings
x11()
plot(pp$x[,1], pp$x[,2], type="p", pch=20, cex=0.7, xlab="", ylab="",
     axes=FALSE, xlim=c(-0.1,1.1), ylim=c(-0.1,1.1), asp=1)
points(pp$c[,1], pp$c[,2], type="p", pch=1, cex=0.5, col="gray70")
lines(c(0,0), c(0,1), lwd=2)
lines(c(0,1), c(0,0), lwd=2)
lines(c(1,1), c(0,1), lwd=2)
lines(c(0,1), c(1,1), lwd=2)
## plot distribution of number of offsprings
x11()
M <- max(round(pp$r)) + 1
hist(pp$r, col="gray60", freq=FALSE, main="", breaks=c(0:M),
     xlab="Expected number of offsprings") 

