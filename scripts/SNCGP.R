#A script which run inference with MCMC for a SNCGP
rm(list=ls())
library(SNCPbayes)

set.seed(123)
#Simulate data-------------------------------------
par = list()
par$theta = 0.0005;par$beta = 0.05; par$alpha = 1.2; par$lambda = 1;w.delta = 0.1;r_x =1
GMT = simFiniteGmeasure(par,w.delta = w.delta, r_x = r_x,minGamma = 5)
#--------------------------------------------------

#Plot simulation-----------------------------------
x11()
par(mfrow = c(ceiling(sqrt(r_x)),ceiling(sqrt(r_x))))
for(i in 1:r_x){
  plot(GMT$x[[i]],xlim = c(-w.delta, 1+w.delta), ylim = c(-w.delta,1 + w.delta),cex = 0.5,
       xlab = "x",ylab = "y")
  points(GMT$c[[i]][GMT$r[[i]]>5,],cex = 2, col = 'red', pch = 20) 
  points(GMT$c[[i]][GMT$r[[i]]<5,],cex = 1, col = 'green', pch = 20) 
  rect(0, 0, 1, 1) 
  rect(-w.delta, -w.delta, 1+w.delta, 1+w.delta, lty = 2)
}
#--------------------------------------------------


#run MCMC------------------------------------------
x11()
x = GMT$x
conf = setConf(R = 2000,kappaVonMises=0,starAtTruth = FALSE, lBext = 1.2,jumpFactor = 0.1,splitFactor = 1,truncGamma = 5)
run = run.MCMC.SNGCP(x,conf)
#--------------------------------------------------

#Plot some results---------------------------------
x11()
reportRun(run = run,path = getwd(), par = par,truthKnown = TRUE)
x11()
plot(run[[4]][,1], main = "loglikelihood")
abline(h =run[[4]][1,1])
#--------------------------------------------------





#Run with real data
data(list = "treesPP", package = "SNCPbayes")
xReal = list()
treesMod = list()
i = 1
indeks = 1:8
for(j in indeks){
  treesMod[[i]] = trees[[j]][which(trees[[j]]$x<=6 & trees[[j]]$MY ==(-1)),]
  treesMod[[i]] = treesMod[[i]][which(treesMod[[i]]$y<=6),]
  treesMod[[i]] = treesMod[[i]][which(treesMod[[i]]$x>=0),]
  treesMod[[i]] = treesMod[[i]][which(treesMod[[i]]$y>=0),]
  treesMod[[i]]$x = treesMod[[i]]$x/6
  treesMod[[i]]$y = treesMod[[i]]$y/6
  xReal[[i]] = as.matrix(data.frame(treesMod[[i]]$x,treesMod[[i]]$y))
  i = i+1
}

conf = setConf(R = 5000,lBext = 1.2,truncGamma = 5)
x11()
run = run.MCMC.SNGCP(x = xReal,conf)
reportRun(run = run)
