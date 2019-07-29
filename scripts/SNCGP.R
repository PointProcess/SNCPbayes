#A script which run inference with MCMC for a SNCGP
rm(list=ls())
library(SNCPbayes)

set.seed(1311242)
#Simulate data-------------------------------------
par = list()
par$theta = 0.0005;par$beta = 0.05; par$alpha = 1.2; par$lambda = 1;w.delta = 0;r_x =1
#par$theta = 0.0005;par$beta = 0.02; par$alpha = 0.2; par$lambda = 3;w.delta = 0;r_x =1
#par$theta = run[[3]]$theta;par$beta = run[[3]]$beta; par$alpha = run[[3]]$alpha; par$lambda = run[[3]]$lambda;w.delta = 0;r_x = 1
GMT = simFiniteGmeasure(par,w.delta = w.delta, r_x = r_x,minGamma = 0)
#--------------------------------------------------

#Plot simulation-----------------------------------
x11()
par(mfrow = c(ceiling(sqrt(r_x)),ceiling(sqrt(r_x))))
for(i in 1:r_x){
  plot(GMT$x[[i]],xlim = c(-w.delta, 1+w.delta), ylim = c(-w.delta,1 + w.delta),cex = 0.5,
       xlab = "x",ylab = "y")
  points(GMT$c[[i]][GMT$r[[i]]>5,],cex = 2, col = 'red', pch = 20) 
  points(GMT$c[[i]][GMT$r[[i]]<5,],cex = 1, col = 'green', pch = 20) 
  
  j =9
  points(GMT$x[[i]][GMT$m[[i]]==j,],cex = 1, col = 'blue', pch = 20)
  points(GMT$c[[i]][j,1],GMT$c[[i]][j,2],cex = 2, col = 'blue', pch = 20)
  
  rect(0, 0, 1, 1) 
  rect(-w.delta, -w.delta, 1+w.delta, 1+w.delta, lty = 2)
}
#--------------------------------------------------


#run MCMC------------------------------------------
x11()
x = GMT$x
Rprof()
conf = setConf(R = 10000,starAtTruth = F, lBext = 1,jumpFactor = 0.1,splitFactor = 1,truncGamma = 0)
#conf$map = c("alpha", "beta","theta","lambda")
conf$kappaVonMises=0
#conf$map = c("mark", "center","theta","fixedChilderen")
run = run.MCMC.SNGCP(x,conf)
Rprof(NULL)
#summaryRprof()
#--------------------------------------------------

#Plot some results---------------------------------
x11()
reportRun(run = run,path = getwd(), par = par,truthKnown = TRUE)

plot(run[[1]][,5])
abline(h = dim(GMT$c[[1]])[1])

x11()
plot(run[[4]][,1], main = "loglikelihood")
abline(h =run[[4]][1,1])
#--------------------------------------------------







#Run with real data

data(list = "treesPP", package = "ShotSliceSampler")

xReal = list()
treesMod = list()
i = 1
indeks = 1:8
#indeks = c(1,8)
#indeks = c(1,2,4,5,6)
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

conf = setConf()
conf$lBext = 2
run = run.MCMC.SNGCP(x = xReal,conf)
reportRun(run = run)
#Plot the locations of trees in the first dataframe in the list
par(mfrow = c(3,3))
r_x = length(xReal)
for(i in 1:r_x){
  plot(treesMod[[i]]$x,treesMod[[i]]$y)
}





#Simulation experiment--------------------
rm(list=ls())
library(ShotSliceSampler)
set.seed(9991)#set.seed(298885)
nSim = 10
R = 2000
misesList = list()
normalList = list()

par = list()
par$theta = 0.0005;par$beta = 0.05; par$alpha = 1.2; par$lambda = 1;w.delta = 0.2;r_x = 4

for(simulation in 1:nSim){
  #Simulate data
  GMT = simFiniteGmeasure(par,w.delta = w.delta, r_x = r_x)
  
  #run MCMC------------------------------------------
  x = GMT$x
  
  conf = setConf(R = R,useVonMises = TRUE,starAtTruth = TRUE)
  runMises = run.MCMC.SNGCP(x,conf)
  conf = setConf(R = R,useVonMises = FALSE,starAtTruth = TRUE)
  run = run.MCMC.SNGCP(x,conf)
  #--------------------------------------------------

  misesList[[simulation]]=runMises
  normalList[[simulation]]=run
}


splitProp = rep(NA,nSim)
mergeProp = rep(NA,nSim)

for(i in 1:nSim)
{
  splitProp[i] = misesList[[i]][[2]]$split/normalList[[i]][[2]]$split
  mergeProp[i] = misesList[[i]][[2]]$merge/normalList[[i]][[2]]$merge
}

plot(splitProp)



plot(rvonmises(100,1,kappa = 0.1))



