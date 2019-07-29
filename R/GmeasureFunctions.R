#' t
#' @description Simulates a spatial Poisson process with constant intensity mu on a sqarea area
#' @param mu The constant intensity
#' @param w.min Lower bound of the sqare area
#' @param w.max Upper bound of the sqare area
#' @export
#' @return An array with the locations of the samples
#' @examples
simPoisson <- function(mu,w.min,w.max)
{
    ## simulate a spatial Poisson process
    N <- rpois(n=1,lambda=mu*(w.max-w.min)^2)
    p <- array(data=runif(2*N,min=w.min,max=w.max),dim=c(N,2))
    return(p)
}


#' Sim a finite Gmeasure 
#' @description  simulate a cluster process on [0,1]x[0,1] with cluster centers
#' having a G measure with alpha < 0 and iid Gaussian dispersion
#' densities alpha < 0, beta > 0 and lambda > 0 are parameters of
#' the G measure theta is the spread of the isotropic Gaussian
#' dispersion density w.delta is the edge effect factor for
#' simulating the cluster centers
#' @param 
#' @export
#' @return 
#' @examples
#' GMT = simFiniteGmeasureThomas(alpha=1.5,beta = 0.1, lambda =5,theta = 0.03,givenLaten = TRUE)
#' plot(GMT$x,xlim = c(-0.2,1.2), ylim = c(-0.2,1.2),cex = 0.5)
#' points(GMT$c,cex = 0.7, col = 'red')
simFiniteGmeasure <- function(par, w.delta=0.1, r_x = 1,minGamma = 0)
{
 
  mu <- par$lambda * abs(par$alpha)^(-1) * par$beta^(-par$alpha)
  pp <- list()
  pp$x = list()
  pp$m = list()
  pp$c = list()
  pp$r = list()
  for(realisation in 1:r_x){
    c <- simPoisson(mu=mu, w.min=-w.delta, w.max=1+w.delta)
    
    n.c <- dim(c)[1]
    r <- rgamma(n=n.c, shape=abs(par$alpha), rate=par$beta)
    done = FALSE
    while(!done){
      if(min(r)<minGamma) {
        r[r<minGamma] = rgamma(n=sum(r<minGamma), shape=abs(par$alpha), rate=par$beta)
      }else{
        done = TRUE
      }
    }
    n.offspring <- rpois(n=n.c, lambda=r)
    x <- NULL
    m <- NULL
    for(j in 1:n.c)
    {
      if(n.offspring[j]>0)
      {
        n.tmp <- n.offspring[j] 

        x.tmp <- cbind(rnorm(n=n.tmp, mean=c[j,1], sd=sqrt(par$theta)),
                       rnorm(n=n.tmp, mean=c[j,2], sd=sqrt(par$theta)))
        ind.tmp <- apply((x.tmp>0) * (x.tmp<1) * 1, 1, sum)
        ind.tmp <- ind.tmp==2
        if(sum(ind.tmp)>0)
        {
          x.tmp <- x.tmp[which(ind.tmp),,drop=FALSE]
          n.tmp <- dim(x.tmp)[1]
          x <- rbind(x, x.tmp)
          m <- c(m, rep(j, n.tmp))
        }
      }
    }

    pp$x[[realisation]] <- x
    pp$m[[realisation]] <- m
    pp$c[[realisation]] <- c
    pp$r[[realisation]] <- r
  }

  return(pp)
}


#' Joint posterior
#' @description Calculates the joint posterior
#' @export
#' @return 
#' @examples
logJointPost <- function(D,P,B=1,B.ext,conf)
{
  
  truncateAdditional = 0
  for(replicate in 1:D$r_x){
    #Calculate variables in the posterior---
    eOffspring = rep(0,P$n_c[replicate])
    for(j in 1:P$n_c[replicate])
    {
      eOffspring[j] = sum(P$m[[replicate]]==j)
    }
    
    if(conf$truncGamma>0){
      truncateAdditional = truncateAdditional + sum(log(1-pgamma(conf$truncGamma, shape = eOffspring + abs(P$alpha), 
                                                                 rate = P$beta + P$int_kernel[[replicate]])))
    }
  }
  
  jointPost = dgamma(P$alpha,shape = P$priors$a_alpha, rate = P$priors$a_alpha,log = TRUE)+
    dgamma(P$beta,shape = P$priors$a_beta, rate = P$priors$a_beta,log = TRUE)+
    dgamma(P$lambda,shape = P$priors$a_lambda, rate = P$priors$a_lambda,log = TRUE)+
    dgamma(P$theta,shape = P$priors$a_theta, rate = P$priors$a_theta,log = TRUE)
  
  jointPost = jointPost - truncateAdditional
  for(replicate in 1:D$r_x){
    jointPost = jointPost +
      P$n_c[replicate]*B - sum(P$gamma[[replicate]] * P$int_kernel[[replicate]]) + sum(log(P$gamma[[replicate]][P$m[[replicate]]])) + sum(diag(P$log_kernel[[replicate]][1:D$n_x[replicate],P$m[[replicate]]]))-
      (B.ext *P$alpha^(-1)*P$beta^(-P$alpha)*P$lambda)+P$n_c[replicate]*(-log(P$alpha) -P$alpha *log(P$beta) +log(P$lambda) ) +
      sum(dgamma(P$gamma[[replicate]],shape = P$alpha, rate = P$beta,log = TRUE))
  }
  return(jointPost)
}


#' Split: Sample one random center which is proposed for a split. 
#' A split is proposed by adding and subtracting the sampled vector to the original center, the angle of the
#' proposed vector is sampled from Von Mises distribution with mean equal an observed center of a cluster of the children. 
#' The mark to each proposed split is sampled by randomly portion out the mark of the original center, 
#' where one proportion is beta-distributed with parameters corresponding to the number of children in each cluster. 
#' The observed points previously assigned to the center is then divided randomly to new centers accordinlgy to the isotropic Gaussian dispersion density.
#' If the latent point has in the current less than 3 children, a more simple procedure is used.
#' @description  
#' @param D The data
#' @param P The current state of the MCMC
#' @param u Random uniform vector of numbers used in the MCMC for selctin split/merge/move
#' @param r Random uniform vector of numbers used for acceptance
#' @param accept Information regarding acceptance
#' @param lBext Boundary
#' @param kappaVonMises Paramter used in the Von Mises distribtuion for sampling angles for split
#' @export
#' @return 
#' @examples
#UpdateGprocessCenter <- function(D, P, u,r, accept, lBext,kappaVonMises = NULL, shrinkBoundary = FALSE)
UpdateGprocessCenter <- function(D, P, u,r, accept, conf)
{

  #Some constants needed
  numOfChildredToUtilizeCluster = 3
  ext = (conf$lBext -1)/2
  splitSd = sqrt(P$theta)/conf$splitFactor
  jumpSd = sqrt(P$theta)/conf$jumpFactor
  
  
  for(replicate in 1:D$r_x){
    n_c = P$n_c[replicate]
    
    #Sample the center process with split-merge-move----------------------------------------------------
    if(u[replicate] <= 1/3)
    {
        #Sample center to split
        j = sample(1:n_c, prob = P$gamma[[replicate]],1)
        probSelectSplit = P$gamma[[replicate]][j]/sum(P$gamma[[replicate]])
        
        #Sample direction of split and portion of gamma to divide between the splitted center------------
        children = D$x[[replicate]][which(P$m[[replicate]] == j),, drop = FALSE]
        if(sum(P$m[[replicate]] == j)>=numOfChildredToUtilizeCluster& conf$kappaVonMises!=0){

          #Divide the offspring into two clusters
          clusters = find.clusters(children, n.pca = 1, n.clust = 2)$grp
          if(clusters[1] ==2){
            tmpClusters = clusters
            tmpClusters[which(clusters==1)] = 2
            tmpClusters[which(clusters==2)] = 1
            clusters = tmpClusters
          }
          
          clusterCenter1 = c(mean(children[clusters ==1,1]),mean(children[clusters ==1,2])  )
          clusterCenter2 = c(mean(children[clusters ==2,1]),mean(children[clusters ==2,2])  )
             
          #Select the cluster center with highest y-value to make it reversible
          clusterCenterUse = clusterCenter2
          vA = sum(clusters ==2);vB = sum(clusters ==1) #Paramters in the beta distibution for dividing gamma
          if(clusterCenter1[2]>clusterCenter2[2]){
            clusterCenterUse = clusterCenter1
            vA = sum(clusters ==1);vB = sum(clusters ==2)
          } 
        
          #Select proportion of gamma to divide betweene the splitted center
          v <- rbeta(1,vA,vB) 
          
          #Find the angle used as mean in Von Mises
          vectorTmp1 =c(1,0)
          vectorTmp2 =c(clusterCenterUse[1]-P$c[[replicate]][j,1],clusterCenterUse[2]-P$c[[replicate]][j,2])

          if(sum(vectorTmp2)==0){
            vectorTmp2 = c(0,1)
            stop("The direction vector in Von Mises proposal is zero, this should not happen. Stop the MCMC")
          }
          angle = acos( sum(vectorTmp1*vectorTmp2) / ( sqrt(sum(vectorTmp1 * vectorTmp1)) * sqrt(sum(vectorTmp2 * vectorTmp2)) ) )
          if(vectorTmp1[2] >vectorTmp2[2]) angle = -angle #Angle between (-pi,pi)
  
          #Sample angle
          angle = circular(angle) #Avoid warning in vonmises as it is casted to a circular in vonmises()
          angleSampledOriginal = as.numeric(rvonmises(1,mu = angle,kappa = conf$kappaVonMises))
          angleSampled = angleSampledOriginal
          if(angleSampledOriginal>pi){
            angleSampled = angleSampled-2*pi#Angle between (-pi,pi)
          }
          angleSampledOriginal = circular(angleSampledOriginal) #Avoid warning in vonmises as it is casted to a circular in vonmises()
          
          #Sample distance in split
          splitDist = abs(rnorm(1,0,splitSd))
          #splitDist = splitSd*sqrt(rchisq(1,2))

          delta = c(splitDist*cos(angleSampled),splitDist*sin(angleSampled))
        }else{#Sample in a simple way
          vA = 2
          vB = 2
          v <- rbeta(1,vA,vB)
          angle = circular(0) 
          
          angleSampledOriginal = as.numeric(rvonmises(1,mu = angle,kappa = 0)) #uniform sampling
          angleSampled = angleSampledOriginal
          if(angleSampledOriginal>pi){
            angleSampled = angleSampled-2*pi
          }
          angleSampledOriginal = circular(angleSampledOriginal) #Avoid warning in vonmises as it is casted to a circular in vonmises()
          
          #Sample distance in split
          splitDist = abs(rnorm(1,0,splitSd))
          
          delta = c(splitDist*cos(angleSampled),splitDist*sin(angleSampled))
        }
        #-------------------------------------------------------------------------------------
        
        eta1.new <- P$c[[replicate]][j,] + delta
        eta2.new <- P$c[[replicate]][j,] - delta

        #Controll for that we do not move outside the border-------------
        if(eta1.new[1]>(1+ext)){
          eta1.new[1]=(1+ext)
        } 
        if(eta1.new[1]<(-ext)){
          eta1.new[1]=(-ext)
        } 
        if(eta1.new[2]>(1+ext)){
          eta1.new[2]=(1+ext)
        } 
        if(eta1.new[2]<(-ext)){
          eta1.new[2]=(-ext)
        } 
        if(eta2.new[1]>(1+ext)){
          eta2.new[1]=(1+ext)
        } 
        if(eta2.new[1]<(-ext)) {
          eta2.new[1]=(-ext)
        } 
        if(eta2.new[2]>(1+ext)){
          eta2.new[2]=(1+ext)
        } 
        if(eta2.new[2]<(-ext)){
          eta2.new[2]=(-ext)
        }
        #-----------------------------------------------------------------
        
        #Calculate some values in the posterior---------------------------
        r.old = P$gamma[[replicate]][j] 
        
        r1.new <- v * P$gamma[[replicate]][j] 
        r2.new <- (1-v) * P$gamma[[replicate]][j] 

        if(r1.new< conf$truncGamma | r2.new< conf$truncGamma){
          return(list(P=P,accept = accept))
        }
        
        int_kernel1.new <- calcInt(eta1.new, P$theta)
        log_kernel1.new <- calcLogKernel(D$x[[replicate]], eta1.new, P$theta)
        int_kernel2.new <- calcInt(eta2.new, P$theta)
        log_kernel2.new <- calcLogKernel(D$x[[replicate]], eta2.new, P$theta)
        #-----------------------------------------------------------------
        
        #Sample m---------------------------------------------------------
        ind.old <- which(P$m[[replicate]]==j)
        m.new <- rep(j, length(ind.old))
        k <- 1
        new.labels <- c(j, P$n_c[replicate] + 1)
        logProbOfSampled = rep(0,length(ind.old))
        tmp = runif(length(ind.old))

        for(l in ind.old)
        {
          tmpProb = (r1.new*exp(log_kernel1.new[l]))/( (r1.new*exp(log_kernel1.new[l])) + (r2.new * exp(log_kernel2.new[l])))
          if(tmp[k]< tmpProb){
            logProbOfSampled[k] = log(tmpProb)
          }else{
            logProbOfSampled[k] = log(1 - tmpProb)
            m.new[k] = new.labels[2]
          }
          k <- k + 1
        }
        ind1.new <- ind.old[which(m.new==new.labels[1])]
        ind2.new <- ind.old[which(m.new==new.labels[2])]
        #-----------------------------------------------------------------
        
      
        #calculate Hastings ratio------------------------------------------------------- 
        log.post.ratio <- 1-r1.new*int_kernel1.new - r2.new*int_kernel2.new +
          sum(log(r1.new) + log_kernel1.new[ind1.new]) + sum(log(r2.new) + log_kernel2.new[ind2.new]) +
          r.old*P$int_kernel[[replicate]][j] - sum(log(r.old) + P$log_kernel[[replicate]][ind.old,j])+ 
          log(P$lambda/P$alpha) - 
          log(gamma(P$alpha)) + (P$alpha-1)*(log(r1.new) + log(r2.new)- log(r.old)) 
        

        getDist <- as.matrix(dist(rbind(P$c[[replicate]][-j,, drop = FALSE],eta1.new, eta2.new)))[which(lower.tri(as.matrix(dist(rbind(P$c[[replicate]][-j,, drop = FALSE],eta1.new, eta2.new)))))]
        probSelect = (exp(-getDist[length(getDist)]/(sqrt(P$theta)*2)))/ sum(exp(-getDist/(sqrt(P$theta)*2)))

        #Find intervall of v which produces alowed gamma
        v1 = conf$truncGamma/P$gamma[[replicate]][j] 
        v2 = -conf$truncGamma/P$gamma[[replicate]][j]  +1
        addBeta = log(pbeta(max(v1,v2), vA,vB) - pbeta(min(v1,v2), vA,vB))
        
      
         log.prop.ratio = -log(probSelectSplit) - 
           #        log(dchisq((splitDist/splitSd)^2,2)*2*(splitDist/splitSd^2))-
           dnorm(splitDist,0,splitSd,log = TRUE) - log(2)-
           sum(logProbOfSampled) +
           log(probSelect) 
         
         log.prop.ratio = log.prop.ratio - addBeta
         
        
         
         
        if(sum(P$m[[replicate]]==j)>=numOfChildredToUtilizeCluster & conf$kappaVonMises!=0){
          log.prop.ratio = log.prop.ratio - 
            log(dvonmises(angleSampledOriginal,mu =angle,kappa = conf$kappaVonMises) * dbeta(v,vA,vB) + dvonmises(angleSampledOriginal-pi,mu =angle,kappa = conf$kappaVonMises) * dbeta(1-v,vA,vB) )
        }else{
          log.prop.ratio = log.prop.ratio- #            log(dvonmises(angleSampledOriginal,mu =angle,kappa = 0) )
            log(dvonmises(angleSampledOriginal,mu =angle,kappa = 0)* dbeta(v,vA,vB) + dvonmises(angleSampledOriginal-pi,mu =angle,kappa = 0)* dbeta(1-v,vA,vB)) #Pairs of v and theta that results in the split can be sampled in two different ways.
        }
        log.ratio = log.post.ratio + log.prop.ratio
 
        
        #--------------------------------------------------------------------------------

        
        if(log(r[replicate]) < log.ratio)
        {
            ## accept and update object---------------------------------------------------------
            if(n_c>1)P$c[[replicate]][j,] <- eta1.new
            if(n_c==1) P$c[[replicate]] <- eta1.new
            P$c[[replicate]] <- rbind(P$c[[replicate]], eta2.new)
            P$n_c[replicate] <-  P$n_c[replicate] + 1
            P$gamma[[replicate]][j] <- r1.new
            P$gamma[[replicate]] <- c(P$gamma[[replicate]], r2.new)
            P$m[[replicate]][ind2.new] <- P$n_c[replicate]
            P$log_kernel[[replicate]][,j] <- log_kernel1.new
            P$log_kernel[[replicate]] <- cbind(P$log_kernel[[replicate]], log_kernel2.new)
            P$int_kernel[[replicate]][j] <- int_kernel1.new
            P$int_kernel[[replicate]] <- c(P$int_kernel[[replicate]], int_kernel2.new)
            P$pairdist_c[[replicate]] <- as.matrix(dist(P$c[[replicate]]))
            accept$split = accept$split+1
            #--------------------------------------------------------------------------------
        }        
    }
  
    if(u[replicate] > 1/3 & u[replicate] <= 2/3 & n_c>2)
    {
        ## Sample centers to merge------------------------------------------------------------
        N <- n_c
        getDistProb <- (exp(-P$pairdist_c[[replicate]][which(lower.tri(P$pairdist_c[[replicate]]))]/(sqrt(P$theta)*2))) / sum(exp(-P$pairdist_c[[replicate]][which(lower.tri(P$pairdist_c[[replicate]]))]/(sqrt(P$theta)*2)))

        tmp <- array(rep(1:N), dim=c(N, N))
        get.ind1 <- tmp[which(lower.tri(tmp))]
        tmp <- array(rep(1:N, each=N), dim=c(N,N))
        get.ind2 <- tmp[which(lower.tri(tmp))]
        M <- length(getDistProb)

        ## select a pair and create merged proposal 
        d <- sample(1:M, 1, prob=getDistProb)
        j1 <- get.ind1[d]
        j2 <- get.ind2[d]
        
        if(FALSE){
          j2 = j
          j1 = n_c
          d = which(get.ind1==j1 &get.ind2 == j2 )
        }
        
        if(max(P$c[[replicate]][j1,]==P$c[[replicate]][j2,])){
          #Exactly the same point, must be on the corner of the boarder. Do nothing
          break;
        }
        
        eta.new <- (P$c[[replicate]][j1,] + P$c[[replicate]][j2,])/2
        
        r.new <- P$gamma[[replicate]][j1] + P$gamma[[replicate]][j2]  

        probSelectSplit = r.new/sum(P$gamma[[replicate]])
        
        int_kernel.new <- calcInt(eta.new, P$theta)
        log_kernel.new <- calcLogKernel(D$x[[replicate]], eta.new, P$theta)
        v <- P$gamma[[replicate]][j1]/r.new
        delta <- eta.new - P$c[[replicate]][j2,]
        #--------------------------------------------------------------------------------
        
        
        #Calculate what needed in the splitting routine    
        children = D$x[[replicate]][which(P$m[[replicate]] == j1|P$m[[replicate]] == j2),, drop = FALSE]
        mergedLocation = eta.new 
        

        if(sum(P$m[[replicate]]==j1 | P$m[[replicate]]==j2)>=numOfChildredToUtilizeCluster & conf$kappaVonMises!=0){

          clusters = find.clusters(children, n.pca = 1, n.clust = 2)$grp
          if(clusters[1] ==2){
            tmpClusters = clusters
            tmpClusters[which(clusters==1)] = 2
            tmpClusters[which(clusters==2)] = 1
            clusters = tmpClusters
          }
          
          clusterCenter1 = c(mean(children[clusters ==1,1]),mean(children[clusters ==1,2])  )
          clusterCenter2 = c(mean(children[clusters ==2,1]),mean(children[clusters ==2,2])  )
          
          #Select the cluster center with highest y-value to make it reversible
          clusterCenterUse = clusterCenter2

          vA = sum(clusters ==2)
          vB = sum(clusters ==1)
          if(clusterCenter1[2]>clusterCenter2[2]){
            clusterCenterUse = clusterCenter1
            vA = sum(clusters ==1)
            vB = sum(clusters ==2)
          } 
      
          #Find the angles used as mean in Von Mises----------------
          tmp = as.matrix(dist(rbind(P$c[[replicate]][j1,],P$c[[replicate]][j2,],clusterCenterUse)))[3,]

          selected = j1
          if(tmp[2]<tmp[1]){
            selected = j2
          } 
          vectorTmp1 =c(1,0)
          vectorTmp2 =c(P$c[[replicate]][selected,1]-mergedLocation[1],P$c[[replicate]][selected,2]-mergedLocation[2])

          angleSampledOriginalFictive = acos( sum(vectorTmp1*vectorTmp2) / ( sqrt(sum(vectorTmp1 * vectorTmp1)) * sqrt(sum(vectorTmp2 * vectorTmp2)) ) )
          angleSampledOriginalFictive = circular(angleSampledOriginalFictive)
          vectorTmp1 =c(1,0)
          vectorTmp2 =c(clusterCenterUse[1]-mergedLocation[1],clusterCenterUse[2]-mergedLocation[2])
          
          if(sum(vectorTmp2)==0){
            vectorTmp2 = c(0,1)
          }
          
          angle = acos( sum(vectorTmp1*vectorTmp2) / ( sqrt(sum(vectorTmp1 * vectorTmp1)) * sqrt(sum(vectorTmp2 * vectorTmp2)) ) )
          if(vectorTmp1[2] >vectorTmp2[2]) angle = -angle
          angle = circular(angle)
          v = P$gamma[[replicate]][selected]/(P$gamma[[replicate]][j1] + P$gamma[[replicate]][j2])
          #----------------------------------------------------------
          
          #Calculate distance in split
          splitDist = sqrt((mergedLocation[1]-P$c[[replicate]][j1,1])^2 + (mergedLocation[2]-P$c[[replicate]][j1,2])^2)
        }else{
          splitDist = sqrt((mergedLocation[1]-P$c[[replicate]][j1,1])^2 + (mergedLocation[2]-P$c[[replicate]][j1,2])^2)
          vA =2
          vB =2
        }

        #Calculate the probability of sample the children------------------------
        ind <- which(P$m[[replicate]]==j1 | P$m[[replicate]]==j2)
        ind1 <- which(P$m[[replicate]]==j1)
        ind2 <- which(P$m[[replicate]]==j2)
        k = 1
        logProbOfSampled = rep(0,length(ind))
        for(l in ind)
        {
          rmpProb = (P$gamma[[replicate]][j1] * exp(P$log_kernel[[replicate]][l,j1]))/ (P$gamma[[replicate]][j1] * exp(P$log_kernel[[replicate]][l,j1]) + P$gamma[[replicate]][j2] * exp(P$log_kernel[[replicate]][l,j2]))
          if(is.element(l,ind1)){
            logProbOfSampled[k] = log(rmpProb)
          }else{
            logProbOfSampled[k] = log(1-rmpProb)
          }
          k <- k + 1
        }
        #-------------------------------------------------------------------------
        
        
        #Find intervall of v which produces alowed gamma
        v1 = conf$truncGamma/ r.new
        v2 = -conf$truncGamma/r.new  +1
        addBeta = log(pbeta(max(v1,v2), vA,vB) - pbeta(min(v1,v2), vA,vB))
        
        
        #Calculate Hastings ratio -------------------------------------------------Olav: Has compared this to the difference of logJointPost(), The difference between logJointPost() gave same result.
        log.post.ratio <- -(1-P$gamma[[replicate]][j1]*P$int_kernel[[replicate]][j1] - P$gamma[[replicate]][j2]*P$int_kernel[[replicate]][j2]+
          sum(log(P$gamma[[replicate]][j1]) + P$log_kernel[[replicate]][ind1,j1]) + sum(log(P$gamma[[replicate]][j2]) + P$log_kernel[[replicate]][ind2,j2]) +
          r.new*int_kernel.new - sum(log(r.new) + log_kernel.new[ind])+ 
          log(P$lambda/P$alpha)  -  
          log(gamma(P$alpha)) + (P$alpha-1)*(log(P$gamma[[replicate]][j1]) + log(P$gamma[[replicate]][j2])- log(r.new)))

      
        log.prop.ratio = log(probSelectSplit) +# dbeta(v,vA,vB,log=TRUE) + 
  #       log(dchisq((splitDist/splitSd)^2,2)*2*(splitDist/splitSd^2))+
          dnorm(splitDist,0,splitSd,log = TRUE) + log(2)+
         sum(logProbOfSampled) - 
          log(getDistProb[d]) 

        log.prop.ratio = log.prop.ratio + addBeta
        
        
        if(sum(P$m[[replicate]]==j1 | P$m[[replicate]]==j2)>=numOfChildredToUtilizeCluster& conf$kappaVonMises!=0){
          log.prop.ratio = log.prop.ratio + log(dvonmises(angleSampledOriginalFictive,mu =angle,kappa = conf$kappaVonMises)* dbeta(v,vA,vB) +dvonmises(angleSampledOriginalFictive-pi,mu =angle,kappa = conf$kappaVonMises)* dbeta(1-v,vA,vB) )
        }else{
          log.prop.ratio = log.prop.ratio + log(dvonmises(circular(0),mu =circular(0),kappa = 0)* dbeta(v,vA,vB) + dvonmises(circular(0),mu =circular(0),kappa = 0)* dbeta(1-v,vA,vB))
        }
        log.ratio = log.post.ratio + log.prop.ratio

        #-------------------------------------------------------------------------
        if(log(r[replicate]) < log.ratio)
        {
            ## accept and update object-----------------------------------------------
            tmp1 = min(j1,j2)
            tmp2 = max(j1,j2)
            P$c[[replicate]][tmp1,] <- eta.new
            P$c[[replicate]] <- P$c[[replicate]][-tmp2,, drop = FALSE]
            P$m[[replicate]][which(P$m[[replicate]]>tmp2)] = P$m[[replicate]][which(P$m[[replicate]]>tmp2)]-1
            P$n_c[replicate] <-  P$n_c[replicate] - 1
            P$gamma[[replicate]][tmp1] <- r.new
            P$gamma[[replicate]] <- P$gamma[[replicate]][-tmp2]
            P$m[[replicate]][ind] <- tmp1
            P$log_kernel[[replicate]][,tmp1] <- log_kernel.new
            P$log_kernel[[replicate]] <- P$log_kernel[[replicate]][,-tmp2, drop = FALSE]
            P$int_kernel[[replicate]][tmp1] <- int_kernel.new
            P$int_kernel[[replicate]] <- P$int_kernel[[replicate]][-tmp2, drop = FALSE]
            P$pairdist_c[[replicate]] <- as.matrix(dist(P$c[[replicate]]))
            accept$merge = accept$merge+1
            #-------------------------------------------------------------------------
        }        
    }
    if(u[replicate] > 2/3)
    {
      
        ## Sample center to move-------------------------------------------------------
        j <- sample(1:P$n_c[replicate], 1)
        
        eta.new <- rnorm(2, mean=P$c[[replicate]][j,], sd=jumpSd)
        #-----------------------------------------------------------------------------
        
        adjustProp = min(abs(eta.new[1]- (-ext)), abs(eta.new[2]-  (-ext)),  abs((1+ext)-eta.new[1]),abs((1+ext)-eta.new[2]) )+0.01

        #Controll for that we do not move outside the border-------------
        if(eta.new[1]>(1+ext)){
          eta.new[1] = (1+ext)
        } 
        if(eta.new[1]<(-ext)) {
          eta.new[1]=(-ext)
        } 
        if(eta.new[2]>(1+ext)){
          eta.new[2]=(1+ext)
        } 
        if(eta.new[2]<(-ext)) {
          eta.new[2]=(-ext)
        } 
        #-------------------------------------------------------------------------
        
        #Calculates some values in the posterior----------------------------------
        int_kernel.new <- calcInt(eta.new, P$theta)
        ind <- which(P$m[[replicate]]==j)
        log_kernel.new <- calcLogKernel(D$x[[replicate]], eta.new, P$theta)
        #-------------------------------------------------------------------------
        
        
        ## calculate Hastings ratio-----------------------------------------------
        log.post.ratio <- -P$gamma[[replicate]][j] * (int_kernel.new - P$int_kernel[[replicate]][j]) +
            sum(log_kernel.new[ind] - P$log_kernel[[replicate]][ind,j])
        
        prop_ratio = sum(dnorm(P$c[[replicate]][j,], mean=eta.new, sd=jumpSd,log =TRUE)) -
          sum(dnorm(eta.new, mean=P$c[[replicate]][j,], sd=jumpSd,log =TRUE))
        
        log.ratio = log.post.ratio + prop_ratio
        #-------------------------------------------------------------------------
        
        if(log(r[replicate]) < log.ratio)
        {
            ## accept and update object---------------------------------------------------
            if(n_c>1)P$c[[replicate]][j,] <- eta.new
            if(n_c==1)P$c[[replicate]] <- t(as.matrix(eta.new))
            P$log_kernel[[replicate]][,j] <- log_kernel.new
            P$int_kernel[[replicate]][j] <- int_kernel.new
            if(n_c>1){
              new.dist <- sqrt( (P$c[[replicate]][,1]-eta.new[1])^2
                                + (P$c[[replicate]][,2]-eta.new[2])^2 )
              P$pairdist_c[[replicate]][j,] <- new.dist
              P$pairdist_c[[replicate]][,j] <- new.dist
            }
            accept$move = accept$move+1
            #------------------------------------------------------------------------------
        }
    }
  }
  #end split-merge-move--------------------------------------------------------------------------------------
  return(list(P=P,accept = accept))
}


#'
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
UpdateGprocessMark <- function(D,P,conf)
{
  
  if(!conf$NS){ #If NS, a Newman-Scot process is used. Work in progress...
  for(replicate in 1:D$r_x){
    #Calculate variables in the posterior---
      eOffspring = rep(0,P$n_c[replicate])
      for(j in 1:P$n_c[replicate])
      {
        eOffspring[j] = sum(P$m[[replicate]]==j)
      }
      
      kk = P$gamma[[replicate]]
      
      tmp = rgamma(P$n_c[replicate], shape = eOffspring + abs(P$alpha), 
                                     rate = P$beta + P$int_kernel[[replicate]])
      
      P$gamma[[replicate]] = tmp
      
      P$gamma[[replicate]][which(tmp<conf$truncGamma)] <- kk[which(tmp<conf$truncGamma)]
    }
  }else{
    
    helperKernel = 0
    xx = 0
    for(i in 1:D$r_x){
      helperKernel = helperKernel + sum(P$int_kernel[[replicate]])
      xx = xx +  sum(D$n_x[replicate])
    }
    gamma = rgamma(1, shape = xx + P$priors$a_gamma, 
                   rate = P$priors$b_gamma + helperKernel)
    for(i in 1:D$r_x){
      P$gamma[[replicate]] <- rep(gamma,length(P$gamma[[replicate]]))
    }
  }
  return(list(P = P))
}


#' 
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
updateGMeasureHyperparameters <- function(P,accept,a,r_x,conf) {

  #Define some variables ------------------
   Bext = conf$lBext^2 #The area of the area where the center process live.
  #----------------------------------------
   truncateAdditional = 0
   for(replicate in 1:r_x){
     #Calculate variables in the posterior---
     eOffspring = rep(0,P$n_c[replicate])
     for(j in 1:P$n_c[replicate])
     {
       eOffspring[j] = sum(P$m[[replicate]]==j)
     }
     truncateAdditional = truncateAdditional + sum(log(1-pgamma(conf$truncGamma, shape = eOffspring + abs(P$alpha), 
                                                            rate = P$beta + P$int_kernel[[replicate]])))
   }
   
   
   
  #Sample alpha----------------------------
   if(!"alpha" %in% conf$map){
     if(!conf$NS){ #If NS, a Newman-Scot process is used. Work in progress...
       f = P$alpha^(-1)*P$beta^(-P$alpha)
       fDot = -f*(P$alpha^(-1) + log(P$beta))
       fDotDot = -fDot*(P$alpha^(-1) + log(P$beta)) + P$alpha^(-2)*f
         
       g1 =  (P$priors$a_alpha-1)*P$alpha^(-1) - P$priors$b_alpha #prior
       g1 = g1 - length(P$n_c)*P$lambda *Bext*fDot -sum(P$n_c*(P$alpha^(-1) + digamma(P$alpha)))
       for(i in 1:r_x){
         g1 = g1 + sum(log(P$gamma[[i]]))
       }
      
       g2 = -(P$priors$a_alpha-1)*P$alpha^(-2) #prior
       g2 = g2 - length(P$n_c)*P$lambda*Bext*fDotDot + sum(P$n_c*(P$alpha^(-2) - trigamma(P$alpha)))
       
       b = g1-g2*P$alpha
       c = -g2
       
       sdProp = sqrt(1/c)
       muProp = b/c
       
       alphaProp <- rnorm(1, muProp, sdProp)
       
       truncateAdditionalProp=0
       for(replicate in 1:r_x){
         #Calculate variables in the posterior---
         eOffspring = rep(0,P$n_c[replicate])
         for(j in 1:P$n_c[replicate])
         {
           eOffspring[j] = sum(P$m[[replicate]]==j)
         }
         truncateAdditionalProp = truncateAdditionalProp + sum(log(1- pgamma(conf$truncGamma, shape = eOffspring + alphaProp, 
                                                                        rate = P$beta + P$int_kernel[[replicate]])))
       }
       
       
       if(alphaProp >0){
         fBack = alphaProp^(-1)*P$beta^(-alphaProp)
         fDotBack = -fBack *(alphaProp^(-1) + log(P$beta))
         fDotDotBack = -fDotBack *(alphaProp^(-1) + log(P$beta)) + alphaProp^(-2)*fBack 
         
         g1Back =  (P$priors$a_alpha-1)*alphaProp^(-1) - P$priors$b_alpha #prior
         g1Back = g1Back - length(P$n_c)*P$lambda *Bext*fDotBack  -sum(P$n_c*(alphaProp^(-1) + digamma(alphaProp)))
         for(i in 1:r_x){
           g1Back = g1Back + sum(log(P$gamma[[i]]))
         }
         
         g2Back = -(P$priors$a_alpha-1)*alphaProp^(-2) #prior
         g2Back = g2Back - length(P$n_c)*P$lambda*Bext*fDotDotBack  + sum(P$n_c*(alphaProp^(-2) - trigamma(alphaProp)))
         
         b = g1Back-g2Back*alphaProp
         c = -g2Back
         
         sdBack = sqrt(1/c)
         muBack = b/c
         
          logPostAlphaProp = (P$priors$a_alpha-1)*log(alphaProp)- P$priors$b_alpha*alphaProp - truncateAdditionalProp;
          logPostAlphaOld = (P$priors$a_alpha-1)*log(P$alpha)- P$priors$b_alpha*P$alpha - truncateAdditional;
          for(replicate in 1:r_x){
            logPostAlphaProp = logPostAlphaProp -P$n_c[replicate]*(log(alphaProp) + log(gamma(alphaProp))) -Bext* P$lambda*alphaProp^(-1)*P$beta^(-alphaProp)+ 
              sum((alphaProp-1)*log(P$gamma[[replicate]])) 
            logPostAlphaOld = logPostAlphaOld -P$n_c[replicate]*(log(P$alpha) +log(gamma(P$alpha))) -Bext* P$lambda*P$alpha^(-1)*P$beta^(-P$alpha)+ 
              sum((P$alpha-1)*log(P$gamma[[replicate]]))
          }
          
           rAlpha = logPostAlphaProp - logPostAlphaOld+
             dnorm(P$alpha, muBack, sdBack,log = TRUE)-
             dnorm(alphaProp, muProp, sdProp,log = TRUE)
          if( log(a[1]) < rAlpha )
          {
            P$alpha <- alphaProp
            accept$alpha <- accept$alpha + 1
          }
       }
     }else{ 
       gamma = P$gamma[[1]][1]

       sdProp = 0.1
       alphaProp <- exp(rnorm(1, log(P$alpha), sdProp))
       
       alpha.prop.ratio =  dlnorm(P$alpha, meanlog = log(alphaProp), sdProp,log = TRUE)-
         dlnorm(alphaProp, meanlog = log(P$alpha), sdProp,log = TRUE)
       
       logPostAlphaProp = alphaProp*log(P$beta) - log(gamma(alphaProp)) + (alphaProp-1)*log(gamma) + (P$priors$a_alpha-1)*log(alphaProp)- P$priors$b_alpha*alphaProp
       logPostAlphaOld = P$alpha*log(P$beta) - log(gamma(P$alpha)) + (P$alpha-1)*log(gamma) + (P$priors$a_alpha-1)*log(P$alpha)- P$priors$b_alpha*P$alpha
       
       rAlpha = logPostAlphaProp - logPostAlphaOld+alpha.prop.ratio
         
       if( log(a[1]) < rAlpha )
       {
         P$alpha <- alphaProp
         accept$alpha <- accept$alpha + 1
       }
     }
   }
  ##----------------------------------------
     
   #Sample beta-----------------------------
   if(!"beta" %in% conf$map){
     if(!conf$NS){#If NS, a Newman-Scot process is used. Work in progress...
         
       g1 =  (P$priors$a_beta-1)/P$beta - P$priors$b_beta + r_x*P$beta^(-P$alpha-1)*Bext*P$lambda 
       for(replica in 1:r_x){
         g1 = g1 - sum(P$gamma[[replica]])
       }
       g2 = -(P$priors$a_beta-1)/P$beta^2 - r_x*(P$alpha+1)*P$beta^(-P$alpha -2)*Bext*P$lambda 
       
       b = g1-g2*P$beta
       c = -g2
       
       sdProp = sqrt(1/c)
       muProp = b/c
       
       betaProp <- rnorm(1, muProp, sdProp)
       
       truncateAdditionalProp=0
       for(replicate in 1:r_x){
         #Calculate variables in the posterior---
         eOffspring = rep(0,P$n_c[replicate])
         for(j in 1:P$n_c[replicate])
         {
           eOffspring[j] = sum(P$m[[replicate]]==j)
         }
         truncateAdditionalProp = truncateAdditionalProp + sum(log(1-pgamma(conf$truncGamma, shape = eOffspring + P$alpha, 
                                                                   rate = betaProp + P$int_kernel[[replicate]])))
       }
       
       
       if(betaProp >0){
         g1Back =  (P$priors$a_beta-1)/betaProp - P$priors$b_beta + r_x*betaProp^(-P$alpha-1)*Bext*P$lambda 
         for(replica in 1:r_x){
           g1Back = g1Back - sum(P$gamma[[replica]])
         }
         g2Back = -(P$priors$a_beta-1)/betaProp^2 - r_x*(P$alpha+1)*betaProp^(-P$alpha -2)*Bext*P$lambda 
         bBack = g1Back-g2Back*betaProp
         cBack = -g2Back
         sdBack = sqrt(1/cBack)
         muBack = bBack/cBack
         
         logPostBetaProp = (P$priors$a_beta-1)*log(betaProp)- P$priors$b_beta*betaProp -truncateAdditionalProp
         logPostBetaOld = (P$priors$a_beta-1)*log(P$beta)- P$priors$b_beta*P$beta - truncateAdditional
         
         for(replicate in 1:r_x){
           logPostBetaProp = logPostBetaProp -Bext*P$lambda*P$alpha^(-1)*betaProp^(-P$alpha) - betaProp*sum(P$gamma[[replicate]])  
           logPostBetaOld  = logPostBetaOld -Bext*P$lambda*P$alpha^(-1)*P$beta^(-P$alpha) - P$beta*sum(P$gamma[[replicate]])
         }
         
         rBeta = logPostBetaProp- logPostBetaOld +
           dnorm(P$beta, muBack, sdBack,log = TRUE)-
           dnorm(betaProp, muProp, sdProp,log = TRUE)
         if( log(a[2]) < rBeta )
         {
           P$beta <- betaProp
           accept$beta <- accept$beta + 1
         }
       }
     }else{
       
       gamma = P$gamma[[1]][1]
       sdProp = 0.1
       betaProp <- exp(rnorm(1, log(P$beta), sdProp))
       
       beta.prop.ratio =  dlnorm(P$beta, meanlog = log(betaProp), sdProp,log = TRUE)-
         dlnorm(betaProp, meanlog = log(P$beta), sdProp,log = TRUE)
       
       logPostBetaProp = P$alpha*log(betaProp) - betaProp*gamma + (P$priors$a_beta-1)*log(betaProp)- P$priors$b_beta*betaProp
       logPostBetaOld = P$alpha*log(P$beta) - P$beta*gamma + (P$priors$a_beta-1)*log(P$beta)- P$priors$b_beta*P$beta
       
       rAlpha = logPostAlphaProp - logPostAlphaOld+alpha.prop.ratio
       
       if( log(a[1]) < rAlpha )
       {
         P$alpha <- alphaProp
         accept$alpha <- accept$alpha + 1
       }
     }
   }
     
   #----------------------------------------
 

   #Sample lambda---------------------------
   if(!"lambda" %in% conf$map){
     if(!conf$NS){#If NS, a Newman-Scot process is used. Work in progress...
       P$lambda <- rgamma(1, shape = P$priors$a_lambda + sum(P$n_c),
                          rate = P$priors$b_lambda + r_x*(Bext * P$alpha^(-1) * P$beta^(-P$alpha)))
     }else{
       P$lambda <- rgamma(1, shape = P$priors$a_lambda + sum(P$n_c),
                          rate = P$priors$b_lambda + r_x*(Bext))
     }
  }
   #----------------------------------------
   
  return(list(P=P,accept=accept))
}

#'Update the parent label for the observed offsprings.
#'
#' @description Update the parent label for the observed offsprings
#' @param D The data
#' @param P The parameters
#' @export
#' @return 
#' @examples
updateM = function(D,P)
{
  for(replicate in 1:D$r_x){
    for(i in 1:D$n_x[replicate]) 
    {
      P$m[[replicate]][i] <- sample(1:P$n_c[replicate], 1,prob = P$gamma[[replicate]] *exp(P$log_kernel[[replicate]][i,]))
    }
  }
  return(list(P=P))
}

#'t
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
MetropolisHastingsTheta <- function(D,P,v,accept,conf )
{

  Bext = conf$lBext^2 #The area the center process live.
  
  #Sample proposition -------------------------------
  sdProp = 0.1
  thetaProp <- exp(rnorm(1, log(P$theta), sdProp))
  #-----------------------------------------------------------------
  theta.post.ratio =  (P$priors$a_theta-1)*log(thetaProp)- P$priors$b_theta*thetaProp  -
    (P$priors$a_theta-1)*log(P$theta)+ P$priors$b_theta*P$theta
  
  
  intKernelProp = list()
  logKernelProp = list()
  
  for(replicate in 1:D$r_x){
  #Calculates parts of the loglikelihood ratio----------------------
    n_c = P$n_c[replicate]
    n_x = D$n_x[replicate]
    
    logKernelProp[[replicate]] = P$log_kernel[[replicate]]
    intKernelProp[[replicate]] = rep(0,n_c)
    
    for(j in 1:n_c)
    {
      m = which(P$m[[replicate]]==j)
      if(n_c>1){
        if(length(m)==1) logKernelProp[[replicate]][m,j] = calcLogKernel(t(as.matrix(D$x[[replicate]][m,])),P$c[[replicate]][j,], thetaProp)
        if(length(m)>1) logKernelProp[[replicate]][m,j] = calcLogKernel(as.matrix(D$x[[replicate]][m,]),P$c[[replicate]][j,], thetaProp)
        intKernelProp[[replicate]][j]  = calcInt(P$c[[replicate]][j,], thetaProp)
      }else{
        logKernelProp[[replicate]][m,j] = calcLogKernel(as.matrix(D$x[[replicate]][m,]),P$c[[replicate]], thetaProp)
        intKernelProp[[replicate]][j]  = calcInt(P$c[[replicate]], thetaProp)
      }
    }
    
    #Calcalculate the log-kernel of interest for the ratio-----
      logKernelOfInterest = rep(0,n_x)
      logKernelPropOfInterest = rep(0,n_x)
      for(i in 1:n_x)
      {
        logKernelOfInterest[i] = P$log_kernel[[replicate]][i,P$m[[replicate]][i]]
        logKernelPropOfInterest[i] = logKernelProp[[replicate]][i,P$m[[replicate]][i]]
      }
    #----------------------------------------------------------
    
    theta.post.ratio = theta.post.ratio-sum(P$gamma[[replicate]]*intKernelProp[[replicate]]) + sum(logKernelPropOfInterest)+
      sum(P$gamma[[replicate]]*P$int_kernel[[replicate]]) - sum(logKernelOfInterest)
  }
  
  
  theta.prop.ratio =  dlnorm(P$theta, meanlog = log(thetaProp), sdProp,log = TRUE)-
    dlnorm(thetaProp, meanlog = log(P$theta), sdProp,log = TRUE)
  
  theta.ratio = theta.post.ratio + theta.prop.ratio
  #--------------------------------------------------
  

  #Accept or reject----------------------------------
  if(log(v) < theta.ratio)
  {
    ## accept and update object
    P$theta =  thetaProp
    
    for(replicate in 1:D$r_x){
      P$int_kernel[[replicate]] <- intKernelProp[[replicate]]
      
      #Change the whole log_kernel-matrix------------
      P$log_kernel[[replicate]] = matrix(NA,D$n_x[replicate],P$n_c[replicate]) 
      for(j in 1:P$n_c[replicate])
      {
        P$log_kernel[[replicate]][,j] <- calcLogKernel(D$x[[replicate]], P$c[[replicate]][j,], P$theta)
      }
    }

    #-----------------------------------------------
    
    accept$theta = accept$theta+1
  }  
  #--------------------------------------------------

  return(list(P=P,accept=accept))
}


#' Vector of log(k(x_i|nu_j,theta)) for a given nu_j and theta.
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
calcLogKernel <- function(x, eta, theta)
{
    return(dnorm(x[,1], mean=eta[1], sd=sqrt(theta), log=TRUE) +
           dnorm(x[,2], mean=eta[2], sd=sqrt(theta), log=TRUE))
}

#'We define this as h(nu,theta) in our article
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
calcInt <- function(eta,theta)
{
    return((pnorm((1-eta[1])/sqrt(theta))-pnorm(-eta[1]/sqrt(theta)))*
           (pnorm((1-eta[2])/sqrt(theta))-pnorm(-eta[2]/sqrt(theta))))
}

#' 
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
UpdateSNCGP <- function(D,P,accept,u,r,v,a, conf) {
  if(!"center" %in% conf$map){
    #Sample the latent centers------------------------
    tmp <- UpdateGprocessCenter(D=D,P=P,u=u,r = r,accept=accept,conf)
    
  accept <- tmp$accept
  P <- tmp$P
  }
  #-------------------------------------------------
  
  #Sample the mark to each latent center------------
  if(!"mark" %in% conf$map){
  tmp <- UpdateGprocessMark(D = D,P=P,conf = conf)
  P <- tmp$P
  }
  #-------------------------------------------------

  #Sample the parameters in G measure---------------
  #tmp <- updateGMeasureHyperparameters(P=P,accept=accept,lBext=lBext,a=a,r_x = D$r_x,map = map) 
  tmp <- updateGMeasureHyperparameters(P=P,accept=accept,a=a,r_x = D$r_x,conf = conf) 
  P <- tmp$P
  accept <- tmp$accept
  #-------------------------------------------------
  
  #Sample the affiliation of each point to centers--
  if(!"fixedChilderen" %in% conf$map){
    tmp = updateM(D=D,P=P)
    P <- tmp$P
  }
  #-------------------------------------------------


  #Sample the Gaussian disperion parameter----------
  if(!"theta" %in% conf$map){
    tmp <- MetropolisHastingsTheta(D=D,P=P,v=v,accept=accept,conf)
    P <- tmp$P
    accept <- tmp$accept
  }
  #-------------------------------------------------
  
  return(list(P=P,accept=accept))
}


#'
#' @description 
#' @param R number of iterations
#' @export
#' @return 
#' @examples
run.MCMC.SNGCP = function(x,conf)
{
  #Define data object--------
  D = define.data.object(x = x)
  #--------------------------
  
  #Find staring values and priors-----
  if(conf$starAtTruth){
    P  = find.staring.values.truth(D,lBext = conf$lBext)
  }else{
    P = find.staring.values(D,lBext= conf$lBext)
  }
  accept = define.accept.object()
  #-----------------------------------
  
  for(replicate in 1:D$r_x){
    P$gamma[[replicate]][which(P$gamma[[replicate]]<conf$truncGamma)] = conf$truncGamma
  }
  
  # random values needed for proposals---------------------------
  u <- array(runif(conf$R*D$r_x),dim=c(conf$R,D$r_x))#select split, merge, or move
  v <- cbind(runif(conf$R)) # dispersion
  r <- array(runif(conf$R*D$r_x),dim=c(conf$R,D$r_x)) # accept split-merge-move
  a <- cbind(runif(conf$R),runif(conf$R)) # for hyperparameters of the G-measure update
  #--------------------------------------------------------------
  
  #Start MCMC---------------------------------------------------
  par.history <- array(NA,dim=c(conf$R,5)) # to keep parameters
  post = logJointPost(P=P,D=D,B.ext = conf$lBext^2,conf = conf)
  for(i in 1:conf$R)
  { # simulations
    tmp <- UpdateSNCGP(D,P,accept,u[i,],r[i,],v[i],a[i,], conf)
    
    P <- tmp$P
    
    post = rbind(post,logJointPost(P=P,D=D,B.ext = conf$lBext^2,conf = conf))

    accept <- tmp$accept
    par.history[i,] <- c(P$lambda,P$alpha,P$beta,P$theta,P$n_c[1])
    
    
    if(i%%500==0)
    { 
      print(c(i,date()))
      Sys.sleep(0.05)
      print(paste("Alpha is: ", round(P$alpha,6)))
      print(paste("Beta is: ", round(P$beta,6)))
      print(paste("Lambda is: ", round(P$lambda,6)))
      print(paste("Theta is: ", round(P$theta,6)))
      print(paste("n_c is: ", round(P$n_c[1],6)))

      print(paste("Accept split: ", round(accept$split/(i*D$r_x/3),6), " Merge: ", round(accept$merge/(i*D$r_x/3),6), " move:", round(accept$move/(i*D$r_x/3),6), sep = ""))
      print(paste("Accept theta: ", round(accept$theta/(i),6), sep = ""))
      print(paste("Accept alpha: ", round(accept$alpha/(i),6), sep = ""))
      print(paste("Accept beta: ", round(accept$beta/(i),6), sep = ""))

      

      
      ## plotting the current pattern
      extra = (conf$lBext-1)/2
      par(mfrow = c(ceiling(sqrt(D$r_x)),ceiling(sqrt(D$r_x))))
      for(replicate in 1:D$r_x){
       plot(D$x[[replicate]],xlim = c(-extra,1+extra), ylim = c(-extra,1+extra),cex = 0.5,
             xlab = "x",ylab = "y")
        if(P$n_c[[replicate]]>1){ 
          points(P$c[[replicate]],cex = 1, col = 'green', pch = 20)
          if(sum(P$gamma[[replicate]]>5)>1){
            points(P$c[[replicate]][P$gamma[[replicate]]>5,],cex = 2, col = 'red', pch = 20)
          }else if(sum(P$gamma[[replicate]]>5)==1){
            points(P$c[[replicate]][P$gamma[[replicate]]>5,1],P$c[[replicate]][P$gamma[[replicate]]>5,2],cex = 2, col = 'red', pch = 20)
      
          }
        }else points(P$c[[replicate]][1],P$c[[replicate]][2],cex = 0.7, col = 'red', pch = 20)
        rect(0, 0, 1, 1)
        rect(-extra, -extra, 1+extra, 1+extra, lty = 2)
      
        jj = sample(1:P$n_c[replicate],1)
        points(P$c[[replicate]][jj,1],P$c[[replicate]][jj,2], col = 'blue',cex = 2)
        points(D$x[[replicate]][P$m[[replicate]]==jj,1],D$x[[replicate]][P$m[[replicate]]==jj,2],col = "blue",cex = 1.7)
      }
      print(paste("Log joint posterior: ",logJointPost(P=P,D=D,B.ext = conf$lBext^2,conf = conf),sep =""))
      print("--------------------------------------------")
    } 
  } 
  #End MCMC-------------------------------------------------------------
  
  return(list(par.history,accept,P,post))
}


#'
#' @description 
#' @export
#' @return 
#' @examples
nInside = function(c)
{
  toReturn = 0;
  for(i in 1:dim(c)[1]){
    if(max(c[i,])<1 & min(c[i,])>0)toReturn = toReturn+1     
  }
  return(toReturn)
}
