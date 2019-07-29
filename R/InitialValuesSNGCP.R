#'
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
find.staring.values <- function(D,lBext)
{
  #Some helping values
  Bext <- lBext*lBext # size of extended observation window
  
  #Define paramters-----------------------------------
  P             = list()
  P$theta       <- 0.0005
  P$lambda      <- 2
  P$alpha       <- 1
  P$beta        <- 0.1
  P$c = D$x
  P$int_kernel = list()
  P$log_kernel = list()
  P$pairdist_c = list()
  for(replicate in 1:D$r_x){
    P$m[[replicate]] = 1:D$n_x[replicate]
    P$gamma[[replicate]]       <- rep(1,dim(P$c[[replicate]])[1])
  }
  #----------------------------------------------------
  
  for(replicate in 1:D$r_x){
    #Define the latent points. Divide the area into small rectangles,
    #and select the latent point as the average of the observed points in the rectangle.------
    step = 0.2
    ind <- seq((1-lBext)/2,lBext +(1-lBext)/2  ,by = step)
    m = list()
    c = matrix(-1,10000,2)
    neste = 1
    for(i in 1:(length(ind)-1))

    {
      for(j in 1:(length(ind)-1))
      {
        area = cbind(c(ind[i],ind[i+1],ind[i+1], ind[i]), c(ind[j], ind[j], ind[j+1],ind[j+1]))
        area = data.frame(area)
        names(area) = c("x","y")
        if(sum(sp::point.in.polygon(D$x[[replicate]][,1],D$x[[replicate]][,2], area$x, area$y))>0)
        {
          hvilke = which(sp::point.in.polygon(D$x[[replicate]][,1],D$x[[replicate]][,2], area$x, area$y)!=0)
          x = mean(D$x[[replicate]][hvilke,1])
          y = mean(D$x[[replicate]][hvilke,2])
          c[neste,] = c(x,y)
          P$m[[replicate]][hvilke] = neste
          P$gamma[[replicate]][neste] = length(hvilke)
          neste = neste + 1
        }else if(min(area)<0){
          c[neste,] = c(min(area$x),min(area$y))
          P$gamma[[replicate]][neste] = 1
          neste = neste + 1
        }else if(max(area)>1){
          c[neste,] = c(max(area$x),max(area$y))
          P$gamma[[replicate]][neste] = 1
          neste = neste + 1
        }
      }
    }
    c = c[which(c[,1]>(-1)),]
    P$c[[replicate]] = c
    P$gamma[[replicate]] = P$gamma[[replicate]][which(c[,1]>(-1))]
    P$n_c[replicate]         <- dim(P$c[[replicate]])[1]
    #----------------------------------------------------
    
    
    #Define variables needed in the calculation of the posterior----
    P$int_kernel[[replicate]] = rep(NA,P$n_c[replicate]) 
    for(j in 1:P$n_c[replicate])
    {
      P$int_kernel[[replicate]][j] <- calcInt(P$c[[replicate]][j,], P$theta)
    }
    
    P$log_kernel[[replicate]] = matrix(NA,D$n_x[replicate],P$n_c[replicate]) 
    for(j in 1:P$n_c[replicate])
    {
      P$log_kernel[[replicate]][,j] <- calcLogKernel(D$x[[replicate]], P$c[[replicate]][j,], P$theta)
    }
    
    P$pairdist_c[[replicate]] = as.matrix(dist(P$c[[replicate]])) 
    #----------------------------------------------------
    
    
  
  }
  
  
  P$priors = set.priors.SNGCP()
  return(P)
}



#' find.staring.values.truth
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
find.staring.values.truth <- function(D,lBext)
{
  r_x = D$r_x
  P <- list()
  #Some helping values
  for(i in 1:r_x){
    n_x <- D$n_x[i]
    Bext <- lBext*lBext # size of extended observation window
    
    #Define paramters-----------------------------------
    P$theta       <- par$theta
    P$lambda      <- par$lambda
    P$alpha       <- par$alpha
    P$beta        <- par$beta
    P$c[[i]]          <- GMT$c[[i]] 
    P$gamma[[i]]        <- GMT$r[[i]] 
    P$m[[i]]            <- GMT$m[[i]] 
    #----------------------------------------------------
    
    #Define variables often needed in the MCMC-----------
    P$n_c[i]         <- dim(P$c[[i]])[1]
    
    P$int_kernel[[i]]  = rep(NA,P$n_c[i]) 
    for(j in 1:P$n_c[i])
    {
      P$int_kernel[[i]] [j] <- calcInt(P$c[[i]][j,], P$theta)
    }
    
    P$log_kernel[[i]]  = matrix(NA,n_x,P$n_c[i]) 
    for(j in 1:P$n_c[i])
    {
      P$log_kernel[[i]][,j] <- calcLogKernel(D$x[[i]] , P$c[[i]][j,], P$theta)
    }
    
    P$pairdist_c[[i]]  = as.matrix(dist(P$c[[i]] )) 
    #----------------------------------------------------
    
    P$priors = set.priors.SNGCP()
    
  }
  return(P)
}






#'
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
set.priors.SNGCP = function()
{
  priors <- NULL
  
  #Priors for lambda (gamma-distributed)----
  priors$a_lambda = 1#100000000
  priors$b_lambda = 0.1#100000000
  #-----------------------------------------
  
  #Priors for theta (gamma-distributed)-----
  priors$a_theta = 1#5000000
  priors$b_theta = 0.1#0.1#1000000000
  #-----------------------------------------
  
  #Priors for alpha (gamma-distributed)-----
  priors$a_alpha = 1#100000000
  priors$b_alpha = 0.1#0.1#100000000
  #-----------------------------------------
  
  #Priors for beta (gamma-distributed)-----
  priors$a_beta = 1#50000000
  priors$b_beta = 0.1#1000000000
  #-----------------------------------------
  
  return(priors)
}

#'
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
define.data.object = function(x)
{
  r_x = length(x)
  #Define the data object-----------------
  D = list()
  D$x = x
  D$n_x = rep(0,r_x)
  for(i in 1:r_x){
    D$n_x[i] = dim(x[[i]])[1]
  }
  #---------------------------------------
  D$r_x = r_x
return(D)
  
}


#'
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
define.accept.object = function()
{
  accept <- NULL
  accept$split <- 0
  accept$merge <- 0
  accept$move <- 0
  accept$m <- 0
  accept$theta <-0
  accept$alpha <- 0
  accept$beta <- 0

  return(accept)
}


