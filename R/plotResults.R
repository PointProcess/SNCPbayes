#'
#' @description 
#' @export
#' @return 
#' @examples
reportRun = function(run, path =NULL, par = NULL,truthKnown = FALSE)
{
  #save(run,file = paste(path,"/run.Rdata",sep = ""))
  
  par(mfrow = c(2,2))
  #some plots----------------------------------------
  plot(run[[1]][,1],main = "lambda")
  if(truthKnown){
    abline(h = par$lambda,col = 'red')
  }
  plot(run[[1]][,2],main ="alpha")
  if(truthKnown){
    abline(h = par$alpha,col = 'red')
  }
  plot(run[[1]][,3],main ="beta")
  if(truthKnown){
    abline(h = par$beta,col = 'red')
  }

  plot(run[[1]][,4],main ="theta")
  if(truthKnown){
    abline(h = par$theta,col = 'red')
  }
  #--------------------------------------------------
}

