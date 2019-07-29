#'
#' @description 
#' @param 
#' @export
#' @return 
#' @examples
setConf = function(r_x,lBext=1.4,R = 1e4,  kappaVonMises = 20, starAtTruth = FALSE, shrinkBoundary = FALSE,map = NULL,splitFactor = 2, jumpFactor = 2, NS = FALSE,truncGamma = 1, increaseGamma= 1)
{
  conf = list()

  conf$R = R
  conf$lBext = lBext
  conf$kappaVonMises = kappaVonMises
  conf$kappaVonMises = 0
  conf$starAtTruth= starAtTruth
  conf$map = map
  conf$splitFactor = splitFactor
  conf$jumpFactor = jumpFactor
  conf$NS = NS #If true, a Newman-Scot process is used. Work in progress...
  conf$truncGamma = truncGamma
  conf$increaseGamma = increaseGamma
  return(conf)
  
}
