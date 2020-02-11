#' Check input
#'@title
#'@param inout.list
checkInput <- function(input.list){
  ##----step1: process input


  ## check quantification parameters
  algo.control <- input.list[[4]]
  default.algo.control <- list(ncore = detectCores(),
                               method = "two-step",
                               convcontrol = 10^(-3))

if(is.null(algo.control)){
  algo.control <- default.algo.control
}else{
  if(is.null(algo.control[["ncore"]])|(as.numeric(algo.control[["ncore"]])>default.algo.control[["ncore"]])){
    algo.control[["ncore"]] <- default.algo.control[["ncore"]]
  }
  if(is.null(algo.control[["method"]])|(!(algo.control[["method"]] %in% c("one-step","two-step")))){
    algo.control[["method"]] <- default.algo.control[["method"]]
  }
  if(is.null(algo.control[["method"]])){
    algo.control[["convcontrol"]] <- default.algo.control[["convcontrol"]]
  }
}
process.input <- list(dt, algo.control)
return(process.input)
}
