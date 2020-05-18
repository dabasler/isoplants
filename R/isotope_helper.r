
#' Isotopic Delta conversion small d to D
#'
#' @param smpl d value fo sample
#' @param ref d value fo reference
#' @return D notation of sample over reference
#' @export
#'
delta2D <- function(smpl, ref) {
  return((smpl - ref) / (1 + ref/1000))
}

#' Isotopic Delta conversion D to d
#'
#' @param smpl D value fo sample
#' @param ref d value fo reference
#' @return d notation of sample over reference
#' @export
#'
D2delta <- function(smpl, ref) {
  return(  ((1+(ref/1000))*(smpl/1000)+(ref/1000))*1000  )
}
