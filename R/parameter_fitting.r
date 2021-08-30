####PARAMTER FITTING

#' adds fitted parameter to constant data
#'
#' @param constpar data.frame containing static values used to run the model
#' @param fitpar named vector containing fitted values used to run the model
#' @return data.frame containing all values used to run the model
#' @export
get_run_parameters<-function(constpar,fitpar){
  constpar<-constpar[,! names(constpar) %in% names(fitpar)]
  constpar[,names(fitpar)]<-rep(fitpar,each=nrow(constpar))
  return (constpar)
}


#' Calculate RMSE for parameter fitting with isoplants
#'
#' @param fitpar data.frame containing fitted values used to run the model
#' @param fitparnames  names of these values
#' @param constpar static parmeters that are reruired to run the model but are not part of the fitting process
#' @param obsvalue name of the column in fixpar datafrme contoing the observation data
#' @param modvalue name of the output variable to be tested againtst
#' @param leaftemperture name of the output variable to be tested againtst
#' @return RMSE value
#' @export

plant18O_rmse<-function(fitpar,fitparnames,constpar,obsvalue, modvalue ='d18O_c',leaftemperture=FALSE){
  val<-constpar[,obsvalue] # Validation data
  names(fitpar)<-fitparnames
  runpar<-get_run_parameters(constpar,fitpar)
  if (leaftemperture){
    runpar<-calculate_Tleaf(runpar)
  }
  out<-plant18O_model(runpar,output = modvalue)# model data
  if (any(is.na(out))) {
    return(9999)}
  return(sqrt(mean((val - out) ^ 2, na.rm = T)))
}



#' Calculate RMSE for parameter fitting with isoplants
#'
#' @param fitparnames  vecort of value names to be fitted
#' @param data data.frame of plant18O model parameters (except the ones to be fitted)
#' including a column of observation data
#' @param obsvalue name of the column in fixpar datafrme contoing the observation data
#' @param modvalue name of the output variable to be tested againtst
#' @param method optimization methods. Supported methods are all methods in \code{optim} or 'GenSA' , which is more likely to find the glabal optimum in complex parameter space
##@param control optinal customized set of control parameters for the selected optimization method
#' @param pd optional parameter definition data.frame (pd <- get_parameter_definition()) with custom boundaries
#' @param leaftemperture boolean. Calculate leaf temperature using the tealeaves model. Parameters 'swrad', 'wind' and 'leafsize must be provided'
#' @param verbose boolean. Output some information while parameter fitting
#'
#' @return output from optimization routine including best fitting parameters
#' @export

fit_plant18O<-function(fitparnames, data, obsvalue ,modvalue ='d18O_c',  method = 'L-BFGS-B', pd=NULL, leaftemperture=FALSE, verbose =TRUE) {
  # Prepare initial parameters and limits
  if (is.null(pd)) pd<-get_parameter_definition()
  fitparlower     <- pd$lower  [pd$name %in% fitparnames]
  fitparupper     <- pd$upper  [pd$name %in% fitparnames]
  fitpardefault   <- pd$default[pd$name %in% fitparnames]
  fitparnames     <- pd$name[pd$name %in% fitparnames] # ensure order stays the same
  names(fitparlower)  <-fitparnames
  names(fitparupper)  <-fitparnames
  names(fitpardefault)<-fitparnames
  if (leaftemperture){
    if (!sum(c("swrad","wind","leafsize") %in% names(data))==3)
      stop("leaftemperture=TRUE, but not all nescessary parameters are provided")
  }
  if (verbose) for ( i in 1:length(fitparnames)) print(sprintf('parameter: %s  |   boundaries: %0.4f <= %s   >= %0.4f     | initial value: %s = %0.4f' , fitparnames[i], fitparlower[i],fitparnames[i],fitparupper[i], fitparnames[i],fitpardefault[i]))
  # Ensure values to be fitted are not in dataframe
  data<-data[,!names(data) %in% fitparnames]
  # optimization method
  if ( tolower(method) == 'gensa'){
    global.min <- 0
    tol        <- 1e-13
    #if (is.null(control))
    control <- list(threshold.stop=global.min+tol,verbose=TRUE,maxit=1000)
    out<-GenSA::GenSA( par = fitpardefault, fn = plant18O_rmse, lower = fitparlower,  upper = fitparupper, control=control,constpar=data,fitparnames=fitparnames,obsvalue=obsvalue , modvalue = modvalue , leaftemperture = leaftemperture)
  } else {
    out<-stats::optim(par = fitpardefault, fn = plant18O_rmse, lower = fitparlower,  upper = fitparupper, method = method, constpar=data,fitparnames=fitparnames,obsvalue=obsvalue , modvalue = modvalue , leaftemperture = leaftemperture)
  }
  if (verbose) print('done')
  return (out)
}
