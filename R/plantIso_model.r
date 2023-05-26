
# Péclet-modified Craig Gordon and cellulose model
# This is a general implementation of the Peclet-modified Craig Gordon model, offering a bunch of different options

#' (Péclet mofified) Craig-Gordon and cellulose model for calculation of plant tissue 18O enrichment
#' Calculates d18O of leaf water, cellulose and bulk tissue with different models based on the parameter provided
#' @param par a named vector or data.frame of parameter values (see details below)
#' @param output vector of output selected output values or 'all' (default)
#' @param addpar additional parameters, that will my join internally with the main parameter list. To be used for fixed parameters for model fitting/sensitivity analysis
#' @param verbose print messages while running
#' @return a dataframe of calculated values (see details below)
#' @keywords isotope, craig-gordon, leaf cellulose, leaf water
#' @details
#'
#' \bold{Input parameters}\cr A data.frame with the following numeric columns:
#'
#'\tabular{llll}{
#' \bold{Environment}\cr
#' \code{Tair} \tab Air Temperature (C)\cr
#' \code{RH} \tab Relative humidity (%)\cr
#' \tab or\cr
#' \code{ea} \tab atmospheric vapor pressure  (kPa) \cr
#' \code{P} \tab barometric pressure (kPa) \cr
#' \tab or\cr
#' \code{elevation} \tab elevation (m) \tab calculate P(elevation)\cr
#' \code{d18O_sw} \tab d18O soil water (permille)\cr
#' \code{D18O_wv} \tab D18O of water vapor over soil water (permille) \tab optional, if not supplied water vapour is assubed to be in equilibrium with soil water\cr
#' \cr
#' \bold{Leaf}\cr
#' \code{DTleaf} \tab delta leaf temperature (K)\cr
#' \tab or\cr
#' \code{Tleaf} \tab absolute leaf temperature \tab C\cr
#' \code{rb} \tab boundary resistance (m^2 s mol^-1) \cr
#' \code{gs} \tab stomatal conductance (mol m^-2 s^-1)\cr
#' \tab or\cr
#' \code{gsref} \tab ref. stomatal conductance @ 1kPa VPD (mol m^-2 s^-1) \tab calculate gs as gs(vpd)\cr
#' \code{Lm} \tab Peclet-model scaled path length (m) \tab use Peclet-model\cr
#' \tab or\cr
#' \code{phi} \tab Two pool model mixing parameter \tab use Two pool-model\cr
#' \code{px} \tab Proportion of O exchanged during cellulose synthesis\cr
#' \code{pex} \tab Proportion of unenriched xylem water in developing cell\cr
#' \code{ewc} \tab Biosynthetic fracination during cellulose synthesis (assumed to be 27 permil)\cr
#' \code{ecp} \tab offset from cellulose to bulk leaf material\cr
#' \cr
#' \bold{additional leaf temperature model \code{tealeaves} parameter}\cr
#' \code{swrad} \tab shortwave radiation \tab (Wm^-2)\cr
#' \code{wind} \tab windspeed  \tab (ms^-1)\cr
#' \code{ls} \tab leaf size  \tab (m)\cr
#'}
#'
#' @export
#' @examples
#'
#' # The model returns the ....
#' \dontrun{
#' result = plant18O_model(par)
#' }
plant18O_model <- function(par, addpar=NULL, output = "all", verbose = FALSE){ # Wrapper for backward compability
  plantIso_model(par, addpar, element="O", output, verbose)
}


#' (Péclet mofified) Craig-Gordon and cellulose model for calculation of plant tissue 18O enrichment
#' Calculates d18O of leaf water, cellulose and bulk tissue with different models based on the parameter provided
#' @param par a named vector or data.frame of parameter values (see details below)
#' @param output vector of output selected output values or 'all' (default)
#' @param addpar additional parameters, that will my join internally with the main parameter list. To be used for fixed parameters for model fitting/sensitivity analysis
#' @param verbose print messages while running
#' @return a dataframe of calculated values (see details below)
#' @keywords isotope, craig-gordon, leaf cellulose, leaf water
#' @details
#'
#' \bold{Input parameters}\cr A data.frame with the following numeric columns:
#'
#'\tabular{llll}{
#' \bold{Environment}\cr
#' \code{Tair} \tab Air Temperature (C)\cr
#' \code{RH} \tab Relative humidity (%)\cr
#' \tab or\cr
#' \code{ea} \tab atmospheric vapor pressure  (kPa) \cr
#' \code{P} \tab barometric pressure (kPa) \cr
#' \tab or\cr
#' \code{elevation} \tab elevation (m) \tab calculate P(elevation)\cr
#' \code{d18O_sw} \tab d18O soil water (permille)\cr
#' \code{D18O_wv} \tab D18O of water vapor over soil water (permille) \tab optional, if not supplied water vapour is assubed to be in equilibrium with soil water\cr
#' \cr
#' \bold{Leaf}\cr
#' \code{DTleaf} \tab delta leaf temperature (K)\cr
#' \tab or\cr
#' \code{Tleaf} \tab absolute leaf temperature \tab C\cr
#' \code{rb} \tab boundary resistance (m^2 s mol^-1) \cr
#' \code{gs} \tab stomatal conductance (mol m^-2 s^-1)\cr
#' \tab or\cr
#' \code{gsref} \tab ref. stomatal conductance @ 1kPa VPD (mol m^-2 s^-1) \tab calculate gs as gs(vpd)\cr
#' \code{Lm} \tab Peclet-model scaled path length (m) \tab use Peclet-model\cr
#' \tab or\cr
#' \code{phi} \tab Two pool model mixing parameter \tab use Two pool-model\cr
#' \code{px} \tab Proportion of O exchanged during cellulose synthesis\cr
#' \code{fH} \tab Proportion of unenriched xylem water in developing cell\cr
#' \code{epsA} \tab Biosynthetic isotopic fractionation between exchangeable H and medium water during autotrophic metabolism (assumed to be -171 permil)\cr
#' \code{epsH} \tab Biosynthetic isotopic fractionation between exchangeable H and medium water during heterotrophic metabolism (assumed to be 158 permil)\cr
#' \code{ecp} \tab offset from cellulose to bulk leaf material\cr
#' \cr
#' \bold{additional leaf temperature model \code{tealeaves} parameter}\cr
#' \code{swrad} \tab shortwave radiation \tab (Wm^-2)\cr
#' \code{wind} \tab windspeed  \tab (ms^-1)\cr
#' \code{ls} \tab leaf size  \tab (m)\cr
#'}
#'
#' @export
#' @examples
#'
#' # The model returns the ....
#' \dontrun{
#' result = plant2H_model(par)
#' }
plant2H_model <- function(par, addpar=NULL, output = "all", verbose = FALSE){ # Wrapper for backward compability
  plantIso_model(par, addpar, element="H", output, verbose)
}

#' (Péclet mofified) Craig-Gordon and based model for isotopic enrichment in plant tissue (18O and 2H)
#' @param par a named vector or data.frame of parameter values (see details below)
#' @param addpar additional parameters, that will my join internally with the main parameter list. To be used for fixed parameters for model fitting/sensitivity analysis
#' @param element additional parameters, run for oxygen 'O' or hydrogen 'H'
#' @param output vector of output selected output values or 'all' (default)
#' @param verbose print messages while running
#' @return a dataframe of calculated values (see details below)
#' @keywords isotope, craig-gordon, leaf cellulose, leaf water
#' @details
#' Main code for the model, for detailed input parameter descriptions see ?plant18O_model() or ?plant2H_model() functions
#' @export
#'
plantIso_model <- function(par, addpar=NULL, element="O", output = "all", verbose = FALSE) {
  # Define output names
  if (element=="H"){
      header <- c("ea", "ei", "ea_ei", "vpd", "eq", "ek", "gs", "E", "D", "pn", "D2H_e", "d2H_e", "D2H_lw", "d2H_lw", "D2H_c", "d2H_c","D2H_pt","d2H_pt")
      isotope<-"2H"
  } else{
      header <- c("ea", "ei", "ea_ei", "vpd", "eq", "ek", "gs", "E", "D", "pn", "D18O_e", "d18O_e", "D18O_lw", "d18O_lw", "D18O_c", "d18O_c","D18O_pt","d18O_pt")
      isotope<-"18O"
  }
  if (!is.null(addpar)) par <- set_parameters(addpar,par) # Allows to have some fixed parameters (used for optimization, sensitivity analisys,...)
  #Parse parameter
  if (output=='?') return (header) # Get header
  parameter_names <- names(par)
  valueNA <- rep(NA, length(par[, 1]))
  # Check required parameters and set options
  ## METEO   ########################
  # Tair
  if ("Tair" %in% parameter_names) {
    Tair <- par["Tair"]
  } else {stop("Missing required parameter: Tair")}
  # ea
  if ("ea" %in% parameter_names) {
    ea <- par["ea"]
  } else if ("RH" %in% parameter_names) {
    RH <- par["RH"]
    ea <- get_saturation_vapor_pressure (Tair) * (RH / 100)
  } else {
    stop("Could not calculate vapor pressure in air: provide either RH or ea")
  }
  # pressure
  if ("P" %in% parameter_names) {
    P <- par["P"]
  } else if ("elevation" %in% parameter_names) {
    P <- get_pressure_at_elevation(par["elevation"],par["Tair"])
    if (verbose) print("Calculating atmopheric pressure based on elevation")
  } else {
    stop("could not calculate atmopheric pressure: provide either P or elevation")
  }
  ## ISOTOPES ########################
  # d18O_sw /  d2H_sw
  dsw_name<-sprintf("d%s_sw",isotope)
  if (dsw_name %in% parameter_names) {
    dsw <- par[dsw_name]
  } else {
    stop(sprintf("Missing parameter %s",dsw_name))
  }
  # D18O_wv / # D2H_wv
  Dwv_name<-sprintf("D%s_sw",isotope)
  if (Dwv_name %in% parameter_names) {
    Dwv <- par[Dwv_name]
  } else {
    #stop("Missing parameter D18O_wv")
    if (verbose) message (sprintf("assuming equilibrium water vapour. %s = -eeq (tair)",Dwv_name))
    #Dwv <- -get_equilibrium_fractionation(Tair,element = element)
    Dwv <- delta2D(0,get_equilibrium_fractionation(Tair,element = element))
  }

  if ("Tleaf" %in% parameter_names) { #Use custom parametername for tealeaves calculated leaf temperature?
    Tleaf <- par["Tleaf"]
  } else if ("DTleaf" %in% parameter_names) {
    Tleaf <- Tair + par["DTleaf"]
  } else {
    stop("could not calculate Leaf temperature: provide either Tleaf or DTleaf")
  }

  # ei
  ei <- get_saturation_vapor_pressure (Tleaf)

  # vpd
  vpd <- ei - ea

  ## PLANT   ########################

  # Stomatal conductance gs
  if ("gs" %in% parameter_names) {
    gs <- par["gs"]
  } else if ("gsref" %in% parameter_names) {
    gsref <- par["gsref"]
    gs <- get_stomatal_conductance(gsref, vpd)
    if (verbose) print("Calculating stomatal conductance from gsref according to [Oren 1999]")
  } else {
    stop("could not calculate stomatal conductance (provide gs or gsref)")
  }

  # boundary layer resistance rb
  if ("rb" %in% parameter_names) {
      rb <- par["rb"]
    } else {
      stop("Missing parameter rb (boundary layer resistance)")
    }
  ## Isotope Fractination
  eq <- get_equilibrium_fractionation(Tleaf,element = element)
  D  <- get_H2O_diffusivity(Tleaf,element = element)

  rs <- valueNA
  E  <- valueNA
  ek <- valueNA

  # VPD > 0, we have transpiration
  rs[vpd > 0] <- 1 / gs[vpd > 0]
  ek[vpd > 0] <- get_kinetic_fractionation(rs[vpd > 0], rb[vpd > 0],element = element)
  E [vpd > 0] <- get_transpiration(vpd[vpd > 0], P[vpd > 0], rs[vpd > 0], rb[vpd > 0])
  # if vpd<=0
  ea[vpd <= 0] <- 1
  ei[vpd <= 0] <- 1
  ek[vpd <= 0] <- 0
  gs[vpd <= 0] <- NA
  E [vpd <= 0] <- 1E-16 # TO DO: What happens with transpiration in saturated air and how does this effect peclet mixing??

  # Evaporative front
  De <- apply_craig_gordon_model(eq, ek, Dwv, ea, ei)
  de <- D2delta(De,dsw)

  ## Mixing the evaporative front with leafwater
  if ("Lm" %in% parameter_names) {
    # PECLET MIXING
    Lm <- par["Lm"]
    if (verbose) print("Peclet mixing")
    pn <- get_peclet_number(Lm, E, D)
    Dlw <- apply_peclet_1D(De, pn)
  } else if ("phi" %in% parameter_names) {
    # TWO POOL
    phi <- par["phi"]
    if (verbose) print("two-pool model")
    pn <- NA
    Dlw <- apply_twopool_mixing(phi, De)
  } else {
    if (verbose) message("Could not calculate leaf water mixing. Provide Lm for Peclet model or phi for two-pool model")
    pn <- NA
    Dlw <- valueNA
  }

  # Leaf water
  dlw = D2delta(Dlw,dsw)

  # Cellulose
  if (element=="H"){
    if (!"epsA" %in% parameter_names) { epsA <- -171.0 }
    else  {epsA <- par["epsA"]}
    if (!"epsH" %in% parameter_names) { epsH <- 158.0 }
    else  {epsH <- par["epsH"]}
    if (!"px" %in% parameter_names | !"fH" %in% parameter_names) {
      Dc <- valueNA
      dc <- valueNA
      if (verbose) message("cannot calculate d2H cellulose. Provide values for px and fH")
    } else {
      px <- par["px"]
      fH <- par["fH"]
      dc <- get_cellulose2H(dlw,dsw,px,fH, epsA, epsH)
    }
  }else{ # Oxygen
    if (!"ewc" %in% parameter_names) { ewc<-27.0 }
    else  {ewc <- par["ewc"]}
    if (!"px" %in% parameter_names | !"pex" %in% parameter_names) {
      Dc <- valueNA
      dc <- valueNA
      if (verbose) message("cannot calculate d18O cellulose. Provide values for px and pex")
    } else {
      px <- par["px"]
      pex <- par["pex"]
      dc <- get_cellulose18O(dlw,dsw,px,pex,ewc)
    }
  }
  Dc <- delta2D(dc,dsw)
  # Bulk Materials
  if (!"ecp" %in% parameter_names) {
    Dpt <- valueNA
    dpt <- valueNA
    if (verbose) message("Parameter ecp missing. Cannot calculate values for bulk tissue")
  } else {
    ecp <- par["ecp"]
    Dpt <- Dc + ecp
    dpt <- D2delta(Dpt,dsw)
  }
  # Prepare Output
  #header <- c("ea", "ei", "ea_ei", "vpd", "eq", "ek", "gs", "E", "D", "pn", "D18O_e", "d18O_e", "D18O_lw", "d18O_lw", "D18O_c","d18O_c","D18O_pt","d18O_pt")
  values <- data.frame(ea, ei, ea / ei, vpd, eq, ek, gs, E, D, pn, De, de, Dlw, dlw, Dc, dc,Dpt,dpt)
  names(values) <- header
  if (output == "all") {
    return(values)
  }
  return(values[,output])
}




#' Basic description of input parameters for the leaf18O-model
#' This data.frame object returned by this function dataframe contains basic information about model inpiut parameters and is used for the check_parameters function
#' @param element additional parameters, run for oxygen 'O' or hydrogen 'H'
#' @return a dataframe containing vaild parameter names,  ranges, descriptors and default values
#' @keywords isotope, craig-gordon, leaf cellulose, leaf water, parameters
#' @export
#' @examples
#'
#' \dontrun{
#' parameters = get_parameter_definition ()
#' }

get_parameter_definition<-function(element="O"){
  parameter_definition<-data.frame(group=NA,name=NA,lower=NA,upper=NA,default=NA,pargroup=NA,unit=NA,description=NA)
  # items with the same pargroup variable are options to be set by providing one or the other parameter

  # Environment
  parameter_definition[ 1,]<-c('environment', 'Tair',      -30,  50,       20,  1,'[deg C]','Air Temperature')
  parameter_definition[ 2,]<-c('environment', 'RH',          0, 100,       70,  3,'[%]','relative humidity')
  parameter_definition[ 3,]<-c('environment', 'ea',        0.5,   5, 1.642629,  3,'[kPa]','atmospheric vapor pressure') # Default value calculated
  parameter_definition[ 4,]<-c('environment', 'P',         40, 110,  101.325,  4,'[kPa]','barometric pressure')
  parameter_definition[ 5,]<-c('environment', 'elevation', -50,8000,        0,  4,'[m]','elevation')
  if (element=="O") parameter_definition[ 6,]<-c('environment', 'd18O_sw',   -30,   0,      -20,  5,'[permil]','d18O soil water')
  if (element=="O") parameter_definition[ 7,]<-c('environment', 'D18O_wv',   -15,   0,-9.793879,  6,'[permil]','D18O of water vapor over soil water') # Default value calculated
  if (element=="H") parameter_definition[ 6,]<-c('environment', 'd2H_sw',   -120,   0,      -80,  5,'[permil]','d2H soil water')
  if (element=="H") parameter_definition[ 7,]<-c('environment', 'D2H_wv',   -100,   0, -85.0313,  6,'[permil]','D2H of water vapor over soil water') # Default value calculated


  # Leaf
  parameter_definition[ 8,]<-c('leaf', 'DTleaf',    -20,  20,        0,  2,'[deg C]','delta leaf temperature')
  parameter_definition[ 9,]<-c('leaf', 'Tleaf',     -30,  50,       20,  2,'[deg C]','absolute Leaf temperature')
  parameter_definition[10,]<-c('leaf', 'rb',        0.4,   6,        1,  7,'[m2 s mol-1]','boundary resistance')
  parameter_definition[11,]<-c('leaf', 'gs',          0,   2,      0.4,  8,'[mol m-2 s-1]','stomatal conductance')
  parameter_definition[12,]<-c('leaf', 'gsref',       0,   2,0.3304147,  8,'[mol m-2 s-1]','ref. stomatal conductance @ 1kPa VPD') # Default value calculated

  # LeafWater
  parameter_definition[13,]<-c('leaf', 'Lm',     0.0001,   2,     0.03,  9,'[m]','Peclet-model scaled path length')
  parameter_definition[14,]<-c('leaf', 'phi',         0,   1,      0.6,  9,'[]','Two pool model mixing parameter')

  # Cellulose
  if (element=="O") parameter_definition[15,]<-c('leaf', 'pex',          0,   1,      0.4, 10,'[]','Proportion of O exchanged during cellulose synthesis')
  if (element=="H") parameter_definition[15,]<-c('leaf', 'fH',          0,   1,      0.4, 10,'[]','Proportion of H exchanged during cellulose synthesis')
  parameter_definition[16,]<-c('leaf', 'px',           0,   1,        1, 11,'[]','Proportion of unenriched xylem water in developing cell')

  # Bulk leaf material
  parameter_definition[17,]<-c('leaf', 'ecp',     -10,  10,        0, 12,'[permil]','offset from cellulose to bulk leaf material')

  # Additional Parameters for use with tealeaves leaf temperature model
  parameter_definition[ 18,]<-c('leafTmodel', 'swrad',   0,   2500,1000,  2,'[Wm^-2]','shortwave (solar) radiation')
  parameter_definition[ 19,]<-c('leafTmodel', 'wind',    0,   20,     2,  2,'[m/s]','windspeed')
  parameter_definition[ 20,]<-c('leafTmodel', 'leafsize',    0.001,   0.3, 0.1,  2,'[m]','leaf size (used for leaf temperature model)')

  # Format
  parameter_definition[ , c(3:6)] <- apply(parameter_definition[ ,c(3:6)], 2,function(x) as.numeric(as.character(x)))

  return(parameter_definition)
}

#' Get a set of default environment and leaf parameters
#' @param n number of rows in the retuned data.frame
#' @param mode 'peclet' (default) or 'twopool', specify the mixing model used
#' @param element get parameters for oxygen 'O' or hydrogen 'H'
#' @param tealeaves add default parameters for leaf temperature calculation with the tealeaves model
#' @return a dataframe containing vaild parameter names,  ranges, descriptors and default values
#' @export
#'
get_default_parameters<-function(n=1,mode='peclet',element="O",tealeaves=FALSE){
  mix<- 13
  if (mode=='twopool') mix<- 14
  lt<-9
  if (tealeaves) lt<-c(18,19,20)
  pd <- get_parameter_definition(element)
  dp <- data.frame(matrix(rep(pd$default[c(1,2,4,6,lt,10,11,mix,15,16,17)],n),nrow=n,byrow=TRUE))
  names(dp)<-pd$name[c(1,2,4,6,lt,10,11,mix,15,16,17)]
  attr(dp,'units')<-pd$unit[c(1,2,4,6,lt,10,11,mix,15,16,17)]
  attr(dp,'element')<-element
  return (dp)
}


#' Check parameters for *plant18O_model()* or *plant2H_model()*
#'
#'This function does some basic range checking for the provided parameters and informs about missing parameters to run the model
#'correctly. It is intended to be used once before running the model with the parameter set.
#'
#' @param par a data.frame containting all required leaf18O parameter values
#' @param element specify whether the check should be performed for oxygen 'O' of hydrogen 'H' parameters
#' @param parameter_check defaults to get_parameter_definition(), but allows to set own valid parameter ranges
#' @return a list consisting of  1.) a general description of the tested parameters,
#'  2.) a data frame including details of the parameter check and  3.) a list with indexes to parameter sets that failed the check
#'
#' @keywords isotope, craig-gordon, leaf cellulose, leaf water, parameter-check
#' @export
#' @examples
#'
#' \dontrun{
#' check = check_parameters(par)
#' }

check_parameters <- function(par,element=NULL,parameter_check=NULL ) {
  if (is.null(element)) { # Try to get the element from the provided parameters
    if (sum(grepl("18O",names(par)))) element<-"O" else if (sum(grepl("2H",names(par)))) element<-"H" else stop("Please provide a valid argument for element ('O' or 'H')")
  }
  #ALL POSSIBLE PARAMETERS (SOME ARE REDUNDAT)
  if (is.null(parameter_check)) parameter_check<-get_parameter_definition(element)
  parnames<-names(par) # Get names of provided parameter table

  # Caluclate optional parmeter
  # pargroup 6 D18O_vp
  if (element=="O"){
    wvstring=''
    if (!('D18O_wv' %in% parnames) & ('Tair' %in% parnames)){
      #par$D18O_wv <- -get_equilibrium_fractionation(par$Tair,"O")
      par$D18O_wv <- delta2D(0,get_equilibrium_fractionation(par$Tair,"O"))
      wvstring<-'D18O of water vapor over soil water assumed to be in equilibrium'
    }
  }else{
    wvstring=''
    if (!('D2H_wv' %in% parnames) & ('Tair' %in% parnames)){
      #par$D2H_wv <- -get_equilibrium_fractionation(par$Tair,"H")
      par$D2H_wv  <- delta2D(0,get_equilibrium_fractionation(par$Tair,"H"))
      wvstring<-'D2H of water vapor over soil water assumed to be in equilibrium'
    }
  }
  parnames<-names(par) # Get names of provided parameter table
  parnames<-intersect(parnames,parameter_check$name) # drop all non-parameters from provided parameter table
  par<-par[,names(par) %in% parnames]
  np  <- length(parnames)
  npc <-nrow(par)
  parameter_check$check<-NA
  parameter_check$valid<-NA
  parameter_check$invalid<-NA
  parameter_check$na<-NA
  parameter_check$range<-NA
  parameter_check$type<-NA
  failed<-list()

  for (i in 1:np){
    pid<-which(parameter_check$name==parnames[i])
    parameter_check$na[pid]<-length(par[,i][is.na(par[,i])])
    if (min(par[,i],na.rm=TRUE)==max(par[,i],na.rm=TRUE)){
      parameter_check$range[pid]<-sprintf('%f',min(par[,i],na.rm=TRUE))
      parameter_check$type[pid]<-'constant'
    } else {
      parameter_check$range[pid]<-sprintf('%4.2f to %4.2f',min(par[,i],na.rm=TRUE),max(par[,i],na.rm=TRUE))
      parameter_check$type[pid]<-'variable'
    }
    check<-(par[,i]>=parameter_check$lower[pid] & par[,i]<=parameter_check$upper[pid])
    parameter_check$valid[pid]   <- sum(check)
    parameter_check$invalid[pid] <- length(check)-sum(check)
    #
    if (all(check)){
      parameter_check$check[pid]   <-'passed'
    }else{
      parameter_check$check[pid]<-'failed'
      failed[[parnames[i]]] <- which(check==FALSE)
    }
  }

  if (element=="O") if (length(wvstring)>0) parameter_check$type[which(parameter_check$name=='D18O_wv')]<-'calculated'
  if (element=="H") if (length(wvstring)>0) parameter_check$type[which(parameter_check$name=='D2H_wv')]<-'calculated'

  ## CHECK IF ALL REQUIRED PARAMETERS ARE SET
  setpar<-unique(parameter_check$pargroup[!is.na(parameter_check$check)])
  # Remove redundant parameters if one of the is set
  parameter_check<-parameter_check[!(is.na(parameter_check$check) & parameter_check$pargroup %in% setpar),]
  parameter_check$check[is.na(parameter_check$check)]<-'missing'

  # pargroup 2 specifies leaf temperature value/model
  if ((2 %in% setpar)){
    ltpar<-parameter_check$name[which(parameter_check$pargroup==2)]
    if (length(ltpar)==1) {
      if (ltpar == "Tleaf"| ltpar == "Dleaf") tstring<- 'using provided leaf temperature'
    } else if (length(ltpar)>=3) {if ("swrad" %in% ltpar & "wind" %in% ltpar &  "leafsize" %in% ltpar) tstring<- 'using leaf temperature model (tealeaves)'}
    else {  tstring<-'leaf temperature assumed to be equal to air temperature\n'}
  }

  # pargroup 9 specifies mixing model
  if ((9 %in% setpar)){
    mixpar<-parameter_check$name[which(parameter_check$pargroup==9)]
    if (length(mixpar)==2) {
      mixstring<-'redundant mixing parameter set, using peclet model\n'
    }else if (mixpar=='Lm'){
      mixstring<-'peclet mixing\n'
    }else{
      mixstring<-'two pool mixing\n'
    }
  }else{
    mixstring<-'no mixing parameter set\n'
  }
  parameter_check<-parameter_check[,-which(names(parameter_check)=='pargroup')] # this column is for internal use only
  header=c(
    sprintf('Provided %i set(s) for %i parameters:\t %s\n',npc,np,paste(parnames,collapse = ', ')),
    sprintf('Leaf temperature:\t %s\n',tstring),
    sprintf('Mixing:\t %s\n',mixstring)
    )
    if (wvstring!='') header=c(header,sprintf('Other:\t %s\n',wvstring))
  header=data.frame(model_options=header)
  output<-list(header,parameter_check,failed)
  names(output)<-c('model_options','model_parameters','errors')
  return(output)
}


################## Specific functions


#' atmospheric pressure at elevation
#'
#' calculate mean pressure based on mean sea level pressure (1013.25 mb) and elevation
#' this approximation doesn't take changing atmospheric pressure into account, but these changes have a very minor impact on calculated
#'  The main thing is to get the pressure changes due to elevation.
#'
#' @param elevation elevation in m
#' @param Tair Air temperature in C
#' @return a dataframe containing vaild parameter names,  ranges, descriptors and default values
#' @export
## ---- get_pressure_at_elevation
get_pressure_at_elevation <- function(elevation,Tair) { #[m],[deg Celsius]
  #atp <- 101.325 * exp(-1 * (elevation / 7990)) / 10 # [kPa]
  M<-0.02896 #[kg/mol] Molar mass of earth
  g<-9.807 # ms-2 gravitation
  Tair<-Tair+273.15
  R<- 8.3143 #[Nm mol-1 K-1] # universal Gas constant
  atp <- 101.325 * exp (-((M*g)/(R*Tair))*elevation)
  return(atp)
}
## ----


#' Saturation Vapor pressure
#'
#' @param T Temperature in Celsius
#' @return saturation Vapor pressure in kPa
#' \deqn{svp = h_\mathrm{vap} g_\mathrm{tw} d_\mathrm{wv}}{L = h_vap g_tw d_wv}
#' @export
# ? svp<-0.611121*exp( (18.678-(T/234.5))*(T/(257.14+T))) Alden eq.
## ---- get_saturation_vapor_pressure
get_saturation_vapor_pressure <- function(T) {
  svp <- 0.61365 * exp((17.502 * T / (240.97 + T))) # kPa
  return(svp)
}
## ----


#' Transpiration
#'
#' @param VPD Vapor pressure deficit
#' @param P Atmospheric pressure
#' @param rs stomatal resistance
#' @param rb boundary layer resistance
#' @return Leaf transpiration in mol m^-2 s^-2
#' @export
## ---- get_transpiration
get_transpiration <- function(VPD, P, rs, rb) {
  E <- (1 / (rs + rb)) * (VPD / P) # [mol m^-2 s^-2]
  return(E)
}
## ----


#' Craig-Gordon model for evaporative isotope enrichment
#'
#' Models the heavy isotope enchichment due to evaporation
#' @param eq equilibrium fractination (permil)
#' @param ek kinetic fractination (permil)
#' @param Dwv Water vapor isotope signal over source water (oxygen or hydrogen, depending on application)
#' @param ea atmospheric water pressure
#' @param ei leaf internal water pressure (assumed to be saturated)
#' @return D18O at evaporative front over source water
## ---- craig_gordon_model
apply_craig_gordon_model <- function(eq, ek, Dwv, ea, ei) {
  eaei <- ea / ei
  De <- ((1 + eq / 1000) * ((1 + ek / 1000) * (1 - eaei) + eaei * (1 + Dwv / 1000)) - 1) * 1000
  return(De)
}
## ----


#' Equilibrium fractionation
#'
#' Temperature dependence of equilibrium fractinations see Majoube 1971.
#' implemented as described in review by Cernusak et al. 2016
#' @param Tleaf Leaf temperature in Celsius
#' @param element specify whether it should be calculated for oxygen 'O' of hyrdogen 'H'
#' @return equilibrium_fractionation in permil
#' @references
#'
#' Cernusak, L.A., Barbour, M.M., Arndt, S.K., Cheesman, A.W., English, N.B., Feild, T.S., Helliker, B.R., Holloway-Phillips, M.M., Holtum, J.A., Kahmen, A. and McInerney, F.A., 2016. Stable isotopes in leaf water of terrestrial plants. Plant, Cell & Environment, 39(5), pp.1087-1102.
#' @export
## ---- get_equilibrium_fractionation
get_equilibrium_fractionation <- function(Tleaf, element="O") { #[C]
  T <- Tleaf + 273.15
  if (element=="H"){
    es <- (exp((24.844 / T^2) * 10^3 - (76.248 / T) + 52.612 * 10^(-3)) - 1) * 1000 # [per mil] for Hydrogen
  } else{
    es <- (exp((1.137 / T^2) * 10^3 - (0.4156 / T) -2.0667 * 10^(-3)) - 1) * 1000 # [per mil] for Oxygen
  }
  return(es)
}
## ----
#es<- 2.644-3.206(10^3/T)+1.534*(10^6/T^2) #[per mil] [as described in review by Barbour 2007]


#' Kinetic fractionation
#'
#' Resistance (stomatal/boundary layer) dependence of kinetic fractinations.
#' implemented as described in review by Cernusak et al. 2016:
#' ek <- (28 x rs + 19 x rb) / (rs + rb)
#' other suggest:
#' ek<-(32 x rs+22 x rb)/(rs+rb) #described in review by Barbour 2007
#' ek<-(32 x rs+21 x rb)/(rs+rb) #as described in paper by Song 2013
#' @param rs stomatal resistance
#' @param rb boundary layer resistance
#' @param element specify whether it sould be calculated for oxygen 'O' of hyrdogen 'H'
#' @return kinetic_fractionation in permil
#' @references
#'
#' Cernusak, L.A., Barbour, M.M., Arndt, S.K., Cheesman, A.W., English, N.B., Feild, T.S., Helliker, B.R., Holloway-Phillips, M.M., Holtum, J.A., Kahmen, A. and McInerney, F.A., 2016. Stable isotopes in leaf water of terrestrial plants. Plant, Cell & Environment, 39(5), pp.1087-1102.
#' @export
## ---- get_kinetic_fractionation
get_kinetic_fractionation <- function(rs, rb,element="O") { # [m^2s mol^-1]
  if (element=="H"){
    ek <- (25 * rs + 17 * rb) / (rs + rb) # [per mil] for Hydrogen
  } else{
    ek <- (28 * rs + 19 * rb) / (rs + rb) # [per mil] for Oxygen
  }
  return(ek)
}
## ----


#' Stomatal conductance
#'
#' Stomatal conductance as calculated by the formula of Oren 1999
#' @param gsref stomatal conductance reference, measured at 1 kPa VPD
#' @param vpd vapour pressure deficit
#' @return stomatal_conductance in mol m-2 s-1
#' @export
## ---- get_stomatal_conductance
get_stomatal_conductance <- function(gsref, vpd) { #gs@1kPA VPD [mol m-2 s-1], [kPa]
  vpd[vpd<=0]<-NA # Check for invalid VPD values
  gs <- gsref - 0.6 * gsref * log(vpd) #[mol m-2 s-1]
  return(gs)
}
## ----

#' Peclet number
#'
#' @param Lm scaled path lenth, i m
#' @param E  transpiration in mol m-2 s-1
#' @param D  Diffusivity
#' @return peclet number
#' @export
## ---- get_peclet_number
get_peclet_number <- function(Lm, E, D) { #[m],[mol m-2 s-1],[]
  C <- 5.55 * 10^4 # molar concentration of water [mol m^-3]
  pn <- (Lm * E) / (C * D) #[]
  return(pn)
}
## ----


#' 1D Peclet mixing
#'
#' reduces the d18O at evaporative front based on the peclet value (Farquhar & Lloyd 1993)
#' @param De D18O or D2H value at evaporative front
#' @param pn  peclet number
#' @return D18O of total leaf water from d18O value at evaporative front
#' @export
## ---- apply_peclet_1D
apply_peclet_1D <- function(De, pn) {  #[permil],[]
  Dlw <- (De * (1 - exp(-pn))) / pn #[permil]
  return(Dlw)
}
## ----


#' Diffusivity
#'
#' Temperature (C) dependence of diffusivity of heavy isotopologue in water (H2 and
#' 18O) (Cuntz et al. 2007) as described in review by Cernusak et al. 2016
#'
#' @param T Temperature in Celsius
#' @param element specify whether it sould be calculated for oxygen 'O' of hyrdogen 'H'
#' @return diffusivity of heavy isotopologue  in water (H2 18O)
#' @export
## ---- get_H2O_diffusivity
get_H2O_diffusivity <- function(T,element="O") {
  if (element=="H"){ # Hydrogen
    D <- (98.5 * 10^-9) * exp(-577 / (T + 128)) # [m^2 s^-1]
  }else{ # Oxygen
    D <- (97.5 * 10^-9) * exp(-577 / (T + 128)) # [m^2 s^-1]}
  }
  return(D)
}

## ----
# D<- (97.5*10^-8)*exp(-577/(T+128))  #[m^2 s^-1] [Typo corrected version suggested by Keel 2016]
# D<- 119*10^-9*exp(-637/(136.15+T))  #[m^2 s^-1] [as described in paper by Song 2013]
# get_cellulose18O <- function(pex, px, d18O_sw, d18O_lw) { # [], [], [permil], [permil]
#   ewc <- 27.0 # [permil] fractination between source water and primary products of photosyntesis
#   d18Oc <- pex * px * (d18O_sw + ewc) + (1 - pex * px) * (d18O_lw + ewc) # [permil]
#   return(d18Oc)
#}

#' Leaf water to cellulose 18O
#'
#' @param d18O_lw  d180 of leaf water
#' @param d18O_sw  d180 of source (xylem) water
#' @param px fraction of unenriched water in this water pool used for cellulose synthesis
#' @param pex exchange of leaf water signal in carbohydrates with local water before incorporation into the cellulose polymer
#' @param ewc  Equilibrium fractionation between carbonyloxygen and medium water (default 27 permil)
#' @return Returns Oxygen isotope composition of cellulose (permil)
#' @export
## ---- get_cellulose18O

get_cellulose18O <- function(dlw,dsw, px, pex ,ewc = 27.0) { #[permil],[permil] [], [], [permil]
  # Same Formulation as for 2H, except epsA and espH are both equal to ewc
  f<-pex
  epsA<-ewc
  epsH<-ewc
  dc <- (1-f)*(dlw*(1+epsA/1000 )+epsA)+f*((px*dsw +(1-px)*dlw )*(1+epsH/1000 )+epsH) # Full formulation including second order terms
  return(dc)
}
## ----


#' Leaf water to cellulose 2HO
#'
#' @param d2H_lw  d2H of leaf water
#' @param d2H_sw  d2H of source (xylem) water
#' @param px fraction of unenriched water in this water pool used for cellulose synthesis
#' @param fH Proportion of H exchange during heterotrophic reactions
#' @param epsA Biosynthetic fractionation between exchangeable H and medium water during autotrophic metabolism
#' @param epsB  Biosynthetic fractionation between exchangeable H and medium water during heterotrophic metabolism
#' @return Returns Hydrogen isotope composition of cellulose (permil)
#' @export
## ---- get_cellulose2H
get_cellulose2H <- function(dlw, dsw, px, fH, epsA= -171 ,epsH= 158) { #[permil],[permil] [], [], [permil],[permil]
  f<-fH
  dc<-(1-f)*(dlw*(1+epsA/1000 )+epsA)+f*((px*dsw +(1-px)*dlw )*(1+epsH/1000 )+epsH)  # Full formulation including second order terms
  return(dc)
}
## ----


#' Transpiration dependent scaled path length
#'
#' calculate path length after song et al. 2013
#' @param E Transpiration
#' @return path length value in m
#' @export
## ---- get_pathlength
get_pathlength <- function(E) { #[mol m-2 s-1]
  pl <- 2.36 * 10^-5 * E^-1.2 # [m]
  return(pl)
}
## ----


#' Equilibrium fractionation between carbonyl oxygen and medium water
#'
#' calculate ewc after by Sternberg & Ellsworth (2011)
#' @param Tleaf Leaf temperature
#' @return ewc value in permil
#' @export
## ---- get_ewc
get_ewc <- function(Tleaf) { #[mol m-2 s-1]
  ewc=0.0084*(Tleaf^2)-0.51*Tleaf+33.172
  return(ewc)
}
## ----


#' Two pool mixing
#' two pool model phi is the proportion of unenriched water
#' @param phi mixing coefficient
#' @param De D18O or D2H of water at the evaporative front
#' @return D18O or D2H of total leaf water from d18O value at evaporative front
#' @export
## ---- apply_twopool_mixing
apply_twopool_mixing <- function(phi, De) { #[], [permil]
  Dlw <- (1 - phi) * De # [permil]
  return(Dlw)
}
## ----

#' sample parameters from specified distribution for leaf18O model
#'
#' @param n number of parameter sets
#' @param parranges  data.frame (name,pdist,pdm,pdv)
#'
#' @return data.frame containing parameters sampled from the defined distrubution
#' @importFrom  stats rnorm runif
#' @export
#' @examples
#' \dontrun{
#' # RUN 1000 more time on a randomly varied parameters
#' n<-10000
#' rndpar<-get_randomized_parameters(n,parranges)
#' check_parameters(rndpar)[[2]]
#' modelout<-plant18O_model(rndpar,output='d18O_pt')
#' quantile(modelout,c(0.025,0.25,0.5,0.75,0.975))
#' mean(modelout[!is.infinite(modelout)])
#'}
#'
get_randomized_parameters<-function(n,parranges){
  rndpv<-data.frame(do.call(cbind,lapply(1:nrow(parranges),FUN= function (i) {
    if (parranges$pdist[i]=='norm') rnorm(n,mean=parranges$pdm[i],sd=parranges$pdv[i])
    else if (parranges$pdist[i]=='unif') runif(n,min=parranges$pdm[i]-parranges$pdv[i],max=parranges$pdm[i]+parranges$pdv[i])
    else if (parranges$pdist[i]=='minmax') runif(n,min=parranges$pdm[i],max=parranges$pdv[i])
    else rep(parranges$pdm[i],n)
    })))
  names(rndpv)<-parranges$name
  return(rndpv)
}

#' Set parameters by replacing or adding parameters to parameter table
#'
#'if multiple values are assigned to d parameter table with only a single row, the nuber of rows will be extended to match the number of rows in the newparameter table.
#' @param newvalues data.frame with new parameter values and columnnames machting a model parameter
#' @param parameter_table  data.frame of parameters (same format as returned by get_default_parameters())
#'
#'
#' @return updated parameter data.frame
#' @export

set_parameters<-function(newvalues,parameter_table){

  if (nrow(parameter_table)==1) {parameter_table<-parameter_table[rep(seq_len(nrow(parameter_table)), each = nrow(newvalues)), ]}
  if (nrow(newvalues)==1) {newvalues<-newvalues[rep(seq_len(nrow(newvalues)), each = nrow(parameter_table)), ]}

  if (nrow(parameter_table)!=nrow(newvalues)) {stop('data.frames have different numbers of rows, cannot merge')}
  update<-intersect(names(newvalues),names(parameter_table))
  if (length(update)>0) parameter_table[,update]<-newvalues[,update]

  add<-setdiff(names(newvalues),names(parameter_table))
  if (length(add)>0) {
    parameter_table<-cbind(parameter_table,newvalues[,add])
    #names(parameter_table)[ (ncol(parameter_table)-length(add)) : ncol(parameter_table) ]<- names(newvalues[,add])
    }
  return (parameter_table)

}


#' get leaf temperature for
#' Interface with the tealeaves package in order to calculate leaf temperaure and boundary conductance. Sets air temperauture,
#' barometric pressure, relative humidity and stomatal conductance from isoplants parameters
#' @param isoplantspar data.frame of randomized model parameters for the siensitvity analysis
#' @param tealeaves_enviro_par (optional) customized tealeaves parameterset
#' @param tealeaves_leaf_par   (optional) customized tealeaves parameterset
#' @param tealeaves_constants  (optional) customized tealeaves parameterset
#' @return vector containing Tleaf and rb value
#' @export

tealeaves_getValues<-function(isoplantspar,tealeaves_enviro_par=NULL,tealeaves_leaf_par=NULL, tealeaves_constants=NULL ){
  #TEALEAVES default parameters.
  if (is.null(tealeaves_leaf_par))   tealeaves_leaf_par   <- tealeaves::make_leafpar()   # leaf parameters
  if (is.null(tealeaves_enviro_par)) tealeaves_enviro_par <- tealeaves::make_enviropar() # environmental parameters
  if (is.null(tealeaves_constants))  tealeaves_constants  <- tealeaves::make_constants() # physical constants
  tealeaves_enviro_par$T_air<-units::set_units(isoplantspar$Tair+273.15,'K')
  tealeaves_enviro_par$P<-units::set_units(isoplantspar$P,'kPa')
  tealeaves_enviro_par$RH<-units::set_units(isoplantspar$RH/100,'1')
  if ("swrad" %in% names(isoplantspar)) tealeaves_enviro_par$S_sw <- units::set_units(isoplantspar$swrad,'W/m^2')
  if ("wind" %in% names(isoplantspar))  tealeaves_enviro_par$wind <- units::set_units(isoplantspar$wind,'m/s')
  if ("leafsize" %in% names(isoplantspar)) tealeaves_leaf_par$leafsize <- units::set_units(isoplantspar$leafsize,'m')
  tealeaves_leaf_par$g_sw<-tealeaves::convert_conductance(units::set_units(isoplantspar$gs, "mol/m^2/s"),Temp = tealeaves_enviro_par$T_air, P = tealeaves_enviro_par$P)[[2]]
  t<-tealeaves::tleaf(tealeaves_leaf_par, tealeaves_enviro_par, tealeaves_constants, quiet = TRUE)[,1]
  gwb<-tealeaves:::.get_gbw(t, "lower", c(tealeaves_constants, tealeaves_enviro_par, tealeaves_leaf_par), FALSE)
  gwb<-units::set_units(gwb,"m/s") # fix unit form m^2/m/s to m/s
  gbw<-tealeaves::convert_conductance(gwb, Temp = tealeaves_enviro_par$T_air, P = tealeaves_enviro_par$P)[[3]]
  return(c(as.numeric(t)-273.15, 1/as.numeric(gbw)))
}



#' get leaf temperature for
#' Interface with the tealeaves package in order to calculate leaf temperaure and boundary conductance.
#' Wrapper for tealeaves_getValues
#' @param isoplantspar data.frame of randomized model parameters for the siensitvity analysis
#' @return input datafrme with updated Tleaf and boundary resistance values
#' @export

calculate_Tleaf<-function(isoplantspar){
  #BASIC PARAMTER CHECK
  tmp<-sapply(1:nrow(isoplantspar),FUN= function(i) tealeaves_getValues(isoplantspar[i,]))
  isoplantspar$Tleaf<-tmp[1,]
  isoplantspar$rb<-tmp[2,]
  return(isoplantspar)
}



# # SET PARAMETER UNCERTAINTY
# #names(specimen_model_par)[names(specimen_model_par) %in% get_parameter_definition()$name][-1]
# parranges<-data.frame(
#   name = c( "Tair",   "RH",      "P",    "d18O_sw", "D18O_wv", "DTleaf",  "px",     "pex",     "rb",    "Lm",     "gs",      "ecp"),
#   pdist = c( 'norm',   'norm',    'norm',  'norm',    'norm',    'norm',    "unif",   "unif",    "unif",  "unif",   "unif",    "unif" ),
#   pdm   = c(       1,       1,         1,       1,         1,         1,         1,        1,         1,       1,        1,         1 ),
#   pdv   = c(       1,       1,         1,       1,         1,         1,         1,        1,         1,       1,        1,         1 )
#   ,stringsAsFactors = FALSE)
#
# #set means (done for each specimen)
# parranges$pdm<-merge(parranges,get_parameter_definition()[,c('name','default')],by='name')$default[order(order(parranges$name))]
#

