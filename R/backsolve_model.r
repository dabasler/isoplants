
#' Solves the PMCG model for pex
#'
#' This function solves the PMCG model for pex, provided that the nescessary parameters are supplied to the function
#'
#' \bold{Input parameters}\cr A data.frame with the following numeric columns:
#'
#'\tabular{llll}{
#' \bold{Environment}\cr
#' \code{Tair} \tab Air Temperature (C)\cr
#' \code{RH} \tab Relative humidity (%)\cr
#' \code{P} \tab barometric pressure (kPa) \cr
#' \code{d18O_sw} \tab d18O soil water (permille)\cr
#' \cr
#' \bold{Leaf}\cr
#' \code{Tleaf} \tab absolute leaf temperature \tab C\cr
#' \code{rb} \tab boundary resistance (m^2 s mol^-1) \cr
#' \code{gs} \tab stomatal conductance (mol m^-2 s^-1)\cr
#' \code{Lm} \tab Peclet-model scaled path length (m) \tab use Peclet-model\cr
#' \code{px} \tab Proportion of O exchanged during cellulose synthesis\cr
#' \code{ewc} \tab Biosynthetic fracination during cellulose synthesis (assumed to be 27 permil)\cr
#' \cr
#' \bold{Isotope signal}\cr
#' \code{d18O_cellulose} \tab Measured d18O_cellulose values, columnname may be specified as function parameter\cr
#'}
#'
#'
#' @param par a named vector or data.frame of model parameters (see description oablve) and (observed) d18O_cellulose values
#' @param d18Oc_colname name of column in the input parameter vector/dataframe containg the d18O_cellupose values used for calculation  (default 'd18O_c')
#'
#' @return vector of calculated pex values
#' @export
#' @examples
#'
#' \dontrun{
#' #Test
#' par <- get_default_parameters(n=3)
#' par$pex <- seq(0.3,0.5,0.1)                            # set test pex values
#' par$d18O_c <- plant18O_model(par,output = 'd18O_c')    # calculate d18O with test pex values
#' plant18O_solve_pex(par)                                # backsolve model for pex
#' }
#'
#'
plant18O_solve_pex<-function(par, d18Oc_colname = "d18O_c"  ){
  # Get Basic Values: THis is a short form of the functions in plant18O_model()
  if (! d18Oc_colname %in% names(par)) stop ("no value for d18O_cellulose. Provide valid column name for d18Oc_colname")
  required_par <- c("Tair", "Tleaf", "RH",  "P", "px", "gs", "rb", "Lm", "d18O_sw")
  if (sum(!required_par %in% names(par))) stop (sprintf( "Missing parameters: %s",paste(required_par[!required_par %in% names(parameters)],collapse = ",")))

  # Backward compatibility - supply D18O_c directly
  if (substr(d18Oc_colname,1,1)=="D") D18O_c<- par[,d18Oc_colname]
  else D18O_c<- delta2D(par[,d18Oc_colname],par[,"d18O_sw"])

  valueNA <- rep(NA, length(par[, 1]))

  Tair  <- par[,"Tair"]
  Tleaf <- par[,"Tleaf"]
  RH    <- par[,"RH"]
  P     <- par[,"P"]
  #P <- get_pressure_at_elevation(par["elevation"],par["Tair"])
  px    <- par[,"px"]
  gs    <- par[,"gs"]
  rb    <- par[,"rb"]
  gb    <- 1/rb
  Lm    <- par[,"Lm"]

  ea <- get_saturation_vapor_pressure (Tair) * (RH / 100)
  ei <- get_saturation_vapor_pressure (Tleaf)
  D18O_wv <- -get_equilibrium_fractionation(Tair)
  vpd <- ei - ea
  eq <- get_equilibrium_fractionation(Tleaf)
  D <- get_H2O_diffusivity(Tleaf)
  C <- 5.55 * 10^4
  ## Isotope Fractination
  rs <- valueNA
  E  <- valueNA
  ek <- valueNA
  # VPD > 0, we have transpiration
  rs[vpd > 0] <- 1 / gs[vpd > 0]
  ek[vpd > 0] <- get_kinetic_fractionation(rs[vpd > 0], rb[vpd > 0])
  E [vpd > 0] <- get_transpiration(vpd[vpd > 0], P[vpd > 0], rs[vpd > 0], rb[vpd > 0])
  # if vpd<=0
  ea[vpd <= 0] <- 1
  ei[vpd <= 0] <- 1
  ek[vpd <= 0] <- 0
  gs[vpd <= 0] <- NA
  E [vpd <= 0] <- 1E-16 # TO DO: What happens with transpiration in saturated air and how does this effect peclet mixing??
  ewc<-27

  # ## USE SAGE TO SOLVE FOR GS LS AND PEX: https://sagecell.sagemath.org/
  # var('ei','ea','gs', 'gb','P','Lm', 'C', 'D','eq','D18O_wv','D18O_lw','D18O_c','ewc','px','pex')
  # VPD    = ei-ea
  # eaei   = ea / ei
  # E      = 1/(1/gs+1/gb)*(VPD/P)
  # pn     = (Lm*E)/(C*D)
  # ek     = (28*(1/gs)+19*(1/gb))/((1/gs)+(1/gb))
  # D18O_e  = ((1 + eq / 1000) * ((1 + ek / 1000) * (1 - eaei) + eaei * (1 + D18O_wv / 1000)) - 1) * 1000
  # D18O_lw = (D18O_e * (1 - exp(-pn))) / pn
  # print (solve([D18O_c == D18O_lw * (1 - pex * px) +  ewc], pex))

  pex <-(
        ((28000*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
           + 1000*(C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                   - 28*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea
           + (1028*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
              + (C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                 - 28*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea)*eq)*gb
          + (19000*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
             + 1000*(C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                     - 19*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea
             + (1019*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                + (C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                   - 19*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea)*eq
             - 1000*(D18O_c*Lm*ea*ei - D18O_c*Lm*ei^2 - (Lm*ea*ei - Lm*ei^2)*ewc)*gb)*gs
          )/
        (((28000*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
            + 1000*(C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                    - 28*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea
            + (1028*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
               + (C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                  - 28*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea)*eq)*gb
           + (19000*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
              + 1000*(C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                      - 19*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea
              + (1019*C*D*P*ei*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                 + (C*D*D18O_wv*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1)
                    - 19*C*D*(exp(Lm*ea/(C*D*P/gb + C*D*P/gs) - Lm*ei/(C*D*P/gb + C*D*P/gs)) - 1))*P*ea)*eq)*gs)*px
          )
    )

  return (pex)
}

#' Solves the PMCG model for gs
#'
#' This function solves the PMCG model for gs, provided that the necessary parameters are supplied to the function
#'
#' \bold{Input parameters}\cr A data.frame with the following numeric columns:
#'
#'\tabular{llll}{
#' \bold{Environment}\cr
#' \code{Tair} \tab Air Temperature (C)\cr
#' \code{RH} \tab Relative humidity (%)\cr
#' \code{P} \tab barometric pressure (kPa) \cr
#' \code{d18O_sw} \tab d18O soil water (permille)\cr
#' \cr
#' \bold{Leaf}\cr
#' \code{Tleaf} \tab absolute leaf temperature \tab C\cr
#' \code{rb} \tab boundary resistance (m^2 s mol^-1) \cr
#' \code{gs} \tab stomatal conductance (mol m^-2 s^-1)\cr
#' \code{Lm} \tab Peclet-model scaled path length (m) \tab use Peclet-model\cr
#' \code{px} \tab Proportion of O exchanged during cellulose synthesis\cr
#' \code{pex} \tab Proportion of unenriched xylem water in developing cell\cr
#' \code{ewc} \tab Biosynthetic fracination during cellulose synthesis (assumed to be 27 permil)\cr
#' \cr
#' \bold{Isotope signal}\cr
#' \code{d18O_cellulose} \tab Measured d18O_cellulose values, columnname may be specified as function parameter\cr
#'}
#'
#'
#' @param par a named vector or data.frame of model parameters (see description oablve) and (observed) d18O_cellulose values
#' @param d18Oc_colname name of column in the input parameter vector/dataframe containg the d18O_cellupose values used for calculation  (default 'd18O_c')
#'
#' @return vector of calculated gs values
#' @export
#' @examples
#'
#' \dontrun{
#' #Test
#' par <- get_default_parameters(n=3)
#' par$gs <- seq(0.3,0.5,0.1)                             # set test gs values
#' par$d18O_c <- plant18O_model(par,output = 'd18O_c')    # calculate d18O with test gs values
#' plant18O_solve_gs(par)                                 # backsolve model for gs
#' }
#'
#'
plant18O_solve_gs<-function(par, d18Oc_colname = "d18O_c"  ){
  # Get Basic Values: THis is a short form of the functions in plant18O_model()
  if (! d18Oc_colname %in% names(par)) stop ("no value for d18O_cellulose. Provide valid column name for d18Oc_colname")
  required_par <- c("Tair", "Tleaf", "RH",  "P", "px", "gs", "rb", "Lm", "d18O_sw")
  if (sum(!required_par %in% names(par))) stop (sprintf( "Missing parameters: %s",paste(required_par[!required_par %in% names(parameters)],collapse = ",")))

  # Backward compatibility - supply D18O_c directly
  if (substr(d18Oc_colname,1,1)=="D") {
    D18O_c<- par[,d18Oc_colname]
  } else {
    D18O_c<- delta2D(par[,d18Oc_colname],par[,"d18O_sw"])
    }

  valueNA <- rep(NA, length(par[, 1]))
  Tair  <- par[,"Tair"]
  Tleaf <- par[,"Tleaf"]
  RH    <- par[,"RH"]
  P     <- par[,"P"]
  #P <- get_pressure_at_elevation(par["elevation"],par["Tair"])
  px    <- par[,"px"]
  pex   <- par[,"pex"]
  gs    <- par[,"gs"]
  rb    <- par[,"rb"]
  gb    <- 1/rb
  Lm    <- par[,"Lm"]

  ea <- get_saturation_vapor_pressure (Tair) * (RH / 100)
  ei <- get_saturation_vapor_pressure (Tleaf)
  D18O_wv <- -get_equilibrium_fractionation(Tair)
  vpd <- ei - ea
  eq <- get_equilibrium_fractionation(Tleaf)
  D <- get_H2O_diffusivity(Tleaf)
  C <- 5.55 * 10^4
  ## Isotope Fractination
  rs <- valueNA
  E  <- valueNA
  ek <- valueNA
  # VPD > 0, we have transpiration
  rs[vpd > 0] <- 1 / gs[vpd > 0]
  ek[vpd > 0] <- get_kinetic_fractionation(rs[vpd > 0], rb[vpd > 0])
  E [vpd > 0] <- get_transpiration(vpd[vpd > 0], P[vpd > 0], rs[vpd > 0], rb[vpd > 0])
  # if vpd<=0
  ea[vpd <= 0] <- 1
  ei[vpd <= 0] <- 1
  ek[vpd <= 0] <- 0
  gs[vpd <= 0] <- NA
  E [vpd <= 0] <- 1E-16 # TO DO: What happens with transpiration in saturated air and how does this effect peclet mixing??

  if ("ewc" %in% names(par)) ewc <- par[,"ewc"]
  else ewc<-27

  # ## USE SAGE TO SOLVE FOR GS LS AND PEX: https://sagecell.sagemath.org/
  # var('ei','ea','gs', 'gb','P','Lm', 'C', 'D','eq','D18O_wv','D18O_lw','D18O_c','ewc','px','pex')
  # VPD    = ei-ea
  # eaei   = ea / ei
  # E      = 1/(1/gs+1/gb)*(VPD/P)
  # pn     = (Lm*E)/(C*D)
  # ek     = (28*(1/gs)+19*(1/gb))/((1/gs)+(1/gb))
  # D18O_e  = ((1 + eq / 1000) * ((1 + ek / 1000) * (1 - eaei) + eaei * (1 + D18O_wv / 1000)) - 1) * 1000
  # D18O_lw = (D18O_e * (1 - exp(-pn))) / pn
  # print (solve([D18O_c == D18O_lw * (1 - pex * px) +  ewc], gs))
  gs <-
    -((((pex*px - 1)*eq*gb + 1000*(pex*px - 1)*gb)*D18O_wv*P*ea
       - 4*(7*((pex*px - 1)*eq*gb + 1000*(pex*px - 1)*gb)*ea - (257*(pex*px - 1)*eq*gb + 7000*(pex*px - 1)*gb)*ei)*P
       )*C*D*exp((Lm*ea - Lm*ei)*gb*gs/(C*D*P*gb + C*D*P*gs))
      - (((pex*px - 1)*eq*gb + 1000*(pex*px - 1)*gb)*D18O_wv*P*ea
         - 4*(7*((pex*px - 1)*eq*gb + 1000*(pex*px - 1)*gb)*ea - (257*(pex*px - 1)*eq*gb + 7000*(pex*px - 1)*gb)*ei)*P
         )*C*D
      )/(
        (((pex*px - 1)*eq + 1000*pex*px - 1000)*D18O_wv*P*ea
         - (19*((pex*px - 1)*eq + 1000*pex*px - 1000)*ea - (1019*(pex*px - 1)*eq + 19000*pex*px - 19000)*ei)*P
         )*C*D*exp((Lm*ea - Lm*ei)*gb*gs/(C*D*P*gb + C*D*P*gs))
        - (((pex*px - 1)*eq + 1000*pex*px - 1000)*D18O_wv*P*ea
           - (19*((pex*px - 1)*eq + 1000*pex*px - 1000)*ea - (1019*(pex*px - 1)*eq + 19000*pex*px - 19000)*ei)*P
           )*C*D
        + 1000*(ea*ei*gb - ei^2*gb)*D18O_c*Lm - 1000*(ea*ei*ewc*gb - ei^2*ewc*gb)*Lm
        )

  return(gs)
}




