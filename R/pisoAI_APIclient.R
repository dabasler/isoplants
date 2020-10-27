#'  Get multi site data through the Piso.AI API
#'
#' [Piso.AI](https://isotope.bot.unibas.ch/PisoAI/) is tool for predicting monthly timeseries of oxygen and hydrogen isotope values of precipitation that uses a machine learning model trained on geographic and climate data.
#'
#' @param location data.frame object with colums 'site','latitude','longitude'. Colums 'elevation' and 'year' are optional
#' @param years single year or vector with range of years (default c(1950,2019))
#' @param months single month or vector with range of months (default c(1,12))
#' @param storelocal optional path to store downloaded data (default NA, files are not stored), use pisoai_readlocal() to read local data
#' @return data.frame with PISO.AI prediction for precipitation isotopic composition
#' @export
#' @examples
#'
#' # Piso.AI data may be requested as
#' \dontrun{
#'  # Simple request for Piso.AI data of two locations and the full time span of Piso.AI
#'  location<-data.frame(site= c("TEST1","TEST2") ,latitude=c(45.7,46.5),longitude=c(7.6,8.0))
#'  pisoai_data <- pisoai_get_data(location)
#'
#'
#'  # Request for Piso.AI data of two locations and specific timespan
#'  location<-data.frame(site= c("TEST1","TEST2") ,latitude=c(45.7,46.5),longitude=c(7.6,8.0))
#'  pisoai_data <- pisoai_get_data(location, years=c(1990,2000),months=c(4,6))
#'
#'  # Request for Piso.AI data of two locations individual years
#'  location<-data.frame(site= c("TEST1","TEST2"),
#'                       latitude=c(45.7,46.5),longitude=c(7.6,8.0),
#'                       years=c(1979,2010))
#'  pisoai_data <- pisoai_get_data(location)
#'
#'  # Similar request, additionally storing data locally
#'  pisoai_data <- pisoai_get_data(location,storelocal='~/PisoAI_data')
#'  # read back local data
#'  pisoai_data <- pisoai_readlocal('~/PisoAI_data')
#' }

pisoai_get_data<-function(location,years=c(1950,2019),months=c(1,12),storelocal=NA){
  i=0
  if (is.data.frame(location) & length(intersect (c('site','latitude','longitude'), names(location)))==3){
    if (!'elevation' %in% names(location)) location$elevation<-NA
    if ('year' %in% names(location)) {
      pisoai_data<-lapply(1:nrow(location),FUN = function (i) pisoai_request (
        location$site[i],location$latitude[i], location$longitude[i],location$elevation[i]),location$year[i],
        months,storelocal)
    } else {
      pisoai_data<-lapply(1:nrow(location),FUN = function (i) pisoai_request (
        location$site[i],location$latitude[i], location$longitude[i],location$elevation[i],years,months,storelocal))
    }
    pisoai_data<-do.call("rbind",pisoai_data)
    return(pisoai_data)
  }else{
    stop("The first parameter needs to a data.frame with the colums 'site','latitude','longitude'. Colums 'elevation' and 'year' are optional")
  }
}

#'simple Piso.AI API request
#'
#' [Piso.AI](https://isotope.bot.unibas.ch/PisoAI/) is tool for predicting monthly timeseries of oxygen and hydrogen isotope values of precipitation that uses a machine learning model trained on geographic and climate data.
#'
#' @param site site name
#' @param latitude latitude in degrees
#' @param longitude latitude in degrees
#' @param elevation optional custom elevation (otherwise it will be extracted from DEM)
#' @param years single year or vector with range of years (default c(1950,2019))
#' @param months single month or vector with range of months (default c(1,12))
#' @param storelocal optional path to store downloaded data (default NA, files are not stored), use pisoai_readlocal() to read local data
#' @return data.frame with PISO.AI prediction for precipitation isotopic composition
#' @export
#' @examples
#'
#' # single site Piso.AI data may be requested as
#' \dontrun{
#'  pisoai_reqest('BSL',45.7,7.6)
#' }


pisoai_request <- function(site,latitude,longitude,elevation=NA,years=c(1950,2019),months=c(1,12),storelocal=NA){
  print(sprintf("downloading data for %s..."),site)
  request_url<-sprintf('https://isotope.bot.unibas.ch/PisoAI/api?site=%s&lat=%s&lon=%s&from=%s&to=%s',site,latitude,longitude,sprintf('%4g-%02g-01',min(years),min(months)),sprintf('%4g-%02g-31',max(years),max(months)))
  if (!is.na(elevation))  request_url<- sprintf("%s&elevation=%s", request_url,elevation)
  res<-httr::GET(request_url)
  if (res$status_code==200) { # valid data
    data<-utils::read.csv(text=httr::content(res,type='text',encoding = 'UTF-8'),sep=',',skip=2,header=TRUE,stringsAsFactors = FALSE)
    if (!is.na(storelocal)) {
      if (dir.exists(storelocal)){
      fileConn<-file(file.path(storelocal,sprintf("pisoAI_%s.csv",site)))
      write(httr::content(res,type='text',encoding = 'UTF-8'), fileConn)
      close(fileConn)
      } else {
        stop (sprintf("Invalid storelocal parameter: directory %s does not exists",storelocal))
      }
    }
    return(data)
  } else if (res$status_code==400) { # PISO AI API ERROR
      warning(utils::read.csv(text=httr::content(res,type='text',encoding = 'UTF-8'),stringsAsFactors = FALSE)[1,2])
      return(data.frame(Site=site,Date=NA,Latitude=latitude,Longitude=longitude,Elevation=if (!is.na(elevation)) elevation else NA, d18O.Piso.AI =NA,  d2H.Piso.AI = NA ,stringsAsFactors = FALSE))
  } else{ # Other ERROR
    warning(sprintf ('Error: HTTP statuscode %i',httr::status_code(res)))
    return(data.frame(Site=site,Date=NA,Latitude=latitude,Longitude=longitude,Elevation=if (!is.na(elevation)) elevation else NA, d18O.Piso.AI =NA,  d2H.Piso.AI = NA ,stringsAsFactors = FALSE))
  }
}

#' load locally stored Piso_AI data downloaded with through this API
#'
#' @param path local path where the downloads Piso.AI csv files are stored
#' @return data.frame with PISO.AI prediction for precipitation isotopic composition
#' @export
#' @examples
#'
#' # Load local or download if local is not available (no checks are performed testing if local data is complete)
#' \dontrun{
#'  PATH= getwd()
#'  pisoai_data <- pisoai_readlocal(PATH)
#'  if (is.null(pisoai_data)) pisoai_data <- pisoai_get_data(sites, years=c(1990,2018),storelocal=PATH)
#' }
#'
pisoai_readlocal<-function(path){
  pisoai_data<-lapply (list.files(path,'pisoAI_.*.csv',full.names=TRUE), FUN = function (f) utils::read.csv(f,sep=',',skip=2,header=TRUE,stringsAsFactors = FALSE))
  pisoai_data<-do.call("rbind",pisoai_data)
  return(pisoai_data)
}

