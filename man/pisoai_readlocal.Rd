% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pisoAI_APIclient.r
\name{pisoai_readlocal}
\alias{pisoai_readlocal}
\title{load locally stored Piso_AI data downloaded with through this API}
\usage{
pisoai_readlocal(path)
}
\arguments{
\item{path}{local path where the downloads Piso.AI csv files are stored}
}
\value{
data.frame with PISO.AI prediction for precipitation isotopic composition
}
\description{
load locally stored Piso_AI data downloaded with through this API
}
\examples{

# Load local or download if local is not available (no checks are performed testing if local data is complete)
\dontrun{
 PATH= getwd()
 pisoai_data <- pisoai_readlocal(PATH)
 if (is.null(pisoai_data)) pisoai_data <- pisoai_get_data(sites, years=c(1990,2018),storelocal=PATH)
}

}
