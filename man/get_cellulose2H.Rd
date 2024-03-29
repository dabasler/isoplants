% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plantIso_model.r
\name{get_cellulose2H}
\alias{get_cellulose2H}
\title{Leaf water to cellulose 2HO}
\usage{
get_cellulose2H(dlw, dsw, px, fH, epsA = -171, epsH = 158)
}
\arguments{
\item{px}{fraction of unenriched water in this water pool used for cellulose synthesis}

\item{fH}{Proportion of H exchange during heterotrophic reactions}

\item{epsA}{Biosynthetic fractionation between exchangeable H and medium water during autotrophic metabolism}

\item{d2H_lw}{d2H of leaf water}

\item{d2H_sw}{d2H of source (xylem) water}

\item{epsB}{Biosynthetic fractionation between exchangeable H and medium water during heterotrophic metabolism}
}
\value{
Returns Hydrogen isotope composition of cellulose (permil)
}
\description{
Leaf water to cellulose 2HO
}
