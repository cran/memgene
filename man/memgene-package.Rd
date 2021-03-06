\name{memgene-package}
\alias{memgene-package}
\alias{memgene}
\docType{package}
\title{
Spatial pattern detection in genetic distance data using Moran's Eigenvector Maps
}
\description{
Memgene can detect relatively weak spatial genetic patterns by using Moran's Eigenvector Maps (MEM) to extract only the spatial component of genetic variation.  Memgene has applications in landscape genetics where the movement and dispersal of organisms are studied using neutral genetic variation.
}
\details{
\tabular{ll}{
Package: \tab memgene\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-07\cr
License: \tab GPL (>=2)\cr
}

}
\author{
Paul Galpern (pgalpern@ucalgary.ca)\cr
Pedro Peres-Neto (peres-neto.pedro@uqam.ca)

Maintainer: Paul Galpern (pgalpern@ucalgary.ca)
}
\references{
Galpern, P., Peres-Neto, P., Polfus, J., and Manseau, M.  2014.  MEMGENE: Spatial pattern detection in genetic distance data.  Submitted.
}
\keyword{ package }

\examples{

## The basic interface to MEMGENE is mgQuick()
?mgQuick

## For landscape genetic analysis with MEMGENE see mgLandscape()
?mgLandscape
}
