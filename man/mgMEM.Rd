\name{mgMEM}
\alias{mgMEM}
\title{
Extraction of MEM eigenvectors given distances among sampling locations
}
\description{
Extract Moran's eigenvectors given a distance matrix
among sampling locations of genetic material.  This matrix could be Euclidean
or otherwise.  If truncation and/or transformation parameters
are provided these operations occur. Truncation implies that distances
that exceed a threshold amount are assigned to 4 * threshold.  Minimum
spanning tree truncation is the recommended default.
Transformation performs an exponential or gaussian transformation of the
distance matrix after truncation.
}
\usage{
mgMEM(locD, truncation = NULL, transformation = NULL)
}
\arguments{
  \item{locD}{
  
A symmetric distance matrix giving the distances (typically Euclidean)
among the sampling locations of genetic material (i.e. of genotyped individuals).


}
  \item{truncation}{

See details (EXPERIMENTAL)
  
  }
  \item{transformation}{
Can be character "exponential" or "gaussian" or NULL for no transformation (EXPERIMENTAL)
}
}
\details{
If sampling locations are in longitude/latitude and are far apart, be sure to
supply the geodesic distance as \code{locD}.  (Note that \code{\link{mgQuick}} implements
geodesic distances using the \code{longlat=TRUE} parameter when provided with sampling coordinates)

\code{truncation}\cr

1.  Can be numeric from 0 to 1 specifying the
proportion of the maximum distance in locD to truncate
following this a minimum spanning tree is used to further
truncate as in PCNM (aka dbMEM or classical MEM) \cr

2.  Can be NULL (default) indicating only the minimum spanning tree (MST) truncation where
links that exceed the longest link in the MST (dMST) are replaced with 4 * dMST\cr

3.  Can be FALSE indicating that nothing is done to the distance matrix
which is only suitable when locD is non-euclidean (i.e. will have
negative eigenvectors
}
\value{
A list\cr
\code{$valuesMEM} gives the eigenvalues all Moran's eigenvectors\cr
\code{$vectorsMEM} gives the Moran's eigenvectors in columns
}
\references{
Legendre, P., and Legendre L.  2012. Numerical Ecology, 3rd. ed.  Elsevier, Amsterdam.
}
\author{
Pedro Peres-Neto (peres-neto.pedro@uqam.ca)\cr
Paul Galpern (pgalpern@ucalgary.ca)
}

\examples{

## Prepare the radial data for analysis
radialData <- read.csv(system.file("extdata/radial.csv", package="memgene"))
radialXY <- radialData[, 1:2]

## Find Moran's eigenvectors given sampling locations
## by first finding the Euclidean distance matrix
radialEuclid <- dist(radialXY)
radialMEM <- mgMEM(radialEuclid)

}
