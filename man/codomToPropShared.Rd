\name{codomToPropShared}
\alias{codomToPropShared}
\title{
Find proportion of shared alleles genetic distances from a codominant alleles matrix
}
\description{
Given a matrix with two adjacent columns for each locus (e.g. LOCUS1a, LOCUS1b, LOCUS2a, LOCUS2b, ...)
containing codominant alelles, where individual genotypes are in rows, find the proportion of shared alleles
(Bowcock et al., 1994) among individuals using functions in the adegenet package.\cr

 This is a convenience function that wraps adegenet routines.\cr
 
 Note that any type of genetic distance matrix can be used in MEMGENE, and the proportion of shared alleles metric is not a requirement. 
}
\usage{
codomToPropShared(alleles, missingData = c(-98,-99))
}
\arguments{
  \item{alleles}{
A matrix with two adjacent columns for each locus containing codominant alleles

}
  \item{missingData}{
 A vector of any length giving the values in the alleles matrix representing missing data (\code{NA} also represents missing)
}
}
\details{
First prepares the alleles matrix into a format that can be converted using functions in the adegenet package to the \code{genind} format.  \code{propShared} is then run on this object.
}
\value{
Returns a genetic distance matrix using the proportion of shared alleles metric (Bowcock et al., 1994)
}
\references{
Bowcock AM, Ruizlinares A, Tomfohrde J, et al. 1994 High resolution of human evolutionary trees with polymorphic microsatellites. Nature, 368, 455-457.
}
\author{
Pedro Peres-Neto (peres-neto.pedro@uqam.ca)\cr
Paul Galpern (pgalpern@ucalgary.ca)
}

\examples{

radialData <- read.csv(system.file("extdata/radial.csv", package="memgene"))
radialGen <- radialData[, -c(1,2)]
radialDM <- codomToPropShared(radialGen)

}
