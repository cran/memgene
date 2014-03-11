\name{mgVarPart}
\alias{mgVarPart}
\title{
Variation partitioning of the genetic distance matrix 
}
\description{
This function performs a variation partitioning of the genetic distance matrix
using the supplied MEM eigenvectors and spatial coordinates.  Randomization tests
are conducted to determine the significance of the [a] fraction representing
the MEM eigenvectors, the [c] fraction representing the spatial coordinates and
the [abc] fraction representing the spatial genetic variation.  It is called
by \code{\link{mgLandscape}}.
}
\usage{
mgVarPart(genD, vectorsMEM, coords, perm=1000)
}
\arguments{
   \item{genD}{
A symmetrical distance matrix giving the genetic distances among
individual genotypes or populations
}
  
\item{vectorsMEM}{
A matrix giving a set of any number of MEM eigenvectors 
}
  \item{coords}{
A two column \code{matrix} or \code{data.frame} of x and y coordinates of sampling locations of individual genotypes or populations. 
}
\item{perm}{
The number of permutations to use when testing the significance of the [a], [c] and [abc] fractions.
}
}
\details{
See \code{\link{mgLandscape}} for explanation of the fractions.
}

\author{
Pedro Peres-Neto (peres-neto.pedro@uqam.ca)\cr
Paul Galpern (pgalpern@ucalgary.ca)
}

\examples{

\dontrun{
## Prepare the radial data for analysis
radialData <- read.csv(system.file("extdata/radial.csv", package="memgene"))
radialGen <- radialData[, -c(1,2)]
radialXY <- radialData[, 1:2]

if (require(adegenet)) {
  radialDM <- codomToPropShared(radialGen)
} else {
  stop("adegenent package required to produce genetic distance matrix in example.")
}


## Find MEM eigenvectors given sampling locations
## by first finding the Euclidean distance matrix
radialEuclid <- dist(radialXY)
radialMEM <- mgMEM(radialEuclid)

## Forward select significant MEM eigenvectors using RDA
## Positive MEM eigenvectors (positive spatial autocorrelation) first
radialPositive <- mgForward(radialDM,
    radialMEM$vectorsMEM[ , radialMEM$valuesMEM > 0])
## Negative MEM eigenvectors (negative spatial autocorrelation) second
radialNegative <- mgForward(radialDM,
    radialMEM$vectorsMEM[ , radialMEM$valuesMEM < 0])


## Summarize the selected MEM eigenvectors
allSelected <- cbind(radialMEM$vectorsMEM[, radialMEM$valuesMEM > 0][
                    , na.omit(radialPositive$selectedMEM)],
                 radialMEM$vectorsMEM[, radialMEM$valuesMEM < 0][
                    , na.omit(radialNegative$selectedMEM)])

## Use the selected MEM eigenvectors and coordinates in
## variation partitioning
radialVarPart <- mgVarPart(radialDM, allSelected, radialXY)
}
}