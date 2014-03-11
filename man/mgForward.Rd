\name{mgForward}
\alias{mgForward}
\title{
Forward selection of MEM eigenvectors against genetic distance data
}
\description{
This function calls \code{mgRDA} repeatedly in order to identify a reduced
set of all MEM eigenvectors (i.e. spatial patterns).
}
\usage{
mgForward(genD, vectorsMEM, perm = 100, alpha = 0.05)
}
\arguments{
  \item{genD}{
  A symmetrical distance matrix giving the genetic distances among
individual genotypes
}
  \item{vectorsMEM}{
  A matrix giving a set of Moran's eigenvectors 
}
  \item{perm}{
The number of permutations in a randomization test
}
  \item{alpha}{
The 1-alpha level for forward selection
}
}
\details{
A wrapper for \code{\link{mgRDA}} designed for forward selection
}
\value{
A list\cr
\code{$selectedMEM} gives the indices of the input \code{vectorsMEM} that
were selected and can then be used in a call to \code{mgRDA(..., full=TRUE)}
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

## Find Moran's eigenvectors given sampling locations
## by first finding the Euclidean distance matrix
radialEuclid <- dist(radialXY)
radialMEM <- mgMEM(radialEuclid)

## Forward select significant Moran's eigenvectors using RDA
## Positive Moran's eigenvectors (positive spatial autocorrelation) first
radialPositive <- mgForward(radialDM,
    radialMEM$vectorsMEM[ , radialMEM$valuesMEM > 0])
## Negative Moran's eigenvectors (negative spatial autocorrelation) second
radialNegative <- mgForward(radialDM,
    radialMEM$vectorsMEM[ , radialMEM$valuesMEM < 0])


## Summarize the selected MEM eigenvectors
allSelected <- cbind(radialMEM$vectorsMEM[, radialMEM$valuesMEM > 0][
                    , na.omit(radialPositive$selectedMEM)],
                 radialMEM$vectorsMEM[, radialMEM$valuesMEM < 0][
                    , na.omit(radialNegative$selectedMEM)])

## Use the selected Moran's eigenvectors in a final model
radialAnalysis <- mgRDA(radialDM, allSelected, full=TRUE)

}

}
