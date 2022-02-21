\name{mgRDA}
\alias{mgRDA}

\title{
Extraction of MEMGENE variables using redundancy analysis (RDA)
}
\description{
Performs a redundancy analysis (RDA) given MEM eigenvectors and
a genetic distance matrix.  Optionally performs a permutation test
for the RDA. Returns the MEMGENE variables, which are the product of a PCA
conducted on the fitted values of this RDA. 
}
\usage{
mgRDA(genD, vectorsMEM, perm = NULL, full = TRUE)
}

\arguments{
  \item{genD}{
A symmetrical distance matrix giving the genetic distances among
individual genotypes or populations
}
  \item{vectorsMEM}{
A matrix giving a set of any number of MEM eigenvectors 
}
  \item{perm}{
The number of permutations in a randomization test
}
  \item{full}{
If \code{TRUE} returns the MEMGENE variables.  \code{FALSE} is used
primarily by \code{\link{mgForward}} which calls this function.
}
}
\details{
Any type of genetic distance matrix \code{genD} giving pairwise
distances among individual genotypes could be used.  Population genetic distances (e.g. pairwise
Fst among populations) could also be used in principle, in which case the sampling centroids
of populations should be used to develop the MEM eigenvectors.
}
\value{
A list:\cr
\code{$RsqAdj} is the adjusted R2 for the RDA, understood as the proportion of
all genetic variation that is explicable by spatial pattern (i.e. spatial genetic
signal)\cr
\code{$memgene} gives the MEMGENE variables ordered according to the eigenvalues
which are given in \code{$sdev}
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

## Use the selected MEM eigenvectors in a final model
radialAnalysis <- mgRDA(radialDM, allSelected, full=TRUE)
}
}

