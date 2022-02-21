\name{mgQuick}
\alias{mgQuick}
\title{
Memgene analysis of genetic distance data (main interface for package)
}
\description{
 Performs multiple--typical--steps in a memgene
 analysis of genetic distance data.  Gracefully handles
 potential errors.  Steps are as follows:

 1.   Find MEM eigenvectors given 
      coordinates (\code{coords})\cr

 2.   Perform separate forward selections of positive and
      negative MEM eigenvectors against genetic distance (\code{genD}),
      to identify a significant
      subset, using parameters \code{forwardPerm} as the number of
      permutations and \code{forwardAlpha} as the alpha level
      for a significant eigenvector.  NOTE: The number of permutations
      \code{forwardPerm} is set at 100 by default to reduce analysis time for
      exploratory analyses.  This number should be increased for final analyses.\cr

 3.   Find the fit of the selected eigenvectors to the
      genetic distance data (using RDA).\cr

 4.   Optionally run a permutation test (\code{finalPerm}) for
      the fit of the selected eigenvectors to the genetic distance
      data.\cr
      
 5.   Produce MEMGENE variables using the fitted values from the RDA
      analysis.  MEMGENE variables are the eigenvectors from a PCA of the fitted
      values.  These are the product of memgene and can be used for visualization
      and subsequent analyses.

 6.   Optionally produce plots of the scores for the
      first \code{n} MEMGENE variables if \code{doPlot = n}.
}
\usage{
mgQuick(genD, coords, longlat = FALSE, truncation = NULL,
    transformation = NULL, forwardPerm = 100, forwardAlpha = 0.05,
    finalPerm = NULL, doPlot = NULL, verbose = TRUE)
}
\arguments{
  \item{genD}{
A symmetrical distance matrix giving the genetic distances among
individual genotypes or populations
}
  \item{coords}{
A two column \code{matrix} or \code{data.frame} of x and y coordinates of sampling locations of individual genotypes
}
  \item{longlat}{
If \code{TRUE} then \code{coords} are longitude and latitude, so
find the geodesic distances among sampling locations using the \code{geosphere} package
}
  \item{truncation}{
\code{NULL} under typical usage.  See \code{\link{mgMEM}} for experimental options.
}
  \item{transformation}{
\code{NULL} under typical usage.  See \code{\link{mgMEM}} for experimental options.  
}
  \item{forwardPerm}{
The number of permutations in the randomization test for the forward selection of MEM eigenvectors. The default \code{forwardPerm=100} is sufficient for exploratory purposes, however this should be increased for final analyses.
}
  \item{forwardAlpha}{
The 1-alpha level for the forward selection process
}
  \item{finalPerm}{
The number of permutations for the final randomization test of the reduced model.  \code{NULL} by
default does not perform a final randomization test.

}
  \item{doPlot}{
Plot \code{doPlot = n} MEMGENE variables
}
  \item{verbose}{
If \code{TRUE} then report progress to the console
}
}
\value{
A list\cr
\code{$P} gives the probability of the null hypothesis for the RDA on the final model\cr
\code{$RSqAdj} is the adjusted R2 for the RDA, understood as the proportion of
all genetic variation that is explicable by spatial pattern (i.e. spatial genetic
signal)\cr
\code{$memgene} contains a matrix with the MEMGENE variables in columns\cr
\code{$memSelected} gives a matrix containing the selected MEM eigenvectors in columns\cr
\code{$whichSelectPos} gives the indices of the selected MEM eigenvectors with positive eigenvalues (i.e. from \code{$mem})\cr
\code{$whichSelectNeg} gives the indices of the selected MEM eigenvectors with negative eigenvalues (i.e. from \code{$mem})\cr
\code{$mem} the output of \code{\link{mgMEM}} given \code{coords}
}
\references{
Galpern, P., Peres-Neto, P., Polfus, J., and Manseau, M.  2014.  MEMGENE: Spatial pattern detection in genetic distance data.  Submitted.
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

## Run the MEMGENE analysis
radialAnalysis <- mgQuick(radialDM, radialXY)

## Extract the scores on the first 3 MEMGENE variables
## for subsequent analysis
radialMEMGENE1 <- radialAnalysis$memgene[, 1]
radialMEMGENE2 <- radialAnalysis$memgene[, 2]
radialMEMGENE3 <- radialAnalysis$memgene[, 3]

## Find the proportion of variation explained by all MEMGENE variables
propVariation <- radialAnalysis$sdev/sum(radialAnalysis$sdev)
}
}
