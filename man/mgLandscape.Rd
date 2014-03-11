\name{mgLandscape}
\alias{mgLandscape}
\alias{print.mgLandscape}
\title{
Landscape genetic analysis using MEMGENE
}
\description{
Use least-cost path distances among sampling locations on a resistance surface, rather than Euclidean distances (as in \code{\link{mgQuick}}), to extract MEM eigenvectors.  The goal is to compare multiple resistance surfaces (i.e. representing alternative hypotheses about landscape resistance) in terms of the proportion of variation in genetic distance they explain.  This is often a goal in landscape genetic analysis.  By default Euclidean distances (i.e. representing a surface with no landscape resistance) are also analyzed unless \code{euclid=FALSE}.\cr

The analysis steps are as follows:\cr

 1.   Find MEM eigenvectors given a distance matrix extracted from the
      coordinates (\code{coords}).  In the case of a resistance surface the
      distances are least-cost paths among sampling locations found using the
      function \code{gdistance::costDistance}.  In the Euclidean case
      Euclidean distances are used.  For all distance matrices a minimum spanning tree
      of the locations is found, followed by truncation of the tree (see \code{\link{mgMEM}})\cr

 2.   Perform separate forward selections of positive and
      negative MEM eigenvectors against genetic distance (\code{genD}),
      to identify a significant
      subset, using parameters \code{forwardPerm} as the number of
      permutations and \code{forwardAlpha} as the alpha level
      for a significant eigenvector. NOTE: The number of permutations
      \code{forwardPerm} is set at 100 by default to reduce analysis time for
      exploratory analyses.  This number should be increased for final analyses.\cr

 3.   Use variation partitioning against the genetic distance matrix
      to find the proportion of variation in genetic
      distance explained by the selected positive and negative MEM eigenvectors (i.e.
      fraction [a] representing spatial genetic variation explained by the resistance surface hypothesis)
      and the matrix of coordinates (i.e. fraction [c] representing spatial genetic variation
      not explained by the resistance hypothesis).  These [a] and [c] fractions
      can be used to inform model selection (see below). 
}
\usage{
mgLandscape(resistance, genD, coords, euclid=TRUE, forwardPerm=100,
forwardAlpha=0.05, finalPerm=1000, verbose=TRUE)
\method{print}{mgLandscape}(x, ...)
}
\arguments{
 \item{resistance}{
 A \code{RasterLayer} produced by the \code{\link{raster}} package in a planar (i.e. not longitude/latitude) projection giving the hypothesized resistance to movement of landscape features (all cells must be either missing as \code{NA} or \code{>=1}) To test multiple resistance hypotheses provide a \code{RasterStack} or \code{RasterBrick} also produced by the \code{\link{raster}} package.
 }
  \item{genD}{
A symmetrical distance matrix giving the genetic distances among
individual genotypes or populations
}
  \item{coords}{
A two column \code{matrix} or \code{data.frame} of x and y coordinates of sampling locations of individual genotypes or populations.  Must be in the same planar projection as the resistance surface.  Geographic coordinates (i.e. longitude/latitude) must be projected before use.
}

  \item{euclid}{
If \code{TRUE} will test the Euclidean distances among sampling locations in addition to the resistance
surface(s) supplied.  Including a Euclidean surface is recommended as a null model.
}
  \item{forwardPerm}{
The number of permutations in the randomization test for the forward selection of MEM eigenvectors.
The default \code{forwardPerm=100} is sufficient for exploratory purposes, however this should be increased for final analyses.
}
  \item{forwardAlpha}{
The 1-alpha level for the forward selection process
}
  \item{finalPerm}{
The number of permutations to test the significance of the [a], [c] and [abc] fractions.
}
  \item{verbose}{
If \code{TRUE} then report progress to the console
}

\item{x}{
  An object of class \code{mgLandscape} produced by the \code{mgLandscape} function
}
\item{...}{
  Additional parameters passed to \code{print}
}
}
\value{
A code{$summary} table giving the results of the variation partitioning.  The following table
provides an interpretation of each of the fractions returned:\cr

Proportion of variation in genetic distance that is... \code{(RsqAdj)}\cr
\code{[abc]   explained by spatial predictors}\cr
\code{[a]     spatial and explained by selected patterns in the model}\cr
\code{[c]     spatial and explained by coordinates not patterns in the model}\cr
\code{[b]     spatial and confounded between the model and coordinates}\cr
\code{[d]     residual not explained by the spatial predictors}\cr

A good model will have a relatively high [a] fraction and relatively low [c] fraction indicating
that the selected patterns in the landscape model have captured a large proportion of the spatial variation in genetic distance.\cr

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
## Compare data generated using the radial data against three landscape models
##
## Prepare two resistance surfaces to test (the true radial, and the false river)
## These are produced as a RasterStack object
if (require(raster)) {
    resistanceMaps <- stack(
           raster(system.file("extdata/radial.asc", package="memgene")),
           raster(system.file("extdata/river.asc", package="memgene")))
} else {
  stop("raster package required for mgLandscape.")
}

## Prepare the radial data for analysis
radialData <- read.csv(system.file("extdata/radial.csv", package="memgene"))
radialGen <- radialData[, -c(1,2)]
radialXY <- radialData[, 1:2]
if (require(adegenet)) {
  radialDM <- codomToPropShared(radialGen)
} else {
  stop("adegenent package required to produce genetic distance matrix in example.")
}

## Analyse the two resistance surfaces and a Euclidean model
## and produce a table comparing the three
## Set permutations at low values for a faster (though less accurate) run
compareThree <- mgLandscape(resistanceMaps, radialDM, radialXY, euclid=TRUE,
   forwardPerm=100, finalPerm=100)
   
print(compareThree)
## Results can vary between runs because selected MEM eigenvectors may vary.
## Setting forwardPerm higher will increase consistency in this regard.
##
## We see that the true radial surface has the highest [a] fraction and
## the lowest [c] fraction indicating that it does well at capturing
## the spatial genetic variation that we expect in this simulated genetic data
}
}
