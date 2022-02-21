\name{mgMap}
\alias{mgMap}
\title{
Visualization of MEMGENE variables
}
\description{
A high-level plotting interface for the bubble plot visualization of MEMGENE variables.\cr
If there are exactly two columns in \code{memgene} and therefore two MEMGENE variables to be plotted, then a single plotting window is created with the two plots side by side.  Otherwise each MEMGENE variable is plotted in its own window unless.\cr
}
\usage{
mgMap(coords, memgene, wid = NULL, hei = NULL, dev.open = FALSE,
    add.plot = FALSE, legend = FALSE, ...)
}
\arguments{
  \item{coords}{
  A two column \code{matrix} or \code{data.frame} of x and y coordinates of sampling locations of individual genotypes
}
  \item{memgene}{
A matrix giving as columns the MEMGENE variables to be plotted (e.g. can be subsetted from the \code{$memgene} element produced by \code{\link{mgQuick}})
}
  \item{wid}{
The width of the plotting device to be created.  If \code{NULL} the decision is made by the function.
}
  \item{hei}{
The width of the plotting device to be created.  If \code{NULL} the decision is made by the function.  
}
  \item{dev.open}{
If \code{TRUE} do not open a new plotting device.  
}
  \item{add.plot}{
If \code{TRUE} superimposing bubbles on an existing plot or map.  
}
  \item{legend}{
If \code{TRUE} add a legend to the plot
}
  \item{\dots}{
Additional parameters passed to the \code{sr.value} function modified from Borcard et al (2011).
}
}
\details{
This function embeds slightly modified versions of \code{sr.value}, \code{scatterutil.legend.bw.circle}, and \code{scatterutil.legend.circle.grey} distributed with Borcard et al. (2012) which are themselves modified from similar functions distributed with the \code{\link{ade4}} package under a GPL-2 license. 
}
\value{
Side effect.  A plot is produced.
}
\references{
Borcard, D., Gillet, F., and Legendre. P.  2011.  Numerical Ecology with R.  Springer, New York.
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

## Visualize the first two MEMGENE variables side-by-side
mgMap(radialXY, radialAnalysis$memgene[, 1:2])

## Visualize the first MEMGENE variable superimposed over a raster map
## with the same coordinate system, AND include a legend
if (require(raster)) {
    resistanceMap <- raster(system.file("extdata/radial.asc", package="memgene"))
    plot(resistanceMap, legend=FALSE)
    mgMap(radialXY, radialAnalysis$memgene[, 1], add.plot=TRUE, legend=TRUE)
} else {
    mgMap(radialXY, radialAnalysis$memgene[, 1], legend=TRUE)
}

}
}
