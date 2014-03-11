### R code from vignette source 'memgeneTutorial.Rnw'

###################################################
### code chunk number 1: memgeneTutorial.Rnw:43-45
###################################################
library(memgene)
options(width=65)


###################################################
### code chunk number 2: memgeneTutorial.Rnw:69-72
###################################################
library(raster)
radialRas <- raster(system.file("extdata/radial.asc", package="memgene"))
plot(radialRas, legend=FALSE)


###################################################
### code chunk number 3: memgeneTutorial.Rnw:91-93
###################################################
## Load the radial genetic data
radialData <- read.csv(system.file("extdata/radial.csv", package="memgene"))


###################################################
### code chunk number 4: memgeneTutorial.Rnw:95-98
###################################################
## Create objects for positional information and genotypes
radialXY <- radialData[ ,1:2]
radialGen <- radialData[, 3:ncol(radialData)]


###################################################
### code chunk number 5: memgeneTutorial.Rnw:100-103
###################################################
## Produce a proportion of shared alleles genetic distance matrix
## using the convenience wrapper function provided with the package
radialDM <- codomToPropShared(radialGen)


###################################################
### code chunk number 6: memgeneTutorial.Rnw:113-116
###################################################
## Run the MEMGENE analysis
## May take several minutes
if (!exists("radialAnalysis")) radialAnalysis <- mgQuick(radialDM, radialXY)


###################################################
### code chunk number 7: memgeneTutorial.Rnw:125-128
###################################################
## Visualize the first two MEMGENE variables
## by providing only the first two columns of the $memgene matrix
mgMap(radialXY, radialAnalysis$memgene[, 1:2])


###################################################
### code chunk number 8: memgeneTutorial.Rnw:136-140
###################################################
library(raster)
radialRas <- raster(system.file("extdata/radial.asc", package="memgene"))
plot(radialRas, legend=FALSE)
mgMap(radialXY, radialAnalysis$memgene[, 1], add.plot=TRUE, legend=TRUE)


###################################################
### code chunk number 9: memgeneTutorial.Rnw:143-147
###################################################
library(raster)
radialRas <- raster(system.file("extdata/radial.asc", package="memgene"))
plot(radialRas, legend=FALSE)
mgMap(radialXY, radialAnalysis$memgene[, 1], add.plot=TRUE, legend=TRUE)


###################################################
### code chunk number 10: memgeneTutorial.Rnw:176-178
###################################################
## Load the caribou genetic data
caribouData <- read.csv(system.file("extdata/caribou.csv", package="memgene"))


###################################################
### code chunk number 11: memgeneTutorial.Rnw:180-183
###################################################
## Create objects for positional information and genotypes
caribouXY <- caribouData[ ,1:2]
caribouGen <- caribouData[, 3:ncol(caribouData)]


###################################################
### code chunk number 12: memgeneTutorial.Rnw:185-188
###################################################
## Produce a proportion of shared alleles genetic distance matrix
## using the convenience wrapper function provided with the package
caribouDM <- codomToPropShared(caribouGen)


###################################################
### code chunk number 13: memgeneTutorial.Rnw:194-197
###################################################
## Run the MEMGENE analysis
## May take several minutes
if (!exists("caribouAnalysis")) caribouAnalysis <- mgQuick(caribouDM, caribouXY)


###################################################
### code chunk number 14: memgeneTutorial.Rnw:209-212
###################################################
plot(caribouXY, type="n", xlab="", ylab="", axes=FALSE)
mgMap(caribouXY, caribouAnalysis$memgene[, 1], add.plot=TRUE, legend=TRUE)
box()


###################################################
### code chunk number 15: memgeneTutorial.Rnw:225-226
###################################################
caribouAnalysis$RsqAdj


###################################################
### code chunk number 16: memgeneTutorial.Rnw:233-235
###################################################
## Find the proportional variation explained by each MEMGENE variable
caribouMEMGENEProp <- caribouAnalysis$sdev/sum(caribouAnalysis$sdev)


###################################################
### code chunk number 17: memgeneTutorial.Rnw:238-240
###################################################
## Neatly print proportions for the first three MEMGENE variables
format(signif(caribouMEMGENEProp, 3)[1:3], scientific=FALSE)


