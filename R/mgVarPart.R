mgVarPart <-
function(mg, X, categorical=NULL, categoryN=NULL, returnX=FALSE) {
    
    if (!require(vegan)) {
        stop("memgene: vegan package must be installed for this function.", call.=FALSE)
    }
    
    if (class(mg) == "mgQuick") {
        Y <- mg$memgene[, which(mg$sdev > 0.01)]
    }
    else {
        if (!is.matrix(mg)) {
            stop("memgene: mg is not a matrix or an object returned by mgQuick()", call.=FALSE)
        }
        Y <- mg
    }
    if (!is.list(X)) {
        stop("memgene:  X must be a list() of length 2, 3, or 4 with matrices as elements", call.=FALSE)
    }
    
    if (!(length(X) %in% c(2,3,4))) {
        stop("memgene:  X must be a list() of length 2, 3, or 4 with matrices as elements", call.=FALSE)
    }
    
    if (any(sapply(1:length(X), function(x) NROW(Y) != NROW(X[[x]])))) {
        stop("memgene:  Y matrix and X matrices must all have the same number of rows", call.=FALSE)
    }
    
    if (!is.null(categorical)) {
        if (length(categorical) != length(X)) {
            stop("memgene:  categorical must be a logical vector of the same length as X", call.=FALSE)
        }
        if (any(sapply(1:length(X), function(x) categorical[x] && (NCOL(X[[x]]) > 1)))) {
            stop("memgene:  categorical variables in X must be a single column", call.=FALSE)
        }
    }
    
    if (!is.null(categoryN) && length(categoryN) != length(X)) {
        stop("memgene:  categoryN must be a vector of the same length as X", call.=FALSE)
    }

    ind <- vector("list", length(X))
    testable <- rep(TRUE, nrow(Y))
    for (iX in 1:length(X)) {
        if (categorical[iX]) {
            if (!is.null(categoryN)) {
                untestable <- names(table(X[[iX]])[table(X[[iX]]) < categoryN[iX]])
                untestable <- (as.character(X[[iX]]) %in% untestable)
                testable <- !untestable & testable
            }
        }
    }
    
    if (sum(testable)==0) {
        stop("memgene:  With the given restrictions on categorical variable sample size (categoryN) there are no observations to test", call.=FALSE)
    }
    for (iX in 1:length(X)) {
        if (categorical[iX]) {
            ind[[iX]] <- model.matrix(~as.factor(X[[iX]]))
            ind[[iX]] <- as.matrix(ind[[iX]][, 2:ncol(ind[[iX]])])
        }
        else {
            ind[[iX]] <- as.matrix(X[[iX]])
        }
        ind[[iX]] <- as.matrix(ind[[iX]][testable, ])
        if (categorical[iX]) ind[[iX]] <- ind[[iX]][, apply(ind[[iX]], 2, sum) > 0]
    }
    
    vp <- list()
    vp$summary <- NA
    Y <- Y[testable, ]
    if (length(X) == 2) {
        X1 <- ind[[1]]
        X2 <- ind[[2]]
        vp$vp <- tryCatch(suppressWarnings(varpart(Y, X1, X2)), error=function(e) return(e))
        if (class(vp$vp)[1] != "simpleError") {
            vp$summary <- data.frame(vp$vp$part$n, t(vp$vp$part$indfrac$Adj.R.square))
            names(vp$summary) <- c("N",  paste("[", letters[1:4], "]", sep=""))
        }
    }
    if (length(X) == 3) {
        X1 <- ind[[1]]
        X2 <- ind[[2]]
        X3 <- ind[[3]]        
        vp$vp <- tryCatch(suppressWarnings(varpart(Y, X1, X2, X3)), error=function(e) return(e))
        if (class(vp$vp)[1] != "simpleError") {
            vp$summary <- data.frame(vp$vp$part$n, t(vp$vp$part$indfrac$Adj.R.square))
            names(vp$summary) <- c("N",  paste("[", letters[1:8], "]", sep=""))
        }
    }
    if (length(X) == 4) {
        X1 <- ind[[1]]
        X2 <- ind[[2]]
        X3 <- ind[[3]]
        X4 <- ind[[4]]
        vp$vp <- tryCatch(suppressWarnings(varpart(Y, X1, X2, X3, X4)), error=function(e) return(e))
        if (class(vp$vp)[1] != "simpleError") {
            vp$summary <- data.frame(vp$vp$part$n, t(vp$vp$part$indfrac$Adj.R.square))
            names(vp$summary) <- c("N",  paste("[", letters[1:16], "]", sep=""))
        }
    }
    
    if (!is.null(categorical) && !is.null(categoryN)) vp$obsUsed <- which(testable)
    
    if (returnX) vp$X <- ind
    
    return(vp)
    
}
