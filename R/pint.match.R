pint.match <- function(X, Y, max.dist = 1e7, chrs = NULL, useSegmentedData = FALSE, impute = TRUE, replace.inf = TRUE, remove.duplicates = TRUE){

  # X <- ge; Y <- cn; useSegmentedData = FALSE; impute = TRUE; replace.inf = TRUE
  # X <- geneExp; Y <- geneCopyNum; max.dist <- 1e7; impute = TRUE; replace.inf = TRUE

  # Match genes/probes/clones in X/Y data sets based on location information
  
  if (all(is.na(X$data))) { stop("X data is empty/NA.") }
  if (all(is.na(Y$data))) { stop("Y data is empty/NA.") }
  
  # Impute missing values, replace infinite values etc.
  X <- pint.data(X$data, X$info, impute, replace.inf, remove.duplicates = !useSegmentedData)
  Y <- pint.data(Y$data, Y$info, impute, replace.inf, remove.duplicates = !useSegmentedData)
  
  # Find unique chromosome names present in the data
  # First order chromosomes 1...24, then chromosomes with other names
  if ( is.null(chrs) ) { chrs <- 1:24 }

  message("Matching probes between the data sets..")
  xindices <- yindices <- vector()
  for ( chr in chrs ){
    # Note: chromosome arm information is used in the matching if it is available
    tmp <- get.neighboring.probes(X, Y, chr, max.dist, remove.duplicates = remove.duplicates)
    xindices <- c(xindices, tmp$xinds)
    yindices <- c(yindices, tmp$yinds)
  }

  xdat <- as.matrix(X$data[xindices,], length(xindices))
  ydat <- as.matrix(Y$data[yindices,], length(yindices))

  # Sometimes file reading (with csv at least) leads to situation where the last column is NA.
  # To avoid this and other cases, remove 'NA samples'.
  nainds <- (colMeans(is.na(xdat)) == 1 | colMeans(is.na(ydat)) == 1)
  if (sum(nainds) > 0) {
    xdat <- xdat[, !nainds]
    ydat <- ydat[, !nainds]
    warning(paste("Samples ", colnames(X$data)[nainds], " contained exclusively NA's; removed."))
  }

  newX <- list(data = xdat, info = X$info[xindices,])
  newY <- list(data = ydat, info = Y$info[yindices,])

  # Check if arm info is missing

  list(X = newX, Y = newY)

}

closest <- function(a, vec){which.min(abs(a - vec))}

