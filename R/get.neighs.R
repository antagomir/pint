#' @title Get Neighbors
#' @description Get neighbors.
#' @return A list.
get.neighs <- function (X, Y, xchrinds, ychrinds, max.dist, remove.duplicates = TRUE) {

  xinds <- yinds <- NULL

  if ( length(xchrinds) > 0 && length(ychrinds) > 0 ){
    
    #Find indices of closest probe from Y for each from X
    xi <- 1:length(xchrinds)
    yi <- sapply(as.numeric(as.character(X$info$loc[xchrinds])), closest, vec = as.numeric(as.character(Y$info$loc[ychrinds])))

    # Remove duplicates
    if (remove.duplicates) {
      keep <- !duplicated(yi)
      xi <- xi[keep]
      yi <- yi[keep]
    }
    
    # Corresponding indices between X and Y
    xinds <- xchrinds[xi]
    yinds <- ychrinds[yi]

    # delete indices which are further from each other than threshold
    near <- (abs(X$info$loc[xinds] - Y$info$loc[yinds]) < max.dist)
    xinds <- xinds[near]
    yinds <- yinds[near]
  
    # calculate mean location for each pair for ordering of the pairs
    xy.loc <- (X$info$loc[xinds] + Y$info$loc[yinds])/2

    # ensure probes are ordered by location
    ord <- order(xy.loc)
    xinds <- xinds[ord]
    yinds <- yinds[ord]

  }

  list(xinds = xinds, yinds = yinds)
  
}
