#' @title Get Neighboring Probes
#' @description Get neighboring probes.
#' @return A list.
get.neighboring.probes <- function (X, Y, chr, max.dist, control.arms = TRUE, remove.duplicates = TRUE) {

  xinds <- yinds <- c()
  
  # Use arm information if it is available and not blocked
  if (("arm" %in% names(X$info)) && ("arm" %in% names(Y$info)) && control.arms) {    
    for (arm in c('p', 'q')){      
      # Investigate specified arm
      xchrinds <- which(as.character(X$info$chr) == chr & X$info$arm == arm)
      ychrinds <- which(as.character(Y$info$chr) == chr & Y$info$arm == arm)
      inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist, remove.duplicates)
      xinds <- c(xinds, inds$xinds)
      yinds <- c(yinds, inds$yinds)
    }
  } else {
    # Investigate the whole chromosome
    xchrinds <- which(as.character(X$info$chr) == chr)
    ychrinds <- which(as.character(Y$info$chr) == chr)
    inds <- get.neighs(X, Y, xchrinds, ychrinds, max.dist, remove.duplicates)
    xinds <- c(xinds, inds$xinds)
    yinds <- c(yinds, inds$yinds)
  }

  # return the indices
  list(xinds = xinds, yinds = yinds)
  
}

