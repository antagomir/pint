#' @title Order Features
#' @description Order feature info.
#' @return A list
order.feature.info <- function (feature.info) {

  feature.info$chr <- as.character(feature.info$chr)
  feature.info$arm <- as.character(feature.info$arm)
  
  # Remove genes with no location information from the annotations
  nainds <- is.na(feature.info$chr) | (is.na(feature.info$chr) | is.na(feature.info$arm))
  if (sum(nainds) > 0) {
    feature.info <- feature.info[!nainds, ]
  }
  
  # Order by chromosomal locations
  if ("X" %in% feature.info$chr) {
    feature.info$chr[feature.info$chr == "X"] <- "23"
  }
  if ("Y" %in% feature.info$chr) {
    feature.info$chr[feature.info$chr == "Y"] <- "24"
  }  

  feature.info$chr <- as.numeric(feature.info$chr)
  
  feature.info.ordered <- NULL
  chrs <- sort(unique(feature.info$chr))
  arms <- sort(unique(feature.info$arm))
  for (chr in chrs) {
    if (!is.null(arms)) {
      for (arm in arms) {
        #arm.info <- subset(feature.info, chr == chr & arm == arm)
        inds <- (feature.info$chr == chr & feature.info$arm == arm)
        arm.info <- feature.info[inds,,drop=FALSE]
        if (nrow(arm.info) > 0) {
          o <- order(arm.info$loc)
          feature.info.ordered <- rbind(feature.info.ordered, arm.info[o, ])
        }
      }
    } else {
      arm.info <- feature.info[feature.info$chr == chr,,drop=FALSE]      
      o <- order(arm.info$loc)
      feature.info.ordered <- rbind(feature.info.ordered, arm.info[o, ])
    }
  }
  feature.info.ordered
}

