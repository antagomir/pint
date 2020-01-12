#' @title Join Top
#' @description Join top regions.
#' @return A list.
join.top.regions <- function (model, feature.info, quantile.th = 0.95, augment = FALSE) {

  # Pick all top regions (based on the given threshold) and join overlapping windows
  # to get easily interpretable summaries

  # model: pint model object
  # feature.info: for instance ge$info; location info for genes
  # ntop: number top windows to check

  feature.info <- order.feature.info(feature.info)
  df <- as.data.frame(model)
  th <- quantile(df[["dependencyScore"]], quantile.th)
  #topg <- subset(df, dependencyScore > th)$geneName
  topg <- df$geneName[df$dependencyScore > th]

  # Get window around each top gene
  mygenes <- lapply(topg, function(gn) { rownames(findModel(model, gn)@W$X) })

  # Add other known genes in each region 
  if (augment && !is.null(feature.info)) {
    mygenes.aug <- lapply(mygenes, function (reg) { augment.region(reg, feature.info) })
  }

  # List all genes residing on the top regions
  topw <- unique(unname(unlist(mygenes.aug)))

  # Match to feature table, which is ordered(!) by chromosomal locations
  ord <- match(topw, rownames(feature.info))
  o <- order(ord) 
  ord <- ord[o]

  # Mark break points
  do <- c(1, diff(ord) == 1) 

  regs <- list()
  i <- 1
  k <- 0
  reg <- c()
  while (i <= length(do)) {
  
    if (do[[i]] == 1) {
      reg <- c(reg, ord[[i]])
    } else {
      k <- k + 1
      regs[[k]] <- reg
      reg <- c(ord[[i]])
    }
    i <- i + 1
  }

  regs <- lapply(regs, function (x) {rownames(feature.info)[x]})
  names(regs) <- as.character(1:length(regs))

  regs

}
