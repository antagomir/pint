#' @title Augment
#' @description Augment region.
#' @return A list.
augment.region <- function (region.genes, gene.info) {

  # Check gene.info for genes that are, according to the annotations
  # in gene.info, within the same region than the listed
  # region.genes. Useful for instance when some genes have been left
  # out from modeling due to required one-to-one match between the two
  # data sets.

  # Pick annotations for this region
  #reg.info <- subset(gene.info, geneName %in% region.genes)
  reg.info <- gene.info[gene.info$geneName %in% region.genes, ]

  chrs <- reg.info$chr
  if ( length(unique(chrs)) > 1 ) { stop("Multiple chromosomes listed for the region!") }

  arms <- reg.info$arm
  if ( length(unique(arms)) > 1 ) { stop("Multiple arms listed for the region!") }
  
  reg.start <- as.numeric(min(reg.info$loc))
  reg.end <- as.numeric(max(reg.info$loc))

   # List all genes that are within the same region
  inds <- (gene.info$chr %in% chrs & gene.info$arm %in% arms & gene.info$loc >= reg.start & gene.info$loc <= reg.end )
  rg <-gene.info$geneName[inds]

  # Ensure that also original region.genes are listed
  # even if they are not in gene.info
  union(rg, region.genes)
  
}
