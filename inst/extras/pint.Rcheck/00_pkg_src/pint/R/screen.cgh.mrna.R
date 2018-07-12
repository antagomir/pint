#' Fit Dependency Model
#' 
#' Fits dependency models for whole chromosomal arm, chromosome or genome
#' depending on arguments with chosen window size between two data sets.
#' 
#' Function \code{screen.cgh.mrna} assumes that data is already paired. This
#' can be done with \code{\link{pint.match}}. It takes sliding gene windows
#' with \code{\link{fixed.window}} and fits dependency models to each window
#' with \code{\link{fit.dependency.model}} function. If the window exceeds
#' start or end location (last probe) in the chromosome in the
#' \code{\link{fixed.window}} function, the last window containing the given
#' probe and not exceeding the chromosomal boundaries is used. In practice,
#' this means that dependency score for the last n/2 probes in each end of the
#' chromosome (arm) will be calculated with an identical window, which gives
#' identical scores for these end position probes. This is necessary since the
#' window size has to be fixed to allow direct comparability of the dependency
#' scores across chromosomal windows. Function \code{screen.cgh.mir} calculates
#' dependencies around a chromosomal window in each sample in \code{X}; only
#' one sample from \code{X} will be used. Data sets do not have to be of the
#' same size and\code{X} can be considerably smaller. This is used with e.g.
#' miRNA data. If method name is specified, this overrides the corresponding
#' model parameters, corresponding to the modeling assumptions of the specified
#' model. Otherwise method for dependency models is determined by parameters.
#' Dependency scores are plotted with \link{dependency score plotting}.
#' 
#' @param X Gene/mirna expression data
#' @param Y Copy number data. A list with the following items: data: Data in
#' mat rix form. Genes are in rows and samples in columns. e.g. gene copy
#' number; info Data frame which contains following information about genes in
#' data matrix; chr Number indicating the chrosome for the gene: (1 to 24).
#' Characters 'X' or 'Y' can be used;
#' @param windowSize Number of nearest genes to be included in the chromosomal
#' window of the model, and therefore the scale of the investigated chromosomal
#' region. If not specified, using the default ratio of 1/3 between features
#' and samples or 15 if the ratio would be greater than 15
#' @param chromosome Specify the chromosome for model fitting. If missing,
#' whole genome is screened;
#' @param arm Character indicating the chromosomal arm for the gene 'p' or 'q';
#' loc: Location of the gene in base pairs; pint.data can be used to create
#' data sets in this format; Specify chromosomal arm for model fitting. If
#' missing, both arms are modeled.
#' @param method Dependency screening can utilize any of the functions from the
#' package dmt at CRAN. Particular options include pSimCCA probabilistic
#' similarity constrained CCA (Lahti et al. 2009). This is the default method.
#' pCCA probabilistic canonical correlation analysis (Bach & Jordan 2005) pPCA
#' probabilistic principal component analysis (Tipping & Bishop 1999) pFA
#' probabilistic factor analysis (Rubin & Thayer 1982) TPriorpSimCCA
#' probabilistic similarity constrained canonical correlation analysis with
#' possibility to tune T prior (Lahti et al. 2009). If anything else, the model
#' is specified by the given parameters.
#' @param params List of parameters for the dependency model.  sigmas Variance
#' parameter for the matrix normal prior distribution of the transformation
#' matrix T. This describes the deviation of T from H H Mean parameter for the
#' matrix normal prior distribution prior of transformation matrix T zDimension
#' Dimensionality of the latent variable mySeed Random seed.  covLimit
#' Convergence limit. Default depends on the selected method: 1e-3 for pSimCCA
#' with full marginal covariances and 1e-6 for pSimCCA in other cases.
#' @param max.dist Maximum allowed distance between probes. Used in automated
#' matching of the probes between the two data sets based on chromosomal
#' location information.
#' @param outputType Specifies the output type of the function. possible values
#' are "models" and "data.frame"
#' @param useSegmentedData Logical. Determines the useage of the method for
#' segmented data
#' @param match.probes To be used with segmented data, or nonmatched probes in
#' general. Using nonmatched features (probes) between the data sets.
#' Development feature, to be documented later.
#' @param regularized Regularization by nonnegativity constraints on the
#' projections. Development feature, to be documented later.
#' @return Return value is defined by the the function argument
#' \code{outputType}. With the argument \code{outputType = "models"}, the
#' return value depends on the other arguments; returns a
#' \linkS4class{ChromosomeModels} which contains all the models for
#' dependencies in chromosome or a \linkS4class{GenomeModels} which contains
#' all the models for dependencies in genome. With the argument
#' \code{outputType = "data.frame"}, the function returns a data frame with
#' eachs row representing a dependency model for one gene.  The columns are:
#' \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}.
#' @author Olli-Pekka Huovilainen and Leo Lahti \email{leo.lahti#' @@iki.fi}
#' @seealso To fit a dependency model: \code{\link{fit.dependency.model}};
#' \linkS4class{ChromosomeModels} holds dependency models for chromosome,
#' \linkS4class{GenomeModels} holds dependency models for genome. For plotting,
#' see: \link{dependency score plotting}
#' @references Dependency Detection with Similarity Constraints, Lahti et al.,
#' 2009 Proc. MLSP'09 IEEE International Workshop on Machine Learning for
#' Signal Processing, See
#' \url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf} A
#' Probabilistic Interpretation of Canonical Correlation Analysis, Bach Francis
#' R. and Jordan Michael I. 2005 Technical Report 688. Department of
#' Statistics, University of California, Berkley.
#' \url{http://www.di.ens.fr/~fbach/probacca.pdf} Probabilistic Principal
#' Component Analysis, Tipping Michael E. and Bishop Christopher M. 1999.
#' Journal of the Royal Statistical Society, Series B,
#' 61, Part 3, pp. 611--622.
#' \url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}
#' EM Algorithms for ML Factoral Analysis, Rubin D. and Thayer D. 1982.
#' Psychometrika, vol. 47, no. 1.
#' @keywords iteration
#' @examples
#' 
#' \dontrun{
#'   data(chromosome17)
#'   ## pSimCCA model on chromosome 17
#'   models17pSimCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
#'                                      windowSize = 10, chr = 17)
#'                                     
#'   plot(models17pSimCCA)
#' 
#'   ## pCCA model on chromosome 17p with 3-dimensional latent variable z
#'   models17ppCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
#'                                    windowSize = 10,
#'                                    chromosome = 17, arm = 'p',method="pCCA", 
#' 	      	 	           params = list(zDimension = 3))
#'   plot(models17ppCCA)
#' }
#' 
#' @export screen.cgh.mrna
screen.cgh.mrna <- function(X, Y, windowSize = NULL,
                            chromosome,
                            arm,
                            method = "pSimCCA",
                            params = list(),
                            max.dist = 1e7, 
                            outputType = "models",
                            useSegmentedData = TRUE,
                            match.probes = TRUE,
                            regularized = FALSE)
{

  # FIXME: quick hack - later modify genomeModels class
  if (is.null(X$info$arm) && is.null(Y$info$arm)) {
    warning("Arm information missing, artificially adding p arm for all probes.")
    X$info$arm <- rep("p", nrow(X$info))
    Y$info$arm <- rep("p", nrow(Y$info))
  }

  if (is.null(X$info$arm) && !is.null(Y$info$arm) && nrow(X$info) == nrow(Y$info)) {
    warning("Arm information missing from X data, borrowing the arm info from Y data.")
    X$info$arm <- Y$info$arm
  }

  if (!is.null(X$info$arm) && is.null(Y$info$arm) && nrow(X$info) == nrow(Y$info)) {
    warning("Arm information missing from Y data, borrowing the arm info from X data.")
    Y$info$arm <- X$info$arm
  }    

  if (is.null(windowSize)) {
    windowSize <- min(floor(ncol(X$data)/3),15)
    cat(paste("Chromosomal window (windowSize) not specified. Using
default ratio of 1/3 between features and samples (with max window
size 15 probes). Using window size ", windowSize,"\n"))
  }
      
  # FIXME: move all preprocessing stuff to dedicated preprocessing functions
  # pint.data and pint.match
  
  # Check ordering of samples
  if (any(colnames(X$data) != colnames(Y$data))) {
    warning("Samples not in the same order in the two data sets. Using samples that are found in both data sets and reordering the samples.\n")

    commons <- intersect(colnames(X$data), colnames(Y$data))
    if (length(commons) > 1) {
      X$data <- X$data[, commons]
      Y$data <- Y$data[, commons]
    } else {
      stop("Not enough common samples found. Check that the corresponding samples in the two data sets have identical names.\n")
    }
  }

  ## Checks that segmented data is not used when not implicitely
  ## indicated by argument
  ## FIXME: this is independent of 'segmented' option, join these later
  if (!useSegmentedData){
    if (test.segmented(X$data) || test.segmented(Y$data)){
      warning("Segmented data found while method for non-segmented data is selected.\n", immediate. = TRUE)    
    }
  }

  if (match.probes) {   # FIXME: change name to match.probes or something 
    # Match probes
    tmp <- pint.match(X, Y, max.dist, useSegmentedData = useSegmentedData)
    X <- tmp$X
    Y <- tmp$Y
    match.probes <- FALSE # now the probes are matched.
  } else {
    # If user claims that no matching is needed
    # this will implicate that matching has already been
    # performed. Verify: the number of probes should be
    # identical in ge and cn data sets
    if (!nrow(X$data) == nrow(Y$data)) {
      stop("If match.probes == FALSE then the number of probes in ge and cn data sets (X$data, Y$data) need to match!")      
    }  
  }
  
  # Remove probes where observations are not available in either data set
  # TODO

  ############################################################################

  # Set priors

  # Currently imlement only pSimCCA for the screening
  if (method == "pSimCCA") {
        # similarity prior: Wx = Wy identically
	priors <- list(sigma.w = 0) 
  } else if (method == "nonmatched") {
    # no similarity prior for Wx ~ Wy, completely uncoupled 
    priors <- list(sigma.w = Inf)
  } else {
    message("Method not among pSimCCA, nonmatched. Using empty prior.")
    priors <- NULL
  }
	      
  if (regularized) {priors$W <- 1e-3} # W>=0 prior
		  
  ###########################################################################

  if (method == "pSimCCA") {
    if (is.null(params$H)) {
      params$H <- diag(1, windowSize, windowSize)
    }
  } else if (method == "pPCA" || method == "pCCA" || method == "pFA") {
    params$H <- NA
  } else {
    if (is.null(params$H))
      params$H <- diag(1, windowSize, windowSize)
  }

  if (is.null(params$sigmas))
    params$sigmas <- 0
  if (method == "pPCA") {
    params$marginalCovariances <- "identical isotropic"
  } else if (method == "pFA") {
    params$marginalCovariances <- "diagonal"
  } else if (method == "pCCA") {
    if (is.null(params$marginalCovariances)) {
      params$marginalCovariances <- "full"
    }
  } else {
    if (is.null(params$marginalCovariances)) {
      if (params$sigmas == 0) {
        params$marginalCovariances <- "full"
      } else {
        params$marginalCovariances <- "isotropic"
      }
    }
  }
  if (is.null(params$zDimension))
    params$zDimension <- 1
  if (!is.null(params$H) && any(is.na(params$H))) {
    if (params$marginalCovariances == "full")
      method = "pCCA"
    if (params$marginalCovariances == "isotropic")
      method = "pCCA"
    if (params$marginalCovariances == "diagonal")
      method = "pFA"
    if (params$marginalCovariances == "identical isotropic")
      method = "pPCA"
  } else {
    method = "pSimCCA"
  }

  # Convert X and Y chromosomes to numbers
  if (!missing(chromosome)){
    if (chromosome == 23) chromosome <- "X"
    if (chromosome == 24) chromosome <- "Y"
  }

  # give warning if arm param given and no arm data is available
  if (!missing(arm) && any(X$info$arm == "")){
    warning("No arm information in the data. Calculating thw whole chromosome")
    arm <- NULL
  }

  if (missing(chromosome)) {
    models <- calculate.genome(X, Y, windowSize, method, params, match.probes = match.probes, priors = priors, regularized = regularized)
  } else if (missing(arm) || is.null(arm)) {
    models <- calculate.chr(X, Y, windowSize, chromosome, method, params, match.probes = match.probes, priors = priors, regularized = regularized)
  } else {
    models <- calculate.arm(X, Y, windowSize, chromosome, arm, method, params, match.probes = match.probes, priors = priors, regularized = regularized)
  }


  params <- c(models, params, segmentedData = useSegmentedData)

  if(outputType == "data.frame"){
    message("Convert model to data.frame")
    return(as.data.frame(models))
  }
  else {
    return(models)
  }
}

