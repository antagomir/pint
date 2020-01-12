#' The model parameters z and W
#' 
#' Contribution of each sample to a dependency model, and contribution of each
#' variable.
#' 
#' 
#' \code{z.effects} gives the contribution of each sample to the dependency
#' score. This is approximated by projecting original data to first principal
#' component of \code{Wz}. This is possible only when the data window is
#' smaller than half the number of samples.
#' 
#' \code{W.effects} gives the contribution of each variable to the observed
#' dependency. This is approximated with the loadings of the first principal
#' component of \code{Wz}
#' 
#' Original data can be retrieved by locating the row in \code{X} (or \code{Y})
#' which has the same variable (gene) name than \code{model}.
#' 
#' @aliases z.effects W.effects
#' @param model The fitted dependency model.
#' @param X,Y
#' 
#' Data sets used in fitting the dependency modeling functions
#' (\code{\link{screen.cgh.mrna}} or \code{link{fit.dependency.model}}). Note:
#' Arguments must be given in the same order as in
#' \code{\link{fit.dependency.model}} or \code{\link{screen.cgh.mrna}}.  Only
#' \code{X} is needed for dependency model for one data set.
#' @return \code{z.effects} gives a projection vector over the samples and
#' \code{W.effects} gives a projection vector over the variables.
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com} and Leo Lahti
#' \email{leo.lahti@@iki.fi}
#' @seealso \code{\link{DependencyModel-class}}, \code{\link{screen.cgh.mrna}}
#' @references Dependency Detection with Similarity Constraints, Lahti et al.,
#' 2009 Proc. MLSP'09 IEEE International Workshop on Machine Learning for
#' Signal Processing, See
#' \url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}
#' 
#' A Probabilistic Interpretation of Canonical Correlation Analysis, Bach
#' Francis R. and Jordan Michael I. 2005 Technical Report 688. Department of
#' Statistics, University of California, Berkley.
#' \url{http://www.di.ens.fr/~fbach/probacca.pdf}
#' 
#' Probabilistic Principal Component Analysis, Tipping Michael E. and Bishop
#' Christopher M. 1999. \emph{Journal of the Royal Statistical Society}, Series
#' B, \bold{61}, Part 3, pp. 611--622.
#' \url{http://research.microsoft.com/en-us/um/people/cmbishop/downloads/Bishop-PPCA-JRSS.pdf}
#' @keywords math
#' @examples
#' 
#' #data(chromosome17)
#' #window <- fixed.window(geneExp, geneCopyNum, 150, 10)
#' 
#' ### pSimCCA model around one gene
#' #depmodel <- fit.dependency.model(window$X, window$Y)
#' ## Conversion from DependencyModel to GeneDependencyModel so that
#' ## gene name and location can be stored
#' ##depmodel <- as(depmodel,"GeneDependencyModel")
#' #setGeneName(depmodel) <- window$geneName
#' #setLoc(depmodel) <- window$loc
#' #barplot(z.effects(depmodel, geneExp, geneCopyNum))
#' 
#' ## Plot the contribution of each genes to the model.
#' ## Only the X component is plotted
#' ## here since Wx = Wy (in SimCCA) 
#' #barplot(W.effects(depmodel, geneExp, geneCopyNum)$X)
#' 
#' ## plot.DpenendencyModel shows also sample and variable effects
#' #plot(depmodel,geneExp,geneCopyNum)
#' 
z.effects <- function(model, X, Y = NULL){

  W <- getW(model)

  # for models from 2 data sets
  if (!is.null(Y)){
    # Check if whole data is given instead window for this model
    if (class(X) == "list"){
      # Find correct window for this model
      index <- which(rownames(X$data) == getGeneName(model))

      # Check if model has only 1 variable from X data
      if (nrow(getW(model)$X) == 1) {
        window <- sparse.window(X, Y, index, getWindowSize(model))
        Xm <- window$X
        Ym <- window$Y
      } else {
        #window <- fixed.window(X, Y, index, getWindowSize(model))
        #X <- window$X
        #Y <- window$Y
        Xm <- t(centerData(t(X$data[rownames(W$X), ]), rm.na = TRUE))
        Ym <- t(centerData(t(Y$data[rownames(W$Y), ]), rm.na = TRUE))
      }
    }
    
    # Check that data window is smaller than half the sample size
    # FIXME: not required in all models; loosen this where possible
    if (getWindowSize(model) > ncol(X$data)/2)
      stop("The number of samples must be at least two times higher than number of features (probes)")

    W <- W$total
    z <- z.expectation(model, Xm, Ym)
    
    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]

    # Project data to this component
    data <- rbind(Xm, Ym)
    proj <- t(data)%*%projvec
    
    # Make sure the highest value is always positive
    # FIXME: later adjust sign based on observed data to make directy interpretable
    # this is ok fix for now
    if (abs(min(proj)) > max(proj))
      proj <- -proj   

    return(proj)
  } else {   # for models with one data set
    W <- W$total
    z <- z.expectation(model, X$data[rownames(W$X), ])

    # FIXME: for clarity, make own function for this since same operation is
    # used also above
    # Calculate first component of PCA for W*z
    pca <- princomp(t(W%*%z))
    projvec <- pca$loadings[,1]
      
    # Project data to this component
    data <- t(centerData(t(X$data[rownames(W$X), ]), rm.na = TRUE))
    proj <- t(data)%*%projvec
   
    # Make sure the highest value is allways positive
    if (abs(min(proj)) > max(proj))
      proj <- -proj

    return(proj)
  
  }
}




