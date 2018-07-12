calculate.arm <- function(X, Y, windowSize, chromosome, arm = NULL, method =
"pSimCCA", params = list(), match.probes = FALSE, priors = NULL, regularized = FALSE){

  # Get Indices to X and Y for chosen chromosome and arm
  Xm <- pick.chr.arm(X, chromosome, arm)
  Ym <- pick.chr.arm(Y, chromosome, arm)

  # Storage for dependency scores
  scores <- vector()

  # Storage for gene window location                                         
  locs <- vector()

  # Storage for gene names
  genes <- vector()

  # method name
  methodName <- method
  message(paste("Calculating dependency models for ",chromosome,arm," with method ",methodName, 
		", window size:",windowSize,sep=""))
	
  modelList <- list()
	
  # index for modelList
  k <- 1
	
 if (length(Ym$info$loc) > 0) {                           
   for (n in seq_along(Ym$info$loc)) {       
     message(chromosome, arm, "; window ", n, "/", length(Ym$info$loc))                  
     # Assuming Y is the copy number data (important when match.probes= TRUE)                                                              
        # Get window to dependency modeling                                
        if (!match.probes) {                          
	  #message("Using fixed chromosomal window size.")       


#' Form data with a selected window size for the model fitting
#' 
#' 
#' Forms a chosen window of two data matrices to use for
#' \code{fit.dependency.model} either iteratively picking nearest genes or
#' picking same number of genes from both directions. \code{sparse.window}
#' forms a window around one sample in the first data set with a number of
#' samples from the second data set.
#' 
#' Window contains windowSize nearest genes. Warning is given if windowSize
#' genes is not found in the same chromosomal arm. Data of both data sets is
#' normalised so that each genes data has zero mean.
#' 
#' @aliases fixed.window iterative.window sparse.window
#' @param X First data set. In \code{sparse.window} windows will be formed
#' around each sample in this data set.
#' @param Y Second data set.
#' @param middleIndex Index of middle position for window.
#' @param xIndex Index of middle position in \code{X} for window.
#' @param windowSize Number of genes in window. In \code{sparse.window}
#' \code{X} has always one sample in window.
#' @return List of window data: \item{X}{window of the first data set}
#' \item{Y}{window of the second data set} \item{loc}{location of gene}
#' \item{geneName}{name of the gene} \item{edge}{logical; TRUE if iteration to
#' one direction has stopped because edge of data in chromosomal arm has been
#' found.} \item{fail}{logical; TRUE if chromosomal arm contains less than
#' windowSize genes.}
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com}
#' @seealso Dependency model fitting: \code{\link{fit.dependency.model}}
#' @examples
#' 
#' data(chromosome17)
#' window <- iterative.window(geneExp, geneCopyNum, 30, 10)
#' model <- fit.dependency.model(window$X, window$Y)
#' 
#' # Conversion from DependencyModel to GeneDependencyModel so that gene name and location can be stored
#' model <- as(model,"GeneDependencyModel")
#' setGeneName(model) <- window$geneName
#' setLoc(model) <- window$loc
#' model
#' 
#' window <- fixed.window(geneExp, geneCopyNum, 10, 10)
#' model <- fit.dependency.model(window$X, window$Y)
#' model
#' 
	  window <- fixed.window(Xm, Ym, n, windowSize)
	} else {                           
	  message("Selecting the closest expression probes for each copy number segment.")       
          #NOTE: this calculates for each in X, window in Y      
	  # therefore convert as we have copy number in Y        
	  tmp <- sparse.window(Ym, Xm, n, windowSize)  
	  window <- tmp          
	  window$X <- tmp$Y      
	  window$Y <- tmp$X      
	}
                                                   
	# Skip windows that overlaps chromosome arms                    
       if (!window$fail){                     
   
           if (!match.probes && !regularized) {               
	     #print(c("no-segment", "no-regu"))         
	     # still ensure that no regularization used here:                 
             priors$W <- NULL                           
             model <- fit.dependency.model(window$X,         
                                           window$Y,   
                                           zDimension = params$zDimension,
                                           marginalCovariances = params$marginalCovariances,      
                                           priors = list(Nm.wxwy.mean = params$H,          
                                           Nm.wxwy.sigma = params$sigmas),includeData = FALSE, calculateZ = FALSE)                                     
    
  
             } else if (!match.probes && regularized) {                     
	           #print(c("no-segment", "yes-regu"))             
		   # 1-dimensional cca (in general, Wx != Wy) with nonnegative W                 
                   # assuming matched probes                                
               model <- fit.dependency.model(window$X, window$Y, priors = priors)                 
                                                                 
            } else if (match.probes) {            
	      message("Note: regularization is used with segmented data.")                 
              regularized <- TRUE                          
	      # always use positive prior for W here     
	      if (is.null(priors$W)) {priors$W <- 1e-3} # uninformative                                  
	      	 model <- fit.dependency.model(window$X, window$Y,zDimension = params$zDimension, priors = priors)                                   
                                                                      
     	   }
           model <- as(model, "GeneDependencyModel")
           #setGeneName(model) <- rownames(window$X)[[ trunc((nrow(window$X) + 1)/2) ]]
           setGeneName(model) <- window$geneName
          #setLocs(model) <- window$loc        
          setLoc(model) <- window$loc
        setChromosome(model) <- chromosome
        if (!is.null(arm)) setArm(model)  <- arm          
        modelList[[k]] <- model
        k <- k + 1                                 
      }                
    }
  }
 #Change chromosome and arm factors and get levels from X          
 #chromosome <- factor(chromosome, levels = c(1:22,"X","Y"))     
 #arm <- factor(arm, levels = levels(X$info$arm))      
 
 return(new("ChromosomeModels",                                  
            models = modelList,                                   
            chromosome = chromosome,    
            method = method,                                                
            params = params))                           
}                     

