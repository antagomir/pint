DW#' @title Fit dependency model around one gene between two data sets.
#' @description Takes a window from two datasets around chosen gene and fits a selected dependency model between windows.
#' @details See \code{\link{fit.dependency.model}} for details about dependency models and parameters.
#' @aliases fit.byname fit.cgh.mrna.byname fit.cgh.mir.byname
#' @param X,Y Data sets. Lists containing the following items: \describe{
#' \item{list("data")}{ Data in a matrix form. Genes are in columns and samples
#' in rows. e.g. gene copy number.  } \item{list("info")}{ Data frame which
#' contains following information about genes in data matrix.
#' \describe{ \item{list("chr")}{ Factor indicating the chrosome for the gene:
#' (1 to 23, or X or Y} \item{list("arm")}{ Factor indicating the chromosomal
#' arm for the gene ('p' or 'q')} \item{list("loc")}{ Location of the gene in
#' base pairs.} } } } \code{\link{pint.data}} can be used to create data sets
#' in this format.
#' @param geneName The dependency model is calculated around this gene.
#' @param windowSize Size of the data window.
#' @param ... Arguments to be passed to function
#' \code{\link{fit.dependency.model}}
#' @return \linkS4class{DependencyModel}
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com} and Leo Lahti
#' \email{leo.lahti@@iki.fi}
#' @seealso Reults from this function: \linkS4class{DependencyModel}.
#' \code{\link{fit.dependency.model}}. Calculating dependency models to
#' chromosomal arm, chromosome or genome \code{\link{screen.cgh.mrna}}. For
#' calculation of latent variable z: \code{link{z.expectation}}.
#' @references
#' 
#' Dependency Detection with Similarity Constraints, Lahti et al., 2009 Proc.
#' MLSP'09 IEEE International Workshop on Machine Learning for Signal
#' Processing,
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
#' 
#' EM Algorithms for ML Factorial Analysis, Rubin D. and Thayer D. 1982.
#' \emph{Psychometrika}, \bold{vol. 47}, no. 1.
#' @keywords math iteration
#' @examples \dontrun{
#' 
#' data(chromosome17)
#' 
#' model <- fit.cgh.mrna.byname(geneExp,geneCopyNum,"ENSG00000132361",10)
#' ## With different model parameters (pCCA)
#' model2 <- fit.cgh.mrna.byname(geneExp,geneCopyNum,"ENSG00000132361",10,zDimension=5,priors=list(Nm.wxwy.sigma = NULL))
#'}
NULL



#' Fits dependency models to chromosomal arm, chromosome or the whole genome.
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
#' scores across chromosomal windows.
#' 
#' Function \code{screen.cgh.mir} calculates dependencies around a chromosomal
#' window in each sample in \code{X}; only one sample from \code{X} will be
#' used. Data sets do not have to be of the same size and\code{X} can be
#' considerably smaller. This is used with e.g.  miRNA data.
#' 
#' If method name is specified, this overrides the corresponding model
#' parameters, corresponding to the modeling assumptions of the specified
#' model. Otherwise method for dependency models is determined by parameters.
#' 
#' Dependency scores are plotted with \link{dependency score plotting}.
#' 
#' @aliases screen.cgh.mrna screen.cgh.mir
#' @param X,Y Data sets.  It is recommended to place gene/mirna expression data
#' in X and copy number data in Y.  Each is a list with the following items:
#' \describe{ \item{list("data")}{ Data in a matrix form. Genes are in rows and
#' samples in columnss. e.g. gene copy number.}
#' 
#' \item{list("info")}{ Data frame which contains following information about
#' genes in data matrix.
#' 
#' \describe{ \item{list("chr")}{ Number indicating the chrosome for the gene:
#' (1 to 24). Characters 'X' or 'Y' can be used also.} \item{list("arm")}{
#' Character indicating the chromosomal arm for the gene ('p' or 'q')}
#' \item{list("loc")}{ Location of the gene in base pairs.} } } }
#' \code{\link{pint.data}} can be used to create data sets in this format.
#' @param chromosome Specify the chromosome for model fitting. If missing,
#' whole genome is screened.
#' @param arm Specify chromosomal arm for model fitting. If missing, both arms
#' are modeled.
#' @param windowSize Determine the window size. This specifies the number of
#' nearest genes to be included in the chromosomal window of the model, and
#' therefore the scale of the investigated chromosomal region. If not
#' specified, using the default ratio of 1/3 between features and samples or
#' \code{15} if the ratio would be greater than 15
#' @param method Dependency screening can utilize any of the functions from the
#' package dmt (at CRAN). Particular options include
#' 
#' \describe{ \item{'pSimCCA'}{probabilistic similarity constrained canonical
#' correlation analysis \cite{Lahti et al. 2009}. This is the default method.}
#' \item{'pCCA'}{probabilistic canonical correlation analysis \cite{Bach &
#' Jordan 2005}} \item{'pPCA'}{probabilistic principal component analysis
#' \cite{Tipping & Bishop 1999}} \item{'pFA'}{probabilistic factor analysis
#' \cite{Rubin & Thayer 1982}} \item{'TPriorpSimCCA'}{probabilistic similarity
#' constrained canonical correlation analysis with possibility to tune T prior
#' (Lahti et al. 2009)} } If anything else, the model is specified by the given
#' parameters.
#' @param params List of parameters for the dependency model.  \describe{
#' \item{sigmas}{Variance parameter for the matrix normal prior distribution of
#' the transformation matrix T. This describes the deviation of T from H}
#' \item{H}{Mean parameter for the matrix normal prior distribution prior of
#' transformation matrix T} \item{zDimension}{Dimensionality of the latent
#' variable} \item{mySeed}{Random seed.} \item{covLimit}{Convergence limit.
#' Default depends on the selected method: 1e-3 for pSimCCA with full marginal
#' covariances and 1e-6 for pSimCCA in other cases.} }
#' @param max.dist Maximum allowed distance between probes. Used in automated
#' matching of the probes between the two data sets based on chromosomal
#' location information.
#' @param outputType Specifies the output type of the function. possible values
#' are \code{"models"} and \code{"data.frame"}
#' @param useSegmentedData Logical. Determines the useage of the method for
#' segmented data
#' @param match.probes To be used with segmented data, or nonmatched probes in
#' general. Using nonmatched features (probes) between the data sets.
#' Development feature, to be documented later.
#' @param regularized Regularization by nonnegativity constraints on the
#' projections. Development feature, to be documented later.
#' @return The type of the return value is defined by the the function argument
#' \code{outputType}.
#' 
#' With the argument \code{outputType = "models"}, the return value depends on
#' the other arguments; returns a \linkS4class{ChromosomeModels} which contains
#' all the models for dependencies in chromosome or a
#' \linkS4class{GenomeModels} which contains all the models for dependencies in
#' genome.
#' 
#' With the argument \code{outputType = "data.frame"}, the function returns a
#' data frame with eachs row representing a dependency model for one gene. The
#' columns are:
#' \code{geneName},\code{dependencyScore},\code{chr},\code{arm},\code{loc}.
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com} and Leo Lahti
#' \email{leo.lahti@@iki.fi}
#' @seealso To fit a dependency model: \code{\link{fit.dependency.model}}.
#' \linkS4class{ChromosomeModels} holds dependency models for chromosome,
#' \linkS4class{GenomeModels} holds dependency models for genome. For plotting,
#' see: \link{dependency score plotting}
#' @references
#' 
#' Dependency Detection with Similarity Constraints, Lahti et al., 2009 Proc.
#' MLSP'09 IEEE International Workshop on Machine Learning for Signal
#' Processing, See
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
#' 
#' EM Algorithms for ML Factoral Analysis, Rubin D. and Thayer D. 1982.
#' \emph{Psychometrika}, \bold{vol. 47}, no. 1.
#' @keywords math iteration
#' @examples
#' 
#' data(chromosome17)
#' 
#' ## pSimCCA model on chromosome 17
#' 
#' models17pSimCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
#'                                      windowSize = 10, chr = 17)
#'                                     
#' plot(models17pSimCCA)
#' 
#' ## pCCA model on chromosome 17p with 3-dimensional latent variable z
#' models17ppCCA <- screen.cgh.mrna(geneExp, geneCopyNum,
#'                                    windowSize = 10,
#'                                    chromosome = 17, arm = 'p',method="pCCA", 
#' 	      	 	           params = list(zDimension = 3))
#' plot(models17ppCCA)
#' 
NULL



#' Class "ChromosomeModels"
#' 
#' Collection of dependency models fitting two data sets in particular
#' chromosome.
#' 
#' 
#' @name ChromosomeModels-class
#' @aliases ChromosomeModels-class getPArm getQArm [[ [[<- getChromosome getArm
#' getGeneName getModelMethod getParams getWindowSize topGenes topModels
#' getModelNumbers isEmpty orderGenes findModel [[,ChromosomeModels-method
#' [[<-,ChromosomeModels-method getChromosome,ChromosomeModels-method
#' getPArm,ChromosomeModels-method getQArm,ChromosomeModels-method
#' getArm,ChromosomeModels-method getGeneName,ChromosomeModels-method
#' getLoc,ChromosomeModels-method getScore,ChromosomeModels-method
#' getModelMethod,ChromosomeModels-method getParams,ChromosomeModels-method
#' getWindowSize,ChromosomeModels-method topGenes,ChromosomeModels-method
#' topModels,ChromosomeModels-method getModelNumbers,ChromosomeModels-method
#' isEmpty,ChromosomeModels-method orderGenes,ChromosomeModels-method
#' findModel,ChromosomeModels-method isEmpty,ChromosomeModels-method
#' [[,ChromosomeModels-method as.data.frame,ChromosomeModels-method
#' @docType class
#' @section Objects from the Class: Function \code{\link{screen.cgh.mrna}} and
#' \code{\link{screen.cgh.mir}} returns an object of this class.
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com}
#' @seealso For calculation of dependency models for chromosomal arm:
#' \code{\link{screen.cgh.mrna}}. This class holds a number of
#' \linkS4class{GeneDependencyModel} objects. For plotting dependency scores
#' see \link{dependency score plotting}.  Dependency models for whole genome:
#' \linkS4class{GenomeModels}.
#' @keywords classes
#' @examples \dontrun{
#' data(chromosome17)
#' ## calculate dependency models over chromosome 17
#' model17 <- screen.cgh.mrna(geneExp, geneCopyNum, windowSize = 10, chr
#' = 17) 
#' model17
#' ## Information of the dependency model which has the highest dependency score
#' topGenes(model17, 1)
#' ## Finding a dependency model by its name
#' findModel(model17, "ENSG00000129250")
#' ## Information of the first dependency model
#' model17[[1]]
#' #Plotting
#' plot(model17)
#' # genes in p arm with the highest dependency scores
#' topGenes(model17[['p']], 5)
#' } 
#' 
NULL





#' @title Gene copy number data in chromosome 17
#' @description Preprocessed gene copy number (aCGH) data for 51 patients in chromosome 17. 
#' @name geneCopyNum
#' @docType data
#' @format A list which contain the following data: \describe{
#' \item{data}{ gene copy number data in matrix form. Genes are in columns and
#' samples in rows}
#' \item{info}{ Data frame which contains following information about genes in
#' data matrix.
#' \describe{
#' 
#' \item{chr}{ Factor indicating the chrosome for the gene (1 to 23, or X or Y}
#' 
#' \item{arm}{ Factor indicating the chromosomal arm for the gene ('p' or 'q')}
#' 
#' \item{loc}{ Location of the gene in base pairs.}
#' 
#' } } }
#' @source Integrated gene copy number and expression microarray analysis of
#' gastric cancer highlights potential target genes.  Myllykangas et al.,
#' \emph{International Journal of Cancer}, vol. \bold{123}, \bold{no. 4}, pp.
#' 817--25, 2008.
#' @keywords datasets
NULL





#' Class "GeneDependencyModel"
#' 
#' A Genomic Dependency model for two data sets
#' 
#' 
#' @name GeneDependencyModel-class
#' @aliases GeneDependencyModel-class setLoc<- setGeneName<- setChromosome<-
#' setArm<- getScore getLoc setLoc<-,GeneDependencyModel-method
#' setGeneName<-,GeneDependencyModel-method
#' setChromosome<-,GeneDependencyModel-method
#' setArm<-,GeneDependencyModel-method getScore,GeneDependencyModel-method
#' getLoc,GeneDependencyModel-method getGeneName,GeneDependencyModel-method
#' getWindowSize,GeneDependencyModel-method
#' getChromosome,GeneDependencyModel-method getArm,GeneDependencyModel-method
#' getZ,GeneDependencyModel-method
#' @docType class
#' @section Objects from the Class: Used to represent individual dependency
#' models for screening inside \linkS4class{ChromosomeModels}.
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com}
#' @seealso For calculation of dependency models for chromosomal arm,
#' chromosome or genome: \code{\link{screen.cgh.mrna}}. Dependency models for
#' whole chromosome: \linkS4class{ChromosomeModels}.  Dependency models for
#' whole genome: \linkS4class{GenomeModels}.  For plotting dependency scores
#' see \link{dependency score plotting}.
#' @keywords classes
#' @examples
#' 
#' data(chromosome17)
#' 
#' # First genomic dependency model from screening chromosomal arm
#' models <- screen.cgh.mrna(geneExp, geneCopyNum, 10, chr=17, arm='p')
#' model <- models[[1]]
#' 
#' # Printing information of the model
#' model
#' 
#' # Latent variable Z
#' getZ(model, geneExp,geneCopyNum)
#' 
#' # Contributions of samples and variables to model
#' plot(model,geneExp,geneCopyNum)
#' 
NULL





#' Gene expression data in chromosome 17
#' 
#' Preprocessed gene expression levels of 51 patients in chromosome 17.
#' 
#' 
#' @name geneExp
#' @docType data
#' @format A list which contain the following data: \describe{
#' 
#' \item{data}{ gene expression data in matrix form. Genes are in columns and
#' samples in rows}
#' 
#' \item{info}{ Data frame which contains following information about genes in
#' data matrix.
#' 
#' \describe{
#' 
#' \item{chr}{ Factor of chrosome where the gene is. (1 to 23 or X or Y}
#' 
#' \item{arm}{ Factor of arm of the chromosome arm where the gene is. ('p' or
#' 'q')}
#' 
#' \item{loc}{ Location of the gene from centromere in base pairs.}
#' 
#' } } }
#' @source Integrated gene copy number and expression microarray analysis of
#' gastric cancer highlights potential target genes.  Myllykangas et al.,
#' \emph{International Journal of Cancer}, vol. \bold{123}, \bold{no. 4}, pp.
#' 817--25, 2008.
#' @keywords datasets
NULL





#' Class "GenomeModels"
#' 
#' Collection of dependency models fitting two data sets in whole genome. The
#' dependency models are in a list of \linkS4class{ChromosomeModels}s (which
#' represents each chromosome) that have a list of dependency models in that
#' chromosomal arm.
#' 
#' 
#' @name GenomeModels-class
#' @aliases GenomeModels-class [[,GenomeModels-method [[<-,GenomeModels-method
#' getModelMethod,GenomeModels-method getParams,GenomeModels-method
#' getWindowSize,GenomeModels-method topGenes,GenomeModels-method
#' topModels,GenomeModels-method orderGenes,GenomeModels-method
#' findModel,GenomeModels-method getModelNumbers,GenomeModels-method
#' as.data.frame,GenomeModels-method
#' @docType class
#' @section Objects from the Class: Function \code{\link{screen.cgh.mrna}} and
#' \code{\link{screen.cgh.mir}} returns an object of this class.
#' @author Olli-Pekka Huovilainen
#' @seealso For calculation of dependency models for chromosomal arm:
#' \code{\link{screen.cgh.mrna}}. This class holds a number of
#' \linkS4class{GeneDependencyModel} in each \linkS4class{ChromosomeModels}.
#' For plotting dependency scores see \code{\link{dependency score plotting}}.
#' @keywords classes
NULL





#' Dependency score plotting
#' 
#' Plot the contribution of the samples and variables to the dependency model
#' or dependency model fitting scores of chromosome or genome.
#' 
#' Function plots scores of each dependency model of a gene for the whole
#' chromosome or genome according to used method. \code{plot(x, cancerGenes =
#' NULL, showDensity = FALSE, ...)} is also usable and chosen according to
#' class of \code{models}.
#' 
#' @aliases plot.GeneDependencyModel plot.ChromosomeModels plot.GenomeModels
#' dependency score plotting
#' @param x \code{\link{GeneDependencyModel-class}},
#' \code{\link{ChromosomeModels-class}}, \code{\link{GenomeModels-class}};
#' models to be plotted.
#' @param X,Y data sets used in dependency modeling.
#' @param ann.types a factor for annotation types for samples. Each value
#' corresponds one sample in datasets. Colors are used to indicate different
#' types.
#' @param ann.cols colors used to indicate different annotation types. Gray
#' scale is used if 'NULL' given.
#' @param legend.x,legend.y the x and y co-ordinates to be used to position the
#' legend for annotation types.
#' @param legend.xjust,legend.yjust how the legend is to be justified relative
#' to the legend x and y location.  A value of 0 means left or top justified,
#' 0.5 means centered and 1 means right or bottom justified.
#' @param order logical; if 'TRUE', values for sample contributions are ordered
#' according to their values.
#' @param cex.z,cex.WX,cex.WY Text size for variable names.
#' @param hilightGenes vector of strings; Name of genes to be hilighted with
#' dots.
#' @param showDensity logical; if 'TRUE' small vertical lines are drwan in the
#' bottom of the plot under each gene.
#' @param showTop numeric; Number of models with highest dependencies to be
#' hilighted. A horizontal dashed line is drawn to show threshold value. With
#' \code{0} no line is drawn.
#' @param topName logical; If \code{TRUE}, gene names are printed to hilighted
#' models with highest dependecies. Otherwise hilighted models are numbered
#' according to their rank in dependency score.
#' @param type,xlab,ylab,main plot type and labels. See \code{\link{plot}} for
#' details.  A text for chromosome (and arm if only models from one arm is
#' plotted) is used in \code{main} if \code{NULL} is given. In
#' \code{plot.GenomeModels}, \code{ylab} and \code{xlab} affect only if
#' \code{onePlot} is \code{TRUE}.
#' @param onePlot If \code{TRUE}, all dependency scores are plotted in one plot
#' window. Otherwise one plot window is used for each chromosome.
#' @param pch,cex symbol type and size for hilightGenes. See
#' \code{\link{points}} for details.
#' @param tpch,tcex symbol type and size for genes with highest scores. See
#' \code{\link{points}} for details.
#' @param ylim,xlim axis limits. Default values are calculated from data. Lower
#' limit for y is 0 and upper limit is either 1 or maximum score value. X
#' limits are gene location range. See \code{\link{plot}} for details.
#' @param mfrow,mar,ps,mgp chromosome plots' layout, marginals, text size and
#' margin line. See \code{\link{par}} for details.
#' @param ...  optional plotting parameters
#' @author Olli-Pekka Huovilainen \email{ohuovila@@gmail.com}
#' @seealso \code{\link{DependencyModel-class}},
#' \code{\link{ChromosomeModels-class}}, \code{\link{GenomeModels-class}},
#' \code{\link{screen.cgh.mrna}}, \code{\link{screen.cgh.mir}}
#' @references Dependency Detection with Similarity Constraints Lahti et al.,
#' MLSP'09. See
#' \url{http://www.cis.hut.fi/lmlahti/publications/mlsp09_preprint.pdf}
#' @keywords hplot
#' @examples
#' 
#' 
#' data(chromosome17)
#' 
#' ## pSimCCA model on chromosome 17p
#' models17ppSimCCA <- screen.cgh.mrna(geneExp, geneCopyNum, 10, 17, 'p')
#' plot(models17ppSimCCA,
#'      hilightGenes=c("ENSG00000108342", "ENSG00000108298"), showDensity = TRUE)
#' 
#' ## Dependency model around 50th gene
#' model <- models17ppSimCCA[[50]]
#' 
#' ## example annnotation types
#' ann.types <- factor(c(rep("Samples 1 - 10", 10), rep("Samples 11 - 51", 41)))
#' plot(model, geneExp, geneCopyNum, ann.types, legend.x = 40, legend.y = -4,
#'      order = TRUE)
#' 
#' 
#' 
NULL



Collate:
    AllClasses.R
    AllGenerics.R
    GeneDependencyModel-accessors.R
    ChromosomeModels-accessors.R
    GenomeModels-accessors.R
    imputation.R
    show-methods.R
    plot-methods.R
    screen.cgh.mrna.R
    screen.cgh.mir.R
    calculate.genome.R
    calculate.chr.R
    calculate.arm.R
    pick.chr.arm.R
    fixed.window.R
    iterative.window.R
    z.effects.R
    W.effects.R
    fit.cgh.mrna.byname.R
    fit.cgh.mir.byname.R
    report.R
    sparse.window.R
    calculate.arm.sparse.R
    calculate.chr.sparse.R
    calculate.genome.sparse.R
    pint.data.R
    pint.match.R
    firstlib.R
    test.segmented.R
    centerData.R
