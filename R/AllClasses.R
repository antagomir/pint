setClass("GeneDependencyModel", 
	  representation(
			 loc = "numeric",  
			 #locs = "numeric",
			 geneName = "character",
			 #geneNames = "character",
			 chromosome = "numeric", 
			 arm = "character"), 
			 prototype(geneName = "",  arm = ""),
			 contains="DependencyModel")

setClass("ChromosomeModels", representation(models = "list", chromosome = "numeric",
			method = "character", params = "list"))

setClass("GenomeModels", representation(chromosomeModels = "list", method = "character", params = "list"))
