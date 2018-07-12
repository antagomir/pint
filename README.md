pint
====

Probabilistic dependency analysis for functional genomics.
In particular, tools for integrative screening
of genome-wide RNA and DNA profiling data for cancer gene discovery 
and functional analysis of chromosomal aberrations (Lahti et al. MLSP 2009).

The package implements probabilistic models for integrative analysis
of mRNA expression levels with DNA copy number (aCGH) measurements to
discover functionally active chromosomal alterations. The algorithms
can be used to discover functionally altered chromosomal regions and
to visualize the affected genes and samples. The algorithms can be
applied also to other types of biomedical data, including epigenetic
modifications, SNPs, alternative splicing and transcription factor
binding, or in other application fields. By investigating dependencies
between different functional layers of the genome it is possible to
discover mechanisms and interactions that are not seen in the
individual measurement sources. For instance, integration of gene
expression and DNA copy number can reveal cancer-associated
chromosomal regions and associated genes with potential diagnostic,
prognostic and clinical impact.

The methods are based on latent variable models including
probabilistic canonical correlation analysis and related extensions,
implemented in the
[dmt package](http://cran.fhcrc.org/web/packages/dmt/index.html).
Probabilistic formulation deals rigorously with uncertainty associated
with small sample sizes common in biomedical studies and provides
tools to guide dependency modeling through Bayesian priors.

This is the development version of the package. For a stable release
version, see
[Bioconductor](http://bioconductor.org/packages/release/bioc/html/pint.html).
For further details and examples, see the package
[vignette](vignettes/depsearch.Rnw).




