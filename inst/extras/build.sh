#~/bin/R-3.6.2/bin/R CMD BATCH document.R
~/bin/R-3.6.2/bin/R CMD build ../../ --no-build-vignettes
~/bin/R-3.6.2/bin/R CMD check --as-cran pint_1.37.2.tar.gz --no-build-vignettes
~/bin/R-3.6.2/bin/R CMD INSTALL pint_1.37.2.tar.gz
#~/bin/R-3.6.2/bin/R CMD BiocCheck pint_1.17.13.tar.gz 
