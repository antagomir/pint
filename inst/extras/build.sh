#/usr/bin/R CMD BATCH document.R
/usr/bin/R CMD build ../../ --no-build-vignettes
/usr/bin/R CMD check --as-cran pint_1.25.11.tar.gz --no-build-vignettes
/usr/bin/R CMD INSTALL pint_1.25.11.tar.gz
#/usr/bin/R CMD BiocCheck pint_1.17.13.tar.gz 
