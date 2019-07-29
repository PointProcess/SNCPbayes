#!/bin/bash
pkg=$1
Rscript -e "Rcpp::compileAttributes('./$pkg')"
Rscript -e "devtools::document('./$pkg')"
rm $pkg/src/*.o
rm $pkg/src/*.so
R CMD INSTALL --no-lock $pkg
