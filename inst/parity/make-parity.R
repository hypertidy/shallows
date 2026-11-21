## Regenerate inst/parity/parity.json. Run from the package root with the
## package installed at the version being recorded:
##
##   R CMD INSTALL . && Rscript inst/parity/make-parity.R
##
## Do this only when results are meant to change (a fix or a change to the
## scheme), and say so in the commit. The JavaScript port in
## hypertidy/playproj then needs the same change and a fresh copy of the
## file.
shallows:::.parity_write("inst/parity/parity.json")
