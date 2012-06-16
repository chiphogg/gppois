#!/usr/bin/R

require("roxygen2")
require("devtools")
require("R.oo")

source("setup_env.R")

# This function taken from the editrules package (and slightly reformatted).
# The roxygen2 package uses the list2env function, which the Debian squeeze
# version of R (2.11) lacks.  Adding this function enables me to
# auto-generate documentation of my non-R.oo code.
if (!exists("list2env")){
  list2env <- function(x, envir=NULL, parent=parent.frame()){
    if (is.null(envir)) {
      envir <- new.env(parent=parent)
    }
    for (i in names(x)){
      envir[[i]] <- x[[i]]
    }
    return (envir)
  }
}

non.r.oo.files <- paste(sep="", "R/", c("gppois.R",
    "geometry.R", "layout.R", "poisson.R", "rgl.R", "utils.R"))
r.oo.files <- paste(sep="", "R/", c("Dataset.R", "LazyMatrix.R",
    "Covariance.R", "CovarianceNoise.R", "CovarianceSE.R",
    "CovarianceSELocalized.R", "CovarianceSEVaryingEll.R",
    "CovarianceSEAniso2D"))
sourcefiles <- c(non.r.oo.files, r.oo.files)

# Document non-R.oo code (utility functions and the like):
roclet <- rd_roclet()
for (filename in sourcefiles) {
  roc_out(roclet=roclet, input=filename, base_path=".")
}

# Document R.oo object code; steps listed here:
# http://www.aroma-project.org/developers
install("gppois")

# KEYWORDS.db may be in a different location; if so, the following hack makes
# it work for R.oo<=1.9.3:
rhome <- Sys.getenv("R_HOME")
sane.install <- file.exists(file.path(rhome, "/doc/KEYWORDS.db"))
if (!sane.install) {
  old.home <- rhome
  doc.dir <- sub(x=Sys.getenv("R_DOC_DIR"), pattern="/doc", replacement="")
  Sys.setenv("R_HOME"=doc.dir)
}

author <- "Charles R. Hogg III"

for (filename in r.oo.files) {
  Rdoc$compile(filename=filename, destPath="man/", source=TRUE, verbose=TRUE)
}

if (!sane.install) {
  Sys.setenv("R_HOME"=old.home)
}

# Doing this again adds the new changes to the R.oo parts of the documentation
install("gppois")
