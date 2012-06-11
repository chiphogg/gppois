#!/usr/bin/R

require("roxygen2")
require("R.oo")

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

# Document non-R.oo code (utility functions and the like):
roclet <- rd_roclet()
roc_out(roclet=roclet, input="R/utils.R", base_path=".")
roc_out(roclet=roclet, input="R/poisson.R", base_path=".")
roc_out(roclet=roclet, input="R/geometry.R", base_path=".")
roc_out(roclet=roclet, input="R/rgl.R", base_path=".")

# Document R.oo object code
old.home <- Sys.getenv("R_HOME")
Sys.setenv("R_HOME"="/usr/share/R")
author <- "Charles R. Hogg III"

Rdoc$compile(filename="R/Dataset.R", destPath="man/", source=TRUE)

Sys.setenv("R_HOME"=old.home)
