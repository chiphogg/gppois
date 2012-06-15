#!/usr/bin/R

require("devtools")

# find.package was internal-only prior to R 2.13, but hadley uses it so I need
# to wrap it if I want to install with devtools.
#
# See:
# http://stat.ethz.ch/R-manual/R-patched/library/base/html/find.package.html
if (!exists("find.package")) {
  find.package <- .find.package
}
