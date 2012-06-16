# GPPois.R: Gaussian Process-based inference for Poisson-noised data.
# Copyright (c) 2011 Charles R. Hogg III
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Author: Charles R. Hogg III (2011) <charles.r.hogg@gmail.com>

# FILE DESCRIPTION:
# Provides classes and functions to analyze data using Gaussian Processes.  Key
# examples include:
#  - Dataset: the data to be analyzed
#  - Covariance: (virtual) superclass for the different possible kinds of
#      covariance structures (e.g. Squared-Exponential, periodic, Matern, ...)
#  - Model: an interface to the total covariance structure (different
#      Covariance terms contribute *additively*)
#
# These objects *need* to be mutable, so I am using the ''R.oo'' library to
# enable pass-by-reference.
#
# A lot of computations result in really big matrices.  These should be easy to
# save and restore to disk.
#
# Also, I am trying to adhere to Google's R Style Guide as much as possible:
# http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html

# Object-orientation (with pass-by-reference!!)
library("R.oo")

################################################################################
# CLASS:                             Model
#
# A collection of covariance structures which can be trained on datasets and
# make predictions.
#
# Helper Functions:
#   LogML:  The logarithm of the (M)arginal (L)ikelihood for this model, given
#      a Dataset.
#   GradLogML:  The gradient (w.r.t. the parameter values) of the logarithm of
#      the (M)arginal (L)ikelihood for this model, given a Dataset.
#
# Virtual Fields (R = read only):
# (R) contributionIds:  The id strings of the Covariance objects which
#        contribute to this Model.
# ( ) id:  A character string used to identify this model.
# ( ) lower:  A named vector of lower bounds on parameters governing this
#     Model.
# ( ) params:  A named vector of parameters governing this Model.
# (R) signalIds:  The ids of contributions which are considered to be "signal",
#        i.e., not noise.
# ( ) upper:  A named vector of upper bounds on parameters governing this
#     Model.
# (R) varyingParamNames:  Names of parameters with nonzero range between upper
#     and lower bounds.
#
# Methods:
#   AddCovariance:  Add a new Covariance object to this Model.
#   clone:  Deep clone of this Model object (clones Covariance's as well).
#   Forget:  Clear all large matrices to save memory.
#   Freeze:  Make all parameter values constant.
#   L:  Lower-triangular Cholesky decomposition of the covariance matrix
#      (useful for generating random draws).
#   NamedCovariance:  Retrieve a clone of the first contributing Covariance
#      object with the given id.
#   PlotBubblingSurfaces2D:  Plot smoothly varying random surfaces to visualize
#      the uncertainty.
#   PlotCovariance2D:  Plot the covariance of points in 2D with respect to
#      another point; also plot the derivatives w.r.t. each parameter.
#   PosteriorInterval:  bounds on the uncertainty (in terms of the standard
#      deviation) in the prediction at a given point.
#   PosteriorMean:  The optimal prediction of the underlying function at a set
#      of points.
#   PosteriorStandardDeviation:  pointwise sigma (for *transformed* values) at
#      each prediction point.
#   PredictionMatrix:  A matrix relating function values at the output points
#      to function values at the input points.
#   print:  Prettied-up summary of this Model object.
#   SetNoiseBounds:  Add a 'noise' contribution with the given bounds.
#   Summary:  A data.frame summarizing predictions for multiple model subsets.
#   Train:  Optimize this Model's parameters to describe the given data.
#
# This class should work just fine, as long as\cr
#   a) we are training on all the datapoints together (i.e., not breaking them
#      up into subregions to divide-and-conquer), and\cr
#   b) this$params returns a vector which is amenable to simple optimization
#      routines (i.e., none of the parameters require special treatment).\cr
# If either of these conditions fail, a specialized subclass should be created
# instead.  Note that *both* these conditions fail for *both* scenarios
# considered in our Journal of Applied Crystallography paper, despite the fact
# that I wrote this software to perform the analysis for that paper.  That's
# OK; having these classes still makes it very much easier to build the
# specialized functions I need.

setMethodS3("PlotCovariance2D", "Model", conflict="quiet",
  function(this, d, file.name=NULL, file.type='png', i=1, ...) {
    # Plot the covariance of points in 2D with respect to another point; also
    # plot the derivatives w.r.t. each parameter.
    #
    # Args:
    #   d:  The Dataset to evaluate the Model on.
    #   file.name:  Character, the name of the file to plot to (file extension
    #      should be automatically appended if necessary)
    #   file.type:  Character, one of ("pdf", "png")
    #   i:  numeric; index of datapoint to use for calculating the covariance
    #
    # Returns:
    #   Used for its side-effects of plotting.
    #
    # Notes:
    #   Lay them out in a grid of two rows.  First column: top is covariance
    #   matrix, bottom is for text (parameter values).  Subsequent columns show
    #   derivative matrices.

    # Set the layout:
    n.rows <- 2
    n.cols <- 1 + ceiling(length(this$params) / n.rows)
    Layout <- grid.layout(nrow = n.rows, ncol = n.cols, 
      widths  = unit(rep(1, n.cols), rep("null", n.cols)), 
      heights = unit(rep(1, n.rows), rep("null", n.rows)))

    # Open the file:
    my.file <- SetupFileInfo(name=file.name, type=file.type)
    if (my.file$type == 'png') {
      if (!require("Cairo")) {
        png(filename=my.file$name, ...)
      } else {
        CairoPNG(filename=my.file$name, ...)
      }
    } else if (my.file$type == 'pdf') {
      cairo_pdf(filename=my.file$name, ...)
    } else {
      stop("Sorry, we can only do pdf or png output for now.")
    }
    LayoutNewGridPage(Layout=Layout)

    # Construct a data.frame for all the matrices
    point.out <- matrix(d$X[i, ], nrow=1)
    K.data <<- data.frame(X.1=d$X[, 1], X.2=d$X[, 2],
      Cov=this$KTotal(d=d)[i, ],
      M=t(this$PredictionMatrix(d=d, X.out=point.out,
          contributions=this$contributionIds[
            which(this$contributionIds != 'noise')]))
      )

    # Generic ggplot setup
    ps.out <- 1.3  # outer point size
    ps.in <- 1.0   # inner point size
    p.base <- (ggplot(data=K.data, aes(x=X.1, y=X.2))
      + scale_colour_gradientn(colours=c('white','red','blue'))
      + geom_point(colour='black', size=ps.out)
      + geom_vline(xintercept=K.data[i, 1])
      + geom_hline(yintercept=K.data[i, 2])
      )

    # Print the covariance matrix in the top left
    p.cov <- (p.base
      + geom_point(aes(colour=Cov), size=ps.in)
      )
    print(p.cov, vp=Subplot(1, 1))
    # Weight matrix in bottom left
    p.weight <- (p.base
      + geom_point(aes(colour=M), size=ps.in)
      + scale_colour_gradient2(colours=c('red', 'white', 'blue'))
      )
    print(p.weight, vp=Subplot(2, 1))

    dev.off()
    return (invisible(this))
  })

