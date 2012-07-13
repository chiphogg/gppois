\name{PlotBubblingSurfaces2D.Model}
\alias{PlotBubblingSurfaces2D}
\alias{PlotBubblingSurfaces2D.Model}
\title{Animated uncertainty in a surface}
\arguments{
  \item{this}{The Model object.}

  \item{d}{The Dataset to evaluate the Model on.}

  \item{X.out}{(matrix) The X-points where we want to
  predict the function.}

  \item{contributions}{(character vector) Id's of the
  contributing Covariances we want to predict (other
  contributions are considered noise).}

  \item{n.indep}{The number of independent draws to take
  for each datapoint.}

  \item{n.times}{The final number of animation frames.}

  \item{file.name}{(character) The basename of the file to
  plot to.}

  \item{...}{Not used.}
}
\description{
  This function helps visualize uncertainty using
  animations.  (See Details.)
}
\details{
  So, we have a 2D dataset, and a posterior distribution on
  the underlying true function.  We want to visualize the
  \emph{uncertainty} in that true function. How?

  One way to visualize a posterior distribution is to take
  a number of draws from it.  This is tricky for surfaces;
  they will tend to overlap and intersect and the plot will
  be too cluttered.  So instead of plotting the random
  draws simultaneously, why not show them one at a time?
  This does fix the clutter problem, but consecutive draws
  are unrelated to each other, so now we have a "jumpiness"
  problem.

  Both problems can be solved simultaneously, if the draws
  satisfy these conditions: \enumerate{ \item Consecutive
  draws are \emph{not} independent; they are correlated,
  and the correlation approaches 1 as the time difference
  approaches 0.  \item However, every draw does have the
  same \emph{marginal} distribution, which is equal to the
  posterior distribution we're trying to visualize. } The
  first condition makes the motion continuous.  The second
  makes sure it represents the uncertainty we're actually
  trying to represent.

  This software implements a novel solution which satisfies
  a stronger version of the first condition: specifically,
  the path of every point on the surface is differentiable
  infinitely many times.  Thus, we have the smoothest
  possible animations which actually represent the
  uncertainty we're trying to show.  I came up with this
  solution in May 2012; I am currently planning to write it
  up into a paper after I get back from Japan at the end of
  June.
}
\seealso{
  \code{\link{BubblingRandomMatrix}}

  \code{\link{Model}}
}

