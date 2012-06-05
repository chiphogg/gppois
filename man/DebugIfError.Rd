\name{DebugIfError}
\alias{DebugIfError}
\title{Runs a function with just-in-time debugging}
\usage{
  DebugIfError(FUN, ...)
}
\arguments{
  \item{FUN:}{The function to call (either a character
  string or just the name of the function).}

  \item{...:}{The argument list for \code{FUN.}}
}
\value{
  The result of calling \code{FUN(...)}, or else the
  exception.
}
\description{
  This function tries running \code{FUN(...)} normally.  If
  it's fine, the user should see no difference.  If not, it
  catches the exception, turns on debugging, and re-runs to
  assist debugging.  Debugging is turned off for this
  function at the end.
}

