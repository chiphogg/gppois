#!/usr/bin/R

#-------------------------------------------------------------------------------
# grid/ggplot2 layout functions

LayoutNewGridPage <- function(Layout, ...) {
  # Setup a new grid page with the specified Layout.  
  # (Adapted from http://bit.ly/9m4zyD)
  #
  # Args:
  #   Layout:  Result of calling grid.layout function with the desired
  #      arrangement.
  #
  # Returns:
  #   Used for its side-effect.
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = Layout))
}

Subplot <- function(x, y) {
  # Terse wrapper to return a viewport into the specified row (x) and 
  # column (y) of a viewport.
  # (Adapted from http://bit.ly/9m4zyD)
  #
  # Args:
  #   x:  Numeric; the row(s) of the layout.
  #   y:  Numeric; the column(s) of the layout.
  #
  # Returns:
  #   Used for its side-effect.
  grid::viewport(layout.pos.row=x, layout.pos.col=y)
}
