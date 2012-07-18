# gppois

`gppois` is all about *quantifying uncertainty* in *continuous functions*.

  * **quantifying uncertainty...**
    It's not good enough to say what we think the right answer is.  We also
    have to say how wrong we think we could be.  And to do science, we have to
    be *quantitative* -- which means we need to choose a *mathematical
    language* for uncertainty.

    We choose the language of *probability*.  For every answer which *might* be
    correct, we compute a *probability* that it is.
    * Quantifying uncertainty using **probabilities** is known as [***Bayesian Analysis***](http://bayesian.org/).
  * **continuous functions...**
    * They are *everywhere.*  Examples:
      * ...

# Installation

## Using CRAN (recommended)

`gppois` is [on CRAN](http://cran.r-project.org/web/packages/gppois/).  This makes it very easy to install.

  1. Open an `R` session
  2. Execute one of the following commands, depending on your preferences:
     * **full install** (recommended):
       ```r
       # If you have R 2.15 or newer:
       install.packages("gppois", dependencies=TRUE)
       # If you have R 2.14.xx, uncomment and use this instead:
       #install.packages("gppois", dependencies=c("Depends", "Imports", "LinkingTo", "Suggests")
       ```
     * bare-bones install:
       ```r
       install.packages("gppois")
       ```

And you're done!   

# Package etymology

gppois stands for **G**aussian **P**rocesses for **Pois**son noised data.

(elaborate...)
