# gppois

This package ([installation instructions below](#install)) is for people who want to **quantify uncertainty** about
**continuous functions**.  (See [a few examples from the physical sciences](https://github.com/chiphogg/gppois/wiki/Examples-of-continuous-functions-with-uncertainty).)

We can learn a lot about our function by *measuring its values* at various locations.
But we still have uncertainty, because...

  1. ...a *continuous* function has **(uncountably) infinitely many** values  
     (but we can only perform *finitely many measurements*)
  2. ...even the values we *do* measure are often subject to **noise**.

The best we can hope for is to compute the *probability* that any given
*candidate function* is the true curve.

Using probabilities to quantify uncertainty is known as [***Bayesian
Analysis***](http://bayesian.org/).  So, `gppois` performs Bayesian uncertainty
analysis for continuous functions.

## The Bayesian approach

To find the probability that a given function is the true curve, we ask two questions.

  1. How well does it fit the data we measured?
  2. How plausible is it in the first place?

These questions have standard names in Bayesian analysis:

### 1. Data fit: "Likelihood"

The likelihood depends on the noise model.

The easiest case is
[i.i.d.](http://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables)
Gaussian noise.  The computation becomes extremely simple: we only have to
invert a matrix to get probabilities for all our curves.

Poisson noise (a.k.a. "*counting* statistics") is harder: not only is the noise
level different at different places, but it depends on the value of the true
function!  (And if we knew *that*...)  The probabilities are still perfectly
computable, but it requires
[MCMC](http://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) simulation (and
all its associated headaches).

Fortunately, there's a trick which transforms Poisson-noised data into i.i.d.
Gaussian data.  It's called the [Anscombe
Transform](http://en.wikipedia.org/wiki/Anscombe_transform), and we use it
extensively in this package.

### 2. Plausibility: "Priors"

To define what makes a function "plausible", we need to state our assumptions.

We prefer *robust* assumptions -- e.g., the function is **continuous**, and (in some
sense) **smooth**.  But we don't want to tie ourselves to a particular
functional form without good reason.  How can we assign probabilities to
arbitrary curves, without assuming some functional form?

The answer is **Gaussian Processes** (see [our
paper](http://journals.iucr.org/j/issues/2012/03/00/to5012/index.html), or
[this excellent freely-available text](http://www.gaussianprocess.org/gpml/)).
We break a function into pieces: a function is simply a **collection of
values**, indexed by a **continuous variable**.  Each of these values has
uncertainty, and therefore a probability distribution.

But the individual values only tell half the story. *How they relate* to each other is the real key.

  - Values which are **very close** (in x) are **strongly correlated**
  - Values which are **far apart** are **practically independent**

This gives us smoothness and continuity, without tying us to a particular functional form.

### The package name

The `gp` stands for *Gaussian Processes*, and the `pois` stands for *Poisson noise*.

# <a id="install"></a> Installation

There are several ways to install `gppois`.  In any case, after installation, it can be used by typing

```r
library("gppois")
```
at the prompt in an `R` session.

## Using CRAN (recommended)

`gppois` is [on CRAN](http://cran.r-project.org/web/packages/gppois/).  This makes it very easy to install.

  1. Open an `R` session

  2. Execute one of the following commands, depending on your preferences:
     * **full install (recommended)**:

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

## Manual install

Alternatively, you can download the package file and install from the command line.
Note that you will *have* to take this approach if you want to install an old version.

## [devtools](https://github.com/hadley/devtools) / install_github

The `devtools` package includes a function, `install_github`, which will grab
and install the latest version from github.  I have not tested this myself!
However, it should be fairly straightforward.

Obviously, you will first need to install and load the `devtools` library.
Type or copy-paste the following commands into an `R` session:

```r
# Install the 'devtools' package from CRAN
install.packages("devtools")
# Load the 'devtools' package for this session
library("devtools")
```

Now you can use `devtools` to install `gppois`.  If you want to use the latest
development version, enter the following command from inside an `R` session:

```r
install_github(repo="gppois", username="chiphogg", branch="develop")
```

If you just want the latest stable version, type the following instead:

```r
install_github(repo="gppois", username="chiphogg")
```

Note that the latter will almost always be equivalent to the version on CRAN.

