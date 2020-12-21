
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmcperturb

<!-- badges: start -->

<!-- badges: end -->

The goal of lmcperturb is to provide code to construct a uniformly valid confidence set for a linear function of partially identified parameters defined by linear moment constraints using the procedure described in Cho and Russell (2020, <arXiv:1810.03180>).

## Installation

You can install this GitHub version of lmcperturb with:

``` r
devtools::install_github("cho-joonh/lmcperturb")
```

The package requires R CRAN package `lpSolve` and alternatively suggests the package `Rmosek`. `Rmosek` is suggested so that the user utilizes its large-scale linear programming capabilities that enhances computational performance when the dimension of the parameter space and/or the number of moment (in)equalities are large. `lpSolve` can be installed from CRAN to your library with the usual command `install.packages()`. However, `Rmosek` on CRAN is only a meta-package designed to help with the installation of `Rmosek`. To install and run `Rmosek`, you must request a (free personal academic) licence and follow the installation instructions in Rmosek document on
<https://www.mosek.com/documentation/>.

## Description

The main function `perturbCS` has required inputs:

  + `data`: a data matrix;
  + `psi_fn`: a function that maps data to coefficients of linear functional of interest (the first element must be the intercept);
  + `mom_constr`: a list containing two functions as elements, where the first function maps data to matrix of coefficients of linear moment constraints and the second function maps data to lower and upper bounds of the moment constraints;

and has optional inputs with default values when not user-specified:

  + `par_constr`: a list containing three matrices as elements, where the first matrix has two rows that define the parameter space constraints for each coordinate, the second matrix is a linear map of the parameter vector that together with the third matrix (for lower and upper bounds) determine the boundaries of the parameter space (default `par_constr = NULL`);
  + `max_perturb`: a list containing two vectors as elements, where the first vector of length two determines the support of the perturbation distribution for the linear functional of interest, and the second vector of length two determines that for the moment constraints (default `max_perturb = list(c(-1e-6,1e-6), c(0,1e-6))`);
  + `alpha`: a numeric vector (possibly of length one) for desired significance level of the confidence set to be constructed (default `alpha = c(.10,.05,.01)`);
  + `tol`: a nonnegative number for the tolerance level of finite-sample thresholding (default `tol = 1`);
  + `boot`: a positive integer for bootstrap sample size (default `boot = 999`);
  + `method`: a string; either `"lpSolve"` or `"Rmosek"`.

The lower or upper bound of a constraint in `mom_constr` and `par_constr` may be `-Inf` or `Inf` if no finite bound. When the model has no natural parameter space constraint and the constraint is set to `par_constr = NULL`, default parameter space constraints `[-1e+7,1e+7]` for individual coordinates are imposed, i.e., a large hypercube, to compactify the search space. If a larger search space is desired, `par_constr` should be specified.

  + If the parameter space is defined only by bounds on individual coordinates, i.e., the parameter space is a hyperrectangle `lb <= x <= ub` where `x` is a parameter vector, then define the list `par_constr` with a matrix of bounds `par_constr[[1]] = rbind(lb, ub)` and set `par_constr[[2]] = NULL` and `par_constr[[3]] = NULL`.
  + Similarly, if the parameter space is defined only by jointly defined constraints `lb <= Ax <= ub` for some matrix `A`, set `par_constr[[1]] = NULL`, `par_constr[[2]] = A` and `par_constr[[3]] = rbind(lb, ub)`.

For additional examples, see Examples 1 and 2 below. The package currently supports only two solvers, one from each of the packages `lpSolve` and `Rmosek`.

The output of `perturbCS` contains:
  + `Psihat_lb`: values of the two sample minimization programs for lower bound of identified set;
  + `Psihat_ub`: values of the two sample maximization programs for upper bound of identified set;
  + `PsiCS`: lower and upper confidence bounds for size(s) `alpha`;
  + `threshold`: `TRUE` if finite-sample thresholding activated, and `FALSE` otherwise;
  + `nfeasible`: number of bootstrap samples with nonempty estimated bounds.

## Examples

### Example 1: Missing Data Example

See Appendix B.1 of Cho and Russell (2020, <arXiv:1810.03180>). The DGP and the model can be declared as follows:

``` r
# data generating process
gen_data <- function(n, c) {
  yd = cbind(rep(c(1,2,3,4,5),2), rep(c(0,1), each=5))
  prob_d0 = rep(0.2*c/sqrt(n), 5)
  prob_d1 = rep(0.2*(1-c/sqrt(n)), 5)
  prob = c(prob_d0, prob_d1)
  ind = sample(1:10, n, replace=TRUE, prob=prob)
  yd = yd[ind,]
  data = cbind(yd[,1]*yd[,2], yd[,2])
  return(data)
}

# psi_fn
psi_fn <- function(data) {
  psi = c(0,1,2,3,4,5,1,2,3,4,5)
  return(psi)
}

# mom_constr
mom_constr = list("matfn"=NULL, "bfn"=NULL)

mom_constr$matfn <- function(data) {
  moment_mat = rep(c(1,0), each=5)
  id_mat = diag(5)
  for (y in 1:5) {
    moment_mat = rbind(moment_mat, c(rep(0,5),id_mat[y,]))
  }
  return(moment_mat)
}

mom_constr$bfn <- function(data) {
  moment_b = mean(data[,1]==0 & data[,2]==0)
  for (y in 1:5) {
    moment_b = c(moment_b, mean(data[,1]==y & data[,2]==1))
  }
  moment_b = rbind(moment_b,moment_b)
  return(moment_b)
}

# par_constr
par_constr = list("indivb"=NULL, "jointmat"=NULL, "jointb"=NULL)
par_constr$indivb = rbind(rep(0,10), rep(1,10))
par_constr$jointmat = matrix(rep(1,10), 1, 10)
par_constr$jointb = rbind(1, 1)
```

To construct a confidence set for the average outcome with a random sample of size `n=100` from population `c=2`:

``` r
library(lmcperturb)
n <- 100
c <- 2
data <- gen_data(n, c)
perturbCS(data, psi_fn, mom_constr, par_constr, method="lpSolve") # using lpSolve
perturbCS(data, psi_fn, mom_constr, par_constr, method="Rmosek") # using Rmosek
```

### Example 2: Counterfactual Policy Example

See Appendix B.3 of Cho and Russell (2020, <arXiv:1810.03180>). The DGP and the model can be declared as follows:

``` r
# preset vectors (used to define objective function and moment constraints)
y0order = rep(rep(1:5,each=10),8)
y1order = rep(rep(rep(1:5,each=2),5),8)
dorder = rep(0:1,200)
xorder = rep(1:4,each=100)
zorder = rep(rep(0:1,each=50),4)
horder = rep(c(-0.5,0.25,0.25,0.5), each=100)

# data generating process
gen_data <- function(n, c) {
  z = sample(0:1, n, replace=TRUE)
  x_old = sample(1:4, n, prob=c(0.25,0.25,0.25,0.25), replace=TRUE)
  x_new = sample(1:4, n, prob=c(0.30,0.30,0.20,0.20), replace=TRUE)
  y0prob = rbind(c(0.20,0.20,0.20,0.20,0.20),
                 c(0.30,0.25,0.25,0.10,0.10),
                 c(0.30,0.25,0.25,0.10,0.10),
                 c(0.40,0.35,0.25,0.00,0.00))
  y1prob = y0prob[,5:1]
  y0 = c()
  y1 = c()
  for (x in 1:4) {
    ind = x_old==x
    nx = sum(ind)
    y0[ind] = sample(1:5, nx, prob=y0prob[x,], replace=TRUE)
    y1[ind] = sample(1:5, nx, prob=y1prob[x,], replace=TRUE)
  }
  d = as.numeric(2*z-1 > c/sqrt(n)*rnorm(n))
  y = y0*(1-d) + y1*d
  data = cbind(y,d,x_old,z,x_new)
  return(data)
}

# psi_fn
psi_fn <- function(data) {
  pxorder = rep(sapply(1:4, function(x) mean(x==data[,5])), each=100)
  pzorder = rep(rep(sapply(0:1, function(z) mean(z==data[,4])), each=50), 4)
  psi = c(0, (y1order-y0order)*horder*pxorder*pzorder)
  return(psi)
}

# mom_constr
mom_constr = list("matfn"=NULL, "bfn"=NULL)

mom_constr$matfn <- function(data) {
  moment_mat = c()
  for (x in 1:4) {
    for (z in 0:1) {
      dataxz = (data[,3]==x)*(data[,4]==z)
      xz = (xorder==x)*(zorder==z)*mean(dataxz)
      for (y in 1:5) {
        moment_mat = rbind(moment_mat, (y0order==y)*(dorder==0)*xz, (y1order==y)*(dorder==1)*xz)
      }
    }
  }
  return(moment_mat)
}

mom_constr$bfn <- function(data) {
  moment_b = c()
  for (x in 1:4) {
    for (z in 0:1) {
      dataxz = (data[,3]==x)*(data[,4]==z)
      for (y in 1:5) {
        xzy = (data[,1]==y)*dataxz
        moment_b = c(moment_b, mean(xzy*(data[,2]==0)), mean(xzy*(data[,2]==1)))
      }
    }
  }
  moment_b = rbind(moment_b, moment_b)
  return(moment_b)
}

# par_constr
par_constr = list("indivb"=NULL, "jointmat"=NULL, "jointb"=NULL)
par_constr$indivb = rbind(rep(0,400), rep(1,400))
par_constr$jointmat = t(kronecker(diag(1,8),rep(1,50)))
for (x in 1:4) {
  for (y0 in 1:5) {
    for (y1 in 1:5) {
      indz = (y0order==y0)*(y1order==y1)*(xorder==x)
      indz0 = indz*(zorder==0)
      indz1 = indz*(zorder==1)
      par_constr$jointmat = rbind(par_constr$jointmat, indz0-indz1) # exogeneity constraints P(Y0=y0,Y1=y1|X=x,Z=0) = P(Y0=y0,Y1=y1|X=x,Z=1)
    }
  }
}
par_constr$jointb = rbind(c(rep(1,8), rep(0,100)), c(rep(1,8), rep(0,100)))
```

To construct a confidence set for the average outcome difference between competing policies A and B with a random sample of size `n=1000` from population `c=1`:

``` r
library(lmcperturb)
n <- 1000
c <- 1
data <- gen_data(n, c)
perturbCS(data, psi_fn, mom_constr, par_constr, method="Rmosek")
```
