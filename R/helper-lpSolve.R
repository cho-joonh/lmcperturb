#' Create a list defining the parameter space tailored to lpSolve solver
#'
#' @param data A data matrix.
#' @param mom_constr A list with two functions.
#' @param par_constr A list with three matrices.
#'
#' @return Returns linear constraint matrix, constraint upper and lower bounds
#'     and parameter space constraints reshaped for use in lpSolve.
par_constr_lpSolve <- function(data, mom_constr, par_constr) {
  dtheta = ncol(mom_constr$matfn(data))
  par_mat = rbind(diag(dtheta), -diag(dtheta))
  if (is.null(par_constr)) {
    par_ub = rep(1e+7, 2*dtheta)
  } else {
    if (is.null(par_constr$indivb)) {
      par_ub = rep(1e+7, 2*dtheta)
    } else {
      par_ub = c(par_constr$indivb[2,], -par_constr$indivb[1,])
      par_ub[par_ub == Inf] = 1e+7
    }
    if (!is.null(par_constr$jointmat)) {
      ind_ub = par_constr$jointb[2,] != Inf
      ind_lb = par_constr$jointb[1,] != -Inf
      par_mat = rbind(par_mat, par_constr$jointmat[ind_ub,], -par_constr$jointmat[ind_lb,])
      par_ub = c(par_ub, par_constr$jointb[2,ind_ub], -par_constr$jointb[1,ind_lb])
    }
  }
  return(list("lhs"=par_mat, "rhs"=par_ub))
}

#' Create a list defining the moment constraints tailored to lpSolve solver
#'
#' @param data A data matrix.
#' @param mom_constr A list with two functions.
#'
#' @return Returns linear constraint matrix and  constraint upper and lower
#'     bounds defined by moment inequalities reshaped for use in lpSolve.
mom_constr_lpSolve <- function(data, mom_constr) {
  moment_mat = mom_constr$matfn(data)
  moment_ub = mom_constr$bfn(data)
  ind = moment_ub[1,] != -Inf
  moment_mat = rbind(moment_mat, -moment_mat[ind,])
  moment_ub = c(moment_ub[2,], -moment_ub[1,ind])
  return(list("lhs"=moment_mat, "rhs"=moment_ub))
}

#' Create a list of perturbation draws tailored to lpSolve solver
#'
#' @param data A data matrix.
#' @param mom_constr A list with two functions.
#' @param pconstr A list with three matrices.
#' @param max_perturb A list with two vectors.
#'
#' @return Returns perturbation draws for functional of interest (objective
#'     function of linear programs) and moment constraints (feasible region).
gen_perturb_lpSolve <- function(data, mom_constr, pconstr, max_perturb) {
  mconstr = mom_constr_lpSolve(data, mom_constr)
  dtheta = ncol(pconstr$lhs)
  k = length(pconstr$rhs) + length(mconstr$rhs)
  perturb_psi = stats::runif(dtheta, max_perturb$psi[1], max_perturb$psi[2])
  perturb_constr = stats::runif(k, max_perturb$constr[1], max_perturb$constr[2])
  return(list("psi"=perturb_psi, "constr"=perturb_constr))
}

#' Solve perturbed linear programs using lpSolve solver
#'
#' @param data A data matrix.
#' @param psi_fn A function.
#' @param mom_constr A list with two functions.
#' @param pconstr A list with three matrices.
#' @param perturb A list with a vector and two matrices.
#'
#' @return Returns values of four perturbed linear programs for a given
#'     perturbation draw \code{perturb}. If feasible region empty, a vector
#'     of \code{NA}s is returned.
perturbLP_lpSolve <- function(data, psi_fn, mom_constr, pconstr, perturb) {
  psi = psi_fn(data)
  psi0 = psi[1]
  psim = psi[-1] - perturb$psi
  psip = psi[-1] + perturb$psi
  mconstr = mom_constr_lpSolve(data, mom_constr)
  lhs = rbind(mconstr$lhs, pconstr$lhs)
  rhs = c(mconstr$rhs, pconstr$rhs) + perturb$constr
  minm = lpSolve::lp("min", c(psim,-psim), cbind(lhs,-lhs), rep("<=",length(rhs)), rhs)
  minp = lpSolve::lp("min", c(psip,-psip), cbind(lhs,-lhs), rep("<=",length(rhs)), rhs)
  maxm = lpSolve::lp("max", c(psim,-psim), cbind(lhs,-lhs), rep("<=",length(rhs)), rhs)
  maxp = lpSolve::lp("max", c(psip,-psip), cbind(lhs,-lhs), rep("<=",length(rhs)), rhs)
  if (minm$status != 0) {
    objval = rep(NA, 4)
    return(objval)
  } else {
    objval = psi0 + c(minm$objval, minp$objval, maxm$objval, maxp$objval)
    return(objval)
  }
}
