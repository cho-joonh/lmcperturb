#' Create a list defining the parameter space tailored to Rmosek solver
#'
#' @param data A data matrix.
#' @param mom_constr A list with two functions.
#' @param par_constr A list with three matrices.
#'
#' @return Returns linear constraint matrix, constraint upper and lower bounds
#'     and parameter space constraints reshaped for use in Rmosek.
par_constr_Rmosek <- function(data, mom_constr, par_constr) {
  dtheta = ncol(mom_constr$matfn(data))
  if (is.null(par_constr)) {
    par_b = rbind(rep(-1e+7,dtheta), rep(1e+7,dtheta))
  } else {
    if (is.null(par_constr$indivb)) {
      par_b = rbind(rep(-1e+7,dtheta), rep(1e+7,dtheta))
    } else {
      par_b = par_constr$indivb
      par_b[par_b == -Inf] = -1e+7
      par_b[par_b == Inf] = 1e+7
    }
    if (is.null(par_constr$jointmat)) {
      joint_mat = NULL
      joint_b = NULL
    } else {
      joint_mat = par_constr$jointmat
      joint_b = par_constr$jointb
    }
  }
  return(list("A"=joint_mat, "bc"=joint_b, "bx"=par_b))
}

#' Create a list defining the moment constraints tailored to Rmosek solver
#'
#' @param data A data matrix.
#' @param mom_constr A list with two functions.
#'
#' @return Returns linear constraint matrix and  constraint upper and lower
#'     bounds defined by moment inequalities reshaped for use in Rmosek.
mom_constr_Rmosek <- function(data, mom_constr) {
  moment_mat = mom_constr$matfn(data)
  moment_b = mom_constr$bfn(data)
  return(list("A"=moment_mat, "bc"=moment_b))
}

#' Create a list of perturbation draws tailored to Rmosek solver
#'
#' @param data A data matrix.
#' @param mom_constr A list with two functions.
#' @param pconstr A list with three matrices.
#' @param max_perturb A list with two vectors.
#'
#' @return Returns perturbation draws for functional of interest (objective
#'     function of linear programs) and moment constraints (feasible region).
gen_perturb_Rmosek <- function(data, mom_constr, pconstr, max_perturb) {
  mconstr = mom_constr_Rmosek(data, mom_constr)
  dtheta = ncol(mconstr$A)
  k = ncol(mconstr$bc) + ifelse(is.null(pconstr$bc),0,ncol(pconstr$bc))
  perturb_psi = stats::runif(dtheta, max_perturb$psi[1], max_perturb$psi[2])
  perturb_bc = stats::runif(2*k, max_perturb$constr[1], max_perturb$constr[2])
  perturb_bc = matrix(perturb_bc, nrow=2)*c(-1,1)
  perturb_bx = stats::runif(2*dtheta, max_perturb$constr[1], max_perturb$constr[2])
  perturb_bx = matrix(perturb_bx, nrow=2)*c(-1,1)
  return(list("psi"=perturb_psi, "bc"=perturb_bc, "bx"=perturb_bx))
}

#' Solve perturbed linear programs using Rmosek solver
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
perturbLP_Rmosek <- function(data, psi_fn, mom_constr, pconstr, perturb) {
  psi = psi_fn(data)
  psi0 = psi[1]
  psim = psi[-1] - perturb$psi
  psip = psi[-1] + perturb$psi
  mconstr = mom_constr_Rmosek(data, mom_constr)
  A = rbind(mconstr$A, pconstr$A)
  bc = cbind(mconstr$bc, pconstr$bc) + perturb$bc
  bx = pconstr$bx + perturb$bx
  minm = Rmosek::mosek(list(sense="min", c=psim, A=A, bc=bc, bx=bx), list(verbose=0, soldetail=1))
  minp = Rmosek::mosek(list(sense="min", c=psip, A=A, bc=bc, bx=bx), list(verbose=0, soldetail=1))
  maxm = Rmosek::mosek(list(sense="max", c=psim, A=A, bc=bc, bx=bx), list(verbose=0, soldetail=1))
  maxp = Rmosek::mosek(list(sense="max", c=psip, A=A, bc=bc, bx=bx), list(verbose=0, soldetail=1))
  if (minm$sol$bas$prosta == "PRIMAL_INFEASIBLE") {
    objval = rep(NA, 4)
    return(objval)
  } else {
    objval = psi0 + c(minm$sol$bas$pobjval, minp$sol$bas$pobjval, maxm$sol$bas$pobjval, maxp$sol$bas$pobjval)
    return(objval)
  }
}
