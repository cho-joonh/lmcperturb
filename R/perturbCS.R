#' Construct confidence set for functional of partially identified
#' parameter
#'
#' Construct a confidence set for linear functional of partially identified
#' parameter defined by linear moment inequalities. The confidence set is
#' constructed by bootstrapping perturbed linear programs that define the lower
#' and upper bounds of the identified set for the functional of interest.
#'
#' @param data A data matrix.
#' @param psi_fn A function that returns a numeric vector of intercept and
#'     coefficients of the linear functional.
#' @param mom_constr A list with two elements, first element is a function that
#'     returns the matrix of coefficients for linear moment constraints and the
#'     second element is a function that returns the lower and upper bounds of
#'     the constraints as a matrix with two rows; if no lower or upper bound
#'     for a moment constraint, set \code{-Inf} or \code{Inf} for the
#'     appropriate element of the bound matrix.
#' @param par_constr A list with three elements, first element is a matrix with
#'     two rows for the lower and upper bounds on the constraints for individual
#'     coordinates of the parameter space, second element is a matrix with
#'     coefficients that define the planes that restrict the parameter space,
#'     third element is a matrix with two rows for the lower and upper bounds
#'     of the linear constraints defined by the planes; if no lower or upper
#'     bound for a moment constraint, set \code{-Inf} or \code{Inf} for the
#'     appropriate element of the bound matrix.
#' @param max_perturb A list with two elements, first element is a numeric
#'     vector with two values corresponding to the lower and upper bounds for
#'     perturbations to the functional of interest, second element is a numeric
#'     vector with two values corresponding to the lower and upper bounds for
#'     perturbations to the linear moment constraints; default set to
#'     \code{max_perturb = list("psi"=c(-1e-6,1e-6),"constr"=c(0,1e-6))}.
#' @param alpha A numeric vector of values for desired sizes of the confidence
#'     set; default set to \code{alpha = c(.10,.05,.01)}.
#' @param tol A numeric value for tolerance level to threshold quantile
#'     selection of bootstrap critical values in small samples; default set to
#'     \code{tol = 1}.
#' @param boot A numeric value for the bootstrap sample size; default set to
#'     \code{boot = 999}.
#' @param method A string for the solver to be used; currently only
#'     \code{method = "lpSolve"} and \code{"Rmosek"} are supported; default set
#'     to \code{method = "lpSolve"}.
#'
#' @return A bootstrap confidence set for functional of interest with specified
#'     size. Additionally, estimated lower and upper bounds of the identified
#'     set are reported. Logical output for whether the finite-sample threshold
#'     is used and the number of valid bootstrap samples with nonempty feasible
#'     regions for the perturbed bootstrap linear programs are also reported.
#' @export
perturbCS <- function(data, psi_fn, mom_constr, par_constr=NULL, max_perturb=list("psi"=c(-1e-6,1e-6),"constr"=c(0,1e-6)), alpha=c(.10,.05,.01), tol=1, boot=999, method="lpSolve") {
  check_type(data, psi_fn, mom_constr, par_constr, max_perturb, alpha, tol, boot, method)
  check_dimension(data, psi_fn, mom_constr, par_constr)
  names(mom_constr) <- c("matfn", "bfn")
  if (!is.null(par_constr)) {
    names(par_constr) <- c("indivb", "jointmat", "jointb")
  }
  names(max_perturb) <- c("psi", "constr")
  if (method=="Rmosek") {
    if (!requireNamespace("Rmosek", quietly=TRUE)) {
      stop("Package 'Rmosek' needed to use method 'Rmosek'. Follow installation instructions in Rmosek documentation provided on MOSEK homepage or use method 'lpSolve'.", call.=FALSE)
    }
    pconstr = par_constr_Rmosek(data, mom_constr, par_constr)
    perturb = gen_perturb_Rmosek(data, mom_constr, pconstr, max_perturb)
    perturbLP = perturbLP_Rmosek
  } else {
    if (method != "lpSolve") {
      message(paste0("Method '", method, "' not supported; using method 'lpSolve' instead."))
    }
    pconstr = par_constr_lpSolve(data, mom_constr, par_constr)
    perturb = gen_perturb_lpSolve(data, mom_constr, pconstr, max_perturb)
    perturbLP = perturbLP_lpSolve
  }
  n = nrow(data)
  Psihat = perturbLP(data, psi_fn, mom_constr, pconstr, perturb)
  Psihat_lb = min(Psihat[1:2])
  Psihat_ub = max(Psihat[3:4])
  if (!is.na(Psihat_lb)) {
    Psiboot = c()
    for (b in 1:boot) {
      ind = sample(1:n, n, replace=TRUE)
      Psiboot = cbind(Psiboot, perturbLP(data[ind,], psi_fn, mom_constr, pconstr, perturb))
    }
    feasind = !is.na(Psiboot[1,])
    nfeasible = sum(feasind)
    if (nfeasible/boot >= 0.01) {
      Psiboot = Psiboot[,feasind]
      threshold = tol/sqrt(log(n)) > Psihat_ub-Psihat_lb
      talpha = ifelse(rep(threshold,length(alpha)), alpha/2, alpha)
      PsiCS = c()
      for (a in talpha) {
        meboot = apply(c(1,1,-1,-1)*(Psiboot-Psihat), 1, sort)[ceiling((1-a)*nfeasible),1:4]
        PsiCS = rbind(PsiCS, c(Psihat_lb-min(meboot[1:2]), Psihat_ub-max(-meboot[3:4])))
      }
      PsiCS = cbind(alpha, PsiCS)
    } else {
      threshold = NA
      PsiCS = cbind(alpha, NA, NA)
      warning(paste("More than 99% of the bootstrap samples result in empty estimated set. Only", nfeasible, "of", boot, "bootstrap samples feasible."))
    }
    return(list("Psihat_lb"=sort(Psihat[1:2]), "Psihat_ub"=sort(Psihat[3:4]), "PsiCS"=PsiCS, "threshold"=threshold, "nfeasible"=nfeasible))
  } else {
    warning("Estimated identified set is empty.")
    return(list("Psihat_lb"=c(NA,NA), "Psihat_ub"=c(NA,NA), "PsiCS"=cbind(alpha,NA,NA), "threshold"=NA, "nfeasible"=NA))
  }
}
