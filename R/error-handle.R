#' Checks classes of input objects for perturbCS
#'
#' @param data A data matrix.
#' @param psi_fn A function.
#' @param mom_constr A list with two functions.
#' @param par_constr A list with three matrices where first and third matrices
#'     must have two rows.
#' @param max_perturb A list with two vectors where the second vector must be
#'     nonnegative.
#' @param alpha A numeric vector with element(s) between 0 and 1 (exclusive).
#' @param tol A nonnegative number.
#' @param boot A positive integer.
#' @param method A string.
#'
#' @return Returns an error message is an input object is of invalid class.
check_type <- function(data, psi_fn, mom_constr, par_constr, max_perturb, alpha, tol, boot, method) {
  if (!is.matrix(data)) {
    stop("data must be a matrix.", call.=FALSE)
  }
  if (!is.function(psi_fn)) {
    stop("psi_fn must be a function.", call.=FALSE)
  }
  if (!is.list(mom_constr)) {
    stop("mom_constr must be a list with two elements.", call.=FALSE)
  } else {
    if (length(mom_constr) != 2) {
      stop("mom_constr must be a list with two elements.", call.=FALSE)
    } else if (!is.function(mom_constr[[1]]) || !is.function(mom_constr[[2]])) {
      stop("Both elements of mom_constr must be functions.", call.=FALSE)
    }
  }
  if (!is.null(par_constr) && !is.list(par_constr)) {
    stop("par_constr must either be NULL or a list.", call.=FALSE)
  } else if (is.list(par_constr)) {
    if (length(par_constr) != 3) {
      stop("par_constr must be a list with three elements.", call.=FALSE)
    } else {
      if (!is.null(par_constr[[1]]) && !is.matrix(par_constr[[1]])) {
        stop("par_constr[[1]] must be either NULL or a numeric matrix with two rows.", call.=FALSE)
      } else if (is.matrix(par_constr[[1]])) {
        if (nrow(par_constr[[1]]) != 2 || !is.numeric(par_constr[[1]])) {
          stop("par_constr[[1]] must be a numeric matrix with two rows.", call.=FALSE)
        } else if (any(par_constr[[1]][1,] > par_constr[[1]][2,])) {
          stop("Elements in the first row of par_constr[[1]] cannot exceed those in the second row for each column.", call.=FALSE)
        }
      }
      if (is.null(par_constr[[2]])) {
        if (!is.null(par_constr[[3]])) {
          stop("Both par_constr[[2]] and par_constr[[3]] must be NULL or numeric matrices.", call.=FALSE)
        }
      } else {
        if (is.null(par_constr[[3]])) {
          stop("par_constr[[3]] must be a numeric matrix with two rows.", call.=FALSE)
        } else if (!is.matrix(par_constr[[2]]) || !is.matrix(par_constr[[3]])) {
          stop("Both par_constr[[2]] and par_constr[[3]] must be NULL or numeric matrices.", call.=FALSE)
        } else if (!is.numeric(par_constr[[2]]) || !is.numeric(par_constr[[3]])) {
          stop("Both par_constr[[2]] and par_constr[[3]] must be NULL or numeric matrices.", call.=FALSE)
        } else if (nrow(par_constr[[2]]) != ncol(par_constr[[3]])) {
          stop("Number of columns of par_constr[[2]] and number of rows of par_constr[[3]] must match.", call.=FALSE)
        } else if (nrow(par_constr[[3]]) != 2) {
          stop("par_constr[[3]] must be a numeric matrix with two rows.", call.=FALSE)
        } else if (any(par_constr[[3]][1,] > par_constr[[3]][2,])) {
          stop("Elements in the first row of par_constr[[3]] cannot exceed those in the second row for each column.", call.=FALSE)
        }
      }
      if (is.matrix(par_constr[[1]]) && is.matrix(par_constr[[2]])) {
        if (ncol(par_constr[[1]]) != ncol(par_constr[[2]])) {
          stop("Number of columns of par_constr[[1]] and par_constr[[2]] must match.", call.=FALSE)
        }
      }
    }
  }
  if (!is.list(max_perturb)) {
    stop("max_perturb must be a list with two elements.", call.=FALSE)
  } else if (length(max_perturb) != 2) {
    stop("max_perturb must be a list with two elements.", call.=FALSE)
  } else if (length(max_perturb[[1]]) != 2 || length(max_perturb[[2]]) != 2) {
    stop("Both elements of max_perturb must be numeric vector with length 2.", call.=FALSE)
  } else if (!is.numeric(max_perturb[[1]]) || !is.numeric(max_perturb[[2]])) {
    stop("Both elements of max_perturb must be numeric vector with length 2.", call.=FALSE)
  } else if (max_perturb[[1]][1] >= max_perturb[[1]][2]) {
    stop("max_perturb[[1]][1] must be smaller than max_perturb[[1]][2].", call.=FALSE)
  } else if (max_perturb[[2]][2] <= max_perturb[[2]][1] || max_perturb[[2]][1] < 0) {
    stop("max_perturb[[2]][1] must be smaller than max_perturb[[2]][2]. Both must be nonnegative.", call.=FALSE)
  }
  if (!is.numeric(alpha)) {
    stop("alpha must be numeric.", call.=FALSE)
  } else if (any(alpha <= 0 | alpha >= 1)) {
    stop("alpha must be between 0 and 1 (exclusive).", call.=FALSE)
  }
  if (!is.numeric(tol)) {
    stop("tol must be numeric.", call.=FALSE)
  } else if (length(tol) != 1 || tol < 0) {
    stop("tol must be a nonnegative number.", call.=FALSE)
  }
  if (!is.numeric(boot)) {
    stop("boot must be numeric.", call.=FALSE)
  } else if (length(boot) != 1 || boot < 0 || abs(boot-round(boot)) > 1e-16) {
    stop("boot must be a nonnegative integer.", call.=FALSE)
  }
  if (!is.character(method)) {
    stop("method must be a string; either 'lpSolve' or 'Rmosek'.", call.=FALSE)
  } else if (length(method) != 1) {
    stop("method must be a string.", call.=FALSE)
  }
}

#' Checks the dimensions of input objects for perturbCS
#'
#' @param data A data matrix.
#' @param psi_fn A function.
#' @param mom_constr A list with two functions.
#' @param par_constr A list with three matrices where first and third matrices
#'     must have two rows.
#'
#' @return Returns an error message if dimensions of inputs do not align.
check_dimension <- function(data, psi_fn, mom_constr, par_constr) {
  if(class(try(psi_fn(data), silent=TRUE))[1] == "try-error") {
    stop("Error in psi_fn(data). Check function psi_fn.", call.=FALSE)
  } else if (class(try(mom_constr[[1]](data), silent=TRUE))[1] == "try-error") {
    stop("Error in mom_constr[[1]](data). Check function mom_constr[[1]].", call.=FALSE)
  } else if (class(try(mom_constr[[2]](data), silent=TRUE))[1] == "try-error") {
    stop("Error in mom_constr[[2]](data). Check function mom_constr[[2]].", call.=FALSE)
  } else {
    psi <- psi_fn(data)
    moment_mat <- mom_constr[[1]](data)
    moment_b <- mom_constr[[2]](data)
    if (!is.null(dim(psi)) || !is.numeric(psi)) {
      stop("psi_fn(data) must be a numeric vector with length equal to the length of the parameter vector plus one.", call.=FALSE)
    } else if (!is.matrix(moment_mat) || !is.numeric(moment_mat)) {
      stop("mom_constr[[1]](data) must be a matrix with numeric elements.", call.=FALSE)
    } else if (!is.matrix(moment_b) || !is.numeric(moment_b)) {
      stop("mom_constr[[2]](data) must be a matrix with numeric elements.", call.=FALSE)
    } else if (nrow(moment_b) != 2) {
      stop("mom_constr[[2]](data) must be a matrix with two rows.", call.=FALSE)
    } else if (length(psi) != ncol(moment_mat) + 1) {
      stop("psi_fn(data) must have one more element than the number of columns of mom_constr[[1]](data).", call.=FALSE)
    } else if (nrow(moment_mat) != ncol(moment_b)) {
      stop("Number of columns of mom_constr[[1]](data) and number of rows of mom_constr[[2]](data) must match.", call.=FALSE)
    } else if (!is.null(par_constr)) {
      if (length(psi) != ncol(par_constr[[1]]) + 1) {
        stop("Number of columns of mom_constr[[1]](data) and number of columns of par_constr[[1]] must match.", call.=FALSE)
      }
    }
  }
}
