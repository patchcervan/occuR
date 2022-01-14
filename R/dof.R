#' Compute (estimated) degrees of freedom for fitted occupancy model
#'
#' @param obj fitted model object
#' @param each default is FALSE, if TRUE then list of degrees of freedom is returned, one fixed parameters, smooth
#' for psi, and smooth for p. When FALSE, these are summed and total estimated degrees of freedom returned.
#'
#' @return degrees of freedom (see "each" above)
#' @export
dof.occuR <- function(obj, each = FALSE) {
    fix <- obj$res$par.fixed
    df <- length(fix[!grepl("log_lambda", names(fix))])
    psi_lambda <- exp(fix[names(fix) == "log_lambda_psi"])
    psi_edf <- edf.occuR(obj$mats$X_psi[, -(1:obj$mats$nfix_psi)], obj$mats$S_psi$S, psi_lambda, obj$mats$S_psi$Scols)
    p_lambda <- exp(fix[names(fix) == "log_lambda_p"])
    p_edf <- edf.occuR(obj$mats$X_p[, -(1:obj$mats$nfix_p)], obj$mats$S_p$S, p_lambda, obj$mats$S_p$Scols)
    if (each) return(list(fix = df, psi = psi_edf, p = p_edf))
    return(df + psi_edf + p_edf)
}

#' Compute estimated degrees of freedom for a spline
#'
#' @param X design matrix for spline
#' @param S smoothing matrix for spline
#' @param lambda smoothiing parameters
#' @param Sn number of parameters per spline
#'
#' @description This code was written by David L Miller as part of the CTMCdive R package
#' and is based on Wood et al. (2017) p 212.
#'
#' @return estimated degrees of freedom
edf.occuR <- function(X, S, lambda, Sn){
    if (length(lambda) == 0) return(0)

    # duplicate lambda enough times
    lambda <- rep(lambda, Sn)

    # calculate lambda*S
    Sbig <- S * lambda

    # calculate the hat matrix
    XtX <- t(X) %*% X
    Fi <- solve(XtX + Sbig)
    F <- Fi %*% XtX

    # return the trace
    return(sum(diag(F)))
}


