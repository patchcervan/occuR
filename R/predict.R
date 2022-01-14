#' Predict from fitted occupancy model
#'
#' @param obj fitted model object from fit_occu
#' @param visit_data see fit_occu
#' @param site_data see fit_occu
#' @param nboot number of parametric bootstrap resamples to produce from fitted model
#'
#' @return if nboot = 0 (default), the list of fitted psi and p values; otherwise, list also
#' contains matrix for psi and p where each row is a bootstrap resample
#' @export
#' @importFrom mgcv rmvn
predict.occuR <- function(obj, visit_data, site_data, nboot = 0) {
    mats <- obj$mats
    site_data$psi <- 1:nrow(site_data)
    mats$X_psi <- predict(mats$gam_psi, newdata = site_data, type = "lpmatrix")
    site_data <- site_data[, psi := NULL]
    visit_data$p <- 1:nrow(visit_data)
    mats$X_p <- predict(mats$gam_p, newdata = visit_data, type = "lpmatrix")
    visit_data <- visit_data[, p := NULL]
    fix <- obj$res$par.fixed
    ran <- obj$res$par.random
    pred <- get_predicted_values(fix, ran, mats)
    if (nboot > 0.5) {
        Q <- obj$res$jointPrecision
        Q <- Q[!grepl("log_lambda_", colnames(Q)),
               !grepl("log_lambda_", colnames(Q)), drop = FALSE]
        V <- solve(Q)
        param <- c(fix, ran)
        param <- param[!grepl("log_lambda", names(param))]
        boots <- rmvn(nboot, param, V)
        colnames(boots) <- names(param)
        nfix <- length(fix[!grepl("log_lambda", names(fix))])
        boots_psi <- matrix(0, nr = nboot, nc = length(pred$psi))
        boots_p <- matrix(0, nr = nboot, nc = length(pred$p))
        for (b in 1:nboot) {
            boot_res <- get_predicted_values(boots[b, 1:nfix], boots[b, -(1:nfix)], mats)
            boots_psi[b,] <- boot_res$psi
            boots_p[b,] <- boot_res$p
        }
        pred$psiboot <- boots_psi
        pred$pboot <- boots_p
    }
    return(pred)
}

#' Get predicted values for psi and p from parameters and matrices
#'
#' @param fix fixed parameter vector
#' @param ran random effects parameter vector
#' @param mats matrices from make_matrices
#'
#' @return list of two vectors: predicted values for psi and p
get_predicted_values <- function(fix, ran, mats) {
    nms_fix <- names(fix)
    nms_ran <- names(ran)
    psi_pred <- (mats$X_psi %*% c(fix[nms_fix == "beta_psi"], ran[nms_ran == "z_psi"]))
    p_pred <- (mats$X_p %*% c(fix[nms_fix == "beta_p"], ran[nms_ran == "z_p"]))
    psi_pred <- plogis(psi_pred)
    p_pred <- plogis(p_pred)
    return(list(psi = psi_pred, p = p_pred))
}
