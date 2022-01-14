#' Compute design and smoothing matrices from formulae
#'
#' @param forms list of formulae for psi and p
#' @param site_data see fit_occu function
#' @param visit_data see fit_occu function
#'
#' @return list of X = fixed effects design matrix, Z = random effects design matrix,
#'   and S = smoothing matrix for both psi and p (e.g. X_psi).
#'   make_smoothing_matrix() returns a list of S (block diagonal smoothing matrix),
#'   Scols (number of columns for each smoothing matrix),
#'   bdim (k argument used for each smooth), Snames (names of smoothing parameters).
#'   Returns NULL if no smooths in model.
#'   Based on code written by David L Miller in CTMCdive R Package (https://github.com/r-glennie/CTMCdive).
#' @export
#'
#' @examples
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom methods as
#' @importFrom Matrix bdiag
make_matrices <- function(forms, visit_data, site_data) {
    ## occupancy model
    site_data$psi <- 1:nrow(site_data)
    gam_psi <- gam(forms[["psi"]], data = site_data, method = "REML")
    site_data <- site_data[, psi := NULL]
    # get design matrix
    Xfull <- predict(gam_psi, newdata = site_data, type = "lpmatrix")
    X_psi <- uniquecombs(Xfull)
    # get smoothing matrix
    S_psi <- make_smoothing_matrix(gam_psi)

    ## detection model
    visit_data$p <- 1:nrow(visit_data)
    gam_p <- gam(forms[["p"]], data = visit_data, method = "REML")
    visit_data <- visit_data[, p := NULL]
    # get design matrix
    Xfull <- predict(gam_p, newdata = visit_data, type = "lpmatrix")
    X_p <- uniquecombs(Xfull)
    # get smoothing matrix
    S_p <- make_smoothing_matrix(gam_p)

    ## results
    res <- list(X_psi = X_psi,
                S_psi = S_psi,
                nfix_psi = gam_psi$nsdf,
                X_p = X_p,
                S_p = S_p,
                nfix_p = gam_p$nsdf,
                gam_p = gam_p,
                gam_psi = gam_psi)

    return(res)
}

# Compute block-diagonal smoothing matrix from value returned by
# mgcv::gam with fit = FALSE (gm)

make_smoothing_matrix <- function(gm) {
    if (length(gm$smooth) < 0.5) return(NULL)
    S <- list()
    Sncols <- list()
    bdim <- list()
    smnames <- list()
    for (i in seq_along(gm$smooth)) {
        sm <- gm$smooth[[i]]
        for (j in seq_along(sm$S)) {
            S <- c(S, as(sm$S[[j]], "sparseMatrix"))
            Sncols <- c(Sncols, ncol(sm$S[[j]]))
            bdim <- c(bdim, sm$bs.dim)
            smnames <- c(smnames, attr(sm$sp, "names"))
        }
    }
    Smat <- bdiag(S)
    res <- list(S = Smat, Scols = Sncols, bdim = bdim, Snames = smnames)
    return(res)
}
