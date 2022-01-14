#' Fit occupancy model
#'
#' @param forms list of formulae for psi (occupancy) and p (detection)
#' @param visit_data data.table with row for each visit with a "obs" column for detection record
#' @param site_data data.table with row for each site x occasion
#' @param start named list of starting values for beta_psi and beta_p - fixed effects parameters
#' @param print if TRUE then possibly useful info is printed out
#' @return fitted occupied model object
#' @export
#' @importFrom data.table uniqueN
#' @importFrom TMB MakeADFun sdreport normalize
#' @useDynLib occu_tmb
fit_occu <- function(forms, visit_data, site_data, start = NULL, print = TRUE) {
    ## DATA
    # order data
    site_data <- site_data[order(site_data$site, site_data$occasion),]
    visit_data <- visit_data[order(visit_data$site, visit_data$occasion, visit_data$visit),]
    # check data
    check_inp(forms, visit_data, site_data, start, print)
    # number of sites
    nsites <- uniqueN(site_data$site)
    # number of occasions
    nocc <- site_data[, .(n = .N), .(site)]$n
    # total number of detections per site x occasion
    nvis <- visit_data[, .(totsite = sum(obs), nvisit = .N), .(site, occasion)]
    totsite <- nvis$totsite
    nvisit <- nvis$nvisit

    ## MODEL MATRICES
    # name formulae
    names(forms) <- c(as.character(forms[[1]][[2]]),
                      as.character(forms[[2]][[2]]))
    # get model matrices and smoothing matrices
    mats <- make_matrices(forms, visit_data, site_data)

    ## SETUP DATA FOR TMB
    tmb_dat <- list(flag = 1L,
                    nsites = nsites,
                    nocc = nocc,
                    y = visit_data$obs,
                    totsite = totsite,
                    nvisit = nvisit,
                    X_psi = mats$X_psi,
                    psi_ind = attr(mats$X_psi, "index") - 1,
                    S_psi = mats$S_psi$S,
                    S_psi_n = as.integer(mats$S_psi$Scols),
                    X_p = mats$X_p,
                    p_ind = attr(mats$X_p, "index") - 1,
                    S_p = mats$S_p$S,
                    S_p_n = as.integer(mats$S_p$Scols))

    ## PARAMETERS
    map <- list()
    random <- NULL
    # fixed effects
    beta_psi <- rep(0, mats$nfix_psi)
    beta_p <- rep(0, mats$nfix_p)
    if (!is.null(start)) {
        len <- min(length(start$beta_psi), beta_psi)
        beta_psi[1:len] <- start$beta_psi[1:len]
        len <- min(length(start$beta_p), beta_p)
        beta_p[1:len] <- start$beta_p[1:len]
    }

    # random effects
    # occupancy
    if (is.null(mats$S_psi)) {
        z_psi <- 0
        log_lambda_psi <- 0
        tmb_dat$S_psi <- as(matrix(0, 1, 1), "sparseMatrix")
        tmb_dat$S_psi_n <- -1
        map <- c(map, list(z_psi = as.factor(NA), log_lambda_psi = as.factor(NA)))
    } else {
        z_psi <- rep(0, ncol(mats$X_psi) - mats$nfix_psi)
        log_lambda_psi <- rep(0, length(tmb_dat$S_psi_n))
        random <- c(random, "z_psi")
    }
    # detection
    if (is.null(mats$S_p)) {
        z_p <- 0
        log_lambda_p <- 0
        tmb_dat$S_p <- as(matrix(0, 1, 1), "sparseMatrix")
        tmb_dat$S_p_n <- -1
        map <- c(map, list(z_p = as.factor(NA), log_lambda_p = as.factor(NA)))
    } else {
        z_p <- rep(0, ncol(mats$X_p) - mats$nfix_p)
        log_lambda_p <- rep(0, length(tmb_dat$S_p_n))
        random <- c(random, "z_p")
    }
    if (length(map) < 1) map <- NULL

    ## SETUP PARAMETERS FOR TMB
    tmb_par <- list(beta_psi = beta_psi,
                    beta_p = beta_p,
                    z_psi = z_psi,
                    z_p = z_p,
                    log_lambda_psi = log_lambda_psi,
                    log_lambda_p = log_lambda_p)

    ## CREATE MODEL OBJECT
    oo <- MakeADFun(data = tmb_dat,
                    parameters = tmb_par,
                    map = map,
                    random = random,
                    DLL = "occu_tmb",
                    silent = !print)

    ## NORMALIZE
    oo <- normalize(oo, flag = "flag")

    ## FIT MODEL
    fit <- nlminb(start = oo$par, objective = oo$fn, gradient = oo$gr)

    ## GET INFERENCE
    res <- sdreport(oo, getJointPrecision = TRUE)

    ## RESULTS
    val <- list(res = res, fit = fit, forms = forms, mats = mats)
    class(val) <- "occuR"

    return(val)

}


# Simple input checks for fit_occu
check_inp <- function(forms, visit_data, site_data, start, print) {

    # check visit_data
    if (!("data.table" %in% class(visit_data))) stop("visit_data should be a data.table")
    if (!("obs" %in% colnames(visit_data))) stop("visit_data does not have a column named 'obs'")
    if (!("site" %in% colnames(visit_data))) stop("visit_data does not have a column named 'site'")
    if (!("occasion" %in% colnames(visit_data))) stop("visit_data does not have a column named 'occasion'")
    if (!("visit" %in% colnames(visit_data))) stop("visit_data does not have a column named 'visit'")
    if ("p" %in% colnames(visit_data)) stop("visit_data cannot have a column named 'p'")
    if (any(abs(visit_data$obs) > 1e-10 & abs(visit_data$obs - 1) > 1e-10)) stop("visit_data obs has entries that are not zero or one")

    # check site data
    num <- site_data[, .(max = max(occasion), n = uniqueN(occasion))]
    if (any(num$max != num$n)) stop("site_data has missing occasions or occasions are mis-numered")
    if ("psi" %in% colnames(site_data)) stop("site_data cannot have a column named 'psi'")

    # consistency between visit_data and site_data
    if (!all(visit_data$site %in% site_data$site)) stop("visit_data has sites not included in site_data")
    siteocc <- paste0(visit_data$site, visit_data$occasion)
    siteocc2 <- paste0(site_data$site, site_data$occasion)
    if (!all(siteocc %in% siteocc2)) stop("visit_data has sites with occasions not included in site_data")
    if (!all(siteocc2 %in% siteocc)) stop("site_data has sites with occasions where no visits are recorded in visit_data")

    # check start
    if (!is.null(start)) {
        if (!("beta_p" %in% names(start))) stop("start must have list entry named 'beta_p'")
        if (!("beta_psi" %in% names(start))) stop("start must have a list entry named 'beta_psi'")
        if (!is.numeric(start$beta_p) | !is.numeric(start$beta_psi)) stop("start vectors are invalid")
    }

    if (!is.logical(print)) stop("print must be TRUE or FALSE")

}

# Summary of occuR object
#' @export
summary.occuR <- function(obj, conf.level = 0.95) {
    est <- obj$res$par.fixed
    nms <- names(est)
    sd <- sqrt(diag(obj$res$cov.fixed))
    alpha <- 1 - (1 - conf.level) / 2
    lcl <- est - qnorm(alpha) * sd
    ucl <- est + qnorm(alpha) * sd
    val <- data.frame(estimate = est, sd = sd, lcl = lcl, ucl = ucl)
    val <- val[1:(obj$mats$nfix_psi + obj$mats$nfix_p),]
    nms <- nms[1:(obj$mats$nfix_psi + obj$mats$nfix_p)]
    val <- signif(val, 4)
    nms[nms == "beta_psi"] <- paste("psi:", colnames(obj$mats$X_psi[,1:obj$mats$nfix_psi, drop=FALSE]))
    nms[nms == "beta_p"] <- paste("p:", colnames(obj$mats$X_p[,1:obj$mats$nfix_p, drop = FALSE]))
    rownames(val) <- nms
    print(val)
    dof <- round(as.numeric(dof.occuR(obj, each = TRUE)), 1)
    cat(paste0("Total edf: ", sum(dof), "\t smooth_psi: ", dof[2], "\t smooth_p: ", dof[3]))
    invisible(obj)
}

# Log-likelihood for fitted occupancy model
#' @export
logLik.occuR <- function(object, ...) {
    val <- -object$fit$objective
    attributes(val)$df <- dof.occuR(object)
    return(val)
}

# Print of occuR object
#' @export
print.occuR <- function(obj) {
    summary(obj)
}
