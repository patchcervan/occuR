#' Compute block-diagonal smoothing matrix
#'
#' @param gm value returned by mgcv::gam with fit = FALSE
#'
#' @return list of S (block diagonal smoothing matrix), Scols (number of columns for each smoothing matrix),
#' bdim (k argument used for each smooth), Snames (names of smoothing parameters). Returns NULL if no smooths in
#' model.
#'
#' @description Based on code written by David L Miller in CTMCdive R Package (https://github.com/r-glennie/CTMCdive).
#'
#' @export
#' @importFrom Matrix bdiag
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

#' Compute design and smoothing matrices from formulae
#'
#' @param forms list of formulae for psi and p
#' @param re list with character vectors of random effects variables for psi and p
#' @param site_data see fit_occu function
#' @param visit_data see fit_occu function
#'
#' @return list of X = fixed effects design matrix, Z = random effects design matrix,
#'   and S = smoothing matrix for both psi and p (e.g. X_psi)
#' @export
#'
#' @examples
#' @importFrom mgcv gam uniquecombs
#' @importFrom mgcv predict.gam
#' @importFrom methods as
#' @importFrom stats terms reformulate model.matrix contrasts
make_matrices <- function(forms, re, visit_data, site_data) {

    ## occupancy model
  site_data$psi <- 1:nrow(site_data)
  gam_psi <- gam(forms[["psi"]], data = site_data, method = "REML")

  # get design matrix
  Xfull <- predict(gam_psi, newdata = site_data, type = "lpmatrix")

  if(length(re[[1]]) > 0){
      # Add random effects
      Xre <- model.matrix(as.formula(paste("psi ~", paste(re[[1]], collapse = "+"))),
                          data = site_data,
                          contrasts.arg = lapply(site_data[, re[[1]], with = FALSE],
                                                 contrasts, contrasts = FALSE))
      Xfull <- cbind(Xfull, Xre[,-1])
      U_psi_n <- sapply(site_data[, re[[1]], with = FALSE], nlevels)
      X_psi <- uniquecombs(Xfull)
      U_psi <- X_psi[,tail(seq_len(ncol(Xfull)), sum(U_psi_n))]
      attr(U_psi, "index") <- attr(X_psi, "index")
      X_psi <- X_psi[,-tail(seq_len(ncol(Xfull)), sum(U_psi_n))]
      attr(X_psi, "index") <- attr(U_psi, "index")
      U_psi <- as(U_psi, "sparseMatrix")
  } else {
      X_psi <- uniquecombs(Xfull)
      U_psi_n <- 0
      U_psi <- NULL
  }

  site_data <- site_data[, psi := NULL]

  # get smoothing matrix
  S_psi <- make_smoothing_matrix(gam_psi)

  ## detection model
  visit_data$p <- 1:nrow(visit_data)
  gam_p <- gam(forms[["p"]], data = visit_data, method = "REML")

  # get design matrix
  Xfull <- predict(gam_p, newdata = visit_data, type = "lpmatrix")

  # Add random effects
  if(length(re[[2]]) > 0){
      Xre <- model.matrix(as.formula(paste("p ~", paste(re[[2]], collapse = "+"))),
                          data = visit_data,
                          contrasts.arg = lapply(visit_data[, re[[2]], with = FALSE],
                                                 contrasts, contrasts = FALSE))
      Xfull <- cbind(Xfull, Xre[,-1])
      U_p_n <- sapply(visit_data[, re[[2]], with = FALSE], nlevels)
      X_p <- uniquecombs(Xfull)
      U_p <- X_p[, tail(seq_len(ncol(Xfull)), sum(U_p_n))]
      attr(U_p, "index") <- attr(X_p, "index")
      X_p <- X_p[, -tail(seq_len(ncol(Xfull)), sum(U_p_n))]
      attr(X_p, "index") <- attr(U_p, "index")
      U_p <- as(U_p, "sparseMatrix")
  } else {
      X_p <- uniquecombs(Xfull)
      U_p_n <- 0
      U_p <- NULL
  }

  visit_data <- visit_data[, p := NULL]

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
              gam_psi = gam_psi,
              U_psi = U_psi,
              U_p = U_p,
              U_psi_n = U_psi_n,
              U_p_n = U_p_n,
              re = re)

  return(res)
}

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
#' @importFrom stats terms reformulate
#' @useDynLib occu_tmb_re
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
  nobs <- nrow(visit_data)

  # name formulae
  names(forms) <- c(as.character(forms[[1]][[2]]),
                    as.character(forms[[2]][[2]]))

  ## Extract random effects from model formulas
  re <- extract_re_terms(forms, site_data, visit_data)

  # remove random effects from forms
  forms_fix <- forms
  for(i in seq_along(forms)){
      if(length(re[[i]]) > 0){
          tt <- labels(terms(forms[[i]]))
          tt <- tt[-grep(paste(re[[i]], collapse = "|"), tt)]
          forms_fix[[i]] <- reformulate(tt,
                                        response = c("psi", "p")[i],
                                        intercept = attr(terms(forms[[i]]), "intercept"))
      }
  }

  ## MODEL MATRICES

  # get model matrices and smoothing matrices
  mats <- make_matrices(forms_fix, re, visit_data, site_data)

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
                  S_p_n = as.integer(mats$S_p$Scols),
                  U_psi = mats$U_psi,
                  U_p = mats$U_p,
                  U_psi_n = mats$U_psi_n,
                  U_p_n = mats$U_p_n)

  ## PARAMETERS
  map <- list()
  random <- NULL

  # Initial values for fixed effects
  beta_psi <- rep(0, mats$nfix_psi)
  beta_p <- rep(0, mats$nfix_p)

  # Initial values for occupancy random effects
  if(tmb_dat$U_psi_n[1] == 0){
      tmb_dat$U_psi <- as(matrix(0, ncol = 1, nrow = 1), "sparseMatrix") # this could be any matrix of one column
      gamma_psi <- rep(0, ncol(tmb_dat$U_psi))
      lsig_gamma_psi <- rep(0, length(tmb_dat$U_psi_n))
      map <- c(map, list(gamma_psi = as.factor(NA), lsig_gamma_psi = as.factor(NA)))
  } else {
      gamma_psi <- rep(0, ncol(tmb_dat$U_psi))
      lsig_gamma_psi <- rep(0, length(tmb_dat$U_psi_n))
      random <- c(random, "gamma_psi")
  }

  # Initial values for detection random effects
  if(tmb_dat$U_p_n[1] == 0){
      tmb_dat$U_p <- as(matrix(0, ncol = 1, nrow = 1), "sparseMatrix") # this could be any matrix of one column
      gamma_p <- rep(0, ncol(tmb_dat$U_p))
      lsig_gamma_p <- rep(0, length(tmb_dat$U_p_n))
      map <- c(map, list(gamma_p = as.factor(NA), lsig_gamma_p = as.factor(NA)))
  } else {
      gamma_p <- rep(0, ncol(tmb_dat$U_p))
      lsig_gamma_p <- rep(0, length(tmb_dat$U_p_n))
      random <- c(random, "gamma_p")
  }


  # If initial values are provided by the user
  if (!is.null(start)) {
    len <- min(length(start$beta_psi), beta_psi)
    beta_psi[1:len] <- start$beta_psi[1:len]
    len <- min(length(start$beta_p), beta_p)
    beta_p[1:len] <- start$beta_p[1:len]
    if(!is.null(start$gamma_p)) gamma_p <- start$gamma_p
    if(!is.null(start$gamma_psi)) gamma_psi <- start$gamma_psi
    if(!is.null(start$lsig_gamma_psi)) lsig_gamma_psi <- start$lsig_gamma_psi
    if(!is.null(start$lsig_gamma_p)) lsig_gamma_p <- start$lsig_gamma_p
  }

  # Initial values for (splines) random effects
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
                  log_lambda_p = log_lambda_p,
                  gamma_psi = gamma_psi,
                  gamma_p = gamma_p,
                  lsig_gamma_psi = lsig_gamma_psi,
                  lsig_gamma_p = lsig_gamma_p)

  ## CREATE MODEL OBJECT
  # compile("src/occu_tmb_re.cpp")
  # dyn.load(dynlib("../occuR/src/occu_tmb_re"))
  oo <- MakeADFun(data = tmb_dat,
                  parameters = tmb_par,
                  map = map,
                  random = random,
                  DLL = "occu_tmb_re",
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

#' Simple input checks for fit_occu
#'
#' @inheritParams fit_occu
#'
#' @return stops execution if input is invalid
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

#' Summary of occuR object
#'
#' @param obj output from fit_occu
#'
#' @return prints a summary
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

#' Print of occuR object
#'
#' @param obj output from fit_occu
#'
#' @return prints a summary
#' @export
print.occuR <- function(obj) {
  summary(obj)
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
  psi_pred <- (cbind(mats$X_psi, mats$U_psi) %*% c(fix[nms_fix == "beta_psi"], ran[nms_ran == "z_psi"], ran[nms_ran == "gamma_psi"]))
  p_pred <- (cbind(mats$X_p, mats$U_p) %*% c(fix[nms_fix == "beta_p"], ran[nms_ran == "z_p"], ran[nms_ran == "gamma_p"]))
  psi_pred <- plogis(psi_pred)
  p_pred <- plogis(p_pred)
  return(list(psi = psi_pred, p = p_pred))
}

#' Predict from fitted occupancy model
#'
#' @param obj fitted model object from fit_occu
#' @param visit_data see fit_occu
#' @param site_data see fit_occu
#' @param include_re Logical indicating whether random effects contributions should
#' be included in the predictions. If TRUE (default) they are included, if FALSE
#' they are not (producing "population-level" predictions).
#' @param new_levels If TRUE, then new levels in the random effects are allowed, and
#' they will be sampled from the corresponding distribution. If there is a mix of
#' old and new levels, model estimates will be used for the old levels. If FALSE
#' (default), and `include_re == TRUE` then the function will exit with an error
#' if new levels are present.
#' @param nboot number of parametric bootstrap resamples to produce from fitted model
#'
#' @return if nboot = 0 (default), the list of fitted psi and p values; otherwise, list also
#' contains matrix for psi and p where each row is a bootstrap resample
#' @export
#' @importFrom mgcv rmvn
#' @importFrom Matrix solve
#' @importFrom stats terms reformulate model.matrix contrasts
predict.occuR <- function(obj, visit_data, site_data, include_re = TRUE,
                          new_levels = FALSE, nboot = 0) {
    mats <- obj$mats
    re <- mats$re

    # Occupancy fixed effect design matrix
    site_data$psi <- 1:nrow(site_data)
    mats$X_psi <- predict(mats$gam_psi, newdata = site_data, type = "lpmatrix")

    # Add random effects design matrix
    if(length(re[[1]]) > 0){
        mats$U_psi <- model.matrix(as.formula(paste("psi ~", paste(re[[1]], collapse = "+"))),
                                   data = site_data,
                                   contrasts.arg = lapply(site_data[, re[[1]], with = FALSE],
                                                          contrasts, contrasts = FALSE))
        # Remove intercept
        mats$U_psi <- mats$U_psi[,-1]
        mats$U_psi_n <- sapply(site_data[, re[[1]], with = FALSE], nlevels)
    }
    site_data <- site_data[, psi := NULL]

    # Detection
    visit_data$p <- 1:nrow(visit_data)
    mats$X_p <- predict(mats$gam_p, newdata = visit_data, type = "lpmatrix")

    # Add random effects
    if(length(re[[2]]) > 0){
        mats$U_p <- model.matrix(as.formula(paste("p ~", paste(re[[2]], collapse = "+"))),
                                 data = visit_data,
                                 contrasts.arg = lapply(visit_data[, re[[2]], with = FALSE],
                                                        contrasts, contrasts = FALSE))
        # Remove intercept
        mats$U_p <- mats$U_p[,-1]
        mats$U_p_n <- sapply(visit_data[, re[[2]], with = FALSE], nlevels)
    }
    visit_data <- visit_data[, p := NULL]

    if(include_re && new_levels){
        preds <- pred_mix_re(re, site_data, visit_data, mats, nboot)
    } else if(include_re && !new_levels){
        preds <- pred_re(obj, mats, nboot)
    } else if(!include_re){
        preds <- pred_no_re(obj, mats, nboot)
    }

    return(preds)

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
#' @importFrom Matrix solve diag
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

#' Log-likelihood for fitted occupancy model
#'
#' @param object fitted occupancy model from fit_occu
#' @param ...
#'
#' @return log-likelihood with estimated degrees of freedom as attribute "df"
#' @export
logLik.occuR <- function(object, ...) {
  val <- -object$fit$objective
  attributes(val)$df <- dof.occuR(object)
  return(val)
}


#' Extract random effects from formulas
#'
#' @param forms list of formulae for psi (occupancy) and p (detection)
#' @param visit_data data.table with row for each visit with a "obs" column for detection record
#' @param site_data data.table with row for each site x occasion
#' @return A list with names of random effects variables as per `forms`. The
#' levels of the random effects variables present in the data are returned as
#' attributes.
#' @export
#' @keywords internal
#' @noRd
extract_re_terms <- function(forms, site_data, visit_data){

    re <- lapply(lapply(forms, terms), labels)
    re <- lapply(re, function(x) x[grep("\\|", x)])
    re <- lapply(re, function(x) gsub(" ", "", x)) # remove spaces
    re <- lapply(re, function(x) gsub("*.\\|", "", x))

    # Extract levels for occupancy re
    for(i in seq_along(re[[1]])){
        attr(re[[1]], re[[1]][i]) <- levels(site_data[, get(re[[1]][i])])
    }
    # Extract levels for detection re
    for(i in seq_along(re[[2]])){
        attr(re[[2]], re[[2]][i]) <- levels(visit_data[, get(re[[2]][i])])
    }

    return(re)
}


#' Population-level predictions from fitted mixed-effects occupancy model
#'
#' @description This applies for population level predictions. No known random effects for the present levels.
#' @param obj fitted model object from fit_occu
#' @param mats Model matrices
#' @param nboot number of parametric bootstrap resamples to produce from fitted model
#'
#' @return if nboot = 0 (default), the list of fitted psi and p values; otherwise, list also
#' contains matrix for psi and p where each row is a bootstrap resample
#' @export
#' @keywords internal
#' @noRd
#' @importFrom mgcv rmvn
#' @importFrom Matrix solve
#' @importFrom stats terms reformulate model.matrix contrasts
pred_no_re <- function(obj, mats, nboot){

    fix <- obj$res$par.fixed
    ran <- NULL
    if(sum(unlist(mats$U_psi_n)) > 0){
        ran <- c(ran, rep(0, sum(unlist(mats$U_psi_n))))
    }
    if(sum(unlist(mats$U_p_n)) > 0){
        ran <- c(ran, rep(0, sum(unlist(mats$U_p_n))))
    }

    names(ran) <- c(rep(c("gamma_psi"), sum(unlist(mats$U_psi_n))),
                    rep(c("gamma_p"), sum(unlist(mats$U_p_n))))

    # ran <- obj$res$par.random
    pred <- get_predicted_values(fix, ran, mats)

    if (nboot > 0.5) {
        Q <- obj$res$jointPrecision
        if (!is.null(Q)) {
            Q <- Q[!grepl("log_lambda_|^gamma_", colnames(Q)),
                   !grepl("log_lambda_|^gamma_", colnames(Q)), drop = FALSE]
            V <- solve(Q)
        } else {
            V <- obj$res$cov.fixed
        }
        param <- c(fix, ran)
        param <- param[!grepl("log_lambda|^gamma_", names(param))]
        boots <- rmvn(nboot, param, V)
        colnames(boots) <- names(param)
        nfix <- length(fix[!grepl("log_lambda|gamma_", names(fix))])

        boots_psi <- matrix(0, nr = nboot, nc = length(pred$psi))
        boots_p <- matrix(0, nr = nboot, nc = length(pred$p))


        for (b in 1:nboot) {
            ran <- NULL
            gamma_psi <- NULL
            for(i in seq_along(mats$U_psi_n)){
                sig_psi <- exp(boots[b, colnames(boots) == "lsig_gamma_psi"][i]) + 1e-10 # to avoid zero sd
                gamma_psi <- c(gamma_psi, rnorm(mats$U_psi_n[i], 0, sig_psi))
            }
            ran <- c(ran, gamma_psi)
            gamma_p <- NULL
            for(i in seq_along(mats$U_p_n)){
                sig_p <- exp(boots[b, colnames(boots) == "lsig_gamma_p"][i]) + 1e-10 # to avoid zero sd
                gamma_p <- c(gamma_p, rnorm(mats$U_p_n[i], 0, sig_p))
            }
            ran <- c(ran, gamma_p)
            names(ran) <- c(rep(c("gamma_psi"), sum(unlist(mats$U_psi_n))),
                            rep(c("gamma_p"), sum(unlist(mats$U_p_n))))
            boot_res <- get_predicted_values(boots[b, 1:nfix], ran, mats)
            boots_psi[b,] <- boot_res$psi
            boots_p[b,] <- boot_res$p
        }
        pred$psiboot <- boots_psi
        pred$pboot <- boots_p
    }
    return(pred)
}


#' Factor-level predictions from fitted mixed-effects occupancy model
#'
#' @description This applies for known random effects levels.
#' @param obj fitted model object from fit_occu
#' @param mats Model matrices
#' @param nboot number of parametric bootstrap resamples to produce from fitted model
#'
#' @return if nboot = 0 (default), the list of fitted psi and p values; otherwise, list also
#' contains matrix for psi and p where each row is a bootstrap resample
#' @export
#' @keywords internal
#' @noRd
#' @importFrom mgcv rmvn
#' @importFrom Matrix solve
#' @importFrom stats terms reformulate model.matrix contrasts
pred_re <- function(obj, mats, nboot){

    fix <- obj$res$par.fixed
    ran <- obj$res$par.random
    pred <- get_predicted_values(fix, ran, mats)

    if (nboot > 0.5) {
        Q <- obj$res$jointPrecision
        if (!is.null(Q)) {
            Q <- Q[!grepl("log_lambda_|lsig_gamma_", colnames(Q)),
                   !grepl("log_lambda_|lsig_gamma_", colnames(Q)), drop = FALSE]
            V <- solve(Q)
        } else {
            V <- obj$res$cov.fixed
        }
        param <- c(fix, ran)
        param <- param[!grepl("log_lambda|lsig_gamma_", names(param))]
        boots <- rmvn(nboot, param, V)
        colnames(boots) <- names(param)
        nfix <- length(fix[!grepl("log_lambda|lsig_gamma_", names(fix))])
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

#' Factor-level and population predictions from fitted mixed-effects occupancy model
#'
#' @description This applies when we want factor-level predictions, but there
#' are new levels that we need to sample from the corresponding distribution of random effects.
#' @param re A list with two character vectors with the names of the variables with random
#' effects for occupancy and detection respectively.
#' @param visit_data data.table with row for each visit with a "obs" column for detection record
#' @param site_data data.table with row for each site x occasion
#' @param mats Model matrices
#' @param nboot number of parametric bootstrap resamples to produce from fitted model
#'
#' @return if nboot = 0 (default), the list of fitted psi and p values; otherwise, list also
#' contains matrix for psi and p where each row is a bootstrap resample
#' @export
#' @keywords internal
#' @noRd
#' @importFrom mgcv rmvn
#' @importFrom Matrix solve
#' @importFrom stats terms reformulate model.matrix contrasts
pred_mix_re <- function(re, site_data, visit_data, mats, nboot){

    # Separate random effects levels for occupancy
    if(length(re[[1]]) > 0){
        olds <- unlist(lapply(re[[1]], function(x) paste0(x, attr(re[[1]], x))))
        news <- unlist(lapply(re[[1]], function(x) paste0(x, levels(site_data[, get(x)]))))
        tot <- unique(c(olds, news))
        re_psi <- data.frame(levels = tot,
                             old = tot %in% olds,
                             new = tot %in% news,
                             param = "gamma_psi")
        re_psi$var <- stringr::str_extract(re_psi$levels, paste(unlist(re), collapse = "|"))
        re_psi <- re_psi[order(nchar(re_psi$levels), re_psi$levels),]
    }

    # Separate random effects levels for detection
    if(length(re[[2]]) > 0){
        olds <- unlist(lapply(re[[2]], function(x) paste0(x, attr(re[[2]], x))))
        news <- unlist(lapply(re[[2]], function(x) paste0(x, levels(visit_data[, get(x)]))))
        tot <- unique(c(olds, news))
        re_p <- data.frame(levels = tot,
                           old = tot %in% olds,
                           new = tot %in% news,
                           param = "gamma_p")
        re_p$var <- stringr::str_extract(re_p$levels, paste(unlist(re), collapse = "|"))
        re_p <- re_p[order(nchar(re_p$levels), re_p$levels),]
    }

    re_df <- rbind(re_psi, re_p)
    re_df$id <- paste(re_df$levels, re_df$param, sep = ".")
    re_df$var.proc <- paste(re_df$var, re_df$param, sep = ".")

    # model parameters
    fix <- obj$res$par.fixed
    ran <- obj$res$par.random
    attr(ran, "old_names") <- names(ran)
    names(ran) <- re_df[re_df$old == TRUE, "id"]


    # I need a vector of random param
    rdm <- rep(NA, nrow(re_df[re_df$new == TRUE, ]))
    names(rdm) <-  re_df[re_df$new == TRUE, "id"]
    common <- re_df[re_df$new == TRUE & re_df$old == TRUE, "id"]
    rdm[common] <- ran[common]
    attr(rdm, "missing") <- which(is.na(rdm))

    names(rdm) <- re_df[re_df$id %in% names(rdm), "param"]
    rdm[is.na(rdm)] <- 0
    pred <- get_predicted_values(fix, rdm, mats)

    nfix <- mats$nfix_psi + mats$nfix_p
    keep <- c(rep(TRUE, nfix),   # fixed effects to keep
              re_df[re_df$old == TRUE, "id"] %in% common,  # random effect levels to keep
              unique(re_df$var.proc) %in% (re_df$var.proc[re_df$old == FALSE & re_df$new == TRUE])) # random effects sds to keep
    names(ran) <- attr(ran, "old_names")

    if (nboot > 0.5) {
        Q <- obj$res$jointPrecision
        # For fixed effects and random effects SD
        if (!is.null(Q)) {
            Qf <- Q[!grepl("log_lambda_", colnames(Q)),
                    !grepl("log_lambda_", colnames(Q)), drop = FALSE]
            Qf <- Qf[keep, keep]
            Vf <- solve(Qf)
        } else {
            Vf <- obj$res$cov.fixed
        }

        # Combine fixed and random effects and remove not needed sds
        param <- c(fix[c(rep(TRUE, nfix),
                         unique(re_df$var.proc) %in% (re_df$var.proc[re_df$old == FALSE & re_df$new == TRUE]))],
                   ran[re_df[re_df$old == TRUE, "id"] %in% common])
        # Reorder
        param <- c(param[!grepl("lsig_", names(param))],
                   param[grepl("lsig_", names(param))])
        param <- param[!grepl("log_lambda", names(param))]
        boots <- rmvn(nboot, param, Vf)
        colnames(boots) <- names(param)

        boots_psi <- matrix(0, nr = nboot, nc = length(pred$psi))
        boots_p <- matrix(0, nr = nboot, nc = length(pred$p))

        n_samples <- unique(paste0(re_df$var, ".*", re_df$param, "$"))
        n_samples <- sapply(lapply(n_samples, grep, names(attr(rdm, "missing"))), length)
        to_sample <- colnames(boots)[grep("lsig", colnames(boots))]

        for (b in 1:nboot) {
            gamma_new <- NULL
            for(i in seq_along(to_sample)){
                sig_psi <- exp(boots[b, colnames(boots) == to_sample[i]]) + 1e-10 # to avoid zero sd
                gamma_new <- c(gamma_new, rnorm(n_samples[i], 0, sig_psi))
            }
            rdm[attr(rdm, "missing")] <- gamma_new
            boot_res <- get_predicted_values(boots[b, 1:nfix], rdm, mats)
            boots_psi[b,] <- boot_res$psi
            boots_p[b,] <- boot_res$p
        }
        pred$psiboot <- boots_psi
        pred$pboot <- boots_p
    }

    return(pred)

}


