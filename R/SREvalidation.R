.response_name <- function(object) all.vars(object@f)[1]

#' @rdname SRE
#' @export
simulate <- function(object, newdata = NULL, nsim = 400, conditional_fs = FALSE, ...) {
  
  if (is.null(newdata)) {
    if (length(object@data) > 1) warning("The SRE object contains multiple data sets; the fitted values will only be computed for the first data set.")
    newdata <- object@data[[1]]
  }
  
  if ("obs_fs" %in% names(list(...))) stop("obs_fs should not be provided when calling simulate")
  
  if (object@response %in% c("binomial", "negative-binomial")) {
    
    # the user must provide k_Z with newdata (only extract it later after the predict step):
    if(is.null(newdata@data$k_Z)) stop("The data model is binomial or negative-binomial; please provide the observation-level size parameters in newdata@data$k_Z") 
    
    # Fix the size parameter over the BAUs to 1 (we will use k_Z)
    k <- 1
    k <- rep(k, length(object@BAUs)) # repeat k to prevent message from predict()
    if (!is.null(list(...)$k)) stop("k should not be provided when calling simulate")
  }
  
  # measurement-error variance parameters
  if (object@response == "gaussian") sigma2e <- .measurementerrorvariance(object, newdata)
  
  # fine-scale variance parameters
  obs_fs <- !conditional_fs
  if (obs_fs) {
    # Compute the fine-scale variance parameter for each location. This is 
    # complicated since the constant of proportionality can depend on 
    # object@BAUs$fs, and this complication is compounded by the argument newdata. 
    # The following lines deal with this by constructing the incidence matrix 
    # for newdata, which can then be used to map the fine-scale variance parameters
    # to the observations. 
    Z <- map_data_to_BAUs(newdata, object@BAUs, average_in_BAU = FALSE)
    C <- BuildC(Z, object@BAUs)
    fs <- object@BAUs$fs[C$j_idx]
    sigma2fs <- object@sigma2fshat * fs
  }
  
  ## Simulate the data
  if (object@method == "EM") {
    
    # Simulate the response variable
    pred <- predict(object, newdata = newdata, obs_fs = obs_fs, ...)
    Zhat <- pred$mu
  
    if (obs_fs) {
      sdZ <- sqrt(pred$sd^2 + sigma2fs + sigma2e)
    } else {
      sdZ <- sqrt(pred$sd^2 + sigma2e)
    }
    samples <- t(sapply(seq_along(Zhat), function(i) rnorm(nsim, Zhat[i], sdZ[i]))) 
    
  } else {
    
    ## Simulate the latent process Y()
    if (conditional_fs) {
      
      pred <- predict(object, newdata = newdata, type = "link", percentiles = NULL, nsim = nsim, obs_fs = obs_fs, k = k, ...)
      Y_samples <- pred$MC$Y_samples
      
    } else {
      
      # Obtain samples of the smooth version (i.e., excluding fine-scale variation) 
      # of the latent process, Y(.)
      pred <- predict(object, newdata = newdata, type = "link", percentiles = NULL, nsim = nsim, obs_fs = obs_fs, k = k, ...)
      smooth_Y_samples <- pred$MC$Y_samples
    
      # Add marginal fine-scale variation to the samples of Y(.). Let m 
      # denote the number of locations, where m == nrow(Y_samples) == length(sigma2fs). 
      # For each location, we have nsim samples of Y(.) stored in the columns of
      # Y_samples; hence, before sampling the fine-scale variation xi(.), we 
      # repeat each element of sigm2fs nsim times.  
      m <- nrow(smooth_Y_samples)
      sigma2fs_long <- rep(sigma2fs, each = nsim)
      xi_samples <- rnorm(m * nsim, sd = sqrt(sigma2fs_long))
      xi_samples <- matrix(xi_samples, byrow = TRUE, ncol = nsim)
      Y_samples  <- smooth_Y_samples + xi_samples
    }
    
    # compute mu(.) by applying the inverse link function to Y(.), and 
    # simulate the response, Z
    if (object@response %in% c("binomial", "negative-binomial")) {
      # The user must provide k_Z with newdata (note that this is checked
      # upon function entry). We need to extract k_Z to account for the case 
      # that not all locations in newdata overlap BAUs (in which case pred$newdata
      # will be of a different length to the original newdata object)
      k_Z <- pred$newdata@data$k_Z 
    } else {
      k_Z <- NULL
    }
    mu_samples <- .inverselink(Y_samples, object, k = k_Z)$mu_samples 
    Z_samples  <- .sampleZ(mu_samples = mu_samples, object = object, sigma2e = sigma2e, k = k_Z)
    
    samples <- Z_samples
  }
  
  return(samples)
}

#' @rdname SRE
#' @export
setMethod("fitted", signature="SRE", function(object, ...) {
  
  if (length(object@data) > 1) warning("The SRE object contains multiple data sets; the fitted values will only be computed for the first data set.")
  newdata <- object@data[[1]]
  
  if (object@method == "EM") {
    Zhat <- predict(object, newdata = newdata, ...)$mu
  } else {
    Zhat <- predict(object, newdata = newdata, type = "mean", percentiles = NULL, ...)$newdata$p_mu
  }
  
  return(Zhat)
})


#' @rdname SRE
#' @export
setMethod("residuals", signature="SRE", function(object, type = "pearson") { 

  # Extract the data
  if (length(object@data) > 1) warning("The SRE object contains multiple data sets; the residuals will only be computed for the first data set.")
  newdata <- object@data[[1]]
  Z       <- newdata@data[, .response_name(object)]
  
  # Fitted values. Use simulate() so that we can estimate the variance when 
  # computing pearson residuals. 
  Zsim <- simulate(object, newdata, conditional_fs = TRUE)
  Zhat <- rowMeans(Zsim)
  
  # Compute the "response" residuals
  res <- Zhat - Z
  
  if (type == "pearson") {
    Vhat <- apply(Zsim, 1, var)
    res <- res/sqrt(Vhat)
  } 
  
  return(res)
})

#' @rdname SRE
#' @export
setMethod("AIC", signature="SRE", function(object, k = 2) {
  
  ## Initialise the number of covariance parameters to zero
  p <- 0
  
  ## Fine-scale variance parameters
  if (object@fs_by_spatial_BAU) {
    ## If we are estimating a unique fine-scale variance parameter at each 
    ## spatial BAU, we will have as many fine-scale variance parameters as 
    ## spatial BAUs
    p <- p + dim(object@BAUs)[1]
  } else {
    ## Otherwise, we just have a single fine-scale variance parameter
    p <- p + 1
  }
  
  ## Dispersion parameter
  ## The following distributions do not have a dispersion parameter
  if (!(object@response %in% c("poisson", "binomial", "negative-binomial"))) {
    p <- p + 1
  }
  
  ## Basis-function covariance parameters
  if (is(object@basis,"TensorP_Basis")) {
    # We model the temporal basis-function coefficients as an AR1 process, which 
    # has a marginal variance parameter and a correlation parameter.
    temporal_cov_parameters <- 2 
    spatial_basis <- object@basis@Basis1
  } else {
    temporal_cov_parameters <- 0
    spatial_basis <- object@basis
  }

  R <- nres(spatial_basis)
  if (object@K_type == "block-exponential") {
    spatial_cov_parameters <- 2 * R   
  } else if (object@K_type == "precision") {
    
    if (object@basis@regular) {
      spatial_cov_parameters <- 2 * R 
    } else {
      spatial_cov_parameters <- 3 * R 
    }
    
  } else if (object@K_type == "unstructured") {
    r <- table(spatial_basis@df$res) # number of spatial basis functions at each resolution
    spatial_cov_parameters <- sum(r*(r+1)/2)
  }
  
  p <- p + spatial_cov_parameters + temporal_cov_parameters 
  
  ## Number of fixed effects
  q <- length(coef(object))
  
  ## The criterion
  k * (p + q) - 2 * logLik(object)
})



#' @rdname SRE
#' @export
setMethod("BIC", signature="SRE", function(object) AIC(object, k = log(nobs(object))))