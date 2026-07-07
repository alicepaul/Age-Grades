library(splines)

#' Create a spline basis specification
#'
#' Creates and stores the knot locations, boundary knots, and degree needed to
#' construct a consistent B-spline basis for fitting and prediction.
#'
#' @param x Numeric vector of predictor values.
#' @param knots Numeric vector of internal knot locations.
#' @param degree Degree of the spline basis. Default is 3 for cubic splines.
#'
#' @return A list containing the internal knots, boundary knots, spline degree,
#'   and number of basis functions.
make_spline_spec <- function(x, knots, degree = 3) {
  
  # Construct a B-spline basis so we can extract its attributes
  B <- splines::bs(
    x,
    knots = knots,
    degree = degree,
    intercept = TRUE
  )
  
  list(
    knots = attr(B, "knots"),
    Boundary.knots = attr(B, "Boundary.knots"),
    degree = degree,
    n_basis = ncol(B)
  )
}


fit_asymmetric_spline <- function(x, y, lambda = 1, rho = 10, 
                                  eta = 100, spline_spec, 
                                  max_iter = 100, tol = 1e-8,
                                  convex_after = 40) {
  
  # Sort data by x so second differences are meaningful
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  n <- length(y)
  
  # Construct B-spline design matrix using the fixed spline specification
  B <- splines::bs(
    x,
    knots = spline_spec$knots,
    degree = spline_spec$degree,
    intercept = TRUE,
    Boundary.knots = spline_spec$Boundary.knots
  )
  B <- as.matrix(B)
  
  p <- ncol(B)
  
  # Smoothness penalty on adjacent spline coefficients
  D2_beta <- diff(diag(p), differences = 2)
  Omega <- t(D2_beta) %*% D2_beta
  
  # Second-difference operator for fitted values
  D2_fit <- diff(diag(n), differences = 2)
  
  # Identify which second differences should be subject to convexity.
  # The second difference at row i corresponds to x[i], x[i + 1], x[i + 2],
  # so we use the middle point x[i + 1] as its location.
  convex_region <- x[-c(1, n)] > convex_after
  
  # Initial fit: ordinary penalized spline
  beta <- solve(
    t(B) %*% B + lambda * Omega,
    t(B) %*% y
  )
  
  weights_old <- rep(1, n)
  convex_weights_old <- rep(0, n - 2)
  
  for (iter in seq_len(max_iter)) {
    
    fitted <- as.vector(B %*% beta)
    
    # Asymmetric residual weights
    weights <- ifelse(y < fitted, rho, 1)
    W <- diag(weights)
    
    # Compute second differences of the fitted values
    second_diff <- as.vector(D2_fit %*% fitted)
    
    # Penalize negative second differences only in the region x > convex_after
    convex_weights <- ifelse(second_diff < 0 & convex_region, 1, 0)
    C <- diag(convex_weights)
    
    # Soft convexity penalty
    ConvexPenalty <- t(B) %*% t(D2_fit) %*% C %*% D2_fit %*% B
    
    lhs <- t(B) %*% W %*% B +
      lambda * Omega +
      eta * ConvexPenalty
    
    rhs <- t(B) %*% W %*% y
    
    beta_new <- solve(lhs, rhs)
    
    if (max(abs(beta_new - beta)) < tol &&
        all(weights == weights_old) &&
        all(convex_weights == convex_weights_old)) {
      beta <- beta_new
      break
    }
    
    beta <- beta_new
    weights_old <- weights
    convex_weights_old <- convex_weights
  }
  
  fitted <- as.vector(B %*% beta)
  second_diff <- as.vector(D2_fit %*% fitted)
  
  list(
    x = x,
    y = y,
    fitted = fitted,
    beta = beta,
    lambda = lambda,
    rho = rho,
    eta = eta,
    convex_after = convex_after,
    weights = weights,
    second_diff = second_diff,
    convex_region = convex_region,
    convex_violations = sum(second_diff < 0 & convex_region),
    convex_violations_all = sum(second_diff < 0),
    iterations = iter,
    spline_spec = spline_spec
  )
}


#' Predict from an asymmetric penalized spline
#'
#' Uses a fitted asymmetric spline model to predict fitted values at new
#' predictor values.
#'
#' @param fit A fitted model object returned by `fit_asymmetric_spline()`.
#' @param newx Numeric vector of new predictor values.
#'
#' @return Numeric vector of predicted values.
predict_asymmetric_spline <- function(fit, newx) {
  
  # Construct the same spline basis used during model fitting
  B_new <- splines::bs(
    newx,
    knots = fit$spline_spec$knots,
    degree = fit$spline_spec$degree,
    intercept = TRUE,
    Boundary.knots = fit$spline_spec$Boundary.knots
  )
  
  as.vector(B_new %*% fit$beta)
}


#' Cross-validate lambda for an asymmetric penalized spline
#'
#' Performs K-fold cross-validation over a grid of candidate lambda values.
#' The validation loss uses the same asymmetric squared-error loss used in
#' model fitting.
#'
#' @param x Numeric vector of predictor values.
#' @param y Numeric vector of response values.
#' @param lambdas Numeric vector of candidate lambda values.
#' @param knots Numeric vector of internal knot locations.
#' @param rho Asymmetric residual penalty.
#' @param eta Convexity penalty parameter.
#' @param degree Degree of the spline basis.
#' @param k Number of cross-validation folds.
#' @param seed Random seed for fold assignment.
#'
#' @return A list containing the candidate lambdas, CV errors, selected lambda,
#'   final fitted model, and spline specification.
cv_asymmetric_spline <- function(x, y, lambdas, knots,
                                 rho = 10, eta = 100, degree = 3,
                                 k = 5, seed = 1) {
  
  set.seed(seed)
  n <- length(y)
  
  # Randomly assign observations to folds
  folds <- sample(rep(seq_len(k), length.out = n))
  
  # Use the same spline basis specification across all folds
  spline_spec <- make_spline_spec(x, knots, degree = degree)
  
  cv_error <- numeric(length(lambdas))
  
  for (j in seq_along(lambdas)) {
    
    lambda <- lambdas[j]
    fold_error <- numeric(k)
    
    for (fold in seq_len(k)) {
      
      train <- folds != fold
      test  <- folds == fold
      
      # Fit model on training fold
      fit <- fit_asymmetric_spline(
        x = x[train],
        y = y[train],
        lambda = lambda,
        rho = rho,
        eta = eta,
        spline_spec = spline_spec
      )
      
      # Predict on held-out fold
      yhat <- predict_asymmetric_spline(fit, x[test])
      
      # Asymmetric validation loss:
      # penalize predictions above observed values more heavily
      fold_error[fold] <- mean(ifelse(
        y[test] < yhat,
        rho * (y[test] - yhat)^2,
        (y[test] - yhat)^2
      ))
    }
    
    cv_error[j] <- mean(fold_error)
  }
  
  # Select lambda with lowest CV error
  best_lambda <- lambdas[which.min(cv_error)]
  
  # Refit model on full data using selected lambda
  final_fit <- fit_asymmetric_spline(
    x = x,
    y = y,
    lambda = best_lambda,
    rho = rho,
    eta = eta,
    spline_spec = spline_spec
  )
  
  list(
    lambdas = lambdas,
    cv_error = cv_error,
    best_lambda = best_lambda,
    final_fit = final_fit,
    spline_spec = spline_spec
  )
}
