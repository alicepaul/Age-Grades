library(splines)

#' Create a spline basis specification
#'
#' Creates and stores the knot locations, boundary knots, and degree 
#' needed to construct a consistent B-spline basis for fitting and 
#' prediction.
#'
#' @param x Numeric vector of predictor values.
#' @param knots Numeric vector of internal knot locations.
#' @param degree Degree of the spline basis. Default is 3 for cubic 
#' splines.
#'
#' @return A list containing the internal knots, boundary knots, 
#' spline degree, and number of basis functions.
make_spline_spec <- function(x, knots = seq(13, 80, 4), degree = 3,
                             boundary = range(x), 
                             grid = seq(boundary[1], boundary[2])) {
  
  # check input
  stopifnot(length(x) >= 4, all(is.finite(x)), diff(boundary) > 0)
  
  # keep distinct internal knots inside the domain
  knots <- sort(unique(
    knots[knots > boundary[1] & knots < boundary[2]]
  ))
  
  # construct a list with spline information
  list(knots = knots, Boundary.knots = boundary, degree = degree,
       n_basis = length(knots) + degree + 1L, grid = grid)
}

#' Construct a B-spline design matrix from a stored specification
#'
#' @param x Numeric ages at which to evaluate the basis.
#' @param spec A specification returned by make_spline_spec().
#'
#' @return A numeric matrix with one row per age and one column per
#'   basis function, in the order defined by spec.
spline_basis <- function(x, spec) {
  as.matrix(splines::bs(x, knots = spec$knots, degree = spec$degree,
    intercept = TRUE, Boundary.knots = spec$Boundary.knots
  ))
}

#' Calculate the summed asymmetric squared-error loss
#'
#' Assign weight rho when the fitted standard is slower than the
#' observed time. For positive times, this condition gives an age grade
#' above 100 percent. All other residuals receive weight one.
#'
#' @param y observed times in seconds
#' @param fitted fitted times in seconds, with length equal to y.
#' @param rho penalty weight, rho = 1 gives symmetric squared-error loss
#'
#' @return the sum of weighted squared residuals.
asymmetric_loss <- function(y, fitted, rho) {
  sum(ifelse(y < fitted, rho, 1) * (y - fitted)^2)
}


#' Fit a penalized spline for comparative age-grading standards
#'
#' Minimize asymmetric squared loss plus a second-difference penalty
#' on spline coefficients and a soft penalty on negative curvature.
#'
#' @param x Distinct, finite ages; at least four observations.
#' @param y Finite observed times in seconds, paired with x.
#' @param lambda Positive coefficient smoothness penalty weight.
#' @param rho Positive asymmetric residual weight; default 30.
#' @param eta Nonnegative curvature penalty weight; default 1000.
#' @param spline_spec A stored basis specification. By default, build
#'   one from x. Reuse the full-domain specification for held-out fits.
#' @param max_iter Maximum number of update steps; default 200.
#' @param tol Relative coefficient-change tolerance; default 1e-8.
#' @param convex_after Age threshold for penalized curvature centers.
#'   The default 35 leaves youth and early-adult centers unpenalized.
fit_asymmetric_spline <- function(x, y, lambda = 1, rho = 10, 
                                  eta = 100, spline_spec, 
                                  max_iter = 100, tol = 1e-8,
                                  convex_after = 35) {
  
  # sort data by x 
  o <- order(x)
  x <- x[o]
  y <- y[o]
  n <- length(x)
  
  # initialize basis 
  B <- spline_basis(x, spline_spec)
  p <- ncol(B)
  D <- diff(diag(p), differences = 2)
  
  # smoothness penalty on adjacent spline coefficients
  D2_beta <- diff(diag(p), differences = 2)
  Omega <- t(D2_beta) %*% D2_beta
  
  # evaluate curvature on a fixed annual grid, including held-out ages
  grid <- spline_spec$grid
  n_grid <- length(grid)
  B_grid <- spline_basis(grid, spline_spec)
  D2_fit <- diff(diag(n_grid), differences = 2)
  
  # find second differences
  convex_region <- grid[-c(1, n_grid)] > convex_after
  
  # initial fit: ordinary penalized spline
  beta <- solve(
    t(B) %*% B + lambda * Omega,
    t(B) %*% y
  )
  
  weights_old <- rep(1, n)
  convex_weights_old <- rep(0, n_grid - 2)
  
  for (iter in seq_len(max_iter)) {
    
    fitted <- as.vector(B %*% beta)
    
    # asymmetric residual weights
    weights <- ifelse(y < fitted, rho, 1)
    W <- diag(weights)
    
    # compute second differences on the grid
    second_diff <- as.vector(D2_fit %*% B_grid %*% beta)
    
    # penalize negative curvature at grid centers above convex_after
    convex_weights <- ifelse(second_diff < 0 & convex_region, 1, 0)
    C <- diag(convex_weights)
    
    # soft convexity penalty
    ConvexPenalty <- t(B_grid) %*% t(D2_fit) %*%
      C %*% D2_fit %*% B_grid
    
    lhs <- t(B) %*% W %*% B +
      lambda * Omega +
      eta * ConvexPenalty
    
    rhs <- t(B) %*% W %*% y
    
    beta_new <- solve(lhs, rhs)
    
    # check for convergence
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
  
  # find info on final values
  fitted <- as.vector(B %*% beta)
  second_diff <- as.vector(D2_fit %*% B_grid %*% beta)
  
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
#' @param fit A fitted model object returned by 
#' `fit_asymmetric_spline()`.
#' @param newx Numeric vector of new predictor values.
#'
#' @return Numeric vector of predicted values.
predict_asymmetric_spline <- function(fit, newx) {

  # only predict within age range
  if (any(!is.finite(newx)) ||any( 
    newx < fit$spline_spec$Boundary.knots[1] |
    newx > fit$spline_spec$Boundary.knots[2])) {
    stop("Prediction outside fitted age range")
  }
  
  # find predicted values
  drop(spline_basis(newx, fit$spline_spec) %*% fit$beta)
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
#' @param convex_after Age threshold for penalized curvature centers.
#'
#' @return A list containing the candidate lambdas, CV errors, selected lambda,
#'   final fitted model, and spline specification.
cv_asymmetric_spline <- function(x, y, lambdas, knots,
                                 rho = 10, eta = 100, degree = 3,
                                 k = 5, seed = 1, convex_after = 35) {
  # check input 
  stopifnot(k >= 2, k <= length(y), all(lambdas > 0))
  
  # Reuse the full-domain basis and annual grid in every fold.
  spline_spec <- make_spline_spec(x, knots, degree)
  
  # set seed and folds
  set.seed(seed)
  folds <- sample(rep(seq_len(k), length.out = length(y)))
  errors <- matrix(NA_real_, length(lambdas), k)
  
  # iterate through lambdas and folds
  for (j in seq_along(lambdas)) {
    for (v in c(1:k)) {
      train <- folds != v
      
      # fit on training ages 
      f <- fit_asymmetric_spline(x[train], y[train], lambdas[j],
                                 rho, eta, spline_spec, 
                                 convex_after = convex_after)
      
      # sum held-out losses
      errors[j, v] <- asymmetric_loss(y[!train],
          predict_asymmetric_spline(f, x[!train]), rho)
    }
  }
  
  # find best lambda
  cv_error <- rowSums(errors) / length(y)
  best <- which.min(cv_error)
  
  # refit using best lambda
  final_fit = fit_asymmetric_spline(x, y, lambdas[best], rho, eta,
    spline_spec, convex_after = convex_after)
  
  list(
    lambdas = lambdas,
    cv_error = cv_error,
    best_lambda = lambdas[best],
    final_fit = final_fit,
    spline_spec = spline_spec
  )
}
