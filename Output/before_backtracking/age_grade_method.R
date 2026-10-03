library(splines)

#' Create a spline basis specification
#'
#' Store knots, boundaries, degree, and a fixed evaluation grid for
#' consistent fitting, cross-validation, and prediction.
#'
#' @param x Finite ages; at least four observations.
#' @param knots Candidate internal knots; default integer ages 13--80.
#'   Only distinct knots strictly inside the boundaries are retained.
#' @param degree Spline degree; default 3 (cubic).
#' @param boundary Two increasing boundary ages; default range(x).
#' @param grid Evaluation ages for fitted-value second differences;
#'   defaults to a unit-spaced grid across the boundaries. Use an
#'   equally spaced annual grid for the intended curvature penalty.
#'
#' @return A list of internal knots, boundary knots, degree, number of
#'   basis functions, and the fixed evaluation grid.
make_spline_spec <- function(x, knots = seq(13, 80, 1), degree = 3,
                             boundary = range(x), 
                             grid = seq(boundary[1], boundary[2])) {
  
  # Check the supplied input.
  stopifnot(length(x) >= 4, all(is.finite(x)), diff(boundary) > 0)
  
  # Keep distinct internal knots strictly inside the domain.
  knots <- sort(unique(
    knots[knots > boundary[1] & knots < boundary[2]]
  ))
  
  # Store the basis specification and evaluation grid.
  list(knots = knots, Boundary.knots = boundary, degree = degree,
       n_basis = length(knots) + degree + 1L, grid = grid)
}

#' Construct a B-spline design matrix from a stored specification
#'
#' @param x Numeric ages at which to evaluate the basis.
#' @param spec A specification returned by make_spline_spec().
#'
#' @return A numeric matrix with one row per age and one column per
#'   basis function. The basis includes an intercept.
spline_basis <- function(x, spec) {
  as.matrix(splines::bs(x, knots = spec$knots, degree = spec$degree,
    intercept = TRUE, Boundary.knots = spec$Boundary.knots
  ))
}

#' Calculate the summed asymmetric squared-error loss
#'
#' Assign weight rho when the fitted standard is slower than the
#' observed time. For positive times, this yields a grade above 100
#' percent. All other residuals receive weight one. No smoothness or
#' curvature penalty is included in this loss.
#'
#' @param y Observed times in seconds.
#' @param fitted Fitted times in seconds, with length equal to y.
#' @param rho Positive residual weight; 1 gives symmetric squared loss.
#'
#' @return The sum of weighted squared residuals.
asymmetric_loss <- function(y, fitted, rho) {
  sum(ifelse(y < fitted, rho, 1) * (y - fitted)^2)
}


#' Fit a penalized spline for comparative age-grading standards
#'
#' Use iterative weighted linear solves to target asymmetric squared
#' loss plus a second-difference penalty on spline coefficients and a
#' soft penalty on negative fitted-value second differences. The latter
#' is evaluated on the stored grid, not the observed-age sequence.
#'
#' @param x Finite observed ages corresponding to y.
#' @param y Finite observed times in seconds, paired with x.
#' @param lambda Positive coefficient smoothness weight; default 1.
#' @param rho Positive asymmetric residual weight; default 10.
#' @param eta Nonnegative curvature penalty weight; default 100.
#' @param spline_spec Required specification from make_spline_spec().
#'   Reuse the full-domain basis and grid for held-out fits.
#' @param max_iter Maximum number of update steps; default 100.
#' @param tol Absolute maximum coefficient-change tolerance; default
#'   1e-8. Stopping also requires both sets of weights to be unchanged
#'   from their previous iteration values.
#' @param convex_after Penalize negative second differences only at
#'   grid centers strictly above this age; default 35.
#'
#' @return A list containing sorted observations, coefficients, fitted
#'   times, settings, residual weights from the last update, final grid
#'   second differences, curvature-region indicators, violation counts,
#'   iteration count, and the stored spline specification.
fit_asymmetric_spline <- function(x, y, lambda = 1, rho = 10, 
                                  eta = 100, spline_spec, 
                                  max_iter = 100, tol = 1e-8,
                                  convex_after = 35) {
  
  # Sort ages and their corresponding observed times.
  o <- order(x)
  x <- x[o]
  y <- y[o]
  n <- length(x)
  
  # Construct the observation-level basis.
  B <- spline_basis(x, spline_spec)
  p <- ncol(B)
  D <- diff(diag(p), differences = 2)
  
  # Penalize second differences of adjacent spline coefficients.
  D2_beta <- diff(diag(p), differences = 2)
  Omega <- t(D2_beta) %*% D2_beta
  
  # Evaluate curvature on the stored grid, including held-out ages.
  grid <- spline_spec$grid
  n_grid <- length(grid)
  B_grid <- spline_basis(grid, spline_spec)
  D2_fit <- diff(diag(n_grid), differences = 2)
  
  # Select grid centers strictly above the curvature threshold.
  convex_region <- grid[-c(1, n_grid)] > convex_after
  
  # Initialize with symmetric loss and the coefficient penalty.
  beta <- solve(
    t(B) %*% B + lambda * Omega,
    t(B) %*% y
  )
  
  weights_old <- rep(1, n)
  convex_weights_old <- rep(0, n_grid - 2)
  
  for (iter in seq_len(max_iter)) {
    
    fitted <- as.vector(B %*% beta)
    
    # Weight observations with fitted times above observed times by rho.
    weights <- ifelse(y < fitted, rho, 1)
    W <- diag(weights)
    
    # Compute second differences of fitted times on the fixed grid.
    second_diff <- as.vector(D2_fit %*% B_grid %*% beta)
    
    # Activate the penalty for negative differences above the threshold.
    convex_weights <- ifelse(second_diff < 0 & convex_region, 1, 0)
    C <- diag(convex_weights)
    
    # Construct the active, soft curvature penalty matrix.
    ConvexPenalty <- t(B_grid) %*% t(D2_fit) %*%
      C %*% D2_fit %*% B_grid
    
    lhs <- t(B) %*% W %*% B +
      lambda * Omega +
      eta * ConvexPenalty
    
    rhs <- t(B) %*% W %*% y
    
    beta_new <- solve(lhs, rhs)
    
    # Require a small absolute change and stable weight patterns.
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
  
  # Recompute fitted times and grid differences at the final iterate.
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
#' Evaluate the stored basis and fitted coefficients at new ages.
#'
#' @param fit A fitted model returned by fit_asymmetric_spline().
#' @param newx Finite ages within the stored boundary knots, inclusive.
#'
#' @return Predicted standard times in seconds. Nonfinite ages or ages
#'   outside the fitted domain cause an error.
predict_asymmetric_spline <- function(fit, newx) {

  # Reject nonfinite ages and predictions outside the fitted domain.
  if (any(!is.finite(newx)) ||any( 
    newx < fit$spline_spec$Boundary.knots[1] |
    newx > fit$spline_spec$Boundary.knots[2])) {
    stop("Prediction outside fitted age range")
  }
  
  # Predict with the stored basis and fitted coefficients.
  drop(spline_basis(newx, fit$spline_spec) %*% fit$beta)
}


#' Cross-validate lambda for an asymmetric penalized spline
#'
#' Use the same full-domain basis and annual grid in every fold.
#' Score held-out observations by asymmetric squared loss, without
#' adding coefficient or curvature penalties to the validation score.
#'
#' @param x Numeric observed ages, paired with y.
#' @param y Numeric observed times in seconds.
#' @param lambdas Positive candidate smoothness weights, in tie order.
#' @param knots Candidate internal knots for make_spline_spec().
#' @param rho Positive asymmetric residual weight; default 10.
#' @param eta Nonnegative curvature penalty weight; default 100.
#' @param degree Spline degree; default 3 (cubic).
#' @param k Number of approximately equal-sized folds; default 5.
#' @param seed Random seed for permuting fold labels; default 1.
#' @param convex_after Penalize negative second differences at grid
#'   centers strictly above this age; default 35.
#'
#' @return A list of candidate lambdas, mean held-out losses, selected
#'   lambda, final full-data fit, and the spline specification.
#'
#' @details Each score is the sum of all held-out losses divided by
#'   the total observation count. The first minimum in candidate order
#'   is selected. rho and eta remain fixed throughout cross-validation.
#'   Setting the fold seed resets R's random-number generator state.
cv_asymmetric_spline <- function(x, y, lambdas, knots,
                                 rho = 10, eta = 100, degree = 3,
                                 k = 5, seed = 1, convex_after = 35) {
  # Check the supplied input.
  stopifnot(k >= 2, k <= length(y), all(lambdas > 0))
  
  # Reuse the full-domain basis and annual grid in every fold.
  spline_spec <- make_spline_spec(x, knots, degree)
  
  # Randomly assign approximately equal-sized folds using the seed.
  set.seed(seed)
  folds <- sample(rep(seq_len(k), length.out = length(y)))
  errors <- matrix(NA_real_, length(lambdas), k)
  
  # Use the same folds for every candidate lambda.
  for (j in seq_along(lambdas)) {
    for (v in c(1:k)) {
      train <- folds != v
      
      # Fit training observations with the full-domain specification.
      f <- fit_asymmetric_spline(x[train], y[train], lambdas[j],
                                 rho, eta, spline_spec, 
                                 convex_after = convex_after)
      
      # Sum asymmetric squared losses over the held-out observations.
      errors[j, v] <- asymmetric_loss(y[!train],
          predict_asymmetric_spline(f, x[!train]), rho)
    }
  }
  
  # Average held-out loss per observation; choose the first minimum.
  cv_error <- rowSums(errors) / length(y)
  best <- which.min(cv_error)
  
  # Refit all observations with the selected lambda.
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
