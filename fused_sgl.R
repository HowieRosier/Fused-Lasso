# Fused sparse group lasso on top of sparsegl.
#
# Objective:
#   (1/(2n)) ||y - Xb||^2
#     + (1-alpha)*lambda * sum_g pf_group_g * ||b_g||
#     +    alpha *lambda * sum_j pf_sparse_j * |b_j|
#     +  lambda_fus * sum_m ||delta_{m+1} - delta_m||
#
# Matches sparsegl exactly when lambda_fus = 0 (we just call sparsegl).
# For lambda_fus > 0 we run ADMM with split z = Db; the b-update is a
# small block coordinate descent that mirrors sparsegl's Fortran
# update_step, with an extra rho*D_g'D_g term in the Lipschitz constant
# and rho*D_g'(Db - v) in the gradient.

if (!requireNamespace("sparsegl", quietly = TRUE)) {
  stop("needs the sparsegl package")
}

# group -> list of column indices per group, plus sizes
group_cols_of <- function(group) {
  group <- as.integer(group)
  bs <- as.integer(table(group))
  cum <- cumsum(c(0L, bs))
  list(
    bs = bs,
    bn = length(bs),
    cols = lapply(seq_along(bs), function(g) seq.int(cum[g] + 1L, cum[g] + bs[g]))
  )
}

# difference matrix D with D %*% beta = [delta_2 - delta_1, ..., delta_M - delta_{M-1}]
build_D <- function(p, group, fuse_groups) {
  gi <- group_cols_of(group)
  K <- gi$bs[fuse_groups[1]]
  if (any(gi$bs[fuse_groups] != K)) stop("fuse_groups must all be the same size")
  M <- length(fuse_groups)
  D <- matrix(0, nrow = (M - 1L) * K, ncol = p)
  for (m in seq_len(M - 1L)) {
    rows <- ((m - 1L) * K + 1L):(m * K)
    D[cbind(rows, gi$cols[[fuse_groups[m]]])]     <- -1
    D[cbind(rows, gi$cols[[fuse_groups[m + 1L]]])] <-  1
  }
  list(D = D, K = K, M = M)
}

# index m (within fuse_groups) where ||delta_{m+1} - delta_m|| > tol
change_points <- function(beta, group, fuse_groups, tol = 1e-4) {
  if (is.null(fuse_groups) || length(fuse_groups) < 2L) return(integer(0))
  gi <- group_cols_of(group)
  out <- integer(0)
  for (m in seq_len(length(fuse_groups) - 1L)) {
    d <- beta[gi$cols[[fuse_groups[m + 1L]]]] - beta[gi$cols[[fuse_groups[m]]]]
    if (sqrt(sum(d * d)) > tol) out <- c(out, m)
  }
  out
}

# sparsegl path + BIC-selected lambda
sgl_bic <- function(X, y, group, alpha, intercept, lambda = NULL, ...) {
  fit <- sparsegl::sparsegl(x = X, y = y, group = group, asparse = alpha,
                            intercept = intercept, ...)
  er <- sparsegl::estimate_risk(fit, X, type = "BIC", approx_df = TRUE)
  idx <- if (is.null(lambda)) which.min(er$BIC) else which.min(abs(fit$lambda - lambda))
  co <- as.numeric(stats::coef(fit, s = fit$lambda[idx]))
  if (length(co) == ncol(X) + 1L) {
    b0 <- co[1]; beta <- co[-1]
  } else {
    b0 <- 0;     beta <- co
  }
  list(beta = beta, b0 = b0, lambda = fit$lambda[idx],
       bic = er$BIC[idx], df = er$df[idx])
}

# ---------------------------------------------------------------------------

# Main fit.
#
# X, y, group: as in sparsegl.
# alpha:       sparse/group mix (passed as asparse).
# lambda:      single SGL lambda; if NULL, pick by BIC on the unfused path.
# lambda_fus:  fusion strength. 0 = plain sparsegl.
# fuse_groups: group indices (temporal order) that get the fusion penalty.
# rho, max_iter, tol: ADMM knobs.
fused_sgl <- function(X, y, group,
                      alpha       = 0.85,
                      lambda      = NULL,
                      lambda_fus  = 0,
                      fuse_groups = NULL,
                      intercept   = FALSE,
                      rho         = 1,
                      max_iter    = 500L,
                      tol         = 1e-6,
                      cp_tol      = 1e-4,
                      verbose     = FALSE,
                      ...) {
  n <- nrow(X); p <- ncol(X)
  group <- as.integer(group)
  
  # fast path: no fusion -> just sparsegl
  if (lambda_fus <= 0) {
    r <- sgl_bic(X, y, group, alpha, intercept, lambda = lambda, ...)
    fitted <- as.numeric(X %*% r$beta) + r$b0
    mse <- sum((y - fitted)^2) / n
    return(list(
      beta = r$beta, b0 = r$b0,
      lambda = r$lambda, lambda_fus = 0, alpha = alpha,
      bic = r$bic, mse = mse, df = r$df,
      iters = 0L, converged = TRUE,
      change_points = change_points(r$beta, group, fuse_groups, cp_tol)
    ))
  }
  
  if (intercept) stop("intercept=TRUE not supported with lambda_fus > 0; center y first")
  if (is.null(fuse_groups) || length(fuse_groups) < 2L)
    stop("need fuse_groups (length >= 2) when lambda_fus > 0")
  
  di <- build_D(p, group, fuse_groups)
  D <- di$D; K <- di$K; M <- di$M
  
  if (is.null(lambda)) {
    lambda <- sgl_bic(X, y, group, alpha, intercept = FALSE)$lambda
  }
  
  gi <- group_cols_of(group)
  bs <- gi$bs; bn <- gi$bn
  pf_group  <- sqrt(bs)       # sparsegl default
  pf_sparse <- rep(1, p)      # already sums to nvars
  
  # precompute per-group slices and Lipschitz constants
  Xg <- vector("list", bn)
  Dg <- vector("list", bn)
  Lg <- numeric(bn)
  for (g in seq_len(bn)) {
    cg <- gi$cols[[g]]
    Xg[[g]] <- X[, cg, drop = FALSE]
    Dg[[g]] <- D[, cg, drop = FALSE]
    A <- crossprod(Xg[[g]]) / n + rho * crossprod(Dg[[g]])
    Lg[g] <- if (length(cg) == 1L) as.numeric(A)
    else max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
  }
  Lg <- pmax(Lg, 1e-12)
  
  lama   <- alpha * lambda
  lam1ma <- (1 - alpha) * lambda
  
  # b-update: BCD on (1/(2n))||y-Xb||^2 + (rho/2)||Db - v||^2 + SGL(b).
  # Maintains r = y - Xb and q = Db - v incrementally.
  bcd_step <- function(beta, v, max_inner = 100L, inner_tol = 1e-9) {
    r <- as.numeric(y - X %*% beta)
    q <- as.numeric(D %*% beta - v)
    for (it in seq_len(max_inner)) {
      maxch <- 0
      for (g in seq_len(bn)) {
        cg <- gi$cols[[g]]
        b_old <- beta[cg]
        grad <- -as.numeric(crossprod(Xg[[g]], r)) / n +
          rho * as.numeric(crossprod(Dg[[g]], q))
        t_g <- 1 / Lg[g]
        s <- b_old - t_g * grad
        # L1 soft-threshold
        s <- sign(s) * pmax(abs(s) - t_g * lama * pf_sparse[cg], 0)
        # group soft-threshold
        snorm <- sqrt(sum(s * s))
        thr <- t_g * lam1ma * pf_group[g]
        b_new <- if (snorm > thr) s * (1 - thr / snorm) else numeric(length(cg))
        d <- b_new - b_old
        if (any(d != 0)) {
          beta[cg] <- b_new
          r <- r - as.numeric(Xg[[g]] %*% d)
          q <- q + as.numeric(Dg[[g]] %*% d)
          maxch <- max(maxch, Lg[g] * sum(d * d))
        }
      }
      if (maxch < inner_tol) break
    }
    beta
  }
  
  # ADMM
  beta <- rep(0, p)
  z <- rep(0, nrow(D))
  u <- rep(0, nrow(D))
  thresh <- lambda_fus / rho
  converged <- FALSE
  
  for (it in seq_len(max_iter)) {
    beta <- bcd_step(beta, z + u)
    Db <- as.numeric(D %*% beta)
    
    # z-update: group soft-threshold per K-block
    z_old <- z
    tgt <- Db - u
    z[] <- 0
    for (m in seq_len(M - 1L)) {
      idx <- ((m - 1L) * K + 1L):(m * K)
      nm <- sqrt(sum(tgt[idx]^2))
      if (nm > thresh) z[idx] <- tgt[idx] * (1 - thresh / nm)
    }
    
    u <- u + z - Db
    
    pri <- sqrt(sum((Db - z)^2))
    dua <- rho * sqrt(sum((z - z_old)^2))
    if (verbose && (it == 1 || it %% 10 == 0))
      cat(sprintf("iter %3d  pri=%.2e  dua=%.2e\n", it, pri, dua))
    if (pri < tol * sqrt(length(z) + 1) && dua < tol * sqrt(p + 1)) {
      converged <- TRUE
      break
    }
  }
  
  fitted <- as.numeric(X %*% beta)
  df <- sum(abs(beta) > 1e-8)
  mse <- sum((y - fitted)^2) / n
  
  list(
    beta = beta, b0 = 0,
    lambda = lambda, lambda_fus = lambda_fus, alpha = alpha,
    bic = log(mse) + log(n) * df / n,
    mse = mse, df = df,
    iters = it, converged = converged,
    change_points = change_points(beta, group, fuse_groups, cp_tol)
  )
}

# Grid search over (alpha, lambda, lambda_fus) by BIC.
fused_sgl_grid <- function(X, y, group,
                           alpha_grid      = c(0.5, 0.75, 0.9, 0.95),
                           lambda_fus_grid = c(0, 1e-3, 1e-2, 1e-1),
                           fuse_groups     = NULL,
                           intercept       = FALSE,
                           nlambda         = 30L,
                           verbose         = FALSE,
                           ...) {
  records <- list()
  best <- NULL
  best_bic <- Inf
  
  for (a in alpha_grid) {
    path <- sparsegl::sparsegl(x = X, y = y, group = group, asparse = a,
                               intercept = intercept, nlambda = nlambda)
    lam_grid <- path$lambda
    
    for (lf in lambda_fus_grid) {
      # lambda_fus = 0: let sgl_bic pick lambda directly
      if (lf == 0) {
        r <- fused_sgl(X, y, group, alpha = a, lambda = NULL,
                       lambda_fus = 0, fuse_groups = fuse_groups,
                       intercept = intercept, ...)
        records[[length(records) + 1L]] <-
          data.frame(alpha = a, lambda = r$lambda, lambda_fus = 0,
                     bic = r$bic, df = r$df)
        if (r$bic < best_bic) { best_bic <- r$bic; best <- r }
        next
      }
      for (lam in lam_grid) {
        r <- tryCatch(
          fused_sgl(X, y, group, alpha = a, lambda = lam,
                    lambda_fus = lf, fuse_groups = fuse_groups,
                    intercept = intercept, ...),
          error = function(e) NULL
        )
        if (is.null(r)) next
        records[[length(records) + 1L]] <-
          data.frame(alpha = a, lambda = lam, lambda_fus = lf,
                     bic = r$bic, df = r$df)
        if (verbose)
          cat(sprintf("a=%.2f lf=%.2e lam=%.4g  bic=%.4f\n", a, lf, lam, r$bic))
        if (r$bic < best_bic) { best_bic <- r$bic; best <- r }
      }
    }
  }
  best$grid <- do.call(rbind, records)
  best
}
