# =========================
# 0) parameters
# =========================
parse_char_vector <- function(x, default = character()) {
  if (is.null(x) || !nzchar(x)) return(default)
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  vals[nzchar(vals)]
}

parse_numeric_vector <- function(x, default = numeric()) {
  vals <- parse_char_vector(x, character())
  if (length(vals) == 0) return(default)
  out <- suppressWarnings(as.numeric(vals))
  out[!is.na(out)]
}

parse_integer_vector <- function(x, default = integer()) {
  vals <- parse_char_vector(x, character())
  if (length(vals) == 0) return(default)
  out <- suppressWarnings(as.integer(vals))
  out[!is.na(out)]
}

seed_start <- as.integer(Sys.getenv("SEED_START", unset = "42"))
n_reps <- as.integer(Sys.getenv("N_REPS", unset = "100"))

# simulation settings
n_grid <- parse_integer_vector(Sys.getenv("N_GRID", unset = "50"), default = c(50L))
obs_count_grid <- parse_integer_vector(Sys.getenv("OBS_COUNT_GRID", unset = "10"), default = c(10L))
obs_count_mode_grid <- parse_char_vector(
  Sys.getenv("OBS_COUNT_MODE", unset = "fixed"),
  default = c("fixed")
)
y_mode_grid <- parse_char_vector(
  Sys.getenv("Y_MODE_GRID", unset = "flr_linear,flr_square,flr_cube"),
  default = c("flr_linear", "flr_square", "flr_cube")
)
score_dist_grid <- parse_char_vector(
  Sys.getenv("SCORE_DIST_GRID", unset = "gaussian,laplace,t5,gamma"),
  default = c("gaussian", "laplace", "t5", "gamma")
)
x_noise_dist_grid <- parse_char_vector(
  Sys.getenv("X_NOISE_DIST_GRID", unset = "gaussian,t3,gamma"),
  default = c("gaussian", "t3", "gamma")
)
method_grid <- parse_char_vector(
  Sys.getenv("METHOD_GRID", unset = "pace_rkhs"),
  default = c("pace_rkhs")
)

# training / test sizes
test_n_factor <- 2

# true process settings
K_true <- 50
x_meas_sd <- 0.1
y_noise_sd <- 0.5

# dense grid used for internal numerical implementation
dense_grid_length <- 100

# fpca settings
M_fixed <- 8
fpca_n_reg_grid <- 100
pc_candidates <- seq_len(M_fixed)
pc_cv_folds <- 5
split_repeats <- 2
#split=what

# runtime controls
n_cores <- max(1L, as.integer(Sys.getenv("N_CORES", unset = "1")))

#l_H_grid <- c(0.40, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90,
#              1.00, 1.10, 1.20, 1.30, 1.45, 1.60)

#sigma_h_grid <- c(3.0, 4.0, 5.0, 6.0, 7.5, 9.0, 10.0, 12.5, 15.0, 20.0)

#sigma_m_grid <- c(0.001, 0.0015, 0.002, 0.003, 0.004, 0.006, 0.008, 0.012)

#lambda_grid <- c(5e-6, 7.5e-6, 1e-5, 1.5e-5, 2e-5, 3e-5,
#                 5e-5, 7.5e-5, 1e-4, 2e-4, 4e-4)
#l_H_grid <- c(0.15, 0.20, 0.30, 0.50)
#sigma_h_grid <- c(1.50, 3.00, 5.00)
#sigma_m_grid <- c(0, 0.0001,0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
#lambda_grid <- c(1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
l_H_grid <- c(0.70, 0.90, 1.10, 1.30, 1.60)
sigma_h_grid <- c(5.0, 7.5, 10.0, 15.0)
sigma_m_grid <- c(0.008, 0.012, 0.020, 0.035, 0.060)
lambda_grid <- c(5e-5, 8e-5, 1.2e-4, 2e-4, 4e-4, 8e-4)
# output paths
root_results_dir <- Sys.getenv(
  "RESULT_DIR",
  unset = file.path(getwd(), "results_fdr_all_methods_linear")
)

suppressPackageStartupMessages({
  library(fdapace)
  library(parallel)
  library(VGAM)
  library(splines)
})

# =========================
# 1) helpers
# =========================
sanitize_value <- function(x) {
  x_chr <- format(x, scientific = FALSE, trim = TRUE)
  x_chr <- gsub("\\s+", "", x_chr)
  gsub("\\.", "p", x_chr)
}

build_run_tag <- function(
  n,
  obs_count,
  obs_count_mode,
  y_mode,
  score_dist,
  x_noise_dist,
  method
) {
  paste0(
    "n_", sanitize_value(n),
    "_obs_", sanitize_value(obs_count),
    "_obsmode_", obs_count_mode,
    "_ymode_", y_mode,
    "_score_", score_dist,
    "_xnoise_", x_noise_dist,
    "_method_", method
  )
}

safe_solve_vec <- function(A, b) {
  as.vector(tryCatch(solve(A, b), error = function(e) qr.solve(A, b)))
}

safe_solve_mat <- function(A, B) {
  tryCatch(solve(A, B), error = function(e) qr.solve(A, B))
}

trapz_vec <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

cv_fold_ids <- function(n, k, seed) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n))
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) <= 1) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, na.rm = TRUE, type = 7))
}

serialize_numeric_vector <- function(x, digits = 17) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  paste(format(x, digits = digits, scientific = FALSE, trim = TRUE), collapse = ";")
}

parse_serialized_numeric_vector <- function(x) {
  if (length(x) == 0) return(numeric(0))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(numeric(0))
  vals <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0) return(numeric(0))
  out <- suppressWarnings(as.numeric(vals))
  out[is.finite(out)]
}

centered_gamma <- function(n, shape = 2, scale = 1) {
  g <- rgamma(n, shape = shape, scale = scale)
  (g - shape * scale) / sqrt(shape * scale^2)
}

sample_obs_count <- function(obs_count, obs_count_mode) {
  if (obs_count_mode == "fixed") {
    return(obs_count)
  }
  if (obs_count_mode == "random") {
    return(sample(obs_count:(2L * obs_count), size = 1L))
  }
  stop("Unsupported obs_count_mode.")
}

make_task_grid <- function() {
  expand.grid(
    n = n_grid,
    obs_count = obs_count_grid,
    obs_count_mode = obs_count_mode_grid,
    y_mode = y_mode_grid,
    score_dist = score_dist_grid,
    x_noise_dist = x_noise_dist_grid,
    method = method_grid,
    stringsAsFactors = FALSE
  )
}

# =========================
# 2) live csv append helpers
# =========================
acquire_lock <- function(lock_dir, timeout_sec = 300, sleep_sec = 0.1) {
  start_time <- Sys.time()
  repeat {
    ok <- dir.create(lock_dir, showWarnings = FALSE, recursive = FALSE)
    if (ok) return(TRUE)
    waited <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (waited > timeout_sec) {
      stop(sprintf("Failed to acquire lock: %s", lock_dir))
    }
    Sys.sleep(sleep_sec)
  }
}

release_lock <- function(lock_dir) {
  if (dir.exists(lock_dir)) unlink(lock_dir, recursive = TRUE, force = TRUE)
}

append_row_atomic <- function(df_row, csv_path) {
  lock_dir <- paste0(csv_path, ".lockdir")
  acquire_lock(lock_dir)
  on.exit(release_lock(lock_dir), add = TRUE)
  dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    df_row,
    file = csv_path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(csv_path),
    append = file.exists(csv_path),
    quote = TRUE
  )
}

write_manifest_if_needed <- function(task_grid) {
  manifest_path <- file.path(root_results_dir, "task_manifest.csv")
  if (!file.exists(manifest_path)) {
    write.csv(task_grid, manifest_path, row.names = FALSE)
  }
}

# =========================
# 3) basis functions and true coefficients
# =========================
phi_j <- function(t, j) {
  if (j == 1) {
    rep(1, length(t))
  } else {
    sqrt(2) * cos((j - 1) * pi * t)
  }
}

evaluate_basis_matrix <- function(t, K) {
  out <- vapply(seq_len(K), function(j) phi_j(t, j), numeric(length(t)))
  if (K == 1) out <- matrix(out, ncol = 1)
  out
}

lambda_vec_true <- function(K) {
  (seq_len(K))^(-2)
}

beta_coef_true <- function(K) {
  out <- numeric(K)
  out[1] <- 1
  if (K >= 2) out[2:K] <- 4 * (2:K)^(-3)
  out
}

compute_linear_signal_from_scores <- function(xi_mat, beta_coef) {
  as.vector(xi_mat %*% beta_coef)
}

apply_response_transform <- function(u, y_mode) {
  if (y_mode == "flr_linear") return(u)
  if (y_mode == "flr_square") return(u^2)
  if (y_mode == "flr_cube") return(u + u^2 + 0.001 * u^3)
  stop("Unsupported y_mode.")
}

# =========================
# 4) simulation data generation
# =========================
generate_score_matrix <- function(n, K, score_dist) {
  lam_true <- lambda_vec_true(K)
  xi_mat <- matrix(0, nrow = n, ncol = K)

  for (j in seq_len(K)) {
    if (score_dist == "gaussian") {
      xi_mat[, j] <- rnorm(n, mean = 0, sd = sqrt(lam_true[j]))
    } else if (score_dist == "laplace") {
      target_var <- sqrt(2) * j^(-2)
      b_j <- sqrt(target_var / 2)
      xi_mat[, j] <- VGAM::rlaplace(n, location = 0, scale = b_j)
    } else if (score_dist == "t5") {
      xi_mat[, j] <- j^(-1) * rt(n, df = 5) / sqrt(5 / 3)
    } else if (score_dist == "gamma") {
      xi_mat[, j] <- sqrt(lam_true[j]) * centered_gamma(n)
    } else {
      stop("Unsupported score_dist.")
    }
  }

  xi_mat
}

generate_measurement_noise <- function(m, x_noise_dist, x_meas_sd) {
  if (x_noise_dist == "gaussian") {
    return(rnorm(m, mean = 0, sd = x_meas_sd))
  }
  if (x_noise_dist == "t3") {
    return(rt(m, df = 3))
  }
  if (x_noise_dist == "gamma") {
    return(x_meas_sd * centered_gamma(m))
  }
  stop("Unsupported x_noise_dist.")
}

get_measurement_noise_variance <- function(x_noise_dist, x_meas_sd) {
  if (x_noise_dist == "gaussian") return(x_meas_sd^2)
  if (x_noise_dist == "t3") return(3)
  if (x_noise_dist == "gamma") return(x_meas_sd^2)
  stop("Unsupported x_noise_dist.")
}

generate_curve_sample <- function(
  n,
  obs_count,
  obs_count_mode,
  y_mode,
  score_dist,
  x_noise_dist
) {
  pts_grid <- seq(0, 1, length.out = dense_grid_length)
  beta_true_coef <- beta_coef_true(K_true)

  xi_mat <- generate_score_matrix(n, K_true, score_dist)
  u_true <- compute_linear_signal_from_scores(xi_mat, beta_true_coef)
  Y <- as.numeric(apply_response_transform(u_true, y_mode) + rnorm(n, sd = y_noise_sd))

  Lt <- vector("list", n)
  Ly <- vector("list", n)
  D_list_raw <- vector("list", n)

  for (i in seq_len(n)) {
    m_i <- sample_obs_count(obs_count, obs_count_mode)
    t_i <- sort(runif(m_i, 0, 1))
    x_true_i <- as.vector(evaluate_basis_matrix(t_i, K_true) %*% xi_mat[i, ])
    eps_i <- generate_measurement_noise(m_i, x_noise_dist, x_meas_sd)
    z_i <- x_true_i + eps_i
    Lt[[i]] <- t_i
    Ly[[i]] <- z_i
    D_list_raw[[i]] <- list(T = t_i, Z = z_i)
  }

  list(
    D_list_raw = D_list_raw,
    Lt = Lt,
    Ly = Ly,
    Y = Y,
    xi_mat = xi_mat,
    U_true = u_true,
    pts_grid = pts_grid,
    score_dist = score_dist,
    x_noise_dist = x_noise_dist,
    obs_count = obs_count,
    obs_count_mode = obs_count_mode,
    obs_count_vec = vapply(Lt, length, integer(1))
  )
}

generate_train_test_data <- function(
  n,
  obs_count,
  obs_count_mode,
  y_mode,
  score_dist,
  x_noise_dist
) {
  train <- generate_curve_sample(
    n,
    obs_count,
    obs_count_mode,
    y_mode,
    score_dist,
    x_noise_dist
  )
  test <- generate_curve_sample(
    test_n_factor * n,
    obs_count,
    obs_count_mode,
    y_mode,
    score_dist,
    x_noise_dist
  )
  list(train = train, test = test)
}

# =========================
# 5) embedding + KRR
# =========================
inner_prod_H2_normalized <- function(D1, D2, l_H) {
  T1 <- D1$T
  Z1 <- D1$Z
  T2 <- D2$T
  Z2 <- D2$Z
  M1 <- length(T1)
  M2 <- length(T2)
  all_T <- c(T1, T2)
  K_all <- exp(-(as.matrix(dist(all_T))^2) / (2 * l_H^2))
  K_inter <- K_all[seq_len(M1), M1 + seq_len(M2), drop = FALSE]
  (sum(K_inter) + sum((Z1 %*% t(Z2)) * K_inter)) / (M1 * M2)
}

compute_dist_H2_normalized <- function(D_list, l_H) {
  n0 <- length(D_list)
  self_ip <- vapply(D_list, function(d) inner_prod_H2_normalized(d, d, l_H), numeric(1))
  dist2 <- matrix(0, n0, n0)
  for (i in seq_len(n0)) {
    for (j in i:n0) {
      ip_ij <- inner_prod_H2_normalized(D_list[[i]], D_list[[j]], l_H)
      d2 <- max(0, self_ip[i] + self_ip[j] - 2 * ip_ij)
      dist2[i, j] <- d2
      dist2[j, i] <- d2
    }
  }
  dist2
}

compute_dist_M_inv <- function(D_list) {
  m_inv_vec <- vapply(D_list, function(d) 1 / length(d$T), numeric(1))
  as.matrix(dist(m_inv_vec))^2
}

calculate_gcv <- function(K, Y, lambda) {
  n0 <- length(Y)
  S <- K + n0 * lambda * diag(n0)
  A <- K %*% safe_solve_mat(S, diag(n0))
  resid <- Y - as.vector(A %*% Y)
  mse <- mean(resid^2)
  trA <- sum(diag(A))
  mse / (1 - trA / n0)^2
}

fit_krr_alpha <- function(K, Y, lambda) {
  n0 <- length(Y)
  safe_solve_vec(K + n0 * lambda * diag(n0), Y)
}

score_embedding_for_lh <- function(D_train, Y_train, dist2_M, l_H) {
  best_score <- Inf
  best_local <- NULL
  dist2_H <- compute_dist_H2_normalized(D_train, l_H)
  for (sigma_h in sigma_h_grid) {
    for (sigma_m in sigma_m_grid) {
      K <- exp(-dist2_H / (2 * sigma_h^2) - dist2_M / (2 * sigma_m^2))
      for (lambda in lambda_grid) {
        score <- calculate_gcv(K, Y_train, lambda)
        if (is.finite(score) && score < best_score) {
          best_score <- score
          best_local <- list(
            l_H = l_H,
            sigma_h = sigma_h,
            sigma_m = sigma_m,
            lambda = lambda,
            gcv = score,
            K_train = K
          )
        }
      }
    }
  }
  best_local
}

optimize_embedding_gcv <- function(D_train, Y_train) {
  dist2_M <- compute_dist_M_inv(D_train)
  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    candidate_list <- mclapply(
      l_H_grid,
      function(l_H) score_embedding_for_lh(D_train, Y_train, dist2_M, l_H),
      mc.cores = min(n_cores, length(l_H_grid))
    )
  } else {
    candidate_list <- lapply(
      l_H_grid,
      function(l_H) score_embedding_for_lh(D_train, Y_train, dist2_M, l_H)
    )
  }
  scores <- vapply(candidate_list, function(x) x$gcv, numeric(1))
  candidate_list[[which.min(scores)]]
}

predict_set_regression <- function(D_test, D_train, alpha_hat, params) {
  m_inv_train <- vapply(D_train, function(d) 1 / length(d$T), numeric(1))
  ip_train_self <- vapply(D_train, function(d) inner_prod_H2_normalized(d, d, params$l_H), numeric(1))
  vapply(D_test, function(d_star) {
    m_inv_star <- 1 / length(d_star$T)
    ip_star_self <- inner_prod_H2_normalized(d_star, d_star, params$l_H)
    k_star <- vapply(seq_along(D_train), function(i) {
      ip_cross <- inner_prod_H2_normalized(d_star, D_train[[i]], params$l_H)
      d2_h <- max(0, ip_star_self + ip_train_self[i] - 2 * ip_cross)
      d2_m <- (m_inv_star - m_inv_train[i])^2
      exp(-d2_h / (2 * params$sigma_h^2) - d2_m / (2 * params$sigma_m^2))
    }, numeric(1))
    sum(k_star * alpha_hat)
  }, numeric(1))
}

run_embedding_method <- function(train_dat, test_dat) {
  best <- optimize_embedding_gcv(train_dat$D_list_raw, train_dat$Y)
  alpha_hat <- fit_krr_alpha(best$K_train, train_dat$Y, best$lambda)
  pred_test <- predict_set_regression(
    test_dat$D_list_raw,
    train_dat$D_list_raw,
    alpha_hat,
    best
  )
  list(
    pred_mse = mean((test_dat$Y - pred_test)^2),
    l_H = best$l_H,
    sigma_h = best$sigma_h,
    sigma_m = best$sigma_m,
    lambda = best$lambda,
    gcv = best$gcv
  )
}

# =========================
# 6) fpca and main-aligned beta-method utilities
# =========================
main_sigma2_value <- 0.01

make_zero_mean_object <- function() {
  grid <- seq(0, 1, length.out = dense_grid_length)
  list(t = grid, mu = rep(0, length(grid)))
}

make_fpca_options <- function(method_xi = "IN", num_pc = M_fixed, sigma2_value = main_sigma2_value) {
  list(
    plot = FALSE,
    dataType = "Sparse",
    methodXi = method_xi,
    nRegGrid = fpca_n_reg_grid,
    maxK = M_fixed,
    FVEthreshold = 1,
    methodSelectK = num_pc,
    userMu = make_zero_mean_object(),
    userSigma2 = sigma2_value,
    usergrid = FALSE
  )
}

fit_fpca_model <- function(Ly, Lt, method_xi = "IN", num_pc = M_fixed, sigma2_value = main_sigma2_value) {
  FPCA(Ly, Lt, optns = make_fpca_options(
    method_xi = method_xi,
    num_pc = num_pc,
    sigma2_value = sigma2_value
  ))
}

align_phi_to_reference <- function(phi_mat, ref_phi_mat) {
  k_use <- min(ncol(phi_mat), ncol(ref_phi_mat))
  out <- phi_mat[, seq_len(k_use), drop = FALSE]
  ref <- ref_phi_mat[, seq_len(k_use), drop = FALSE]
  for (k in seq_len(k_use)) {
    if (sum((out[, k] - ref[, k])^2) > sum((out[, k] + ref[, k])^2)) {
      out[, k] <- -out[, k]
    }
  }
  out
}

project_scores_from_phi <- function(Ly_list, Lt_list, phi_grid, work_grid, num_pc) {
  n_sub <- length(Ly_list)
  out <- matrix(0, nrow = n_sub, ncol = num_pc)
  for (i in seq_len(n_sub)) {
    y_i <- as.numeric(Ly_list[[i]])
    t_i <- as.numeric(Lt_list[[i]])
    for (k in seq_len(num_pc)) {
      phi_i <- approx(work_grid, phi_grid[, k], xout = t_i, rule = 2)$y
      out[i, k] <- sum(y_i * phi_i) / length(y_i)
    }
  }
  colnames(out) <- paste0("x", seq_len(num_pc))
  out
}

fit_linear_score_model <- function(score_mat, Y_train) {
  df_train <- data.frame(Y = Y_train, score_mat)
  lm(Y ~ ., data = df_train)
}

beta_from_fpca_lm <- function(fpca_obj, lm_fit, num_pc) {
  coef_full <- as.numeric(coef(lm_fit)[-1])
  if (length(coef_full) < num_pc) {
    coef_full <- c(coef_full, rep(0, num_pc - length(coef_full)))
  }
  coef_vec <- coef_full[seq_len(num_pc)]
  phi_grid <- fpca_obj$phi[, seq_len(num_pc), drop = FALSE]
  as.vector(phi_grid %*% coef_vec)
}

compute_u_from_beta <- function(beta_grid, T_list, Z_list, pts_grid) {
  n_sub <- length(T_list)
  out <- numeric(n_sub)
  for (i in seq_len(n_sub)) {
    beta_i <- approx(pts_grid, beta_grid, xout = T_list[[i]], rule = 2)$y
    out[i] <- sum(beta_i * Z_list[[i]]) / length(Z_list[[i]])
  }
  out
}

nearest_grid_index <- function(grid, x) {
  idx_right <- findInterval(x, grid)
  idx_right[idx_right < 1L] <- 1L
  idx_right[idx_right >= length(grid)] <- length(grid) - 1L
  idx_left <- idx_right
  idx_left[idx_left < 1L] <- 1L
  use_right <- abs(x - grid[idx_right]) > abs(x - grid[idx_right + 1L])
  out <- idx_left
  out[use_right] <- idx_right[use_right] + 1L
  out[x <= grid[1L]] <- 1L
  out[x >= grid[length(grid)]] <- length(grid)
  out
}

GetSmoothedMeanCurve_mainlike <- function(y, t, obsGrid, regGrid, optns){
  userMu = optns$userMu
  methodBwMu = optns$methodBwMu
  npoly = 1
  nder = 0
  userBwMu = optns$userBwMu
  kernel = optns$kernel

  if (is.list(userMu) && (length(userMu$mu) == length(userMu$t))) {
    buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
    rangeUser <- range(optns$userMu$t)
    rangeObs <- range(obsGrid)
    if (rangeUser[1] > rangeObs[1] + buff || rangeUser[2] < rangeObs[2] - buff) {
      stop('The range defined by the user provided mean does not cover the support of the data.')
    }
    mu <- spline(userMu$t, userMu$mu, xout = obsGrid)$y
    muDense <- spline(obsGrid, mu, xout = regGrid)$y
    bw_mu <- NULL
  } else {
    if (userBwMu > 0) {
      bw_mu <- userBwMu
    } else {
      if (any(methodBwMu == c('GCV', 'GMeanAndGCV'))) {
        bw_mu <- unlist(fdapace:::GCVLwls1D1(yy = y, tt = t, kernel = kernel, npoly = npoly, nder = nder, dataType = optns$dataType))[1]
        if (0 == length(bw_mu)) {
          stop('The data is too sparse to estimate a mean function. Get more data!')
        }
        if (methodBwMu == 'GMeanAndGCV') {
          minbw <- fdapace:::Minb(unlist(t), 2)
          bw_mu <- sqrt(minbw * bw_mu)
        }
      } else {
        bw_mu <- fdapace:::CVLwls1D(y, t, kernel = kernel, npoly = npoly, nder = nder, dataType = optns$dataType, kFolds = optns$kFoldMuCov, useBW1SE = optns$useBW1SE)
      }
    }
    xin <- unlist(t)
    yin <- unlist(y)[order(xin)]
    xin <- sort(xin)
    win <- rep(1, length(xin))
    mu <- fdapace:::Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin = yin, xout = obsGrid, win = win)
    muDense <- fdapace:::Lwls1D(bw_mu, kernel_type = kernel, npoly = npoly, nder = nder, xin = xin, yin = yin, xout = regGrid, win = win)
  }

  result <- list(mu = mu, muDense = muDense, bw_mu = bw_mu)
  class(result) <- 'SMC'
  result
}

FPCA_CE_mainlike <- function(Ly, Lt, optns = list()) {
  GetCEScores <- function(y, t, optns, mu, obsGrid, fittedCov, lambda, phi, sigma2) {
    GetMuPhiSig <- function(t, obsGrid, mu, phi, Sigma_Y) {
      lapply(t, function(tvec) {
        if (length(tvec) != 0) {
          muVec <- approx(obsGrid, mu, tvec)$y
          phiMat <- matrix(apply(phi, 2, function(phivec) approx(obsGrid, phivec, tvec)$y), nrow = length(tvec))
          idx <- nearest_grid_index(obsGrid, tvec)
          Sigma_Yi <- matrix(Sigma_Y[idx, idx, drop = FALSE], length(tvec), length(tvec))
          list(muVec = muVec, phiMat = phiMat, Sigma_Yi = Sigma_Yi)
        } else {
          list(muVec = numeric(0), phiMat = numeric(0), Sigma_Yi = numeric(0))
        }
      })
    }

    GetIndCEScores <- function(yVec, muVec, lamVec, phiMat, Sigma_Yi, newyInd = NULL, verbose = FALSE) {
      if (length(yVec) == 0) {
        if (verbose) warning('Empty observation found, possibly due to truncation')
        return(list(xiEst = matrix(NA, length(lamVec)), xiVar = matrix(NA, length(lamVec), length(lamVec)), fittedY = matrix(NA, 0, 0)))
      }
      if (!is.null(newyInd)) {
        if (length(yVec) != 1) {
          newPhi <- phiMat[newyInd, , drop = FALSE]
          newMu <- muVec[newyInd]
          yVec <- yVec[-newyInd]
          muVec <- muVec[-newyInd]
          phiMat <- phiMat[-newyInd, , drop = FALSE]
          Sigma_Yi <- Sigma_Yi[-newyInd, -newyInd, drop = FALSE]
          return(fdapace:::GetIndCEScoresCPPnewInd(yVec, muVec, lamVec, phiMat, Sigma_Yi, newPhi, newMu))
        }
        Lam <- diag(x = lamVec, nrow = length(lamVec))
        LamPhi <- Lam %*% t(phiMat)
        LamPhiSig <- LamPhi %*% solve(Sigma_Yi)
        xiEst <- LamPhiSig %*% matrix(yVec - muVec, ncol = 1)
        xiVar <- Lam - LamPhi %*% t(LamPhiSig)
        return(list(xiEst = xiEst, xiVar = xiVar, fittedY = NA))
      }
      fdapace:::GetIndCEScoresCPP(yVec, muVec, lamVec, phiMat, Sigma_Yi)
    }

    if (length(lambda) != ncol(phi)) stop('No of eigenvalues is not the same as the no of eigenfunctions.')
    if (is.null(sigma2)) sigma2 <- 0
    Sigma_Y <- fittedCov + diag(sigma2, nrow(phi))
    MuPhiSig <- GetMuPhiSig(t, obsGrid, mu, phi, Sigma_Y)
    mapply(function(yVec, muphisig){
      GetIndCEScores(yVec, muphisig$muVec, lambda, muphisig$phiMat, muphisig$Sigma_Yi, verbose = optns$verbose)
    }, y, MuPhiSig)
  }

  firsttsFPCA <- Sys.time()
  fdapace:::CheckData(Ly, Lt)
  inputData <- fdapace:::HandleNumericsAndNAN(Ly, Lt)
  Ly <- inputData$Ly
  Lt <- inputData$Lt
  optns <- fdapace:::SetOptions(Ly, Lt, optns)
  numOfCurves <- length(Ly)
  fdapace:::CheckOptions(Lt, optns, numOfCurves)
  if (optns$usergrid == FALSE & optns$useBinnedData != 'OFF') {
    BinnedDataset <- fdapace:::GetBinnedDataset(Ly, Lt, optns)
    Ly <- BinnedDataset$newy
    Lt <- BinnedDataset$newt
    optns[['nRegGrid']] <- min(optns[['nRegGrid']], BinnedDataset[['numBins']])
    inputData$Ly <- Ly
    inputData$Lt <- Lt
  }
  obsGrid <- sort(unique(c(unlist(Lt))))
  regGrid <- seq(min(obsGrid), max(obsGrid), length.out = optns$nRegGrid)
  outPercent <- optns$outPercent
  buff <- .Machine$double.eps * max(abs(obsGrid)) * 10
  rangeGrid <- range(regGrid)
  minGrid <- rangeGrid[1]
  cutRegGrid <- regGrid[regGrid > minGrid + diff(rangeGrid) * outPercent[1] - buff & regGrid < minGrid + diff(rangeGrid) * outPercent[2] + buff]
  ymat <- fdapace:::List2Mat(Ly, Lt)
  firsttsMu <- Sys.time()
  userMu <- optns$userMu
  if (is.list(userMu) && (length(userMu$mu) == length(userMu$t))) {
    smcObj <- fdapace:::GetUserMeanCurve(optns, obsGrid, regGrid, buff)
    smcObj$muDense <- fdapace:::ConvertSupport(obsGrid, regGrid, mu = smcObj$mu)
  } else if (optns$methodMuCovEst == 'smooth') {
    smcObj <- GetSmoothedMeanCurve_mainlike(Ly, Lt, obsGrid, regGrid, optns)
  } else if (optns$methodMuCovEst == 'cross-sectional') {
    smcObj <- fdapace:::GetMeanDense(ymat, obsGrid, optns)
  }
  mu <- smcObj$mu
  lasttsMu <- Sys.time()
  firsttsCov <- Sys.time()
  if (!is.null(optns$userCov) && optns$methodMuCovEst != 'smooth') {
    scsObj <- fdapace:::GetUserCov(optns, obsGrid, cutRegGrid, buff, ymat)
  } else if (optns$methodMuCovEst == 'smooth') {
    scsObj <- fdapace:::GetSmoothedCovarSurface(Ly, Lt, mu, obsGrid, regGrid, optns, optns$useBinnedCov)
  } else if (optns$methodMuCovEst == 'cross-sectional') {
    scsObj <- fdapace:::GetCovDense(ymat, mu, optns)
    if (length(obsGrid) != length(cutRegGrid) || !identical(obsGrid, cutRegGrid)) {
      scsObj$smoothCov <- fdapace:::ConvertSupport(obsGrid, cutRegGrid, Cov = scsObj$smoothCov)
    }
    scsObj$outGrid <- cutRegGrid
  }
  sigma2 <- scsObj[['sigma2']]
  lasttsCov <- Sys.time()
  firsttsPACE <- Sys.time()
  workGrid <- scsObj$outGrid
  muWork <- fdapace:::ConvertSupport(obsGrid, toGrid = workGrid, mu = smcObj$mu)
  eigObj <- fdapace:::GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, workGrid, optns, muWork = muWork)
  truncObsGrid <- obsGrid
  if (!all(abs(optns$outPercent - c(0, 1)) < .Machine$double.eps * 2)) {
    truncObsGrid <- truncObsGrid[truncObsGrid >= min(workGrid) - buff & truncObsGrid <= max(workGrid) + buff]
    tmp <- fdapace:::TruncateObs(Ly, Lt, truncObsGrid)
    Ly <- tmp$Ly
    Lt <- tmp$Lt
  }
  muObs <- fdapace:::ConvertSupport(obsGrid, truncObsGrid, mu = mu)
  phiObs <- fdapace:::ConvertSupport(workGrid, truncObsGrid, phi = eigObj$phi)
  if (optns$methodXi == 'CE') {
    CovObs <- fdapace:::ConvertSupport(workGrid, truncObsGrid, Cov = eigObj$fittedCov)
  }
  if (optns$methodXi == 'CE') {
    if (optns$methodRho != 'vanilla') {
      if (is.null(optns$userRho)) {
        if (length(Ly) > 2048) {
          randIndx <- sample(length(Ly), 2048)
          rho <- fdapace:::GetRho(Ly[randIndx], Lt[randIndx], optns, muObs, muWork, truncObsGrid, CovObs, eigObj$lambda, phiObs, eigObj$phi, workGrid, sigma2)
        } else {
          rho <- fdapace:::GetRho(Ly, Lt, optns, muObs, muWork, truncObsGrid, CovObs, eigObj$lambda, phiObs, eigObj$phi, workGrid, sigma2)
        }
      } else {
        rho <- optns$userRho
      }
      sigma2 <- rho
    }
    scoresObj <- GetCEScores(Ly, Lt, optns, muObs, truncObsGrid, CovObs, eigObj$lambda, phiObs, sigma2)
  } else if (optns$methodXi == 'IN') {
    scoresObj <- mapply(function(yvec, tvec) fdapace:::GetINScores(yvec, tvec, optns = optns, obsGrid, mu = muObs, lambda = eigObj$lambda, phi = phiObs, sigma2 = sigma2), Ly, Lt)
  }
  if (optns$fitEigenValues) {
    fitLambda <- fdapace:::FitEigenValues(scsObj$rcov, workGrid, eigObj$phi, optns$maxK)
  } else {
    fitLambda <- NULL
  }
  lasttsPACE <- Sys.time()
  ret <- fdapace:::MakeResultFPCA(optns, smcObj, muObs, scsObj, eigObj, inputData = inputData, scoresObj, truncObsGrid, workGrid, rho = if (optns$methodRho != 'vanilla') rho else 0, fitLambda = fitLambda, timestamps = c(lasttsMu, lasttsCov, lasttsPACE, firsttsFPCA, firsttsMu, firsttsCov, firsttsPACE))
  if (optns$plot) plot.FPCA(ret)
  ret
}

smooth_mean_curve_plugin <- function(Ly, Lt, reg_grid, num_pc, sigma2_value = main_sigma2_value) {
  obs_grid <- sort(unique(unlist(Lt, use.names = FALSE)))
  fdapace_fit <- tryCatch({
    opt <- fdapace:::SetOptions(
      Ly,
      Lt,
      list(
        plot = FALSE,
        dataType = "Sparse",
        methodXi = "IN",
        nRegGrid = length(reg_grid),
        maxK = M_fixed,
        FVEthreshold = 1,
        methodSelectK = num_pc,
        userMu = make_zero_mean_object(),
        userSigma2 = sigma2_value,
        usergrid = FALSE
      )
    )
    fdapace:::GetSmoothedMeanCurve(Ly, Lt, obs_grid, reg_grid, opt)
  }, error = function(e) NULL)
  if (!is.null(fdapace_fit) && !is.null(fdapace_fit$muDense)) {
    return(as.numeric(fdapace_fit$muDense))
  }
  x <- unlist(Lt, use.names = FALSE)
  y <- unlist(Ly, use.names = FALSE)
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  fit <- tryCatch(smooth.spline(x = x, y = y, cv = TRUE), error = function(e) smooth.spline(x = x, y = y, df = min(10L, max(4L, floor(length(unique(x)) / 4)))))
  as.numeric(predict(fit, x = reg_grid)$y)
}

integrate_scores_from_curve <- function(curve_grid, phi_grid, work_grid, num_pc) {
  vapply(
    seq_len(num_pc),
    function(k) trapz_vec(work_grid, curve_grid * phi_grid[, k]),
    numeric(1)
  )
}

CVr_main <- function(train_dat, Y, kfold = pc_cv_folds, method = "IN", S = split_repeats) {
  CV <- rep(0, M_fixed)
  CV_In <- rep(0, M_fixed)
  CV_split <- rep(0, M_fixed)
  CV_pi <- rep(0, M_fixed)
  CV_pace <- rep(0, M_fixed)
  n <- length(Y)
  L <- floor(n / kfold)
  if (L < 1) stop("Too few samples for k-fold CV.")
  pts_grid <- seq(0, 1, length.out = dense_grid_length)
  Mu <- make_zero_mean_object()

  for (v in seq_len(kfold)) {
    indout <- ((v - 1) * L + 1):(v * L)
    indin <- sort(setdiff(seq_len(n), indout))
    Y_in <- Y[indin]
    Y_out <- Y[indout]
    obsLt_in <- train_dat$Lt[indin]
    obsLy_in <- train_dat$Ly[indin]
    obsLt_out <- train_dat$Lt[indout]
    obsLy_out <- train_dat$Ly[indout]
    sampX_in <- list(Ly = obsLy_in, Lt = obsLt_in)
    nin <- length(indin)

    for (m in seq_len(M_fixed)) {
      Ly_w <- sampX_in$Ly
      for (i in seq_len(nin)) Ly_w[[i]] <- (Y_in[i] - mean(Y_in)) * sampX_in$Ly[[i]]
      g_dense <- smooth_mean_curve_plugin(Ly_w, sampX_in$Lt, pts_grid, m, main_sigma2_value)

      res_in <- fit_fpca_model(sampX_in$Ly, sampX_in$Lt, method_xi = "IN", num_pc = m, sigma2_value = main_sigma2_value)
      Phi_est <- res_in$phi[, seq_len(m), drop = FALSE]
      G <- integrate_scores_from_curve(g_dense, Phi_est, res_in$workGrid, m)
      b <- as.vector(G) / as.numeric(res_in$lambda[seq_len(m)])
      beta_pi <- as.vector(Phi_est %*% b)

      Xi_1 <- matrix(0, nrow = nin, ncol = m)
      for (i in seq_len(nin)) {
        for (k in seq_len(m)) {
          Xi_1[i, k] <- sum(sampX_in$Ly[[i]] * approx(res_in$workGrid, Phi_est[, k], xout = sampX_in$Lt[[i]], rule = 2)$y) / length(sampX_in$Ly[[i]])
        }
      }
      f1 <- lm(Y_in ~ Xi_1)
      beta_1 <- as.vector(Phi_est %*% as.vector(f1$coef[-1]))

      f_In <- lm(Y_in ~ res_in$xiEst[, seq_len(m), drop = FALSE])
      beta_In <- as.vector(Phi_est %*% as.vector(f_In$coef[-1]))

      res_pace <- FPCA_CE_mainlike(sampX_in$Ly, sampX_in$Lt, optns = list(nRegGrid = fpca_n_reg_grid, methodXi = "CE", maxK = M_fixed, FVEthreshold = 1, userMu = Mu, methodSelectK = m, userSigma2 = main_sigma2_value, usergrid = FALSE))
      Phi_pace <- res_pace$phi[, seq_len(m), drop = FALSE]
      f_pace <- lm(Y_in ~ res_pace$xiEst[, seq_len(m), drop = FALSE])
      beta_pace <- as.vector(Phi_pace %*% as.vector(f_pace$coef[-1]))

      beta_2s <- matrix(0, nrow = S, ncol = length(beta_1))
      a_vec <- numeric(S)
      for (s in seq_len(S)) {
        Index_1 <- sort(sample.int(nin, floor(nin / 2)))
        Index_2 <- sort(setdiff(seq_len(nin), Index_1))
        sampX_1 <- list(Lt = sampX_in$Lt[Index_1], Ly = sampX_in$Ly[Index_1])
        sampX_2 <- list(Lt = sampX_in$Lt[Index_2], Ly = sampX_in$Ly[Index_2])
        Y_1 <- Y_in[Index_1]; Y_2 <- Y_in[Index_2]
        res1 <- fit_fpca_model(sampX_1$Ly, sampX_1$Lt, method_xi = method, num_pc = m, sigma2_value = main_sigma2_value)
        res2 <- fit_fpca_model(sampX_2$Ly, sampX_2$Lt, method_xi = method, num_pc = m, sigma2_value = main_sigma2_value)
        Phi_est21 <- res1$phi[, seq_len(m), drop = FALSE]
        Phi_est22 <- align_phi_to_reference(res2$phi[, seq_len(m), drop = FALSE], Phi_est21)
        Xi_21 <- matrix(0, nrow = length(Index_1), ncol = m)
        Xi_22 <- matrix(0, nrow = length(Index_2), ncol = m)
        for (i in seq_along(Index_1)) {
          for (k in seq_len(m)) {
            Xi_21[i, k] <- sum(sampX_in$Ly[[Index_1[i]]] * approx(res2$workGrid, Phi_est22[, k], xout = sampX_in$Lt[[Index_1[i]]], rule = 2)$y) / length(sampX_in$Ly[[Index_1[i]]])
          }
        }
        for (i in seq_along(Index_2)) {
          for (k in seq_len(m)) {
            Xi_22[i, k] <- sum(sampX_in$Ly[[Index_2[i]]] * approx(res1$workGrid, Phi_est21[, k], xout = sampX_in$Lt[[Index_2[i]]], rule = 2)$y) / length(sampX_in$Ly[[Index_2[i]]])
          }
        }
        Xi_3 <- rbind(Xi_21, Xi_22)
        Y_3 <- c(Y_1, Y_2)
        f3 <- lm(Y_3 ~ Xi_3)
        a_vec[s] <- as.numeric(f3$coefficients[1])
        beta_2s[s, ] <- as.vector(0.5 * (Phi_est21 + Phi_est22) %*% as.vector(f3$coef[-1]))
      }
      a_f <- mean(a_vec)
      beta_2 <- colMeans(beta_2s)

      for (i in seq_along(indout)) {
        temp1 <- approx(res_in$workGrid, beta_1, obsLt_out[[i]], rule = 2)$y
        temp2 <- approx(pts_grid, beta_2, obsLt_out[[i]], rule = 2)$y
        temp3 <- approx(res_in$workGrid, beta_pi, obsLt_out[[i]], rule = 2)$y
        temp4 <- approx(res_pace$workGrid, beta_pace, obsLt_out[[i]], rule = 2)$y
        temp5 <- approx(res_in$workGrid, beta_In, obsLt_out[[i]], rule = 2)$y
        denom_i <- length(obsLy_out[[i]])
        CV[m] <- CV[m] + (Y_out[i] - as.numeric(f1$coefficients[1]) - sum(obsLy_out[[i]] * temp1) / denom_i)^2
        CV_split[m] <- CV_split[m] + (Y_out[i] - a_f - sum(obsLy_out[[i]] * temp2) / denom_i)^2
        CV_pace[m] <- CV_pace[m] + (Y_out[i] - as.numeric(f_pace$coefficients[1]) - sum(obsLy_out[[i]] * temp4) / denom_i)^2
        CV_In[m] <- CV_In[m] + (Y_out[i] - as.numeric(f_In$coefficients[1]) - sum(obsLy_out[[i]] * temp5) / denom_i)^2
        CV_pi[m] <- CV_pi[m] + (Y_out[i] - sum(obsLy_out[[i]] * temp3) / denom_i)^2
      }
    }
  }
  list(
    CV = which.min(CV),
    CV_split = which.min(CV_split),
    CV_In = which.min(CV_In),
    CV_pi = which.min(CV_pi),
    CV_pace = which.min(CV_pace),
    CV_values = CV,
    CV_split_values = CV_split,
    CV_In_values = CV_In,
    CV_pi_values = CV_pi,
    CV_pace_values = CV_pace
  )
}


# =========================
# 7) PACE-RKHS method from Avery et al. (2014)
# =========================
rkhs_rho_multipliers <- c(0.25, 0.50, 0.75, 1.00, 1.25, 1.50)
rkhs_lambda_grid <- c(1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2)
rkhs_cv_folds <- 5

get_valid_score_matrix <- function(score_mat, num_pc = NULL) {
  score_mat <- as.matrix(score_mat)
  if (!is.null(num_pc)) {
    score_mat <- score_mat[, seq_len(num_pc), drop = FALSE]
  }

  finite_cols <- apply(score_mat, 2, function(x) all(is.finite(x)))
  if (!all(finite_cols)) {
    leading_ok <- cumprod(as.integer(finite_cols)) == 1L
    max_leading <- sum(leading_ok)
    if (max_leading < 1L) {
      stop("No finite leading PACE score columns are available for pace_rkhs.")
    }
    score_mat <- score_mat[, seq_len(max_leading), drop = FALSE]
  }

  if (ncol(score_mat) < 1L) {
    stop("No finite PACE score columns are available for pace_rkhs.")
  }
  score_mat
}

get_fpca_mu_on_work_grid <- function(fpca_obj, work_grid) {
  if (!is.null(fpca_obj$mu) && length(fpca_obj$mu) == length(work_grid)) {
    return(as.numeric(fpca_obj$mu))
  }
  if (!is.null(fpca_obj$mu) && !is.null(fpca_obj$obsGrid)) {
    return(as.numeric(approx(
      fpca_obj$obsGrid,
      fpca_obj$mu,
      xout = work_grid,
      rule = 2
    )$y))
  }
  rep(0, length(work_grid))
}

predict_ce_scores_from_fpca <- function(fpca_obj, Ly_list, Lt_list, num_pc) {
  work_grid <- fpca_obj$workGrid
  if (is.null(work_grid)) {
    stop("FPCA object does not contain workGrid.")
  }
  phi <- as.matrix(fpca_obj$phi[, seq_len(num_pc), drop = FALSE])
  lambda <- as.numeric(fpca_obj$lambda[seq_len(num_pc)])
  mu_grid <- get_fpca_mu_on_work_grid(fpca_obj, work_grid)
  sigma2 <- fpca_obj$sigma2
  if (is.null(sigma2) || !is.finite(sigma2)) {
    sigma2 <- main_sigma2_value
  }

  n_sub <- length(Ly_list)
  out <- matrix(NA_real_, nrow = n_sub, ncol = num_pc)
  for (i in seq_len(n_sub)) {
    y_i <- as.numeric(Ly_list[[i]])
    t_i <- as.numeric(Lt_list[[i]])
    mu_i <- as.numeric(approx(work_grid, mu_grid, xout = t_i, rule = 2)$y)
    phi_i <- matrix(
      apply(phi, 2, function(phi_k) {
        approx(work_grid, phi_k, xout = t_i, rule = 2)$y
      }),
      nrow = length(t_i),
      ncol = num_pc
    )
    cov_i <- phi_i %*% diag(lambda, nrow = num_pc) %*% t(phi_i)
    cov_i <- cov_i + sigma2 * diag(length(t_i))
    out[i, ] <- as.vector(
      diag(lambda, nrow = num_pc) %*% t(phi_i) %*%
        safe_solve_vec(cov_i, y_i - mu_i)
    )
  }
  colnames(out) <- paste0("x", seq_len(num_pc))
  out
}

reconstruct_pace_curves <- function(fpca_obj, score_mat, num_pc) {
  work_grid <- fpca_obj$workGrid
  if (is.null(work_grid)) {
    stop("FPCA object does not contain workGrid.")
  }

  score_mat <- as.matrix(score_mat[, seq_len(num_pc), drop = FALSE])
  phi <- as.matrix(fpca_obj$phi[, seq_len(num_pc), drop = FALSE])
  mu_grid <- get_fpca_mu_on_work_grid(fpca_obj, work_grid)

  curve_mat <- score_mat %*% t(phi)
  curve_mat <- sweep(curve_mat, 2, mu_grid, FUN = "+")
  curve_mat
}

trapezoid_weights <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2L) {
    stop("At least two grid points are required for L2 integration.")
  }

  w <- numeric(length(grid))
  dx <- diff(grid)
  w[1L] <- dx[1L] / 2
  w[length(grid)] <- dx[length(dx)] / 2
  if (length(grid) > 2L) {
    w[2L:(length(grid) - 1L)] <- (dx[-length(dx)] + dx[-1L]) / 2
  }
  w
}

pairwise_l2_curve_dist <- function(curve_1, grid, curve_2 = NULL) {
  curve_1 <- as.matrix(curve_1)
  weights <- trapezoid_weights(grid)
  if (ncol(curve_1) != length(weights)) {
    stop("curve_1 and grid have incompatible dimensions.")
  }

  if (is.null(curve_2)) {
    norm_1 <- as.vector((curve_1^2) %*% weights)
    cross <- (curve_1 * matrix(weights, nrow = nrow(curve_1), ncol = ncol(curve_1), byrow = TRUE)) %*% t(curve_1)
    D2 <- outer(norm_1, norm_1, "+") - 2 * cross
  } else {
    curve_2 <- as.matrix(curve_2)
    if (ncol(curve_2) != length(weights)) {
      stop("curve_2 and grid have incompatible dimensions.")
    }
    norm_1 <- as.vector((curve_1^2) %*% weights)
    norm_2 <- as.vector((curve_2^2) %*% weights)
    cross <- (curve_1 * matrix(weights, nrow = nrow(curve_1), ncol = ncol(curve_1), byrow = TRUE)) %*% t(curve_2)
    D2 <- outer(norm_1, norm_2, "+") - 2 * cross
  }

  pmax(D2, 0)
}

gaussian_functional_kernel <- function(D2, rho) {
  exp(-D2 / (rho^2))
}

fit_intercept_krr <- function(K, Y, lambda) {
  n0 <- length(Y)
  block_mat <- rbind(
    c(0, rep(1, n0)),
    cbind(rep(1, n0), K + n0 * lambda * diag(n0))
  )
  rhs <- c(0, Y)
  sol <- safe_solve_vec(block_mat, rhs)

  list(
    intercept = as.numeric(sol[1L]),
    alpha = as.numeric(sol[-1L])
  )
}

predict_intercept_krr <- function(K_test_train, fit) {
  fit$intercept + as.vector(K_test_train %*% fit$alpha)
}

select_pace_rkhs_rho_lambda <- function(curve_train, work_grid, Y_train, seed) {
  n0 <- length(Y_train)
  fold_ids <- cv_fold_ids(n0, rkhs_cv_folds, seed)
  D2_all <- pairwise_l2_curve_dist(curve_train, work_grid)
  dist_vals <- sqrt(D2_all[upper.tri(D2_all)])
  med_dist <- stats::median(
    dist_vals[is.finite(dist_vals) & dist_vals > 0],
    na.rm = TRUE
  )
  if (!is.finite(med_dist) || med_dist <= 0) {
    med_dist <- 1
  }
  rho_grid <- rkhs_rho_multipliers * med_dist

  best_cv <- Inf
  best_rho <- NA_real_
  best_lambda <- NA_real_

  for (rho in rho_grid) {
    K_all <- gaussian_functional_kernel(D2_all, rho)
    for (lambda in rkhs_lambda_grid) {
      fold_mse <- numeric(rkhs_cv_folds)
      for (fold in seq_len(rkhs_cv_folds)) {
        val_idx <- which(fold_ids == fold)
        tr_idx <- setdiff(seq_len(n0), val_idx)
        K_tr <- K_all[tr_idx, tr_idx, drop = FALSE]
        K_val <- K_all[val_idx, tr_idx, drop = FALSE]
        fit <- fit_intercept_krr(K_tr, Y_train[tr_idx], lambda)
        pred_val <- predict_intercept_krr(K_val, fit)
        fold_mse[fold] <- mean((Y_train[val_idx] - pred_val)^2)
      }
      cv_mse <- mean(fold_mse)
      if (is.finite(cv_mse) && cv_mse < best_cv) {
        best_cv <- cv_mse
        best_rho <- rho
        best_lambda <- lambda
      }
    }
  }

  list(
    rho = best_rho,
    lambda = best_lambda,
    cv_mse = best_cv,
    median_pairwise_distance = med_dist
  )
}

run_pace_rkhs_method <- function(train_dat, test_dat, seed) {
  fpca_fit <- FPCA_CE_mainlike(
    train_dat$Ly,
    train_dat$Lt,
    optns = list(
      nRegGrid = fpca_n_reg_grid,
      methodXi = "CE",
      maxK = M_fixed,
      FVEthreshold = 1,
      userMu = make_zero_mean_object(),
      methodSelectK = "BIC",
      userSigma2 = main_sigma2_value,
      usergrid = FALSE
    )
  )

  num_pc <- min(
    M_fixed,
    ncol(as.matrix(fpca_fit$xiEst)),
    ncol(as.matrix(fpca_fit$phi)),
    length(fpca_fit$lambda)
  )
  if (!is.finite(num_pc) || num_pc < 1L) {
    stop("BIC selected no valid FPCA components for pace_rkhs.")
  }

  score_train <- get_valid_score_matrix(fpca_fit$xiEst, num_pc)
  num_pc <- ncol(score_train)
  score_test <- predict_ce_scores_from_fpca(
    fpca_fit,
    test_dat$Ly,
    test_dat$Lt,
    num_pc
  )
  score_test <- score_test[, seq_len(num_pc), drop = FALSE]

  curve_train <- reconstruct_pace_curves(fpca_fit, score_train, num_pc)
  curve_test <- reconstruct_pace_curves(fpca_fit, score_test, num_pc)
  work_grid <- fpca_fit$workGrid

  tune <- select_pace_rkhs_rho_lambda(
    curve_train,
    work_grid,
    train_dat$Y,
    seed
  )

  D2_train <- pairwise_l2_curve_dist(curve_train, work_grid)
  K_train <- gaussian_functional_kernel(D2_train, tune$rho)
  fit <- fit_intercept_krr(K_train, train_dat$Y, tune$lambda)

  D2_test <- pairwise_l2_curve_dist(curve_test, work_grid, curve_train)
  K_test <- gaussian_functional_kernel(D2_test, tune$rho)
  pred_test <- predict_intercept_krr(K_test, fit)

  list(
    pred_mse = mean((test_dat$Y - pred_test)^2),
    num_pc = num_pc,
    sigma = tune$rho,
    lambda = tune$lambda,
    cv_mse = tune$cv_mse,
    median_pairwise_distance = tune$median_pairwise_distance,
    u_test = pred_test,
    cv_choice = "BIC_plus_functional_RKHS_CV"
  )
}

run_beta_method_main <- function(method_name, train_dat, test_dat, seed) {
  set.seed(seed)

  if (method_name == "pace_rkhs") {
    return(run_pace_rkhs_method(train_dat, test_dat, seed))
  }

  cv_res <- CVr_main(train_dat, train_dat$Y, kfold = pc_cv_folds, method = "IN", S = split_repeats)
  pts_grid <- seq(0, 1, length.out = dense_grid_length)
  Mu <- make_zero_mean_object()
  n <- length(train_dat$Y)
  Y <- train_dat$Y
  Y_test <- test_dat$Y

  if (method_name == "plugin") {
    m_use <- cv_res$CV_pi
    res <- fit_fpca_model(train_dat$Ly, train_dat$Lt, method_xi = "IN", num_pc = m_use, sigma2_value = main_sigma2_value)
    Phi_est <- res$phi[, seq_len(m_use), drop = FALSE]
    Ly_w <- train_dat$Ly
    for (i in seq_len(n)) Ly_w[[i]] <- (Y[i] - mean(Y)) * train_dat$Ly[[i]]
    g_dense <- smooth_mean_curve_plugin(Ly_w, train_dat$Lt, res$workGrid, m_use, main_sigma2_value)
    G <- integrate_scores_from_curve(g_dense, Phi_est, res$workGrid, m_use)
    b <- as.vector(G) / as.numeric(res$lambda[seq_len(m_use)])
    beta_hat <- as.vector(Phi_est %*% b)
    preds <- numeric(length(Y_test))
    for (i in seq_along(Y_test)) {
      beta_i <- approx(res$workGrid, beta_hat, xout = test_dat$Lt[[i]], rule = 2)$y
      preds[i] <- sum(beta_i * test_dat$Ly[[i]]) / length(test_dat$Ly[[i]])
    }
    return(list(pred_mse = mean((Y_test - preds)^2), num_pc = m_use, u_test = preds, cv_choice = "CV_pi"))
  }

  if (method_name == "in") {
    m_use <- cv_res$CV_In
    res <- fit_fpca_model(train_dat$Ly, train_dat$Lt, method_xi = "IN", num_pc = m_use, sigma2_value = main_sigma2_value)
    f_In <- lm(Y ~ res$xiEst[, seq_len(m_use), drop = FALSE])
    beta_hat <- as.vector(res$phi[, seq_len(m_use), drop = FALSE] %*% as.vector(f_In$coef[-1]))
    preds <- numeric(length(Y_test))
    for (i in seq_along(Y_test)) {
      beta_i <- approx(res$workGrid, beta_hat, xout = test_dat$Lt[[i]], rule = 2)$y
      preds[i] <- as.numeric(f_In$coefficients[1]) + sum(beta_i * test_dat$Ly[[i]]) / length(test_dat$Ly[[i]])
    }
    return(list(pred_mse = mean((Y_test - preds)^2), num_pc = m_use, u_test = preds, cv_choice = "CV_In"))
  }

  if (method_name == "pace") {
    m_use <- cv_res$CV_pace
    res <- FPCA_CE_mainlike(train_dat$Ly, train_dat$Lt, optns = list(nRegGrid = fpca_n_reg_grid, methodXi = "CE", maxK = M_fixed, FVEthreshold = 1, userMu = make_zero_mean_object(), methodSelectK = m_use, userSigma2 = main_sigma2_value, usergrid = FALSE))
    f_pace <- lm(Y ~ res$xiEst[, seq_len(m_use), drop = FALSE])
    beta_hat <- as.vector(res$phi[, seq_len(m_use), drop = FALSE] %*% as.vector(f_pace$coef[-1]))
    preds <- numeric(length(Y_test))
    for (i in seq_along(Y_test)) {
      beta_i <- approx(res$workGrid, beta_hat, xout = test_dat$Lt[[i]], rule = 2)$y
      preds[i] <- as.numeric(f_pace$coefficients[1]) + sum(beta_i * test_dat$Ly[[i]]) / length(test_dat$Ly[[i]])
    }
    return(list(pred_mse = mean((Y_test - preds)^2), num_pc = m_use, u_test = preds, cv_choice = "CV_pace"))
  }

  if (method_name == "split_s5") {
    m_use <- cv_res$CV_split
    beta_s <- matrix(0, nrow = split_repeats, ncol = dense_grid_length)
    a_vec <- numeric(split_repeats)
    for (s in seq_len(split_repeats)) {
      Index_1 <- sort(sample.int(n, floor(n / 2)))
      Index_2 <- sort(setdiff(seq_len(n), Index_1))
      sampX_1 <- list(Lt = train_dat$Lt[Index_1], Ly = train_dat$Ly[Index_1])
      sampX_2 <- list(Lt = train_dat$Lt[Index_2], Ly = train_dat$Ly[Index_2])
      Y_1 <- Y[Index_1]; Y_2 <- Y[Index_2]
      res1 <- fit_fpca_model(sampX_1$Ly, sampX_1$Lt, method_xi = "IN", num_pc = m_use, sigma2_value = main_sigma2_value)
      res2 <- fit_fpca_model(sampX_2$Ly, sampX_2$Lt, method_xi = "IN", num_pc = m_use, sigma2_value = main_sigma2_value)
      Phi_est21 <- res1$phi[, seq_len(m_use), drop = FALSE]
      Phi_est22 <- align_phi_to_reference(res2$phi[, seq_len(m_use), drop = FALSE], Phi_est21)
      Xi_21 <- matrix(0, nrow = length(Index_1), ncol = m_use)
      Xi_22 <- matrix(0, nrow = length(Index_2), ncol = m_use)
      for (i in seq_along(Index_1)) {
        for (k in seq_len(m_use)) {
          Xi_21[i, k] <- sum(train_dat$Ly[[Index_1[i]]] * approx(res2$workGrid, Phi_est22[, k], xout = train_dat$Lt[[Index_1[i]]], rule = 2)$y) / length(train_dat$Ly[[Index_1[i]]])
        }
      }
      for (i in seq_along(Index_2)) {
        for (k in seq_len(m_use)) {
          Xi_22[i, k] <- sum(train_dat$Ly[[Index_2[i]]] * approx(res1$workGrid, Phi_est21[, k], xout = train_dat$Lt[[Index_2[i]]], rule = 2)$y) / length(train_dat$Ly[[Index_2[i]]])
        }
      }
      Xi_3 <- rbind(Xi_21, Xi_22)
      Y_3 <- c(Y_1, Y_2)
      f3 <- lm(Y_3 ~ Xi_3)
      a_vec[s] <- as.numeric(f3$coefficients[1])
      beta_s[s, ] <- as.vector(0.5 * (Phi_est21 + Phi_est22) %*% as.vector(f3$coef[-1]))
    }
    beta_hat <- colMeans(beta_s)
    preds <- numeric(length(Y_test))
    for (i in seq_along(Y_test)) {
      beta_i <- approx(pts_grid, beta_hat, xout = test_dat$Lt[[i]], rule = 2)$y
      preds[i] <- mean(a_vec) + sum(beta_i * test_dat$Ly[[i]]) / length(test_dat$Ly[[i]])
    }
    return(list(pred_mse = mean((Y_test - preds)^2), num_pc = m_use, u_test = preds, cv_choice = "CV_split"))
  }

  stop("Unsupported method_name.")
}

# =========================
# 8) one replication
# =========================
run_one_replication <- function(cfg, rep_idx) {
  seed <- seed_start + rep_idx - 1L
  run_tag <- build_run_tag(
    cfg$n,
    cfg$obs_count,
    cfg$obs_count_mode,
    cfg$y_mode,
    cfg$score_dist,
    cfg$x_noise_dist,
    cfg$method
  )

  start_time <- proc.time()[["elapsed"]]
  set.seed(seed)
  sim <- generate_train_test_data(
    cfg$n,
    cfg$obs_count,
    cfg$obs_count_mode,
    cfg$y_mode,
    cfg$score_dist,
    cfg$x_noise_dist
  )

  if (cfg$method == "embedding") {
    res <- run_embedding_method(sim$train, sim$test)
    out <- data.frame(
      task_id = cfg$task_id,
      run_tag = run_tag,
      n = cfg$n,
      obs_count = cfg$obs_count,
      obs_count_mode = cfg$obs_count_mode,
      y_mode = cfg$y_mode,
      score_dist = cfg$score_dist,
      x_noise_dist = cfg$x_noise_dist,
      method = cfg$method,
      rep = rep_idx,
      seed = seed,
      pred_mse = res$pred_mse,
      num_pc = NA_integer_,
      l_H = res$l_H,
      sigma_h = res$sigma_h,
      sigma_m = res$sigma_m,
      lambda = res$lambda,
      gcv = res$gcv,
      u_values = NA_character_,
      u_min = NA_real_,
      u_q25 = NA_real_,
      u_median = NA_real_,
      u_q75 = NA_real_,
      u_max = NA_real_,
      runtime_seconds = proc.time()[["elapsed"]] - start_time,
      stringsAsFactors = FALSE
    )
  } else {
    res <- run_beta_method_main(
      cfg$method,
      sim$train,
      sim$test,
      seed = seed
    )
    out <- data.frame(
      task_id = cfg$task_id,
      run_tag = run_tag,
      n = cfg$n,
      obs_count = cfg$obs_count,
      obs_count_mode = cfg$obs_count_mode,
      y_mode = cfg$y_mode,
      score_dist = cfg$score_dist,
      x_noise_dist = cfg$x_noise_dist,
      method = cfg$method,
      rep = rep_idx,
      seed = seed,
      pred_mse = res$pred_mse,
      num_pc = res$num_pc,
      l_H = if ("median_pairwise_distance" %in% names(res)) res$median_pairwise_distance else NA_real_,
      sigma_h = if ("sigma" %in% names(res)) res$sigma else NA_real_,
      sigma_m = NA_real_,
      lambda = if ("lambda" %in% names(res)) res$lambda else NA_real_,
      gcv = if ("cv_mse" %in% names(res)) res$cv_mse else NA_real_,
      u_values = serialize_numeric_vector(res$u_test),
      u_min = if (length(res$u_test) > 0 && any(is.finite(res$u_test))) min(res$u_test, na.rm = TRUE) else NA_real_,
      u_q25 = safe_quantile(res$u_test, 0.25),
      u_median = safe_quantile(res$u_test, 0.50),
      u_q75 = safe_quantile(res$u_test, 0.75),
      u_max = if (length(res$u_test) > 0 && any(is.finite(res$u_test))) max(res$u_test, na.rm = TRUE) else NA_real_,
      runtime_seconds = proc.time()[["elapsed"]] - start_time,
      stringsAsFactors = FALSE
    )
  }

  out
}

# =========================
# 9) task runner
# =========================
run_task <- function(task_id) {
  task_grid <- make_task_grid()
  write_manifest_if_needed(task_grid)
  if (task_id < 1 || task_id > nrow(task_grid)) {
    stop(sprintf("task_id must be between 1 and %d.", nrow(task_grid)))
  }

  cfg <- task_grid[task_id, , drop = FALSE]
  cfg$task_id <- task_id
  run_tag <- build_run_tag(
    cfg$n,
    cfg$obs_count,
    cfg$obs_count_mode,
    cfg$y_mode,
    cfg$score_dist,
    cfg$x_noise_dist,
    cfg$method
  )
  cfg$run_tag <- run_tag
  out_dir <- file.path(root_results_dir, run_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  per_rep_path <- file.path(out_dir, "per_rep_results.csv")
  summary_path <- file.path(out_dir, "method_summary.csv")
  config_path <- file.path(out_dir, "run_config.csv")
  global_seed_live_path <- file.path(root_results_dir, "all_seed_results_live.csv")
  global_summary_path <- file.path(root_results_dir, "all_method_config_summary.csv")

  write.csv(cfg, config_path, row.names = FALSE)

  read_existing_rep_results <- function(csv_path, n_reps_target) {
    if (!file.exists(csv_path)) {
      return(NULL)
    }

    existing <- tryCatch(
      read.csv(csv_path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(existing) || nrow(existing) == 0 || !("rep" %in% names(existing))) {
      return(NULL)
    }

    existing$rep <- suppressWarnings(as.integer(existing$rep))
    existing <- existing[!is.na(existing$rep), , drop = FALSE]
    existing <- existing[existing$rep >= 1L & existing$rep <= n_reps_target, , drop = FALSE]
    if (nrow(existing) == 0) {
      return(NULL)
    }

    # If a previous run appended the same replication more than once, keep the
    # last row for that replication. This preserves resume behavior while
    # avoiding duplicated rows in the final per-replication output.
    existing$.row_order_for_resume <- seq_len(nrow(existing))
    existing <- existing[order(existing$rep, existing$.row_order_for_resume), , drop = FALSE]
    existing <- existing[!duplicated(existing$rep, fromLast = TRUE), , drop = FALSE]
    existing$.row_order_for_resume <- NULL
    existing <- existing[order(existing$rep), , drop = FALSE]
    rownames(existing) <- NULL
    existing
  }

  upsert_summary_atomic <- function(summary_row, csv_path) {
    lock_dir <- paste0(csv_path, ".lockdir")
    acquire_lock(lock_dir)
    on.exit(release_lock(lock_dir), add = TRUE)
    dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)

    if (file.exists(csv_path)) {
      old <- tryCatch(
        read.csv(csv_path, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (!is.null(old) && nrow(old) > 0) {
        if ("task_id" %in% names(old)) {
          old <- old[as.integer(old$task_id) != as.integer(summary_row$task_id), , drop = FALSE]
        } else if ("run_tag" %in% names(old)) {
          old <- old[old$run_tag != summary_row$run_tag, , drop = FALSE]
        }
        out <- rbind(old, summary_row)
      } else {
        out <- summary_row
      }
    } else {
      out <- summary_row
    }

    write.csv(out, csv_path, row.names = FALSE)
  }

  existing_df <- read_existing_rep_results(per_rep_path, n_reps)
  existing_reps <- if (!is.null(existing_df)) sort(unique(existing_df$rep)) else integer(0)
  missing_reps <- setdiff(seq_len(n_reps), existing_reps)

  cat(sprintf(
    "[fdr_2.R] resume check task_id=%s run_tag=%s existing=%d missing=%d total=%d per_rep=%s\n",
    cfg$task_id,
    run_tag,
    length(existing_reps),
    length(missing_reps),
    n_reps,
    per_rep_path
  ))
  flush.console()

  rep_results <- list()
  if (!is.null(existing_df) && nrow(existing_df) > 0) {
    rep_results <- split(existing_df, seq_len(nrow(existing_df)))
  }

  if (length(missing_reps) > 0) {
    for (rep_idx in missing_reps) {
      rep_row <- tryCatch(
        run_one_replication(cfg, rep_idx),
        error = function(e) {
          data.frame(
            task_id = cfg$task_id,
            run_tag = run_tag,
            n = cfg$n,
            obs_count = cfg$obs_count,
            obs_count_mode = cfg$obs_count_mode,
            y_mode = cfg$y_mode,
            score_dist = cfg$score_dist,
            x_noise_dist = cfg$x_noise_dist,
            method = cfg$method,
            rep = rep_idx,
            seed = seed_start + rep_idx - 1L,
            pred_mse = NA_real_,
            num_pc = NA_integer_,
            l_H = NA_real_,
            sigma_h = NA_real_,
            sigma_m = NA_real_,
            lambda = NA_real_,
            gcv = NA_real_,
            u_values = NA_character_,
            u_min = NA_real_,
            u_q25 = NA_real_,
            u_median = NA_real_,
            u_q75 = NA_real_,
            u_max = NA_real_,
            runtime_seconds = NA_real_,
            error_message = e$message,
            stringsAsFactors = FALSE
          )
        }
      )

      if (!("error_message" %in% names(rep_row))) rep_row$error_message <- NA_character_

      rep_results[[length(rep_results) + 1L]] <- rep_row

      # Real-time local write: this file lives inside the folder for this
      # parameter combination and will contain all n_reps rows after completion.
      append_row_atomic(rep_row, per_rep_path)
      append_row_atomic(rep_row, global_seed_live_path)

      pred_mse_msg <- ifelse(
        is.na(rep_row$pred_mse),
        "NA",
        format(rep_row$pred_mse, digits = 8, scientific = TRUE)
      )
      cat(sprintf(
        "[fdr_2.R] task_id=%s rep=%d/%d pred_mse=%s wrote=%s\n",
        cfg$task_id,
        rep_idx,
        n_reps,
        pred_mse_msg,
        per_rep_path
      ))
      flush.console()
    }
  } else {
    cat(sprintf(
      "[fdr_2.R] task_id=%s already complete; rebuilding summary only.\n",
      cfg$task_id
    ))
    flush.console()
  }

  rep_df <- do.call(rbind, rep_results)
  rep_df$rep <- suppressWarnings(as.integer(rep_df$rep))
  rep_df <- rep_df[!is.na(rep_df$rep), , drop = FALSE]
  rep_df <- rep_df[rep_df$rep >= 1L & rep_df$rep <= n_reps, , drop = FALSE]
  if (nrow(rep_df) > 0) {
    rep_df$.row_order_for_resume <- seq_len(nrow(rep_df))
    rep_df <- rep_df[order(rep_df$rep, rep_df$.row_order_for_resume), , drop = FALSE]
    rep_df <- rep_df[!duplicated(rep_df$rep, fromLast = TRUE), , drop = FALSE]
    rep_df$.row_order_for_resume <- NULL
    rep_df <- rep_df[order(rep_df$rep), , drop = FALSE]
    rownames(rep_df) <- NULL
  }

  # Rewrite the per-replication CSV in a clean, sorted, de-duplicated form.
  write.csv(rep_df, per_rep_path, row.names = FALSE)

  u_all <- if ("u_values" %in% names(rep_df)) {
    parse_serialized_numeric_vector(rep_df$u_values)
  } else {
    numeric(0)
  }
  summary_df <- data.frame(
    task_id = cfg$task_id,
    run_tag = run_tag,
    n = cfg$n,
    obs_count = cfg$obs_count,
    obs_count_mode = cfg$obs_count_mode,
    y_mode = cfg$y_mode,
    score_dist = cfg$score_dist,
    x_noise_dist = cfg$x_noise_dist,
    method = cfg$method,
    n_reps = n_reps,
    n_success = sum(!is.na(rep_df$pred_mse)),
    mean_pred_mse = safe_mean(rep_df$pred_mse),
    sd_pred_mse = safe_sd(rep_df$pred_mse),
    u_min = if (length(u_all) > 0) min(u_all, na.rm = TRUE) else NA_real_,
    u_q25 = safe_quantile(u_all, 0.25),
    u_median = safe_quantile(u_all, 0.50),
    u_q75 = safe_quantile(u_all, 0.75),
    u_max = if (length(u_all) > 0) max(u_all, na.rm = TRUE) else NA_real_,
    mean_runtime_seconds = safe_mean(rep_df$runtime_seconds),
    total_runtime_seconds = sum(rep_df$runtime_seconds, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  write.csv(summary_df, summary_path, row.names = FALSE)
  upsert_summary_atomic(summary_df, global_summary_path)
  invisible(summary_df)
}

aggregate_global_summary <- function() {
  task_grid <- make_task_grid()
  files <- vapply(seq_len(nrow(task_grid)), function(i) {
    cfg <- task_grid[i, , drop = FALSE]
    file.path(
      root_results_dir,
      build_run_tag(
        cfg$n,
        cfg$obs_count,
        cfg$obs_count_mode,
        cfg$y_mode,
        cfg$score_dist,
        cfg$x_noise_dist,
        cfg$method
      ),
      "method_summary.csv"
    )
  }, character(1))
  files <- files[file.exists(files)]
  if (length(files) == 0) stop("No method_summary.csv files found.")
  summary_df <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  write.csv(
    summary_df,
    file.path(root_results_dir, "all_method_config_summary_rebuilt.csv"),
    row.names = FALSE
  )
  invisible(summary_df)
}

# =========================
# 10) cli
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript fdr_2.R task <task_id>")
mode <- args[[1]]

if (mode == "task") {
  if (length(args) < 2) stop("Please provide task id.")
  task_id <- as.integer(args[[2]])
  run_task(task_id)
} else if (mode == "manifest") {
  task_grid <- make_task_grid()
  dir.create(root_results_dir, recursive = TRUE, showWarnings = FALSE)
  write_manifest_if_needed(task_grid)
  print(task_grid)
} else if (mode == "aggregate") {
  aggregate_global_summary()
} else {
  stop("Unsupported mode.")
}
