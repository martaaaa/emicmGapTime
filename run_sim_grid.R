########################################################################
##  BIVARIATE INTERVAL-CENSORING · WEIGHTED EM–ICM · SIMULATION SUITE
##  (clean, robust version –  June 2025)
########################################################################

## ────────────────────────────────────────────────────────────────────
##  0 ·  PACKAGES  (install on-the-fly if missing)
## ────────────────────────────────────────────────────────────────────
required <- c("data.table", "Matrix", "copula", "isotone", "ranger",
              "SQUAREM", "future", "future.apply")

for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("→ installing ‘", p, "’ …")
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

## reproducible parallel RNG
plan(multisession, workers = future::availableCores())
options(future.rng.onMisuse = "error")

## data.table uses all logical cores inside each worker
data.table::setDTthreads(0)

## ────────────────────────────────────────────────────────────────────
##  1 ·  GLOBAL CONFIGURATION   (matches manuscript text)
## ────────────────────────────────────────────────────────────────────
CFG <- list(
  tau_max    = 5,
  eps        = 1e-8,        # numerical floor for rectangles/EM
  w_cap      = 25,          # cap for extreme IPCW weights
  min_delta2 = 0.04,        # ≥ 4 % exactly-observed Y
  m_method   = "rf",        # "rf" | "logit" | "kernel"
  ntrees     = 1200,
  ps_clip    = 1e-4,        # bounds for m-hat
  p_exact1   = 0.15,        # probability of exact T1
  p_exact2   = 0.25         # probability of exact Y
)
TRUE_MAX_ITER <- 750        # EM–ICM hard stop

## ────────────────────────────────────────────────────────────────────
##  2 ·  DATA-GENERATION TOOLS
## ────────────────────────────────────────────────────────────────────
## 2·1  (T1,T2)  ~  Exp(1) margins + FGM(θ) dependence
rFGM <- function(n, theta = 0) {
  stopifnot(abs(theta) <= 1)
  U <- copula::rCopula(n, copula::fgmCopula(theta, dim = 2))
  cbind(T1 = -log(U[,1L]), T2 = -log(U[,2L]))
}

## 2·2  Visit schedule – Poisson(λ)+1 uniformly in (0,τ_max)
simVisits <- function(n, lambda, tau_max) {
  J <- rpois(n, lambda) + 1L
  lapply(J, function(k) sort(runif(k, 0, tau_max)))
}

## 2·3  Optionally add exact T1 / exact Y
addExact <- function(vis, T1, Y, p1, p2) {
  Map(function(v, t1, y)
    sort(unique(c(v,
                  if (runif(1) < p1) t1,
                  if (runif(1) < p2) y))),
    vis, T1, Y)
}

## 2·4  Independent administrative right-censoring
rCens <- function(n, c_max) runif(n, 0, c_max)

## ────────────────────────────────────────────────────────────────────
##  3 ·  INTERVAL BUILDER  (guarantees R1 < LY)
## ────────────────────────────────────────────────────────────────────
buildRect <- function(T1, T2, visits, C, tau_max, eps) {

  Y   <- T1 + T2
  out <- vector("list", length(T1))

  for (i in seq_along(T1)) {

    t  <- c(0, visits[[i]], Inf)
    N1 <- as.integer(T1[i] <= t)     # after first event
    NY <- as.integer(Y [i] <= t)     # after second event

    ## first gap
    L1 <- max(t[N1 == 0])
    R1 <- min(t[N1 == 1]); if (!is.finite(R1)) R1 <- tau_max

    ## second event
    LY <- max(t[N1 == 1 & NY == 0], R1)
    RY <- min(t[NY == 1]); if (!is.finite(RY)) RY <- tau_max

    ## exact observations
    if (any(abs(t - T1[i]) < eps)) L1 <- R1 <- T1[i]
    if (any(abs(t - Y [i]) < eps)) LY <- RY <-  Y[i]

    ## independent right-censoring
    R1 <- min(R1, C[i], tau_max)
    RY <- min(RY, C[i], tau_max)

    if (R1 - L1 < eps) R1 <- L1 + eps   # widen degenerate intervals
    if (RY - LY < eps) RY <- LY + eps

    out[[i]] <- list(L1=L1, R1=R1, LY=LY, RY=RY,
                     delta1 = as.integer(abs(R1-L1) < eps),
                     delta2 = as.integer(abs(RY-LY) < eps))
  }
  data.table::rbindlist(out)
}

## ────────────────────────────────────────────────────────────────────
##  4 ·  GEOMETRY  – Turnbull incidence matrix
## ────────────────────────────────────────────────────────────────────
## ──────────────────────────────────────────────────────────────────
##  GEOMETRY · Turnbull incidence matrix (versão final)
## ──────────────────────────────────────────────────────────────────
incidenceMatrix <- function(rect, tau_max)
{
  # 1 · quebras únicas nos eixos
  xs <- sort(unique(c(0, rect$L1, rect$R1, tau_max)))
  zs <- sort(unique(c(0, rect$LY, rect$RY, tau_max)))

  # 2 · todas as faixas (ix, iz)
  cells <- data.table::CJ(ix = seq_len(length(xs) - 1L),
                          iz = seq_len(length(zs) - 1L),
                          sorted = FALSE)

  # 3 · converter índices → limites reais
  cells[, :=(x_lo = xs[ix],
             x_hi = xs[ix + 1L],
             z_lo = zs[iz],
             z_hi = zs[iz + 1L])]

  # 4 · manter apenas células com x < z (estrito)
  cells <- cells[x_hi < z_lo]

  # 5 · descartar colunas auxiliares
  cells[, c("ix", "iz") := NULL]

  # 6 · se nenhuma célula sobreviveu, devolva matriz vazia
  if (nrow(cells) == 0L) {
    A <- Matrix::Matrix(0, nrow(rect), 0L, sparse = TRUE)
    return(list(A = A, cells = cells))
  }

  # 7 · lista de índices: quais células estão 100 % contidas em cada retângulo
  idx_list <- lapply(seq_len(nrow(rect)), function(i)
    which(cells$x_lo >= rect$L1[i] & cells$x_hi <= rect$R1[i] &
            cells$z_lo >= rect$LY[i] & cells$z_hi <= rect$RY[i]))

  # 8 · matriz esparsa (linhas = retângulos, colunas = células)
  A <- Matrix::sparseMatrix(
    i    = rep(seq_along(idx_list), lengths(idx_list)),
    j    = unlist(idx_list),
    x    = 1,
    dims = c(nrow(rect), nrow(cells))
  )

  list(A = A, cells = cells)
}


## ────────────────────────────────────────────────────────────────────
##  5 ·  IPCW WEIGHTS  (KM core + capping)
## ────────────────────────────────────────────────────────────────────
KMcore <- function(RY, v) {

  RY <- as.numeric(RY); v <- as.numeric(v)
  stopifnot(length(RY) == length(v))

  ord   <- order(RY, -v, method = "radix")   # RY ↑  then  v ↓
  n     <- length(ord)
  risk  <- n:1
  surv  <- cumprod(1 - v[ord] / risk)
  W     <- v[ord] / risk * c(1, head(surv, -1))

  out <- numeric(n); out[ord] <- W
  out
}

weightKM <- function(RY, d2, cfg = CFG) {
  w <- KMcore(RY, d2)
  cap <- cfg$w_cap * mean(w)
  pmin(w, cap)
}
weightPS <- function(RY, mhat, cfg = CFG) {
  w <- KMcore(RY, mhat)
  cap <- cfg$w_cap * mean(w)
  pmin(w, cap)
}

## ────────────────────────────────────────────────────────────────────
##  6 ·  m(x,z)  ESTIMATION  (three options)
## ────────────────────────────────────────────────────────────────────
estimate_m <- function(rect, cfg = CFG) {

  train <- rect[delta1 == 1L]                 # exact T1 only
  if (nrow(train) < 25L)
    return(rep(mean(rect$delta2), nrow(rect)))

  if (cfg$m_method == "logit") {

    fit <- glm(delta2 ~ R1 + RY, data = train, family = binomial)
    mp  <- predict(fit, newdata = rect, type = "response")

  } else if (cfg$m_method == "kernel") {

    h  <- 1.06 * stats::sd(train$R1) * nrow(train)^(-1/6)
    K  <- function(x) exp(-0.5*x^2) / sqrt(2*pi)
    num <- rowSums( outer(rect$R1, train$R1, \(x,y) K((x-y)/h)) *
                      outer(rect$RY, train$RY, \(z,w) K((z-w)/h)) *
                      train$delta2 )
    den <- rowSums( outer(rect$R1, train$R1, \(x,y) K((x-y)/h)) *
                      outer(rect$RY, train$RY, \(z,w) K((z-w)/h)) )
    mp  <- num / pmax(den, 1e-9)

  } else {                                    # random-forest default
    train[, f2 := factor(delta2, levels = 0:1)]
    rf  <- ranger::ranger(f2 ~ R1 + RY, data = train,
                          probability     = TRUE,
                          num.trees       = cfg$ntrees,
                          mtry            = 1,
                          min.node.size   = 5,
                          respect.unordered.factors = "partition")
    mp  <- predict(rf, data = rect)$predictions[,"1"]
  }

  pmin(pmax(mp, cfg$ps_clip), 1 - cfg$ps_clip)
}

## ────────────────────────────────────────────────────────────────────
##  7 ·  ISOTONIC PROJECTION Π_{≼}  (robust 2-D IPAVA, no gpava)
## ────────────────────────────────────────────────────────────────────
## ----------  compatibilidade com isotone ≥ 1.3  ----------
pava_safe <- function(x) {
  if (requireNamespace("isotone", quietly = TRUE)) {
    ns <- asNamespace("isotone")
    if (exists("pava", envir = ns, inherits = FALSE)) {
      return(get("pava", envir = ns)(x))    # isotone:::pava
    } else {
      ## fallback: gpava() sem pesos
      return(isotone::gpava(1L:length(x), x)$x)
    }
  } else {
    stop("Pacote 'isotone' não está instalado.")
  }
}
#--------------------------
ICMproj <- function(u, cells, tol = 1e-10, maxit = 1000) {

  ux <- sort(unique(cells$x_hi))
  uz <- sort(unique(cells$z_hi))
  m  <- length(ux); n <- length(uz)

  grid <- matrix(NA_real_, m, n)
  idx  <- matrix(NA_integer_, m, n)
  for (k in seq_len(nrow(cells))) {
    i <- match(cells$x_hi[k], ux)
    j <- match(cells$z_hi[k], uz)
    grid[i, j] <- u[k]; idx[i, j] <- k
  }

  iter <- 0; changed <- TRUE
  while (changed && iter < maxit) {
    changed <- FALSE; iter <- iter + 1

    ## rows : x fixed, pool along z
    for (i in seq_len(m)) {
      good <- which(!is.na(grid[i, ]))
      if (length(good) > 1L) {
        new <- pava_safe(grid[i, good])
        if (max(abs(new - grid[i, good])) > tol) {
          grid[i, good] <- new; changed <- TRUE
        }
      }
    }
    ## cols : z fixed, pool along x
    for (j in seq_len(n)) {
      good <- which(!is.na(grid[, j]))
      if (length(good) > 1L) {
        new <- pava_safe(grid[good, j])
        if (max(abs(new - grid[good, j])) > tol) {
          grid[good, j] <- new; changed <- TRUE
        }
      }
    }
  }

  u_new <- u
  for (i in seq_len(m))
    for (j in seq_len(n))
      if (!is.na(idx[i, j]))
        u_new[idx[i, j]] <- grid[i, j]

  u_new[u_new < 0] <- 0
  u_new / sum(u_new)
}

## ────────────────────────────────────────────────────────────────────
##  8 ·  WEIGHTED EM–ICM (global convergence, SQUAREM-accelerated)
## ────────────────────────────────────────────────────────────────────
emICM <- function(rect, W, cfg = CFG,
                  maxit = TRUE_MAX_ITER, tol = 1e-8) {

  inc   <- incidenceMatrix(rect, cfg$tau_max)
  A     <- inc$A
  cells <- inc$cells
  K     <- ncol(A)

  ## ---------------------------------------------
  ## NOVO: lidar com K == 0  (nenhuma célula)
  ## ---------------------------------------------
  if (K == 0L) {
    return(list(u = numeric(0), iter = 0,
                converged = FALSE, cells = cells))
  }

  W  <- W / pmax(mean(W), .Machine$double.eps)   # μ = 1
  u0 <- rep(1 / K, K)                            # arranque uniforme

  fixpt <- function(u) {
    d  <- pmax(as.numeric(A %*% u), cfg$eps)
    u1 <- u * as.numeric(Matrix::crossprod(A, W / d))
    u1 <- u1 / sum(u1)
    ICMproj(u1, cells)
  }

  ## ---------- SQUAREM ou fallback EM ----------
  if (K > 1 && requireNamespace("SQUAREM", quietly = TRUE)) {
    res <- tryCatch(
      SQUAREM::squarem(par = u0, fixptfn = fixpt,
                       control = list(maxiter = maxit, tol = tol)),
      error = function(e) NULL)
    if (!is.null(res)) {
      return(list(u = res$par, iter = res$iter,
                  converged = (res$convergence == 0), cells = cells))
    }
  }

  ## EM clássico
  u <- u0
  for (it in seq_len(maxit)) {
    u_new <- fixpt(u)
    if (sum(abs(u_new - u)) < tol) break
    u <- u_new
  }
  list(u = u, iter = it, converged = (it < maxit), cells = cells)
}



## ────────────────────────────────────────────────────────────────────
##  9 ·  TRUE  F12(x,y)   (Exp(1) margins + FGM(θ))
## ────────────────────────────────────────────────────────────────────
trueF12 <- function(x, y, theta) {
  u <- 1 - exp(-x); v <- 1 - exp(-y)
  u * v + theta * u * v * (1 - u) * (1 - v)
}

## ────────────────────────────────────────────────────────────────────
## 10 ·  CDF EVALUATION (piecewise-constant estimator)
## ────────────────────────────────────────────────────────────────────
evalF12 <- function(fit, x, y) {
  ## Se não há malha ou u, não há estimativa
  if (is.null(fit) || length(fit$u) == 0L) return(NA_real_)
  ## Caso contrário, use a última iteração, convergida ou não
  idx <- with(fit$cells, which(x_hi <= x & z_hi <= x + y))
  sum(fit$u[idx])
}


## ────────────────────────────────────────────────────────────────────
## 11 ·  ONE REPLICATE
## ────────────────────────────────────────────────────────────────────
oneReplicate <- function(n, lambda, cmax, theta,
                         eval_pt, cfg = CFG, seed = NULL,
                         use_ps = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  ## 1 · gerar dados
  gaps   <- rFGM(n, theta)
  visits <- addExact(simVisits(n, lambda, cfg$tau_max),
                     gaps[, 1L],
                     gaps[, 1L] + gaps[, 2L],
                     cfg$p_exact1, cfg$p_exact2)

  rect <- buildRect(gaps[, 1L], gaps[, 2L], visits,
                    C        = rCens(n, cmax),
                    tau_max  = cfg$tau_max,
                    eps      = cfg$eps)

  ## 2 · pesos IPCW (KM) – sempre disponíveis
  Wkm   <- weightKM(rect$RY, rect$delta2, cfg)
  fitKM <- emICM(rect, Wkm, cfg)

  ## 3 · pesos pré-suavizados (opcional)
  fitPS <- NA
  if (use_ps && mean(rect$delta2) >= cfg$min_delta2) {
    mhat <- estimate_m(rect, cfg)
    Wps  <- weightPS(rect$RY, mhat, cfg)
    fitPS <- emICM(rect, Wps, cfg)
  }

  ## 4 · guardar estatísticas
  x <- eval_pt$x
  y <- eval_pt$y
  list(
    stats = data.table::data.table(
      Fkm      = evalF12(fitKM, x, y),
      Fps      = if (is.list(fitPS)) evalF12(fitPS, x, y) else NA_real_,
      iter_km  = fitKM$iter,
      iter_ps  = if (is.list(fitPS)) fitPS$iter else NA_integer_,
      exactY   = sum(rect$delta2)
    ),
    true = trueF12(x, y, theta)
  )
}

## ────────────────────────────────────────────────────────────────────
## 12 ·  MONTE-CARLO DRIVER  (parallel, future.apply)
## ────────────────────────────────────────────────────────────────────
## ──────────────────────────────────────────────────────────────────
##  RUN ONE SCENARIO  –  Monte-Carlo resumido
## ──────────────────────────────────────────────────────────────────
runScenario <- function(n, B,
                        lambda, cmax, theta,
                        eval_pt,
                        cfg = CFG)
{
  ## 1 · replicates em paralelo
  reps <- future.apply::future_lapply(
    seq_len(B),
    future.seed = TRUE,
    FUN = function(b)
      oneReplicate(n, lambda, cmax, theta,
                   eval_pt = eval_pt, cfg = cfg)
  )

  ## 2 · tabela bruta
  tab  <- data.table::rbindlist(lapply(reps, function(x) x$stats))
  true <- reps[[1]]$true                # valor exato

  tab[, :=(bias_km = Fkm - true,
           bias_ps = Fps - true)]

  ## 3 · função de variância que devolve NA se tudo é NA
  v_ok <- function(v) if (all(is.na(v))) NA_real_ else var(v, na.rm = TRUE)

  ## 4 · agregados
  summary <- tab[, {
    var_km <- v_ok(Fkm)
    var_ps <- v_ok(Fps)

    list(
      mean_km     = mean(Fkm,  na.rm = TRUE),
      var_km      = var_km,
      bias_km     = mean(bias_km, na.rm = TRUE),
      rmse_km     = sqrt(mean(bias_km^2, na.rm = TRUE)),

      mean_ps     = mean(Fps,  na.rm = TRUE),
      var_ps      = var_ps,
      bias_ps     = mean(bias_ps, na.rm = TRUE),
      rmse_ps     = sqrt(mean(bias_ps^2, na.rm = TRUE)),

      RE          = if (is.finite(var_ps)) var_km / var_ps else NA_real_,
      prop_exactY = mean(exactY) / n
    )
  }]

  list(summary = summary, raw = tab)
}

## ────────────────────────────────────────────────────────────────────
## 13 ·  QUICK INTEGRITY CHECK
## ────────────────────────────────────────────────────────────────────

eval_pt <- list(x = stats::qexp(.30), y = stats::qexp(.30))

demo <- runScenario(n = 100, B = 5,
                    lambda = 2, cmax = 3, theta = 0.50,
                    eval_pt = eval_pt)

print(demo$summary)

#----------------------

## ────────────────────────────────────────────────────────────────────
##  SIMULATION MODULE
##  • grid definition
##  • execution
##  • storage (data.frame + csv)
## ────────────────────────────────────────────────────────────────────

library(data.table)

## 1 · utilidade: corre um único cenário e devolve lista(summary, raw)
run_one <- function(n, lambda, cmax, theta,
                    B       = 20,
                    eval_pt = list(x = qexp(0.30), y = qexp(0.30)),
                    cfg     = CFG) {

  runScenario(n       = n,
              B       = B,
              lambda  = lambda,
              cmax    = cmax,
              theta   = theta,
              eval_pt = eval_pt,
              cfg     = cfg)
}

## 2 · grelha de parâmetros + execução paralela
run_sim_grid <- function(ns, lambdas, cmaxs, thetas,
                         B        = 20,
                         eval_pt  = list(x = qexp(0.30), y = qexp(0.30)),
                         cfg      = CFG,
                         out_dir  = "sim_results") {

  ## criar diretório de saída (se necessário)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  ## grelha completa
  grid <- CJ(n = ns,
             lambda = lambdas,
             cmax   = cmaxs,
             theta  = thetas)

  ## contêineres
  res_summary <- vector("list", nrow(grid))
  res_raw     <- vector("list", nrow(grid))

  ## ciclo sobre a grelha
  for (g in seq_len(nrow(grid))) {

    pars <- grid[g]
    cat(sprintf(">>> cenário %d / %d  (n=%d, λ=%g, cmax=%g, θ=%g)\n",
                g, nrow(grid), pars$n, pars$lambda, pars$cmax, pars$theta))

    out <- run_one(n        = pars$n,
                   lambda   = pars$lambda,
                   cmax     = pars$cmax,
                   theta    = pars$theta,
                   B        = B,
                   eval_pt  = eval_pt,
                   cfg      = cfg)

    ## acrescentar identificadores ao summary
    summary_df <- cbind(pars, out$summary)
    res_summary[[g]] <- summary_df

    ## guardar summary em CSV cumulativo (append = TRUE)
    fwrite(summary_df,
           file      = file.path(out_dir, "summary_all.csv"),
           append    = file.exists(file.path(out_dir, "summary_all.csv")))

    ## guardar raw replicados: um ficheiro por cenário
    raw_fname <- sprintf("raw_n%g_l%g_c%g_t%g.csv",
                         pars$n, pars$lambda, pars$cmax, pars$theta)
    fwrite(out$raw, file = file.path(out_dir, raw_fname))

    res_raw[[g]] <- out$raw
  }

  ## devolver resultados em memória
  list(
    grid     = grid,
    summary  = rbindlist(res_summary),
    raw_list = res_raw
  )
}

## ────────────────────────────────────────────────────────────────────
## CHAMADA
##
#
sim_out50_100 <- run_sim_grid(
  ns      = c(100),
  lambdas = c(1, 2),
  cmaxs   = c(3),
  thetas  = c(0, 0.5),
  B       = 50           #
)
sim_out50_100
sim_out50_100$summary

#------------------------

sim_out100 <- run_sim_grid(
  ns      = c(100),
  lambdas = c(1, 2),
  cmaxs   = c(3),
  thetas  = c(0, 0.5),
  B       = 100            #
)
sim_out100$summary
sim_out100

#------------------------

sim_out50_300 <- run_sim_grid(
  ns      = c(300),
  lambdas = c(1, 2),
  cmaxs   = c(3),
  thetas  = c(0, 0.5),
  B       = 50            #
)
sim_out50_300$summary
sim_out50_300
#-------------

sim_out50_300b <- run_sim_grid(
  ns      = c(300),
  lambdas = c(1, 2),
  cmaxs   = c(3),
  thetas  = c(0, 0.5),
  B       = 50            #
)
sim_out50_300b$summary
sim_out50_300b

#------------------------

sim_out50_500 <- run_sim_grid(
  ns      = c(500),
  lambdas = c(1, 2),
  cmaxs   = c(3),
  thetas  = c(0, 0.5),
  B       = 50            #
)
sim_out50_500$summary
sim_out50_500
#-------------

sim_out50_500b <- run_sim_grid(
  ns      = c(500),
  lambdas = c(1, 2),
  cmaxs   = c(3),
  thetas  = c(0, 0.5),
  B       = 50            #
)
sim_out50_500b$summary
sim_out50_500b


attr(sim_objs, "info") <- list(
  date      = Sys.Date(),
  git_commit= "a1b2c3d",
  R_version = R.version.string
)
saveRDS(sim_objs, "sim_objs_v2025-06-08.rds", compress = "xz")
