########################################################################
##  BIVARIATE INTERVAL‑CENSORING · EM–ICM  (MACS / KMsurv::aids demo)
##  ─ Fully self‑contained script – copy‑paste into R ≥ 4.1
########################################################################

## ─────────────────────────────────────────────────────────────────────
##  0 ·  PACKAGES      (install automatically if absent)
## ─────────────────────────────────────────────────────────────────────
pkgs <- c("data.table", "Matrix", "isotone", "KMsurv", "ggplot2")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  message("→ installing ‘", p, "’ …"); install.packages(p, deps = TRUE)
}
invisible(lapply(pkgs, require, character.only = TRUE))

## ─────────────────────────────────────────────────────────────────────
##  1 ·  GLOBAL CONFIGURATION
## ─────────────────────────────────────────────────────────────────────
CFG <- list(
  tau_max = 25,      # administrative horizon (years)
  eps     = 1e-8,    # numerical floor for rectangles / EM
  w_cap   = 25       # capping factor for IPC weights
)

width <- 0.25        # ± quarter‑year artificial uncertainty

## ─────────────────────────────────────────────────────────────────────
##  2 ·  LOAD & BUILD RECTANGLES  (± width around exact times)
## ─────────────────────────────────────────────────────────────────────
data(aids, package = "KMsurv")   # 295 × 3 (infect, induct, adult)
setDT(aids)

aids[, :=(
  L1   = pmax(infect - width, 0),
  R1   = infect + width,
  LY   = fifelse(is.na(induct), infect + width + CFG$eps,
                 infect + pmax(induct - width, 0)),
  RY   = fifelse(is.na(induct), CFG$tau_max,
                 infect + induct + width),
  delta1 = 0L,                           # no exact T1 after widening
  delta2 = as.integer(!is.na(induct))    # exact Y if AIDS observed
)]

## force  R1 < LY  strictly
aids[LY <= R1 + CFG$eps, LY := R1 + CFG$eps]

## ─────────────────────────────────────────────────────────────────────
##  3 ·  GEOMETRY : Turnbull incidence matrix
## ─────────────────────────────────────────────────────────────────────
incidenceMatrix <- function(rect, tau_max) {
  xs <- sort(unique(c(0, rect$L1, rect$R1, tau_max)))
  zs <- sort(unique(c(0, rect$LY, rect$RY, tau_max)))
  cells <- CJ(ix = seq_len(length(xs) - 1L),
              iz = seq_len(length(zs) - 1L), sorted = FALSE)[,
                                                             :=(x_lo = xs[ix], x_hi = xs[ix+1L],
                                                                z_lo = zs[iz], z_hi = zs[iz+1L])][x_hi < z_lo][, c("ix","iz") := NULL]

  idx <- lapply(seq_len(nrow(rect)), \(i)
                which(cells$x_lo >= rect$L1[i] & cells$x_hi <= rect$R1[i] &
                        cells$z_lo >= rect$LY[i] & cells$z_hi <= rect$RY[i]))

  A <- sparseMatrix(i = rep(seq_along(idx), lengths(idx)),
                    j = unlist(idx), x = 1,
                    dims = c(nrow(rect), nrow(cells)))
  list(A = A, cells = cells)
}

inc <- incidenceMatrix(aids, CFG$tau_max)

## ─────────────────────────────────────────────────────────────────────
##  4 ·  WEIGHTS  (simple KM IPCW on δ₂)
## ─────────────────────────────────────────────────────────────────────
KMcore <- function(RY, v) {
  ord <- order(RY, -v); n <- length(ord); risk <- n:1
  surv <- cumprod(1 - v[ord] / risk)
  w    <- v[ord] / risk * c(1, head(surv, -1))
  out  <- numeric(n); out[ord] <- w; out
}
weightKM <- function(RY, d2, cfg = CFG) {
  w <- KMcore(RY, d2); cap <- cfg$w_cap * mean(w); pmin(w, cap)
}
W_KM <- weightKM(aids$RY, aids$delta2, CFG)

## ─────────────────────────────────────────────────────────────────────
##  5 ·  2‑D ISOTONIC PROJECTION  (safe PAVA)
## ─────────────────────────────────────────────────────────────────────
pava_safe <- function(x){
  ns <- asNamespace("isotone")
  if (exists("pava", envir = ns, inherits = FALSE))
    get("pava", envir = ns)(x)
  else isotone::gpava(seq_along(x), x)$x
}
ICMproj <- function(u, cells, tol = 1e-10){
  ux <- sort(unique(cells$x_hi)); uz <- sort(unique(cells$z_hi))
  m <- length(ux); n <- length(uz)
  grid <- matrix(NA_real_, m, n); idx <- matrix(NA_integer_, m, n)
  for (k in seq_len(nrow(cells))){
    i <- match(cells$x_hi[k], ux); j <- match(cells$z_hi[k], uz)
    grid[i,j] <- u[k]; idx[i,j] <- k
  }
  repeat{
    changed <- FALSE
    for (i in seq_len(m)){
      g <- !is.na(grid[i,]); if(sum(g)>1){
        new <- pava_safe(grid[i,g])
        if(max(abs(new-grid[i,g]))>tol){ grid[i,g]<-new; changed<-TRUE }
      }
    }
    for (j in seq_len(n)){
      g <- !is.na(grid[,j]); if(sum(g)>1){
        new <- pava_safe(grid[g,j])
        if(max(abs(new-grid[g,j]))>tol){ grid[g,j]<-new; changed<-TRUE }
      }
    }
    if(!changed) break
  }
  u_new <- u
  for (i in seq_len(m)) for (j in seq_len(n))
    if (!is.na(idx[i,j])) u_new[idx[i,j]] <- grid[i,j]
  u_new[u_new < 0] <- 0; u_new / sum(u_new)
}

## ─────────────────────────────────────────────────────────────────────
##  6 ·  ROBUST EM–ICM
## ─────────────────────────────────────────────────────────────────────
emICM <- function(rect, W, cfg = CFG, maxit = 500, tol = 1e-8){
  inc <- incidenceMatrix(rect, cfg$tau_max); A <- inc$A; cells <- inc$cells
  K <- ncol(A); if(K == 0L) stop("no admissible cells")
  W <- W / mean(W); u <- rep(1/K, K)
  for(it in seq_len(maxit)){
    d  <- pmax(as.numeric(A %*% u), cfg$eps)
    u1 <- u * as.numeric(Matrix::crossprod(A, W/d))
    if(any(!is.finite(u1) | u1<0)) u1[!is.finite(u1)|u1<0] <- cfg$eps
    u1 <- u1 / sum(u1); u1 <- ICMproj(u1, cells)
    if(sum(abs(u1-u)) < tol) break
    u <- u1
  }
  list(u = u, iter = it, converged = (it < maxit), cells = cells)
}

## ─────────────────────────────────────────────────────────────────────
##  7 ·  FITS
## ─────────────────────────────────────────────────────────────────────
fit_U  <- emICM(aids, W = rep(1, nrow(aids)), cfg = CFG)
fit_KM <- emICM(aids, W = W_KM,                cfg = CFG)

cat("\nUnweighted EM–ICM :", fit_U$iter, "iterations – converged",
    fit_U$converged,
    "\nKM‑weighted      :", fit_KM$iter, "iterations – converged",
    fit_KM$converged, "\n")

## ─────────────────────────────────────────────────────────────────────
##  8 ·  QUICK INSPECTION
## ─────────────────────────────────────────────────────────────────────
cat("\n #cells =", length(fit_U$u),
    "\n 5 largest masses (U):\n",  head(sort(fit_U$u,  decreasing = TRUE), 5),
    "\n 5 largest masses (KM):\n", head(sort(fit_KM$u, decreasing = TRUE), 5), "\n")

## ─────────────────────────────────────────────────────────────────────
##  9 ·  OPTIONAL HEAT‑MAP   (requires ggplot2)
## ─────────────────────────────────────────────────────────────────────
plot_heat <- function(fit, title){
  d <- copy(fit$cells)[, u := fit$u]
  d[, :=(cx = (x_lo+x_hi)/2, cz = (z_lo+z_hi)/2)]
  ggplot(d, aes(cx, cz, fill = u)) +
    geom_tile() + scale_fill_viridis_c() +
    coord_equal() + theme_bw() +
    labs(title = title, x = expression(T[1]~"(years)"),
         y = expression(Y~"(years)"), fill = "mass")
}
