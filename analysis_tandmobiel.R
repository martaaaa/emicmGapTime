#-----------
## ────────────────────────────────────────────────────────────────────────────
## 0 · PACKAGES
## ────────────────────────────────────────────────────────────────────────────
libs <- c("icensBKL", "data.table", "ggplot2",
          "viridis", "patchwork", "knitr")
invisible(lapply(libs, require, character.only = TRUE))

## ────────────────────────────────────────────────────────────────────────────
## 1 · LOAD & PRE-PROCESS  –   stack all four first molars
## ────────────────────────────────────────────────────────────────────────────
data(tandmob, package = "icensBKL");   setDT(tandmob)

tooth <- rbind(
  tandmob[, .(tooth="16", L1=L16, R1=R16, LY=CL16, RY=CR16)],
  tandmob[, .(tooth="26", L1=L26, R1=R26, LY=CL26, RY=CR26)],
  tandmob[, .(tooth="36", L1=L36, R1=R36, LY=CL36, RY=CR36)],
  tandmob[, .(tooth="46", L1=L46, R1=R46, LY=CL46, RY=CR46)]
)

tooth <- tooth[complete.cases(tooth) & R1 <= LY]          # logical order

## --- strictly positive rectangles --------------------------------------------
eps <- 1e-6                                               # ≈ 0.02 day
tooth[R1 >= LY     , LY := R1 + eps]
tooth[R1 - L1 < eps, R1 := L1 + eps]
tooth[RY - LY < eps, RY := LY + eps]

## --- flags (no NA allowed) ----------------------------------------------------
tooth[, delta1 := as.integer(abs(R1 - L1) < eps)]
tooth[, delta2 := as.integer(abs(RY - LY) < eps)]
tooth[is.na(delta1), delta1 := 0L]
tooth[is.na(delta2), delta2 := 0L]

## --- administrative horizon ---------------------------------------------------
STUDY_HORIZON <- 120   # months
for(v in c("L1","R1","LY","RY"))
  set(tooth, which(!is.finite(tooth[[v]])), v, STUDY_HORIZON)

cat("clean sample :", nrow(tooth), "tooth records\n")

## ────────────────────────────────────────────────────────────────────────────
## 2 · CONFIGURATION
## ────────────────────────────────────────────────────────────────────────────
cfg         <- CFG
cfg$tau_max <- STUDY_HORIZON
cfg$eps     <- 1e-12

## ────────────────────────────────────────────────────────────────────────────
## 3 · IPCW  (weights)
## ────────────────────────────────────────────────────────────────────────────
n <- nrow(tooth)

if(sum(tooth$delta2, na.rm = TRUE) > 0){
  W_km <- weightKM(tooth$RY, tooth$delta2, cfg)
} else {
  warning("no exact caries times – KM weights = 1")
  W_km <- rep(1, n)
}

if(sum(tooth$delta2, na.rm = TRUE) > 5){
  m_hat <- estimate_m(tooth, cfg)               # random-forest
  W_ps  <- weightPS(tooth$RY, m_hat, cfg)
} else {
  W_ps  <- rep(1, n)
}

## ────────────────────────────────────────────────────────────────────────────
## 4 · SAFE EM–ICM
## ────────────────────────────────────────────────────────────────────────────
emICM_safe <- function(rect, w, cfg, tag=""){
  fit <- tryCatch(icensBKL:::emICM(rect, w, cfg),
                  error=function(e){message(tag," fit failed: ",e$message); NULL})
  if(is.null(fit)) return(NULL)
  stopifnot(all(is.finite(fit$u)), all(fit$u>=0),
            abs(sum(fit$u)-1) < 1e-6)
  fit
}

fitKM <- emICM_safe(tooth, W_km, cfg, "KM")
fitPS <- emICM_safe(tooth, W_ps, cfg, "PS")

## ────────────────────────────────────────────────────────────────────────────
## 5 · DIAGNOSTICS
## ────────────────────────────────────────────────────────────────────────────
if(!is.null(fitKM)){
  cat("\nEM–ICM sweeps:\n",
      "KM :", fitKM$iter, "(converged =", fitKM$converged, ")\n",
      "PS :", fitPS$iter, "(converged =", fitPS$converged, ")\n")
} else {
  stop("Both fits failed – please inspect warnings above.")
}

info_tbl <- tooth[, .(
  n_records    = .N,
  pct_exact_T1 = round(100*mean(delta1),1),
  pct_exact_Y  = round(100*mean(delta2),1)
)]
print(kable(info_tbl, caption = "Cleaned first-molar data (all quadrants)"))

## ────────────────────────────────────────────────────────────────────────────
## 6 · HEAT-MAPS
## ────────────────────────────────────────────────────────────────────────────
plot_heat <- function(fit, ttl){
  d <- cbind(fit$cells, u=fit$u)
  d[, :=(cx=(x_lo+x_hi)/2, cz=(z_lo+z_hi)/2)]
  ggplot(d, aes(cx, cz, fill=u)) +
    geom_tile() +
    scale_fill_viridis() +
    coord_equal() +
    labs(title=ttl, x=expression(T[1]~"(months)"),
         y=expression(Y~"(months)")) +
    theme_bw()
}

(p_km <- plot_heat(fitKM, "KM weights")) |
  (p_ps <- plot_heat(fitPS, "Pre-smoothed weights"))

## ────────────────────────────────────────────────────────────────────────────
## 7 · MARGINAL  S_Y(z)
## ────────────────────────────────────────────────────────────────────────────
FhatY <- function(fit, z){
  idx <- which(fit$cells$z_hi <= z)
  if(!length(idx)) return(0)
  cumsum(c(0,fit$u))[tail(idx,1)+1]
}
z_grid <- seq(0, STUDY_HORIZON, length.out = 300)
S_km <- 1 - vapply(z_grid, FhatY, numeric(1), fit=fitKM)
S_ps <- 1 - vapply(z_grid, FhatY, numeric(1), fit=fitPS)

ggplot() +
  geom_line(aes(z_grid, S_km, colour="KM"), linewidth=0.9) +
  geom_line(aes(z_grid, S_ps, colour="PS"), linewidth=0.9, linetype=2) +
  scale_colour_manual(values=c(KM="black", PS="#D55E00"),
                      breaks=c("KM","PS"), name="Weights") +
  labs(x="months", y=expression(hat(S)[Y](z)),
       title="Marginal survival of first-caries onset") +
  theme_bw()
