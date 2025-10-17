## make_figures.R
## Generate paper-style figures from QUANTINEMO outputs (and a WF schematic).
## Run from the directory that contains the "spatpg/" folder.

## --------- configuration ---------
base_dir <- "spatpg"              # folder that contains simout*/ directories
dir.create("figures", showWarnings = FALSE)

## --------- helpers (no external packages) ---------
rescale01 <- function(x) { (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) }

vec_clean <- function(m) {
  v <- as.numeric(unlist(m[ , -(1:2), drop = FALSE]))
  v[is.finite(v)]
}

make_long <- function(df) {
  nm <- colnames(df)
  locs <- nm[(which(nm == "env")+1):ncol(df)]
  data.frame(
    env  = rep(df$env, times = length(locs)),
    dp   = as.numeric(unlist(df[ , locs, drop = FALSE])),
    locus = rep(locs, each = nrow(df)),
    stringsAsFactors = FALSE
  )
}

## --------- pick a replicate to visualize ---------
sim_dirs <- Sys.glob(file.path(base_dir, "simout*"))
if (length(sim_dirs) == 0) stop("No simout*/ directories found under: ", base_dir)
repdir <- sim_dirs[length(sim_dirs)]  # latest by name
message(sprintf("Using replicate directory: %s", repdir))

## --------- load required files from the replicate ---------
need_files <- c("quantiLocusSelection.txt","ntrlLocusSelection.txt",
                "quantiDp.txt","ntrlDp.txt","fitness.txt")
missing <- need_files[!file.exists(file.path(repdir, need_files))]
if (length(missing) > 0) stop("Missing required files in ", repdir, ": ", paste(missing, collapse=", "))

ql <- read.table(file.path(repdir, "quantiLocusSelection.txt"), header = FALSE)
nl <- read.table(file.path(repdir, "ntrlLocusSelection.txt"),   header = FALSE)
colnames(ql)[1:2] <- c("pop","gen"); colnames(nl)[1:2] <- c("pop","gen")

dq <- read.table(file.path(repdir, "quantiDp.txt"), header = FALSE)
dn <- read.table(file.path(repdir, "ntrlDp.txt"),   header = FALSE)
colnames(dq)[1:3] <- c("pop","gen","env")
colnames(dn)[1:3] <- c("pop","gen","env")

fit <- read.table(file.path(repdir, "fitness.txt"), header = FALSE)
colnames(fit) <- c("pop","gen","fitness")

## =========================================
## Fig. 2: WF schematic fitness curves
## =========================================
env <- seq(-1, 1, length.out = 401)

## Linear model: a1a1 high at negative env; a2a2 high at positive env; hetero constant
lin_a1a1 <- 0.8 + 0.2 * rescale01(-env)
lin_a1a2 <- rep(0.9, length(env))
lin_a2a2 <- 0.8 + 0.2 * rescale01(+env)

## Nonlinear logistic (slope ~10), hetero constant
logi <- function(x, k = 10) 1 / (1 + exp(-k * x))
nl_a2a2 <- 0.8 + 0.2 * logi(env, 10)
nl_a1a1 <- rev(nl_a2a2)
nl_a1a2 <- rep(0.9, length(env))

png("figures/Fig2_fitness_curves.png", width = 1200, height = 700, res = 140)
par(mfrow = c(1,2), mar = c(4,4,1,1))
plot(env, lin_a1a1, type="l", lwd=2, ylim=c(0.78,1.02), xlab="Environmental value", ylab="Fitness")
lines(env, lin_a1a2, lwd=2, lty=2)
lines(env, lin_a2a2, lwd=2)
legend("bottomright", c("a1a1","a1a2","a2a2"), lty=c(1,2,1), lwd=2, bty="n", title="Linear")
plot(env, nl_a1a1, type="l", lwd=2, ylim=c(0.78,1.02), xlab="Environmental value", ylab="Fitness")
lines(env, nl_a1a2, lwd=2, lty=2)
lines(env, nl_a2a2, lwd=2)
legend("bottomright", c("a1a1","a1a2","a2a2"), lty=c(1,2,1), lwd=2, bty="n", title="Nonlinear (logistic, slope≈10)")
dev.off()

## =========================================
## Fig. 3-like panel from current quantiNemo outputs
## (neutral = grey; non-neutral/quanti = blue; dashed line at 0)
## =========================================

# Flatten per-locus slopes
flatten_slopes <- function(df) {
  vv <- as.numeric(unlist(df[ , -(1:2), drop = FALSE]))
  vv[is.finite(vv)]
}
s_quanti <- flatten_slopes(ql)  # non-neutral
s_neutr  <- flatten_slopes(nl)  # neutral

# Build a single index axis: put the highlighted (quanti) first, then all neutrals
y_vals   <- c(s_quanti, s_neutr)
group    <- c(rep("quanti", length(s_quanti)), rep("neutral", length(s_neutr)))
x_idx    <- seq_along(y_vals)

png("figures/Fig3_like_quantiNemo.png", width = 1400, height = 400, res = 140)
par(mar = c(4,4,2,1))
plot(x_idx, y_vals, type = "n", xlab = "Parameter number", ylab = "Selection coefficient",
     main = "Fluctuating selection (quantiNemo)", ylim = range(y_vals, na.rm=TRUE))
abline(h = 0, lty = 2)
# neutrals (grey, semi-transparent)
ne_i <- which(group == "neutral")
points(x_idx[ne_i], y_vals[ne_i], pch = 16, col = rgb(0.5,0.5,0.5,0.5), cex = 0.6)
# quanti (blue/black)
qu_i <- which(group == "quanti")
points(x_idx[qu_i], y_vals[qu_i], pch = 16, col = "blue", cex = 0.7)
legend("topright", c("neutral","quanti (non-neutral)"),
       pch = 16, col = c(rgb(0.5,0.5,0.5,0.8), "blue"), bty = "n")
dev.off()


## =========================================
## Fig. 4(a): distributions of selection coefficients
## =========================================
sel_quanti <- vec_clean(ql)
sel_neutr  <- vec_clean(nl)

message(sprintf("[Fig4a] |s| quanti mean=%.4f sd=%.4f; |s| neutral mean=%.4f sd=%.4f",
                mean(abs(sel_quanti), na.rm=TRUE), sd(abs(sel_quanti), na.rm=TRUE),
                mean(abs(sel_neutr),  na.rm=TRUE), sd(abs(sel_neutr),  na.rm=TRUE)))

png("figures/Fig4a_sel_coeff_distributions.png", width = 1200, height = 800, res = 140)
par(mfrow = c(1,2), mar = c(4,4,2,1))
plot(density(sel_quanti[is.finite(sel_quanti)]), main="Selection coefficient (quanti loci)", xlab="s", lwd=2)
abline(v=0, lty=3)
plot(density(sel_neutr[is.finite(sel_neutr)]),   main="Selection coefficient (neutral loci)", xlab="s", lwd=2)
abline(v=0, lty=3)
dev.off()

## =========================================
## Fig. 4(b): distributions of |Δp|
## =========================================
dp_qu <- abs(as.numeric(unlist(dq[ , -(1:3), drop = FALSE])))
dp_ne <- abs(as.numeric(unlist(dn[ , -(1:3), drop = FALSE])))

message(sprintf("[Fig4b] |Δp| quanti mean=%.4f sd=%.4f; |Δp| neutral mean=%.4f sd=%.4f",
                mean(dp_qu, na.rm=TRUE), sd(dp_qu, na.rm=TRUE),
                mean(dp_ne, na.rm=TRUE), sd(dp_ne, na.rm=TRUE)))

png("figures/Fig4b_delta_p_distributions.png", width = 1200, height = 800, res = 140)
par(mfrow = c(1,2), mar = c(4,4,2,1))
plot(density(dp_qu[is.finite(dp_qu)]), main="Δp distribution (quanti loci)", xlab="|Δp|", lwd=2)
plot(density(dp_ne[is.finite(dp_ne)]), main="Δp distribution (neutral loci)", xlab="|Δp|", lwd=2)
dev.off()

## =========================================
## Fig. 4(c): Δp vs environment (scatter + OLS line)
## =========================================
long_q <- make_long(dq)
long_n <- make_long(dn)

## Subsample to keep file sizes reasonable for huge runs
set.seed(1)
ss <- function(df, n=40000) if (nrow(df) > n) df[sample.int(nrow(df), n), ] else df
plot_q <- ss(long_q); plot_n <- ss(long_n)

## Mean per-locus correlation (descriptive, like in the text)
perloc_cor <- function(df) {
  ## compute cor for each locus; be careful to align env with each slice
  locs <- unique(df$locus)
  cors <- numeric(length(locs))
  for (k in seq_along(locs)) {
    idx <- which(df$locus == locs[k])
    cors[k] <- suppressWarnings(cor(df$env[idx], df$dp[idx], use="complete.obs"))
  }
  mean(cors, na.rm=TRUE)
}
cq <- perloc_cor(long_q)
cn <- perloc_cor(long_n)
message(sprintf("[Fig4c] mean cor(Δp, env): quanti=%.3f neutral=%.3f", cq, cn))

png("figures/Fig4c_delta_p_vs_env.png", width = 1200, height = 600, res = 140)
par(mfrow = c(1,2), mar = c(4,4,2,1))
smoothScatter(plot_q$env, plot_q$dp, xlab="Environment (optimum)", ylab="Δp (quanti)", main="Δp ~ environment (quanti)")
abline(lm(dp ~ env, data = plot_q), lwd=2, lty=2)
smoothScatter(plot_n$env, plot_n$dp, xlab="Environment (optimum)", ylab="Δp (neutral)", main="Δp ~ environment (neutral)")
abline(lm(dp ~ env, data = plot_n), lwd=2, lty=2)
dev.off()

## =========================================
## Fig. 4(d): fitness over time (one population) + rescaled optimum
## =========================================
pop_id <- sort(unique(fit$pop))[1]
fit1 <- fit[fit$pop == pop_id, , drop=FALSE]

## bins
cutfit <- cut(fit1$fitness, breaks = c(-Inf, 0.5, 0.7, 0.9, Inf),
              labels = c("w < 0.5","0.5 < w < 0.7","0.7 < w < 0.9","w > 0.9"))
fit1$bin <- cutfit

## extract the env (optimum) series for that pop from dq
env_series <- aggregate(env ~ gen + pop, data = dq, FUN = function(x) x[1])
env1 <- env_series[env_series$pop == pop_id, , drop=FALSE]

png("figures/Fig4d_fitness_over_time.png", width = 1200, height = 600, res = 140)
par(mar = c(4,4,2,1))
cols <- c("#d73027","#fc8d59","#fee08b","#1a9850")
plot(jitter(fit1$gen, amount=0.15), fit1$fitness,
     pch = 16, col = cols[as.integer(fit1$bin)], xlab="Generation", ylab="Fitness",
     main = sprintf("Fitness over time (pop %s)", pop_id))
legend("topleft", legend = levels(fit1$bin), pch=16, col = cols, bty="n", cex=0.9)

## overlay rescaled env line
rng <- range(fit1$fitness, na.rm=TRUE)
env_scaled <- rescale01(env1$env) * diff(rng) + min(rng)
lines(env1$gen, env_scaled, lwd=3)
mtext("Line = environmental optimum (rescaled)", side=3, line=0.3, cex=0.9)
dev.off()


## =========================================
## Fig.5-like: selection coefficients by "parameter number"
##  - neutral loci grey; quanti (non-neutral) blue; dashed line at 0
## =========================================

flatten_slopes <- function(df) {
  vv <- as.numeric(unlist(df[ , -(1:2), drop = FALSE]))
  vv[is.finite(vv)]
}
s_quanti <- flatten_slopes(ql)  # from quantiLocusSelection.txt
s_neutr  <- flatten_slopes(nl)  # from ntrlLocusSelection.txt

y_vals <- c(s_quanti, s_neutr)
grp    <- c(rep("quanti", length(s_quanti)), rep("neutral", length(s_neutr)))
x_idx  <- seq_along(y_vals)

png("figures/Fig5_like_quantiNemo.png", width = 1600, height = 420, res = 140)
par(mar = c(4,4,2,1))
plot(x_idx, y_vals, type = "n", xlab = "Parameter number", ylab = "Selection coefficient",
     main = "QUANTINEMO (per-locus slopes of fitness ~ genotype)")
abline(h = 0, lty = 2)
ne_i <- which(grp == "neutral")
points(x_idx[ne_i], y_vals[ne_i], pch = 16, col = rgb(0.5,0.5,0.5,0.5), cex = 0.6)
qu_i <- which(grp == "quanti")
points(x_idx[qu_i], y_vals[qu_i], pch = 16, col = "blue", cex = 0.7)
legend("topright", c("neutral loci","quanti (non-neutral) loci"),
       pch = 16, col = c(rgb(0.5,0.5,0.5,0.8), "blue"), bty = "n")
dev.off()

## =========================================
## Fig.6-like: "true" vs "estimated" per-locus coefficients
##  - true  := average slope of fitness ~ genotype per locus across pops/gens
##  - est   := slope of Δp_locus ~ env across pops/gens
## =========================================

# number of loci per class
Lq <- ncol(ql) - 2
Ln <- ncol(nl) - 2

# TRUE per-locus selection: mean slope over pop/gen
true_qu <- colMeans(ql[ , 3:ncol(ql)], na.rm = TRUE)  # length Lq
true_ne <- colMeans(nl[ , 3:ncol(nl)], na.rm = TRUE)  # length Ln

# EST per-locus "b-hat": slope Δp ~ env across all pop/gens
est_from_dp <- function(dp_df) {
  est <- numeric(ncol(dp_df) - 3)
  for (k in 4:ncol(dp_df)) {
    y <- dp_df[[k]]
    x <- dp_df$env
    ok <- is.finite(y) & is.finite(x)
    if (sum(ok) >= 3) {
      est[k-3] <- coef(lm(y[ok] ~ x[ok]))[2]
    } else {
      est[k-3] <- NA_real_
    }
  }
  est
}
est_qu <- est_from_dp(dq)  # length Lq
est_ne <- est_from_dp(dn)  # length Ln

# Combine for plotting and correlations
true_all <- c(true_qu, true_ne)
est_all  <- c(est_qu,  est_ne)
cls      <- c(rep("quanti", Lq), rep("neutral", Ln))

# correlations
r_all <- suppressWarnings(cor(true_all, est_all, use = "complete.obs"))
r_qu  <- suppressWarnings(cor(true_qu,  est_qu,  use = "complete.obs"))

png("figures/Fig6_like_true_vs_est.png", width = 900, height = 900, res = 140)
par(mar = c(4,4,1,1))
plot(true_all, est_all, pch = 16, col = rgb(0.5,0.5,0.5,0.6),
     xlab = "True (mean slope: fitness ~ genotype)",
     ylab = "Estimated (slope: Δp ~ environment)",
     main = "")
points(true_qu, est_qu, pch = 16, col = "blue")
abline(0, 1, lty = 2)  # one-to-one reference
legend("topleft",
       legend = c(sprintf("r (all) = %.2f", r_all), sprintf("r (non-neutral) = %.2f", r_qu)),
       bty = "n")
dev.off()

message(sprintf("[Fig6-like] r(all)=%.3f; r(non-neutral)=%.3f", r_all, r_qu))


message("Figures written to ./figures")
