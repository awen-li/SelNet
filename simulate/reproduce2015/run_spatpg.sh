#!/usr/bin/env bash
set -euo pipefail

BASE_DIR=`pwd`
SPATPG_DIR="${BASE_DIR}/spatpg"

if ! command -v spatpg >/dev/null 2>&1; then
  echo "==> Building spatpg from source"
  cd $SPATPG_DIR/source && make && cp spatpg /usr/bin/ && cd ../../
fi

# =========================
# Config (overridable)
# =========================
RUN_PL="${RUN_PL:-spatpg/runSims.pl}"
QNEMO_BIN_DEFAULT="${QNEMO_BIN_DEFAULT:-quantinemo}"   # quantiNemo executable name
SPATPG_BIN_DEFAULT="${SPATPG_BIN_DEFAULT:-spatpg}"     # spatpg executable name

# Where Dryad-style sims live (if present)
QNEMO_SIMS_DIR="${QNEMO_SIMS_DIR:-spatpg/qnemosims}"
WF_SIMS_DIR="${WF_SIMS_DIR:-spatpg/wfsims}"

# Where to deposit Bayesian outputs
BAYES_OUT_DIR="${BAYES_OUT_DIR:-spatpg/spatpg_out}"

# How to call spatpg.
#   {GENO} -> path to genotype time-series file
#   {ENV}  -> path to environment file
#   {OUT}  -> output file (HDF5 recommended)
SPATPG_ARGS_TEMPLATE='-g {GENO} -e {ENV} -o {OUT} -n 20000 -b 5000 -t 10 -l 20 -u 4000 -p 0.05'
echo "==> Using spatpg args: $SPATPG_ARGS_TEMPLATE"

# Flags (default behavior)
DO_SIM=1            # run quantiNemo/Perl pipeline
DO_CHECK=1          # quick sanity checks
DO_BAYES=1          # run spatpg if available
DO_FIGS=1           # run figure script (point-estimate figs)
DO_BAYES_FIGS=1     # try to extract posteriors (needs R + hdf5r)

OWNER_GROUP="${OWNER_GROUP:-}"   # e.g. "assl:assl" to chown after run

# =========================
# Helpers
# =========================
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -b BIN      quantiNemo executable (default: ${QNEMO_BIN_DEFAULT})
  -p FILE     Perl pipeline (default: ${RUN_PL})
  -s BIN      spatpg executable (default: ${SPATPG_BIN_DEFAULT})
  --qnemo DIR path to qnemosims (default: ${QNEMO_SIMS_DIR})
  --wf DIR    path to wfsims     (default: ${WF_SIMS_DIR})
  --out DIR   spatpg outputs dir (default: ${BAYES_OUT_DIR})
  --args STR  spatpg arg template (overrides SPATPG_ARGS_TEMPLATE env)
  --no-sim    skip quantiNemo/Perl pipeline
  --no-check  skip quick sanity checks
  --no-bayes  skip spatpg runs
  --no-figs   skip point-estimate figures
  --no-bayes-figs  skip posterior-figure extraction
  -o USER:GRP chown -R on outputs (optional)
  -h | --help show this help

Examples:
  $0 -b quantiNemo2 -s /usr/local/bin/spatpg --args '-g {GENO} -e {ENV} -o {OUT} -m env -n 1000'
  QNEMO_BIN=quantiNemo2 SPATPG_ARGS_TEMPLATE='--geno {GENO} --env {ENV} --out {OUT} --model env' $0
EOF
}

# replace placeholders in template
tmpl() {
  local tpl="$1" GENO="$2" ENV="$3" OUT="$4"
  tpl="${tpl//\{GENO\}/${GENO}}"
  tpl="${tpl//\{ENV\}/${ENV}}"
  tpl="${tpl//\{OUT\}/${OUT}}"
  echo "$tpl"
}

have() { command -v "$1" >/dev/null 2>&1; }

# =========================
# Parse args
# =========================
QNEMO_BIN="$QNEMO_BIN_DEFAULT"
SPATPG_BIN="$SPATPG_BIN_DEFAULT"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -b) QNEMO_BIN="$2"; shift 2 ;;
    -p) RUN_PL="$2"; shift 2 ;;
    -s) SPATPG_BIN="$2"; shift 2 ;;
    --qnemo) QNEMO_SIMS_DIR="$2"; shift 2 ;;
    --wf) WF_SIMS_DIR="$2"; shift 2 ;;
    --out) BAYES_OUT_DIR="$2"; shift 2 ;;
    --args) SPATPG_ARGS_TEMPLATE="$2"; shift 2 ;;
    --no-sim) DO_SIM=0; shift ;;
    --no-check) DO_CHECK=0; shift ;;
    --no-bayes) DO_BAYES=0; shift ;;
    --no-figs) DO_FIGS=0; shift ;;
    --no-bayes-figs) DO_BAYES_FIGS=0; shift ;;
    -o) OWNER_GROUP="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

# =========================
# Preflight
# =========================
command -v perl >/dev/null || { echo "ERROR: perl not found"; exit 1; }
command -v R >/dev/null    || { echo "ERROR: R not found"; exit 1; }

# quantiNemo (optional if --no-sim)
if [[ $DO_SIM -eq 1 ]]; then
  have "$QNEMO_BIN" || { echo "ERROR: cannot find '$QNEMO_BIN' (quantiNemo). Use -b."; exit 1; }
fi

# spatpg (optional if --no-bayes)
if [[ $DO_BAYES -eq 1 ]]; then
  if ! have "$SPATPG_BIN"; then
    echo "WARN: spatpg binary '$SPATPG_BIN' not found; skipping Bayesian stage."
    DO_BAYES=0
  fi
fi

echo "==> quantiNemo bin: $(command -v "$QNEMO_BIN" 2>/dev/null || echo '(skipped)')"
echo "==> spatpg bin:     $(command -v "$SPATPG_BIN" 2>/dev/null || echo '(skipped)')"
echo "==> Perl wrapper:    $RUN_PL"
echo "==> Workdir:         $(pwd)"

# =========================
# 1) Run quantiNemo + Perl pipeline
# =========================
if [[ $DO_SIM -eq 1 ]]; then
  export QNEMO_BIN   # your Perl script reads QNEMO_BIN if you added that support
  perl "$RUN_PL"
  echo "==> Perl wrapper finished."
fi

# Fix ownership (optional)
if [[ -n "$OWNER_GROUP" ]]; then
  echo "==> chown -R ${OWNER_GROUP} on spatpg/simout*/ (if present)"
  chown -R "$OWNER_GROUP" spatpg/simout* 2>/dev/null || true
fi

# latest replicate dir (if exists)
LATEST_DIR=""
if ls -d spatpg/simout*/ >/dev/null 2>&1; then
  LATEST_DIR="$(ls -d spatpg/simout*/ | sed 's:/$::' | sort -V | tail -n1)"
  echo "==> Latest replicate: ${LATEST_DIR}"
fi

# =========================
# 2) Quick sanity checks
# =========================
if [[ $DO_CHECK -eq 1 && -n "$LATEST_DIR" ]]; then
  echo "==> Quick sanity check on ${LATEST_DIR}"
  for f in fitness.txt quantiGenotypes.txt ntrlGenotypes.txt quantiDp.txt ntrlDp.txt; do
    if [[ -f "${LATEST_DIR}/${f}" ]]; then
      echo -n "lines ${LATEST_DIR}/${f}: "
      wc -l "${LATEST_DIR}/${f}" | awk '{print $1}'
    else
      echo "MISSING: ${LATEST_DIR}/${f}"
    fi
  done
fi

# =========================
# 3) Run point-estimate figures (your make_figures.R)
# =========================
if [[ $DO_FIGS -eq 1 ]]; then
  if [[ -f make_figures.R ]]; then
    echo "==> Running make_figures.R (point-estimate visuals)"
    Rscript make_figures.R || true
  else
    echo "NOTE: make_figures.R not found; skipping point-estimate figures."
  fi
fi

# =========================
# 4) Run spatpg on Dryad-style sims (if present)
# =========================
run_spatpg_pair() {
  local GENO="$1" ENV="$2" OUT="$3"
  mkdir -p "$(dirname "$OUT")"
  local args; args="$(tmpl "$SPATPG_ARGS_TEMPLATE" "$GENO" "$ENV" "$OUT")"
  echo "   [spatpg] $SPATPG_BIN $args"
  set +e
  $SPATPG_BIN $args
  local rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo "WARN: spatpg failed on GENO=$(basename "$GENO") ENV=$(basename "$ENV") (exit $rc)"
  fi
}

if [[ $DO_BAYES -eq 1 ]]; then
  mkdir -p "$BAYES_OUT_DIR"

  # ---- QUANTINEMO sims (from Dryad) ----
  if [[ -d "$QNEMO_SIMS_DIR" ]]; then
    echo "==> Running spatpg on qnemosims: $QNEMO_SIMS_DIR"
    for i in $(seq 0 9); do
      GENO="${QNEMO_SIMS_DIR}/geno_sim${i}.txt"
      ENVF="${QNEMO_SIMS_DIR}/env_sim${i}.txt"
      [[ -f "$GENO" && -f "$ENVF" ]] || continue
      OUT="${BAYES_OUT_DIR}/qnemo_sim${i}.h5"
      run_spatpg_pair "$GENO" "$ENVF" "$OUT"
    done
  else
    echo "NOTE: qnemosims dir not found at $QNEMO_SIMS_DIR (skip)."
  fi

  # ---- Wright–Fisher sims (from Dryad) ----
  if [[ -d "$WF_SIMS_DIR" ]]; then
    echo "==> Running spatpg on wfsims: $WF_SIMS_DIR"
    # Expect pairs like geno_*.txt + env_*.txt; adjust glob if your names differ
    for GENO in "$WF_SIMS_DIR"/geno_*.txt; do
      [[ -f "$GENO" ]] || continue
      base="$(basename "$GENO" .txt | sed 's/^geno_//')"
      ENVF="$WF_SIMS_DIR/env_${base}.txt"
      [[ -f "$ENVF" ]] || { echo "WARN: no env for $GENO"; continue; }
      OUT="${BAYES_OUT_DIR}/wf_${base}.h5"
      run_spatpg_pair "$GENO" "$ENVF" "$OUT"
    done
  else
    echo "NOTE: wfsims dir not found at $WF_SIMS_DIR (skip)."
  fi
fi

# =========================
# 5) OPTIONAL: Extract posterior medians + 95% intervals from spatpg HDF5
#     and draw Fig.3/5-style panels (requires R + hdf5r)
# =========================
if [[ $DO_BAYES -eq 1 && $DO_BAYES_FIGS -eq 1 ]]; then
  if R -q -e 'quit(status=!requireNamespace("hdf5r", quietly=TRUE))' >/dev/null 2>&1; then
    echo "==> Extracting spatpg posteriors (R + hdf5r) and making Fig.3/5-like panels"
    Rscript - <<'RS'
dir.create("figures", showWarnings = FALSE)
library(hdf5r)

# collect all H5 results
h5s <- Sys.glob(file.path("spatpg/spatpg_out", "*.h5"))
if (length(h5s) == 0) quit(save="no")

# minimal reader (adapt to your file layout if keys differ)
read_summary <- function(h5) {
  f <- H5File$new(h5, mode="r")
  on.exit(f$close_all())
  # The group/dataset names may differ; these are common patterns:
  # a_median, a_q05, a_q95, b_median, b_q05, b_q95, st_median, st_q05, st_q95
  out <- list()
  for (nm in c("a_median","a_q05","a_q95","b_median","b_q05","b_q95")) {
    if (f$exists(nm)) out[[nm]] <- f[[nm]]$read()
  }
  # Per-time selection (vector per locus, or matrix [locus x time])
  for (nm in c("st_median","st_q05","st_q95","s_t_median","s_t_q05","s_t_q95")) {
    if (f$exists(nm)) out[[nm]] <- f[[nm]]$read()
  }
  out
}

summ <- lapply(h5s, read_summary)

# Fig.5-like (many loci): we need b medians and intervals
bmed <- lapply(summ, `[[`, "b_median")
bq05 <- lapply(summ, `[[`, "b_q05")
bq95 <- lapply(summ, `[[`, "b_q95")

# Flatten and plot if we have something
if (any(lengths(bmed) > 0)) {
  b <- unlist(bmed, use.names=FALSE)
  lo <- unlist(bq05, use.names=FALSE)
  hi <- unlist(bq95, use.names=FALSE)
  idx <- seq_along(b)

  png("figures/Fig5_spatpg_posteriors.png", width=1600, height=420, res=140)
  par(mar=c(4,4,2,1))
  plot(idx, b, type="n", xlab="Parameter number", ylab="Selection coefficient (posterior)")
  abline(h=0, lty=2)
  segments(idx, lo, idx, hi, col=gray(0.6))
  points(idx, b, pch=16, cex=0.5)
  dev.off()
}

RS
  else
    echo "NOTE: R package 'hdf5r' not installed; skipping posterior extraction/figures."
  fi
fi

echo "==> Done."
