#!/usr/bin/env bash
set -euo pipefail

# ---------- Paths / Build ----------
BASE_DIR="$(pwd)"
SPATPG_DIR="${BASE_DIR}/spatpg"

# Prefer a local bin dir to avoid /usr/bin installs and path confusion
LOCAL_BIN="${BASE_DIR}/.local/bin"
mkdir -p "${LOCAL_BIN}"
export PATH="${LOCAL_BIN}:$PATH"

if ! command -v spatpg >/dev/null 2>&1; then
  echo "==> Building spatpg from source (local install)"
  ( cd "$SPATPG_DIR/source" && make )
  cp -f "$SPATPG_DIR/source/spatpg" "${LOCAL_BIN}/"
fi

# ---------- Config ----------
RUN_PL="${RUN_PL:-spatpg/runSims.pl}"
QNEMO_BIN_DEFAULT="${QNEMO_BIN_DEFAULT:-quantinemo}"
SPATPG_BIN_DEFAULT="${SPATPG_BIN_DEFAULT:-spatpg}"

QNEMO_SIMS_DIR="${QNEMO_SIMS_DIR:-spatpg/qnemosims}"
WF_SIMS_DIR="${WF_SIMS_DIR:-spatpg/wfsims}"
BAYES_OUT_DIR="${BAYES_OUT_DIR:-spatpg/spatpg_out}"

# Add a model knob (2015 paper uses env / temp models)
SPATPG_MODEL="${SPATPG_MODEL:-env}"  # e.g., env|const|temp
SPATPG_ARGS_TEMPLATE="${SPATPG_ARGS_TEMPLATE:- -g {GENO} -e {ENV} -o {OUT} -m ${SPATPG_MODEL} -n 20000 -b 5000 -t 10 -l 20 -u 4000 -p 0.05 }"
echo "==> Using spatpg args: $SPATPG_ARGS_TEMPLATE"

DO_SIM=1
DO_CHECK=1
DO_BAYES=1
DO_FIGS=1
DO_BAYES_FIGS=1

# If available, auto-disable posterior figs when hdf5r is missing
if ! R -q -e 'quit(status=!requireNamespace("hdf5r", quietly=TRUE))' >/dev/null 2>&1; then
  DO_BAYES_FIGS=0
fi

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
  --model M   spatpg model (env|const|temp), default: ${SPATPG_MODEL}
  --no-sim    skip quantiNemo/Perl pipeline
  --no-check  skip quick sanity checks
  --no-bayes  skip spatpg runs
  --no-figs   skip point-estimate figures
  --no-bayes-figs  skip posterior-figure extraction
  -o USER:GRP chown -R outputs (optional)
  -h | --help show this help
EOF
}

tmpl() {
  local tpl="$1" GENO="$2" ENV="$3" OUT="$4"
  tpl="${tpl//\{GENO\}/${GENO}}"
  tpl="${tpl//\{ENV\}/${ENV}}"
  tpl="${tpl//\{OUT\}/${OUT}}"
  echo "$tpl"
}

have() { command -v "$1" >/dev/null 2>&1; }

# ---------- Parse args ----------
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
    --model) SPATPG_MODEL="$2"; shift 2 ;;
    --no-sim) DO_SIM=0; shift ;;
    --no-check) DO_CHECK=0; shift ;;
    --no-bayes) DO_BAYES=0; shift ;;
    --no-figs) DO_FIGS=0; shift ;;
    --no-bayes-figs) DO_BAYES_FIGS=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

# ---------- Preflight ----------
command -v perl >/dev/null || { echo "ERROR: perl not found"; exit 1; }
command -v R >/dev/null    || { echo "ERROR: R not found"; exit 1; }

if [[ $DO_SIM -eq 1 ]]; then
  have "$QNEMO_BIN" || { echo "ERROR: cannot find '$QNEMO_BIN'"; exit 1; }
fi

if [[ $DO_BAYES -eq 1 ]]; then
  have "$SPATPG_BIN" || { echo "WARN: spatpg not found; skipping Bayesian stage."; DO_BAYES=0; }
fi

echo "==> quantiNemo bin: $(command -v "$QNEMO_BIN" 2>/dev/null || echo '(skipped)')"
echo "==> spatpg bin:     $(command -v "$SPATPG_BIN" 2>/dev/null || echo '(skipped)')"
echo "==> Perl wrapper:    $RUN_PL"
echo "==> Workdir:         $(pwd)"

# Optional: reproducible seeds (if your Perl wrapper supports it)
export QNEMO_SEED="${QNEMO_SEED:-12345}"

# ---------- 1) Run quantiNemo + Perl ----------
if [[ $DO_SIM -eq 1 ]]; then
  export QNEMO_BIN
  echo "==> Running Perl wrapper: $RUN_PL"
  perl "$RUN_PL"
  echo "==> Perl wrapper finished."
fi

# Find newest complete replicate dir by mtime and presence of key files
pick_latest_dir() {
  local candidates
  IFS=$'\n' read -r -d '' -a candidates < <(ls -dt spatpg/simout*/ 2>/dev/null | tr -d '\r' && printf '\0' || true)
  for d in "${candidates[@]:-}"; do
    d="${d%/}"
    for f in fitness.txt quantiGenotypes.txt ntrlGenotypes.txt; do
      [[ -s "${d}/${f}" ]] || { continue 2; }
    done
    echo "$d"; return 0
  done
  return 1
}

LATEST_DIR=""
if LATEST_DIR="$(pick_latest_dir)"; then
  echo "==> Latest complete replicate: ${LATEST_DIR}"
else
  echo "NOTE: No complete simout* directory found."
fi

# ---------- 2) Quick sanity checks ----------
if [[ $DO_CHECK -eq 1 && -n "${LATEST_DIR}" ]]; then
  echo "==> Quick sanity check on ${LATEST_DIR}"
  for f in fitness.txt quantiGenotypes.txt ntrlGenotypes.txt quantiDp.txt ntrlDp.txt; do
    if [[ -f "${LATEST_DIR}/${f}" ]]; then
      printf "lines %-28s: %s\n" "${f}" "$(wc -l < "${LATEST_DIR}/${f}")"
    else
      echo "MISSING: ${f}"
    fi
  done
fi

# ---------- 3) Point-estimate figures (R) ----------
if [[ $DO_FIGS -eq 1 ]]; then
  if [[ -f make_figures.R && -n "${LATEST_DIR}" ]]; then
    echo "==> Running make_figures.R with explicit input dir"
    Rscript make_figures.R --in "${LATEST_DIR}" --out "figures" || true
  else
    echo "NOTE: make_figures.R not found or no LATEST_DIR; skipping."
  fi
fi

# ---------- 4) Run spatpg on Dryad sims ----------
run_spatpg_pair() {
  local GENO="$1" ENV="$2" OUT="$3"
  mkdir -p "$(dirname "$OUT")"
  if [[ ! -s "$GENO" ]]; then echo "WARN: empty/missing GENO: $GENO"; return; fi
  if [[ ! -s "$ENV"  ]]; then echo "WARN: empty/missing ENV:  $ENV";  return; fi
  local args; args="$(tmpl "$SPATPG_ARGS_TEMPLATE" "$GENO" "$ENV" "$OUT")"
  echo "   [spatpg] $SPATPG_BIN $args"
  set +e
  "$SPATPG_BIN" $args 2> "${OUT%.h5}.log"
  local rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo "WARN: spatpg failed on $(basename "$GENO") / $(basename "$ENV") (exit $rc)"
  fi
}

if [[ $DO_BAYES -eq 1 ]]; then
  mkdir -p "$BAYES_OUT_DIR"

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

  if [[ -d "$WF_SIMS_DIR" ]]; then
    echo "==> Running spatpg on wfsims: $WF_SIMS_DIR"
    shopt -s nullglob
    for GENO in "$WF_SIMS_DIR"/geno_*.txt; do
      base="$(basename "$GENO" .txt)"; base="${base#geno_}"
      ENVF="$WF_SIMS_DIR/env_${base}.txt"
      [[ -f "$ENVF" ]] || { echo "WARN: no env for $GENO"; continue; }
      OUT="${BAYES_OUT_DIR}/wf_${base}.h5"
      run_spatpg_pair "$GENO" "$ENVF" "$OUT"
    done
    shopt -u nullglob
  else
    echo "NOTE: wfsims dir not found at $WF_SIMS_DIR (skip)."
  fi
fi

# ---------- 5) Posterior extraction (R + hdf5r) ----------
if [[ $DO_BAYES -eq 1 && $DO_BAYES_FIGS -eq 1 ]]; then
  echo "==> Extracting spatpg posteriors to figures/"
  Rscript - <<'RS'
dir.create("figures", showWarnings = FALSE)
suppressMessages({
  library(hdf5r)
})

h5s <- Sys.glob(file.path("spatpg/spatpg_out", "*.h5"))
if (length(h5s) == 0) quit(save="no")

list_keys <- function(f) {
  walk <- function(g, pref="") {
    kids <- g$ls()$name
    out <- character()
    for (k in kids) {
      obj <- g[[k]]
      nm  <- paste0(pref, k)
      if (inherits(obj, "H5Group")) out <- c(out, walk(obj, paste0(nm, "/")))
      else out <- c(out, nm)
    }
    out
  }
  walk(f)
}

read_first_available <- function(f, names) {
  for (nm in names) {
    if (f$exists(nm)) return(f[[nm]]$read())
  }
  NULL
}

b_all <- list(); lo_all <- list(); hi_all <- list()
for (h in h5s) {
  f <- H5File$new(h, mode="r")
  on.exit(f$close_all(), add=TRUE)

  # Common historical variants for 'b' (fluctuating selection)
  b <- read_first_available(f, c(
    "b_median", "posterior/b_median", "summary/b_median"
  ))
  bq05 <- read_first_available(f, c(
    "b_q05", "posterior/b_q05", "summary/b_q05"
  ))
  bq95 <- read_first_available(f, c(
    "b_q95", "posterior/b_q95", "summary/b_q95"
  ))

  if (is.null(b) || is.null(bq05) || is.null(bq95)) {
    message("NOTE: Could not find b-marginals in: ", h, " (available keys:)")
    print(list_keys(f))
    next
  }

  b_all[[h]]  <- as.numeric(b)
  lo_all[[h]] <- as.numeric(bq05)
  hi_all[[h]] <- as.numeric(bq95)
}

if (length(b_all)) {
  b  <- unlist(b_all, use.names=FALSE)
  lo <- unlist(lo_all, use.names=FALSE)
  hi <- unlist(hi_all, use.names=FALSE)
  idx <- seq_along(b)

  png("figures/Fig_b_posteriors.png", width=1600, height=420, res=140)
  par(mar=c(4,4,2,1))
  plot(idx, b, type="n", xlab="Parameter index", ylab="b (posterior)")
  abline(h=0, lty=2)
  segments(idx, lo, idx, hi, col=gray(0.6))
  points(idx, b, pch=16, cex=0.6)
  dev.off()
}
RS
fi

echo "==> Done."
