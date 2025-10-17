#!/usr/bin/env bash
set -euo pipefail

# --- config (you can override via flags/env) ---
RUN_PL="${RUN_PL:-runSims.pl}"
QNEMO_BIN_DEFAULT="${QNEMO_BIN_DEFAULT:-quantinemo}"  # e.g. quantinemo | quantiNemo | quantiNemo2

usage() {
  cat <<EOF
Usage: $0 [-b /path/to/quantiNemo*] [-p runSims.pl] [-o owner:group] [--no-check]
  -b  Path to quantiNemo executable (or name in PATH). Default: ${QNEMO_BIN_DEFAULT}
      (or set env QNEMO_BIN)
  -p  Perl wrapper file. Default: ${RUN_PL}
  -o  chown -R OWNER:GROUP on output dir (e.g. assl:assl)
      (optional; useful if you ran as root)
  --no-check   Skip quick sanity checks after run
Examples:
  $0 -b /usr/local/bin/quantiNemo2
  QNEMO_BIN=quantiNemo2 $0
EOF
}

OWNER_GROUP=""
DO_CHECK=1

# --- arg parsing ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b) QNEMO_BIN="${2:-}"; shift 2 ;;
    -p) RUN_PL="${2:-}"; shift 2 ;;
    -o) OWNER_GROUP="${2:-}"; shift 2 ;;
    --no-check) DO_CHECK=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

# pick bin
QNEMO_BIN="${QNEMO_BIN:-${QNEMO_BIN_DEFAULT}}"
if ! command -v "${QNEMO_BIN}" >/dev/null 2>&1; then
  echo "ERROR: cannot find '${QNEMO_BIN}'. Try -b /path/to/quantiNemo2 or set QNEMO_BIN." >&2
  exit 1
fi

# deps
command -v perl >/dev/null || { echo "ERROR: perl not found"; exit 1; }
command -v R >/dev/null    || { echo "ERROR: R not found"; exit 1; }

# script present?
[[ -f "${RUN_PL}" ]] || { echo "ERROR: ${RUN_PL} not found in $(pwd)"; exit 1; }

echo "==> Using quantiNemo: $(command -v ${QNEMO_BIN})"
echo "==> Perl wrapper:     ${RUN_PL}"
echo "==> Working dir:      $(pwd)"

# run
export QNEMO_BIN  # used by your Perl if you added the variable
perl "${RUN_PL}"

echo "==> Perl wrapper finished."

# optional ownership fix
if [[ -n "${OWNER_GROUP}" ]]; then
  echo "==> chown -R ${OWNER_GROUP} on simout*/"
  chown -R "${OWNER_GROUP}" simout* || true
fi

# find latest simout (highest index)
LATEST_DIR=""
if ls -d simout*/ >/dev/null 2>&1; then
  LATEST_DIR="$(ls -d simout*/ | sed 's:/$::' | sort -V | tail -n1)"
  echo "==> Latest replicate directory: ${LATEST_DIR}"
else
  echo "WARN: no simout*/ directories found. Nothing to check."
  exit 0
fi

# quick checks
if [[ "${DO_CHECK}" -eq 1 ]]; then
  echo "==> Quick sanity check"
  if [[ -f "${LATEST_DIR}/run.log" ]]; then
    echo "---- ${LATEST_DIR}/run.log (head) ----"
    sed -n '1,60p' "${LATEST_DIR}/run.log" || true
    echo "--------------------------------------"
  fi

  for f in fitness.txt quantiGenotypes.txt ntrlGenotypes.txt quantiDp.txt ntrlDp.txt; do
    if [[ -f "${LATEST_DIR}/${f}" ]]; then
      echo -n "lines ${LATEST_DIR}/${f}: "
      wc -l "${LATEST_DIR}/${f}" | awk '{print $1}'
    else
      echo "MISSING: ${LATEST_DIR}/${f}"
    fi
  done

  # Inline R quick stats (Δp magnitudes + correlation with env)
  Rscript - <<'RS' || true
args <- commandArgs(trailingOnly=TRUE)
repdir <- Sys.glob("simout*")[length(Sys.glob("simout*"))]
if (length(repdir)==0) quit(status=0)
dq <- read.table(file.path(repdir,"quantiDp.txt"))
dn <- read.table(file.path(repdir,"ntrlDp.txt"))
colnames(dq) <- c("pop","gen","env", paste0("q", seq_len(ncol(dq)-3)))
colnames(dn) <- c("pop","gen","env", paste0("n", seq_len(ncol(dn)-3)))
dqv <- abs(unlist(dq[,-(1:3)])); dnv <- abs(unlist(dn[,-(1:3)]))
cat(sprintf("[QUICK] |Δp| quanti mean=%.4f sd=%.4f; |Δp| neutral mean=%.4f sd=%.4f\n",
            mean(dqv,na.rm=TRUE), sd(dqv,na.rm=TRUE), mean(dnv,na.rm=TRUE), sd(dnv,na.rm=TRUE)))
cq <- sapply(colnames(dq)[-(1:3)], function(c) suppressWarnings(cor(dq$env, dq[[c]], use="complete.obs")))
cn <- sapply(colnames(dn)[-(1:3)], function(c) suppressWarnings(cor(dn$env, dn[[c]], use="complete.obs")))
cat(sprintf("[QUICK] mean cor(Δp, env): quanti=%.3f neutral=%.3f\n",
            mean(cq,na.rm=TRUE), mean(cn,na.rm=TRUE)))
RS
fi

echo "==> Done."
