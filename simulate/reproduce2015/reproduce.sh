#!/usr/bin/env bash
set -euo pipefail

SPATPG="spatpg/source/spatpg"
WF_DIR="spatpg/wfsims"

OUT_DIR="results_wfsims"
mkdir -p "$OUT_DIR"

# -------- build spatprg if needed --------
if [[ ! -x "$SPATPG" ]]; then
  echo "Building spatpg..."
  (cd `dirname $SPATPG`  && make )
fi
# -----------------------------------------


run_one () {
  local G="$1" ; local E="$2" ; local tag="$3"
  echo "==> $tag"
  "$SPATPG" -g "$G" -e "$E" -o "${OUT_DIR}/${tag}.hdf5" \
    -n 10000 -b 2000 -t 10 -l 20 -u 4000 -p 0.05
}

# helper: build env filename from genetic filename
env_for () {
  local g="$1" base env
  base="$(basename "$g")"
  # strip sub50/sub100 if present (env files don't have them)
  base="${base#sub50}"
  base="${base#sub100}"
  # core replacements: G -> E, keep scenario
  env="${base//G/E}"
  # genetic files start with 'mod_', env files are 'rsim*'
  env="${env#mod_}"
  # ensure scenario prefixes match:
  #  Sel -> rsimSelE, Fluct -> rsimFluctE, FluctLE -> rsimFluctLE, WeakFluctLE -> rsimWeakFluctLE, WkSel -> rsimWkSelE
  # For plain 'rsimE' case (mod_rsimG#), the line above already produced rsimE#
  echo "${WF_DIR}/${env}"
}

# gather genetic files we know how to map
mapfile -t GFILES < <(ls -1 "$WF_DIR"/mod_rsim*G*.txt "$WF_DIR"/sub*mod_rsim*G*.txt 2>/dev/null)

for G in "${GFILES[@]}"; do
  E="$(env_for "$G")"
  if [[ ! -f "$E" ]]; then
    echo "WARN: No env match for $G -> $E ; skipping" >&2
    continue
  fi
  tag="$(basename "${G%.txt}")"
  run_one "$G" "$E" "$tag"
done

echo "Done. HDF5 outputs at $OUT_DIR/"
