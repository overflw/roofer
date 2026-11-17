#!/bin/bash
# score_all.sh AREA run...
S=${EVAL_DIR:?set EVAL_DIR to a work dir with data/, bin/, runs/ and a venv/}
A=$1; shift
for r in "$@"; do
  $S/venv/bin/python -I $(dirname $0)/../score_lod22.py --city $S/runs/$A/$r/out --pc $S/data/$A.laz --id-attribute identificatie --out-json $S/runs/$A/$r/score.json --out-csv $S/runs/$A/$r/b.csv --out-face-csv $S/runs/$A/$r/f.csv > /dev/null 2>&1
done
$S/venv/bin/python - "$A" "$@" <<'PY'
import json,os,sys
A=sys.argv[1]; runs=sys.argv[2:]
S=os.environ['EVAL_DIR']+'/runs'
keys=['faces','faces_level','faces_sliver','faces_nodata','faces_bad','rms_in_face_p50','rms_in_face_p90','outlier_share','bad_m2','missed_m2','splits','crease_m2','usable_m2','off_m2','usable_frac_eligible','buildings_fallback']
d={r:json.load(open(f'{S}/{A}/{r}/score.json'))['summary'] for r in runs}
print(f"{'metric':22s}"+''.join(f'{r:>12s}' for r in runs))
for k in keys:
    print(f'{k:22s}'+''.join(f'{d[r].get(k,float("nan")):12.4g}' for r in runs))
PY
