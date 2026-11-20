#!/bin/bash
# run.sh BIN NAME AREA [extra toml file appended]
S=${EVAL_DIR:?set EVAL_DIR to a work dir with data/, bin/, runs/ and a venv/}
BIN=$1; NAME=$2; AREA=$3; EXTRA=$4
D=$S/runs/$AREA/$NAME; rm -rf $D; mkdir -p $D
sed -e "s#__PC__#$S/data/$AREA.laz#" -e "s#__FP__#$S/data/${AREA}_fp.gpkg#" $S/runs/base.toml > $D/cfg.toml
if [ -n "$EXTRA" ]; then python3 - "$D/cfg.toml" "$EXTRA" <<'PY'
import sys
cfg=open(sys.argv[1]).read(); extra=open(sys.argv[2]).read()
# merge "[section]" blocks of extra into cfg: append keys under matching headers
sec=None; add={}
for line in extra.splitlines():
    line=line.strip()
    if not line: continue
    if line.startswith('['): sec=line; continue
    add.setdefault(sec,[]).append(line)
out=[];cur=None
lines=cfg.splitlines()
res=[]
for i,l in enumerate(lines):
    res.append(l)
    if l.strip().startswith('['):
        cur=l.strip()
        if cur in add: res+=add.pop(cur)
for s,ls in add.items(): res+=['',s]+ls
out=[];seen={};sec=None
for l in res:
    t=l.strip()
    if t.startswith('['): sec=t
    elif '=' in t:
        k=(sec,t.split('=')[0].strip())
        if k in seen: continue
        seen[k]=len(out)
    out.append(l)
open(sys.argv[1],'w').write('\n'.join(l for l in out if l is not None)+'\n')
PY
fi
cd $D && $S/bin/$BIN --config cfg.toml -j ${JOBS:-2} $D/out > log.txt 2>&1
echo "time $SECONDS s" >> log.txt; tail -3 log.txt
