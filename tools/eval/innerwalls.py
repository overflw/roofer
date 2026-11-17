"""innerwalls.py RUNDIR... : count inner (not on footprint edge) wall surfaces in lod2.2 solids:
n, area, small (<2 m2) inner walls (eave-kink proxy)."""
import sys, json, glob, numpy as np
def area3d(r):
    r=np.asarray(r,float); c=r.mean(0); s=np.zeros(3)
    for i in range(len(r)): s+=np.cross(r[i]-c, r[(i+1)%len(r)]-c)
    return np.linalg.norm(s)/2
for run in sys.argv[1:]:
    n=a=small=tri=0
    for f in glob.glob(run+'/out/*.city.jsonl'):
        L=open(f).read().splitlines(); tr=json.loads(L[0])['transform']; sc=np.array(tr['scale']); tl=np.array(tr['translate'])
        for l in L[1:]:
            ft=json.loads(l); V=np.array(ft['vertices'])*sc+tl
            for o in ft['CityObjects'].values():
                for g in o.get('geometry',[]):
                    if str(g.get('lod'))!='2.2' or g['type']!='Solid': continue
                    sem=g['semantics']; surfs=sem['surfaces']
                    for shell,vals in zip(g['boundaries'],sem['values']):
                        for poly,v in zip(shell,vals):
                            if v is None: continue
                            s=surfs[v]
                            if s['type']=='WallSurface' and s.get('on_footprint_edge') is False:
                                ar=area3d(V[poly[0]]); n+=1; a+=ar
                                if ar<2: small+=1
                                if len(poly[0])==3: tri+=1
    print(f"{run}: inner walls {n}, area {a:.0f} m2, <2m2 {small}, triangles {tri}")
