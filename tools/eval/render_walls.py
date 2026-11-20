"""render_walls.py OUT.png BID RUNDIR... : plan view of roof faces + inner walls (red, width ~ height)."""
import sys, json, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
out,bid,runs=sys.argv[1],sys.argv[2],sys.argv[3:]
fig,ax=plt.subplots(1,len(runs),figsize=(6*len(runs),6),squeeze=False)
for a,run in zip(ax[0],runs):
    for f in glob.glob(run+'/out/*.city.jsonl'):
        L=open(f).read().splitlines(); tr=json.loads(L[0])['transform']; sc=np.array(tr['scale']); tl=np.array(tr['translate'])
        for l in L[1:]:
            ft=json.loads(l)
            if ft['id']!=bid and bid not in json.dumps(ft['CityObjects'].get(ft['id'],{}).get('attributes',{})): continue
            V=np.array(ft['vertices'])*sc+tl; nin=0
            for o in ft['CityObjects'].values():
                for g in o.get('geometry',[]):
                    if str(g.get('lod'))!='2.2' or g['type']!='Solid': continue
                    surfs=g['semantics']['surfaces']
                    for shell,vals in zip(g['boundaries'],g['semantics']['values']):
                        for poly,v in zip(shell,vals):
                            P=V[poly[0]]; t=surfs[v]['type'] if v is not None else ''
                            if t=='RoofSurface': a.fill(P[:,0],P[:,1],fc='#ddd',ec='k',lw=0.5)
                            elif t=='WallSurface' and surfs[v].get('on_footprint_edge') is False:
                                h=P[:,2].max()-P[:,2].min(); nin+=1
                                a.plot(P[:,0],P[:,1],'r-',lw=1+2*min(h,2))
            a.set_title(f"{run.split('/')[-1]}: {nin} inner walls"); a.set_aspect('equal')
fig.tight_layout(); fig.savefig(out,dpi=70)
