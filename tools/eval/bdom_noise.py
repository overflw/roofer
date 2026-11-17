"""bdom_noise.py in.laz out.laz : add image-matching-like noise to z: a smooth random field
(Gaussian blobs ~2-4 m, std 0.22 m) plus 0.05 m white noise. Keeps every 2nd point (density ~bDOM 20 cm raster thinned)."""
import sys, numpy as np, laspy, scipy.ndimage as nd
rng=np.random.default_rng(7)
las=laspy.read(sys.argv[1]); x,y=np.asarray(las.x),np.asarray(las.y)
c=0.25; ix=((x-x.min())/c).astype(int); iy=((y-y.min())/c).astype(int)
f=nd.gaussian_filter(rng.normal(size=(iy.max()+1,ix.max()+1)), sigma=1.2/c)
f*=0.22/f.std()
z=np.asarray(las.z)+f[iy,ix]+rng.normal(0,0.05,len(x))
las.z=z
las.write(sys.argv[2]); print("field std",f.std())
