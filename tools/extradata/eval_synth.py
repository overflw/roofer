import sys
import laspy
import numpy as np
lab = np.load("synth_labels.npy")
c = np.asarray(laspy.read(sys.argv[1]).classification)
for k, name in enumerate(["gable", "green roof", "tree"]):
    s = lab == k
    print(f"  {name:10s} n={s.sum():5d} removed={np.mean(c[s] != 6):.3f}")
