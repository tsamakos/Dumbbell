import numpy as np
import matplotlib.pyplot as plt
import os; from pathlib import Path




N = np.array([5, 10, 20, 40])
D_polymer =   np.zeros(len(N))
D_dum     =   np.zeros(len(N))
D_dum_x   =   np.zeros(len(N))
D_dum_y   =   np.zeros(len(N))
D_dum_z   =   np.zeros(len(N))

for i in range(len(N)):
        
    # Build the paths to the data and trajectory files 
    path = f"/home/atsamo/Documents/dumbbell/3_analysis/c0.85_N{N[i]:02d}/msd/dif.data"
    # read data file in path
    print(path)
    data = np.loadtxt(path)
    for k in range(len(data)):
        D_polymer[i] = data[k]
        D_dum[i]     = data[k]
        D_dum_x[i]   = data[k]
        D_dum_y[i]   = data[k]
        D_dum_z[i]   = data[k]




plt.figure()
plt.plot(N, D_polymer, 'ko-', label="D$_{\\parallel}$")
plt.plot(N, D_dum,     'bo-', label="D$_{\\perp}$")
plt.xlabel("N")
plt.ylabel("D")
plt.legend()
plt.savefig("D.png")

# # plot ratio
# plt.figure()
# plt.plot(N, Dver/Dpar, 'o-')
# plt.xlabel("N")
# plt.ylabel("D$_{\\perp}$/D$_{\\parallel}$")
# plt.savefig("D_ratio.png")

