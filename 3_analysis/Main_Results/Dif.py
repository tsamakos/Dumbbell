import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os; from pathlib import Path

# For Royal Society of Chemistry journals
mpl.rcParams.update({
    "axes.labelsize": 14,
    "figure.autolayout": True,
    "font.size": 12,
    "grid.color": "0",
    "grid.linestyle": (0, (1, 5)),
    "legend.columnspacing": 1,
    "legend.fontsize": 12,
    "legend.handlelength": 1.25,
    "legend.labelspacing": 0.25,
    "savefig.dpi": 1800,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "text.usetex": True
})


N = np.array([5, 10, 20, 40])
D_polymer =   np.zeros(len(N))
D_dum     =   np.zeros(len(N))
D_dum_x   =   np.zeros(len(N))
D_dum_y   =   np.zeros(len(N))
D_dum_z   =   np.zeros(len(N))

time         =  {}
msd_polymer  =  {}
msd_dum      =  {} 
msd_dum_x    =  {} 
msd_dum_y    =  {} 
msd_dum_z    =  {}

for i in range(len(N)):
        
    # Build the paths to the data and trajectory files 
    path     = f"/home/atsamo/Documents/dumbbell/3_analysis/c0.85_N{N[i]:02d}/msd/dif.data"
    path_msd = f"/home/atsamo/Documents/dumbbell/3_analysis/c0.85_N{N[i]:02d}/msd/msd.data"
    # read data file in path
    print(path)
    data  = np.loadtxt(path)

    D_polymer[i] = data[0]
    D_dum[i]     = data[1]
    D_dum_x[i]   = data[2]
    D_dum_y[i]   = data[3]
    D_dum_z[i]   = data[4]

    data2 = np.loadtxt(path_msd)
    time[i]        =  data2[:, 0] 
    msd_polymer[i] =  data2[:, 1]
    msd_dum[i]     =  data2[:, 2] 
    msd_dum_x[i]   =  data2[:, 3] 
    msd_dum_y[i]   =  data2[:, 4] 
    msd_dum_z[i]   =  data2[:, 5]



slope,intercept = np.polyfit(np.log(N), np.log(D_polymer), 1)
D_fit = np.exp(intercept)*N**slope
print(slope, intercept)

plt.figure()
plt.plot(time[0], msd_dum_x[3], 'k-', label=f"x-dimension")
plt.plot(time[0], msd_dum_y[3], 'g-', label=f"y-dimension")
plt.plot(time[0], msd_dum_z[3], 'b--', label=f"z-dimension")
plt.xlim([5e0, 1e5])
plt.xlabel("Time"); plt.ylabel("MSD")
# plt.ylim([1e-1, 1e5])
plt.legend()
plt.xscale("log"); plt.yscale("log")
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/msd_dum.png", dpi=300)


# ----------
plt.figure()
plt.plot(N, D_dum,     'bo-', label="Dumbbell")
plt.plot(N, D_polymer, 'ko-', label="Polymer")
plt.plot(N, D_fit, 'g--', label=f"slope N$^{{{slope:.2f}}}$")
plt.xlabel("N"); plt.ylabel("D")
plt.legend()
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_poly_dum.png", dpi=300)
plt.xscale("log"); plt.yscale("log");
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_poly_dum_log.png", dpi=300)


# ----------
plt.figure()
plt.plot(N, D_dum_x, 'bo-', label="x-direction")
plt.plot(N, D_dum_y, 'ko-', label="y-direction")
plt.plot(N, D_dum_z, 'go-', label="z-direction")
plt.xlabel("N"); plt.ylabel("D")
plt.legend()
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_dum.png", dpi=300)
plt.xscale("log"); plt.yscale("log");
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_dum_log.png", dpi=300)


# ----------
plt.figure()
plt.plot(N, D_dum_x, 'bo-', label="x-direction")
plt.plot(N, (D_dum_y+D_dum_z)/2, 'yo-', label="y,z-direction")
plt.xlabel("N")
plt.ylabel("D")
plt.legend()
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_dum_xyz.png", dpi=300)
plt.xscale("log"); plt.yscale("log");
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_dum_xyz_log.png", dpi=300)

# ----------
plt.figure()
plt.plot(N, (D_dum_y+D_dum_z)/2/D_dum_x, 'ko-')
plt.xlabel("N"); plt.ylabel("D$_{\\perp}$/D$_{\\parallel}$")
plt.savefig("/home/atsamo/Documents/dumbbell/3_analysis/Main_Results/figures/D_dum_ratio.png", dpi=300)
