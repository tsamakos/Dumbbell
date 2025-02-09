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

# Get the name of the current directory, e.g. "c0.85_N05"
script_path = Path(__file__).resolve(); folder = script_path.parent.parent.name 

# Build the paths to the data and trajectory files 
path = f"/home/atsamo/Documents/dumbbell/3_analysis/{folder}/msd/msd.data"

# read data file in path
data = np.loadtxt(path)
# extract time and msd
time = data[:, 0]; msd_polymer = data[:, 1]; msd_dum = data[:, 2]; msd_dum_x = data[:, 3]; msd_dum_y = data[:, 4]; msd_dum_z = data[:, 5]
n_frames = len(time)

# line in log-log scale
time_fit = np.linspace(time[int(n_frames/1000)], time[int(n_frames/100)], 10)
msd_fit  = 0.5e0*time_fit                       
                       

# plot the slope data
plt.figure()
plt.plot(time, msd_dum,     'b-', label="Dumbbell")
plt.plot(time, msd_polymer, 'k-', label="Polymer")
plt.plot(time_fit, msd_fit, '--k', label="slope 1")
plt.legend()
plt.xlabel("Time"); plt.ylabel("MSD")
plt.xscale("log");  plt.yscale("log");
plt.savefig(f"/home/atsamo/Documents/dumbbell/3_analysis/{folder}/msd/msd_poly_dum.png", dpi=300)
plt.close()


def calc_msd_slope(msd, time):  
    slope = np.gradient(np.log(msd[1:]), np.log(time[1:]))
    return slope

slope_msd_polymer  =  calc_msd_slope(msd_polymer, time)
slope_msd_dum      =  calc_msd_slope(msd_dum, time)



# ======================
# Difussion coefficients
# ======================
def calc_diffusion_coeff(t0, tf, time, msd, dimensions):

    idx_t0 = np.searchsorted(time, t0, side='left');     idx_tf = np.searchsorted(time, tf, side='left')
    slope, intercept = np.polyfit(np.log10(time[idx_t0:idx_tf]), np.log10(msd[idx_t0:idx_tf]),1)
    D = 10**intercept/(2*dimensions)
 
    return slope, D


t0 = 1e+3; tf = 3e+3
slope_poly, D_poly = calc_diffusion_coeff(t0, tf, time, msd_polymer, dimensions=3)

t0 = 1e+2; tf = 6.0e+3
slope_dum,  D_dum  = calc_diffusion_coeff(t0, tf, time, msd_dum, dimensions=3)

t0 = 0.5e+2; tf = 1e+3
slope_dum_x,  D_dum_x  = calc_diffusion_coeff(t0, tf, time, msd_dum_x, dimensions=1)

t0 = 0.7e+2; tf = 7e+2
slope_dum_y,  D_dum_y  = calc_diffusion_coeff(t0, tf, time, msd_dum_y, dimensions=1)

t0 = 0.5e+2; tf = 1.8e+2
slope_dum_z,  D_dum_z  = calc_diffusion_coeff(t0, tf, time, msd_dum_z, dimensions=1)


print('--------')
print(slope_poly,slope_dum,slope_dum_x,slope_dum_y,slope_dum_z)
print(D_poly,D_dum,D_dum_x,D_dum_y,D_dum_z)
print('--------')


f = open(f"/home/atsamo/Documents/dumbbell/3_analysis/{folder}/msd/dif.data", "w")
f.write("# D_poly D_dum D_dum_x D_dum_y D_dum_z\n")
f.write(f"{D_poly} {D_dum} {D_dum_x} {D_dum_y} {D_dum_z}\n")
f.close()


# plot the slope of the msd 
plt.figure()
plt.plot(time[1:], calc_msd_slope(msd_dum_y, time),'.')
plt.xlabel("Time ")
plt.ylabel("Slope of MSD")
plt.xscale("log")
plt.ylim([0,3])
plt.show()



# plot the slope of the msd 
plt.figure()
plt.plot(time[1:], calc_msd_slope(msd_dum_z, time),'.')
plt.xlabel("Time ")
plt.ylabel("Slope of MSD")
plt.xscale("log")
plt.ylim([0,3])
plt.show()

