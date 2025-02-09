import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
import matplotlib.pyplot as plt
import os; from pathlib import Path


# Get the name of the current directory, e.g. "c0.85_N05"
script_path = Path(__file__).resolve(); folder = script_path.parent.parent.name 

# Build the paths to the data and trajectory files 
data_path = f"/home/atsamo/Documents/dumbbell/2_runs/{folder}/run02/new.data"
traj_path = f"/home/atsamo/Documents/dumbbell/2_runs/{folder}/run02/dump_unwrapped.dcd"

# ------------------------------------------------------------------
# Load the LAMMPS topology and trajectory

u = mda.Universe(data_path, traj_path, format="LAMMPS")
print("Number of atoms:", len(u.atoms))
print("Number of frames in trajectory:", u.trajectory.n_frames)


# ------------------------------------------------------------------
# Dumbbell 
polymer = u.select_atoms("type 1"); dumbbell = u.select_atoms("type 2"); solvent = u.select_atoms("type 3") 
polymer_positions = [];   polymer_com_positions = []
dum_positions     = [];   dum_com_positions = []
solvent_positions = [];     
 
for ts in u.trajectory:

    # polymer
    temp1 = polymer.positions;        polymer_positions.append(temp1.copy())
    temp2 = polymer.center_of_mass(); polymer_com_positions.append(temp2.copy())

    # dumbbell
    temp1 = dumbbell.positions;        dum_positions.append(temp1.copy())     
    temp2 = dumbbell.center_of_mass(); dum_com_positions.append(temp2.copy())

    # solvent
    temp1 = solvent.positions;          solvent_positions.append(temp1.copy())

# make numpy arrays
polymer_positions  =  np.array(polymer_positions);  polymer_com_positions = np.array(polymer_com_positions)
dum_positions      =  np.array(dum_positions);      dum_com_positions     = np.array(dum_com_positions)
solvent_positions  =  np.array(solvent_positions)


# time array 
dump_freq = 1000; dt = 0.01;
n_frames  = len(dum_positions); time = np.arange(n_frames) * dt*dump_freq
print(n_frames)


# ================
# functions to use
# ================

def compute_msd_1d(positions): 
    n_steps = positions.shape[0]
    # Initialize arrays for MSD in each direction
    msd_x = np.zeros(n_steps); msd_y = np.zeros(n_steps); msd_z = np.zeros(n_steps)
    # For each possible time lag (tau)
    for tau in range(n_steps):
        # Displacements in x, y, z from t to t+tau
        dx = positions[tau:, 0] - positions[:n_steps - tau, 0]
        dy = positions[tau:, 1] - positions[:n_steps - tau, 1]
        dz = positions[tau:, 2] - positions[:n_steps - tau, 2]
        # Mean of the squared displacements in each direction
        msd_x[tau] = np.mean(dx**2); msd_y[tau] = np.mean(dy**2); msd_z[tau] = np.mean(dz**2)
    return msd_x, msd_y, msd_z


import MDAnalysis.analysis.msd as msd

# dumbbell

# msd_x, msd_y, msd_z = compute_msd_1d(dum_com_positions)
MSD_dum   =  msd.EinsteinMSD(u, select='type 2', msd_type='xyz', fft=True); MSD_dum.run()
msd_dum   =  MSD_dum.results.timeseries

MSD_dum_x =  msd.EinsteinMSD(u, select='type 2', msd_type='x', fft=True); MSD_dum_x.run()
msd_dum_x =  MSD_dum_x.results.timeseries

MSD_dum_y =  msd.EinsteinMSD(u, select='type 2', msd_type='y', fft=True); MSD_dum_y.run()
msd_dum_y =  MSD_dum_y.results.timeseries

MSD_dum_z =  msd.EinsteinMSD(u, select='type 2', msd_type='z', fft=True); MSD_dum_z.run()
msd_dum_z =  MSD_dum_z.results.timeseries

# polymer
MSD_polymer =  msd.EinsteinMSD(u, select='type 1', msd_type='xyz', fft=True); MSD_polymer.run()
msd_polymer =  MSD_polymer.results.timeseries


# write in file
f = open(f"/home/atsamo/Documents/dumbbell/3_analysis/{folder}/msd/msd.data", "w")
f.write("# time msd_polymer msd_dumbbell msd_dumbbell_x msd_dumbbell_y msd_dumbbell_z\n")
for i in range(len(time)):
    f.write("%12.4f %12.8f %12.8f %12.8f %12.8f %12.8f\n" % (time[i], msd_polymer[i], msd_dum[i], msd_dum_x[i], msd_dum_y[i], msd_dum_z[i]))
