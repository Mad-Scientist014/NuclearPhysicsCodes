import numpy
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import h5py
import csv
import periodictable
from scipy.optimize import curve_fit

plots =[3,3]
path1 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0017.dat/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0137.dat/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0257.dat/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0377.dat/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0747.dat/"
filename = "WinNet_data.h5"

walletcardPath = "/home/thana1dr/PycharmProjects/WinNetDataGrapher/"

f1 = h5py.File(path1 + filename, 'r')
f2 = h5py.File(path2 + filename, 'r')
f3 = h5py.File(path3 + filename, 'r')
f4 = h5py.File(path4 + filename, 'r')
f5 = h5py.File(path5 + filename, 'r')

def readDecays(path, file):
    decays = {}
    with open(path+file, mode='r') as file:
        csvFile = csv.DictReader(file)
        for lines in csvFile:
            decays[lines['Element']+lines['Atomic Mass (A)']] = {'t1/2': lines['Half-Life'],
                                                                'tUnit': lines['Half-Life (Unit)'],
                                                                'tErr': lines['Half-Life (Error)']
                                                                }
    return decays

data = readDecays(walletcardPath, 'walletcards.csv')
print(data)

np.min([np.sum(f1['tracked_nuclei']['Y'][i, :]*f1['tracked_nuclei']['A'][()]) for i in range(1806)])
plt.plot(f1['tracked_nuclei']['time'][()], [np.sum(f1['tracked_nuclei']['Y'][i, :], axis = 0 ) for i in range(1807)], label="Total Abundance")
plt.plot(f1['tracked_nuclei']['time'][()], [np.sum(f1['tracked_nuclei']['Y'][i, :]*f1['tracked_nuclei']['A'][()], axis = 0) for i in range(1807)], label="Total Mass")

plt.legend()
plt.show()