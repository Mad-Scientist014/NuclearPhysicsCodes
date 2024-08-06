import numpy
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import h5py
import re
from scipy.optimize import curve_fit
from periodictable import elements
from matplotlib.ticker import AutoMinorLocator
from sklearn.preprocessing import MinMaxScaler
import math

plots =[1,1]
path1 = "/home/thana1dr/WinNet/runs/vp_proc_1_flow/traj_0017.dat/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_1_flow/traj_0137.dat/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_1_flow/traj_0257.dat/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_1_flow/traj_0377.dat/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_1_flow/traj_0747.dat/"
filename = "WinNet_data.h5"

trajectorypath = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"

f1 = h5py.File(path1 + filename, 'r')
f2 = h5py.File(path2 + filename, 'r')
f3 = h5py.File(path3 + filename, 'r')
f4 = h5py.File(path4 + filename, 'r')
f5 = h5py.File(path5 + filename, 'r')

def findFlow(dataset, iso):
    p = iso[0]
    n = iso[1]
    flow = {}
    keys = list(dataset['flows'].keys())
    for key in keys:
        filta = np.where(dataset['flows'][key]['p_out'][:] == p, 1, 0)
        filtb = np.where(dataset['flows'][key]['n_out'][:] == p, 1, 0)
        isoFlow = np.where(dataset['flows'][key]['flow']*filta*filtb != 0)[0]
        if len(isoFlow) == 0:
          flow[dataset['flows'][key]['time'][0]] = 0
        else:
            flowTot = 0
            for each in isoFlow:
                flowTot += dataset['flows'][key]['flow'][each]
            flow[dataset['flows'][key]['time'][0]] = flowTot
    return flow

flow = findFlow(f1, (34,46))
print("Max Flow at t = %5.2f" % max(flow, key=flow.get))

i=0
start = 0
while i < len(flow.values()):
    if list(flow.values())[i] == 0 and list(flow.values())[i+1] != 0:
        start = list(flow.keys())[i]
        break
    else:
        i += 1

print("Flow starts at t = %5.2f" % start)
