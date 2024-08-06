import random
from matplotlib import pyplot as plt
import matplotlib.animation as animation
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
from datetime import datetime
from tqdm import tqdm
import time
from sys import getsizeof
import gc
import ffmpeg

path1 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0017.dat/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0137.dat/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0257.dat/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0377.dat/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0747.dat/"
filename = "WinNet_data.h5"

trajectorypath = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"

f1 = h5py.File(path1 + filename, 'r')
f2 = h5py.File(path2 + filename, 'r')
f3 = h5py.File(path3 + filename, 'r')
f4 = h5py.File(path4 + filename, 'r')
f5 = h5py.File(path5 + filename, 'r')

def func(x,a,b):
    return a/(np.power(x, b))

def initGraph(dim):
    if dim[0] != 1 or dim[1] != 1:
        fig, axis = plt.subplots(dim[0],dim[1], layout="constrained", figsize=(20,10))
    else:
        fig, axis = plt.subplots(dim[0], dim[1], layout="constrained")
    return fig, axis

def splitNLast(item):
    return item[0][-1], item[1][-1]

def procNucleiGen(dataset):
    print("Processing Nuclei: ")
    for t in tqdm(range(len(dataset['tracked_nuclei']['time'][()]))):
        maxSize  = dataset['tracked_nuclei']['A'].shape[0]
        data = np.zeros((maxSize, 5))
        data[:, 0] = dataset['tracked_nuclei']['A'][:]
        data[:, 1] = dataset['tracked_nuclei']['N'][:]
        data[:, 2] = dataset['tracked_nuclei']['Z'][:]
        data[:, 3] = dataset['tracked_nuclei']['Y'][t, :]
        data[0, 4] = dataset['tracked_nuclei']['time'][t]
        yield data

def procNuclei(dataset,t):
    maxSize  = dataset['tracked_nuclei']['A'].shape[0]
    data = np.zeros((maxSize, 5))
    data[:, 0] = dataset['tracked_nuclei']['A'][:]
    data[:, 1] = dataset['tracked_nuclei']['N'][:]
    data[:, 2] = dataset['tracked_nuclei']['Z'][:]
    data[:, 3] = dataset['tracked_nuclei']['Y'][t, :]
    data[0, 4] = dataset['tracked_nuclei']['time'][t]
    return data

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def combineLikeRows(data, sort, add):
    i = 1
    data = data[:,data[sort, :].argsort()]
    for col in data[sort, :]:
        if i<len(data[sort, :]) and col == data[sort, i]:
            data[add,i] += data[add, i-1]
            data = np.delete(data, i-1, axis=1)
        else:
            i += 1
    return data

def isotopeFlow(dataset, isotope, time = -1, temp = -1, dens = -1):
    if time == temp == dens == -1:
        print("Invalid Condition")
        return
    else:
        i = 0
        p = isotope[0]
        n = isotope[1]

        timeRange = np.array([dataset['flows'][key]['time'][-1] for key in list(dataset['flows'].keys())])
        tempRange = np.array([dataset['flows'][key]['temp'][-1] for key in list(dataset['flows'].keys())])
        densRange = np.array([dataset['flows'][key]['dens'][-1] for key in list(dataset['flows'].keys())])
        if time != -1:
            key = list(dataset['flows'].keys())[np.where(timeRange == find_nearest(timeRange, time))[0][0]]
        elif temp != -1:
            key = list(dataset['flows'].keys())[np.where(tempRange == find_nearest(tempRange, temp))[0][0]]
        elif dens != -1:
            key = list(dataset['flows'].keys())[np.where(densRange == find_nearest(densRange, dens))[0][0]]

        p_in = np.where(dataset['flows'][key]['p_in'][()] == p, 1, 0)
        n_in = np.where(dataset['flows'][key]['n_in'][()] == n, 1, 0)
        p_out = np.where(dataset['flows'][key]['p_out'][()] == p, 1, 0)
        n_out = np.where(dataset['flows'][key]['n_out'][()] == n, 1, 0)
        mask = p_in*n_in+p_out*n_out
        nums = np.where(mask == 1)[0]
        reactions = np.zeros((len(nums), 10))
        reactions[:, 0] = [dataset['flows'][key]['p_in'][num] for num in nums]
        reactions[:, 1] = [dataset['flows'][key]['n_in'][num] for num in nums]
        reactions[:, 2] = [dataset['flows'][key]['y_in'][num] for num in nums]
        reactions[:, 3] = [dataset['flows'][key]['p_out'][num] for num in nums]
        reactions[:, 4] = [dataset['flows'][key]['n_out'][num] for num in nums]
        reactions[:, 5] = [dataset['flows'][key]['y_out'][num] for num in nums]
        reactions[:, 6] = [dataset['flows'][key]['flow'][num] for num in nums]
        reactions[:, 7] = [dataset['flows'][key]['temp'][-1] for num in nums]
        reactions[:, 8] = [dataset['flows'][key]['time'][-1] for num in nums]
        reactions[:, 9] = [dataset['flows'][key]['dens'][-1] for num in nums]
        return reactions

def graphIsoFlow(dataset, iso, fig, axis, flowScalefactor, minFlowshow, cbar=[], time=-1, temp=-1, dens=-1):
    temp = float(temp)
    xrange = list(range(iso[1] - 6, iso[1] + 6))
    yrange = list(range(iso[0] - 6, iso[0] + 6))
    arrow = isotopeFlow(dataset, iso, time=time, temp=temp, dens=dens)

    adjData = arrow[:, 6]
    adjData = np.where(adjData >= minFlowshow, adjData, 0)
    arrow[:, 6] = adjData[:]

    if yrange[-1]- yrange[1] <= 30 and xrange[-1]- xrange[1] <= 30:
        for i in yrange[1:-1]:
            for j in xrange[1:-1]:
                axis.text(j, i, str(i+j) + "-" + elements[i].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                ha='center', va='center')

    cmap = plt.cm.jet
    if cbar == []:
        cNorm = colors.Normalize(vmin=np.min(arrow[:,6]), vmax=np.max(arrow[:,6]))
    else:
        cNorm = colors.Normalize(vmin=cbar[0], vmax=cbar[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    for arr in arrow:
        if arr[6] != 0:

            colorVal = scalarMap.to_rgba(arr[6])
            axis.arrow(arr[1], arr[0], arr[4]-arr[1], arr[3]-arr[0], color = colorVal, head_width=arr[6]*.4*flowScalefactor, width= arr[6]*.2*flowScalefactor, length_includes_head=True)

    minor_locator1 = AutoMinorLocator(2)
    minor_locator2 = AutoMinorLocator(2)
    axis.xaxis.set_minor_locator(minor_locator1)
    axis.yaxis.set_minor_locator(minor_locator2)
    axis.set_xlim(xrange[0], xrange[-1])
    axis.set_ylim(yrange[0], yrange[-1])
    axis.grid(which='minor', color='black', linewidth=0.5)
    fig.colorbar(scalarMap, ax=axis)

# WORKS JUST TESING GENERATORS EFFICIENCY
plots =[1,1]
fig, axis = initGraph(plots)
minFlowshow =  1e-6
scaleFactor = 3e3

graphIsoFlow(f1, (32, 37), fig, axis, scaleFactor, minFlowshow, time = .3)
plt.show()
