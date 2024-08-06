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
    array = np.sort(array)
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


def dynamicAbundance(dataset, num, xrange, yrange, fig, axis, cbar=[]):
    return dataset

plots =[1,1]
i = 0
source1 = f1
out1 = np.zeros((len(source1['tracked_nuclei']['time'][()]), 2, 2000))
for data in procNucleiGen(source1):
    # updating the data
    abundance = data[:, 3] # Abundance
    mass = data[:, 2] # Z Number
    t = data[0, 4] # Time Stamp
    combined = np.vstack([mass,abundance])
    combined = combineLikeRows(combined, 0, 1)
    out1[i, :, :combined.shape[1]] = combined
    i += 1

i = 0
source2 = f2
out2 = np.zeros((len(source2['tracked_nuclei']['time'][()]), 2, 2000))
for data in procNucleiGen(source2):
    # updating the data
    abundance = data[:, 3] # Abundance
    mass = data[:, 2] # Z Number
    t = data[0, 4] # Time Stamp
    combined = np.vstack([mass,abundance])
    combined = combineLikeRows(combined, 0, 1)
    out2[i, :, :combined.shape[1]] = combined
    i += 1

i = 0
source3 = f3
out3 = np.zeros((len(source3['tracked_nuclei']['time'][()]), 2, 2000))
for data in procNucleiGen(source3):
    # updating the data
    abundance = data[:, 3] # Abundance
    mass = data[:, 2] # Z Number
    t = data[0, 4] # Time Stamp
    combined = np.vstack([mass,abundance])
    combined = combineLikeRows(combined, 0, 1)
    out3[i, :, :combined.shape[1]] = combined
    i += 1

i = 0
source4 = f4
out4 = np.zeros((len(source4['tracked_nuclei']['time'][()]), 2, 2000))
for data in procNucleiGen(source4):
    # updating the data
    abundance = data[:, 3] # Abundance
    mass = data[:, 2] # Z Number
    t = data[0, 4] # Time Stamp
    combined = np.vstack([mass,abundance])
    combined = combineLikeRows(combined, 0, 1)
    out4[i, :, :combined.shape[1]] = combined
    i += 1

i = 0
source5 = f5
out5 = np.zeros((len(source5['tracked_nuclei']['time'][()]), 2, 2000))
for data in procNucleiGen(source5):
    # updating the data
    abundance = data[:, 3] # Abundance
    mass = data[:, 2] # Z Number
    t = data[0, 4] # Time Stamp
    combined = np.vstack([mass,abundance])
    combined = combineLikeRows(combined, 0, 1)
    out5[i, :, :combined.shape[1]] = combined
    i += 1

duration = np.max([len(source['tracked_nuclei']['time'][:]) for source in [f1,f2,f3,f4,f5]])

fig, axis = initGraph(plots)
text = axis.text(85, 1e-2, "Time = 0", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1}, ha='center', va='center')

line1 = axis.plot([0], [0], color="red")[0]
line2 = axis.plot([0], [0], color ="black")[0]
line3 = axis.plot([0], [0], color="green")[0]
line4 = axis.plot([0], [0], color ="blue")[0]
line5 = axis.plot([0], [0], color="orange")[0]


axis.set(xlim=[0, 95], ylim=[1e-30, 1], xlabel='Z', ylabel='Abundance')
axis.set_yscale('log')
axis.grid()

minT = np.min([round(source['tracked_nuclei']['time'][0], 4) for source in [f1,f2,f3,f4,f5]])
maxT = np.max([round(source['tracked_nuclei']['time'][-1], 4) for source in [f1,f2,f3,f4,f5]])


def update_t(frame):
    time = np.linspace(minT, maxT, duration-1)[frame]
    text.set_text("Time = " + str(round(time, 4)))
    if frame < out1.shape[0]:
        t1 = np.where(f1['tracked_nuclei']['time'][()] == find_nearest(f1['tracked_nuclei']['time'][()], time))[0][0]
        x1 = out1[t1, 0, :]
        y1 = out1[t1, 1, :]
        line1.set_xdata(x1)
        line1.set_ydata(y1)

    if frame < out2.shape[0]:
        t2 = np.where(f2['tracked_nuclei']['time'][()] == find_nearest(f2['tracked_nuclei']['time'][()], time))[0][0]
        x2 = out2[t2, 0, :]
        y2 = out2[t2, 1, :]
        line2.set_xdata(x2)
        line2.set_ydata(y2)

    if frame < out3.shape[0]:
        t3 = np.where(f3['tracked_nuclei']['time'][()] == find_nearest(f3['tracked_nuclei']['time'][()], time))[0][0]
        x3 = out3[t3, 0, :]
        y3 = out3[t3, 1, :]
        line3.set_xdata(x3)
        line3.set_ydata(y3)

    if frame < out4.shape[0]:
        t4 = np.where(f4['tracked_nuclei']['time'][()] == find_nearest(f4['tracked_nuclei']['time'][()], time))[0][0]
        x4 = out4[t4, 0, :]
        y4 = out4[t4, 1, :]
        line4.set_xdata(x4)
        line4.set_ydata(y4)

    if frame < out5.shape[0]:
        t5 = np.where(f5['tracked_nuclei']['time'][()] == find_nearest(f5['tracked_nuclei']['time'][()], time))[0][0]
        x5 = out5[t5, 0, :]
        y5 = out5[t5, 1, :]
        line5.set_xdata(x5)
        line5.set_ydata(y5)


    return (text, line1, line2, line3, line4, line5,)

def update_u(frame):
    time = np.linspace(minT, maxT, duration-1)[frame]
    text.set_text("Time1 = " + str(round(time, 4)))
    if frame < out1.shape[0]:
        t1 = math.floor((frame / duration - 1) * out1.shape[0])
        x1 = out1[t1, 0, :]
        y1 = out1[t1, 1, :]
        line1.set_xdata(x1)
        line1.set_ydata(y1)

    if frame < out2.shape[0]:
        t2 = math.floor((frame / duration - 1) * out2.shape[0])
        x2 = out2[t2, 0, :]
        y2 = out2[t2, 1, :]
        line2.set_xdata(x2)
        line2.set_ydata(y2)

    if frame < out3.shape[0]:
        t3 = math.floor((frame / duration - 1) * out3.shape[0])
        x3 = out3[t3, 0, :]
        y3 = out3[t3, 1, :]
        line3.set_xdata(x3)
        line3.set_ydata(y3)

    if frame < out4.shape[0]:
        t4 = math.floor((frame / duration - 1) * out4.shape[0])
        x4 = out4[t4, 0, :]
        y4 = out4[t4, 1, :]
        line4.set_xdata(x4)
        line4.set_ydata(y4)

    if frame < out5.shape[0]:
        t5 = math.floor((frame / duration - 1) * out5.shape[0])
        x5 = out5[t5, 0, :]
        y5 = out5[t5, 1, :]
        line5.set_xdata(x5)
        line5.set_ydata(y5)

    return (text, line1, line2, line3, line4, line5,)

ani = animation.FuncAnimation(fig=fig, func=update_t, frames=duration-1, interval=10)
plt.yscale('log')
plt.show()
if input("Save Animation: (Yes/no)").lower()[0] == "y":
    print("Saving Animation (This might take some time) .  .  .")
    ani.save(filename="TimeAbundancePlot" +str(math.floor(time.time()/1e4))+".mp4", writer="ffmpeg")




"""
    plots =[1,1]
xrange = range(0,40)
yrange = range(0,40)
minFlowshow =  1e-6
scaleFactor = 1e3


fig, axis = initGraph(plots)

plt.ion()
# turning interactive mode on

# the update loop
for data in procNucleiGen(f1):
    # updating the data
    abundance = data[:, 3] # Abundance
    mass = data[:, 0] # Mass Number
    t = data[0, 4] # Time Stamp
    combined = np.vstack([mass,abundance])
    out = combineLikeRows(combined, 0, 1)
    x = out[0, :]
    y = out[1, :]

    # removing the older graph
    plt.cla()

    # plotting newer graph
    graph = axis.plot(x, y, color='r')[0]
    axis.set_xlim(x[0], x[-1])
    plt.title("Time = "+ str(t))
    plt.xlabel('Nuclei Mass')
    plt.ylabel('Nuclei Abundance')
    plt.xscale('linear')
    plt.yscale('log')
    axis.set_ylim(1e-30, 1)
    minor_locator1 = AutoMinorLocator(1)
    minor_locator2 = AutoMinorLocator(1)
    axis.xaxis.set_minor_locator(minor_locator1)
    axis.yaxis.set_minor_locator(minor_locator2)
    axis.grid()
    axis.grid(which='minor', color='black', linewidth=0.5)

    # calling pause function for 0.25 seconds
    plt.pause(.01)
    """


