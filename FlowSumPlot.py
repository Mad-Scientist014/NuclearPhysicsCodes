import periodictable
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import mendeleev
import numpy as np
import h5py
import re
from scipy.optimize import curve_fit
from periodictable import elements
from matplotlib.ticker import AutoMinorLocator
from sklearn.preprocessing import MinMaxScaler
from numpy import log as ln
import math
from datetime import datetime
from tqdm import tqdm
from sys import getsizeof
import time
import csv
import gc

path1 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0017.dat/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0137.dat/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0257.dat/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0377.dat/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0747.dat/"
filename = "WinNet_data.h5"

walletcardPath = "/home/thana1dr/PycharmProjects/WinNetDataGrapher/"

SummedFlowsPath = "/home/thana1dr/Research_Results&Data/Week8/Automatic/ReactionLists/"

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

def procFlow(dataset, key):
    keys = list(dataset['flows'].keys())
    if isinstance(key, int):
        key = keys[key]
    maxSize = np.max([len(dataset['flows'][key]['p_in'][()]) for key in keys])
    data = np.zeros((maxSize, 5))
    columnLen = len(dataset['flows'][key]['p_in'][:])
    data[:columnLen, 0] = dataset['flows'][key]['p_in'][()]
    data[:columnLen, 1] = dataset['flows'][key]['n_in'][()]
    data[:columnLen, 2] = dataset['flows'][key]['p_out'][()]
    data[:columnLen, 3] = dataset['flows'][key]['n_out'][()]
    data[:columnLen, 4] = dataset['flows'][key]['flow'][()]
    return data

def procFlowGen(dataset, key):
    keys = list(dataset['flows'].keys())
    if isinstance(key, int):
        key = keys[key]
    maxSize = len(dataset['flows'][key]['p_in'][()])
    data = np.zeros((6))
    pin = dataset['flows'][key]['p_in'][:]
    nin = dataset['flows'][key]['n_in'][:]
    pout = dataset['flows'][key]['p_out'][:]
    nout = dataset['flows'][key]['n_out'][:]
    fout = dataset['flows'][key]['flow'][:]
    time = dataset['flows'][key]['time'][-1]
    for i in tqdm(range(maxSize)):
        data[0] = pin[i]
        data[1] = nin[i]
        data[2] = pout[i]
        data[3] = nout[i]
        data[4] = fout[i]
        data[5] = time
        yield data

def procFlowGenMasked(dataset, key, mask):
    keys = list(dataset['flows'].keys())
    if isinstance(key, int):
        key = keys[key]
    data = np.zeros((6))
    pin = dataset['flows'][key]['p_in'][mask]
    nin = dataset['flows'][key]['n_in'][mask]
    pout = dataset['flows'][key]['p_out'][mask]
    nout = dataset['flows'][key]['n_out'][mask]
    fout = dataset['flows'][key]['flow'][mask]
    time = dataset['flows'][key]['time'][-1]
    for i in tqdm(range(len(pin))):
        data[0] = pin[i]
        data[1] = nin[i]
        data[2] = pout[i]
        data[3] = nout[i]
        data[4] = fout[i]
        data[5] = time
        yield data
def procFlowRange(dataset, keys):
    maxSize = np.max([dataset['flows'][key]['p_in'].shape[0] for key in keys])
    data = np.zeros((len(keys), maxSize, 5))
    for i in tqdm(range(len(keys))):
        each = keys[i]
        columnLen = len(dataset['flows'][each]['p_in'][:])
        data[i,:columnLen, 0] = dataset['flows'][each]['p_in'][()]
        data[i,:columnLen, 1] = dataset['flows'][each]['n_in'][()]
        data[i,:columnLen, 2] = dataset['flows'][each]['p_out'][()]
        data[i,:columnLen, 3] = dataset['flows'][each]['n_out'][()]
        data[i,:columnLen, 4] = dataset['flows'][each]['flow'][()]
        data[i, :columnLen, 5] = dataset['flows'][each]['y_in'][()]
        data[i, :columnLen, 6] = dataset['flows'][each]['y_out'][()]
    return data

def procFlowRangeGen(dataset, keys):
    maxSize = np.max([dataset['flows'][key]['p_in'].shape[0] for key in keys])
    data = np.zeros((maxSize, 7))
    for i in tqdm(range(len(keys))):
        each = keys[i]
        columnLen = len(dataset['flows'][each]['p_in'][:])
        data[:columnLen, 0] = dataset['flows'][each]['p_in'][()]
        data[:columnLen, 1] = dataset['flows'][each]['n_in'][()]
        data[:columnLen, 2] = dataset['flows'][each]['p_out'][()]
        data[:columnLen, 3] = dataset['flows'][each]['n_out'][()]
        data[:columnLen, 4] = dataset['flows'][each]['flow'][()]
        data[:columnLen, 5] = dataset['flows'][each]['y_in'][()]
        data[:columnLen, 6] = dataset['flows'][each]['y_out'][()]
        yield data

def find_nearest(array,value):
    array = np.sort(array)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def find_nearest_index(array,value):
    array = np.sort(array)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def chunks(xs, n):
    n = max(1, n)
    return (xs[i:i+n] for i in range(0, len(xs), n))

def findCritReactions(dataset, number, split):
    keys = list(dataset['flows'].keys())
    critReacs = {}
    print("Finding Critical Reactions: ")
    chunkedKeys = np.array_split(keys, split)
    for i in range(len(chunkedKeys)):
        data = procFlowRange(f1, chunkedKeys[i])
        for j in range(len(data)):
            ind = np.argpartition(data[j, :, -1], -number)[-number:]
            critReacs[(i+1)*(j+1)] = [data[j, each, :] for each in ind]
        del data
        gc.collect()
    return critReacs

def findCritReactionsGen(dataset, number, split):
    keys = list(dataset['flows'].keys())
    print("Finding Critical Reactions: ")
    i = 1
    for data in procFlowRangeGen(dataset, keys):
        ind = np.argsort(data[:,4]).tolist()
        sort = data[ind.reverse(),:][0][:number]
        i+=1
        yield sort

def plotCritReactions(dataset, fig, axis, flowScalefactor, cbar=[]):
    critReactions = findCritReactionsGen(dataset, 200, 3)
    keys = list(dataset['flows'].keys())
    xr = []
    yr = []
    cmap = plt.cm.jet
    if cbar == []:
        cNorm = colors.Normalize(vmin=dataset['flows'][keys[0]]['time'][-1], vmax=dataset['flows'][keys[-1]]['time'][-1])
    else:
        cNorm = colors.Normalize(vmin=cbar[0], vmax=cbar[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    arrow = np.zeros((7))
    arrow = np.vstack([arrow, np.zeros((7))])
    i=0
    j=0

    for item in critReactions:
        for each in item:
            if each[4] >= 0 and each.tolist()[:4] not in arrow[:,:4].tolist():
                j+=1
                xr.append(np.min([each[1] for each in item]))
                xr.append(np.min([each[3] for each in item]))
                xr.append(np.max([each[1] for each in item]))
                xr.append(np.max([each[3] for each in item]))

                yr.append(np.min([each[0] for each in item]))
                yr.append(np.min([each[2] for each in item]))
                yr.append(np.max([each[0] for each in item]))
                yr.append(np.max([each[2] for each in item]))

                arr = each
                arrow = np.vstack([arrow, arr])
                colorVal = scalarMap.to_rgba(dataset['flows'][keys[i]]['time'][-1])
                axis.arrow(arr[1], arr[0], arr[3]-arr[1], arr[2]-arr[0], color = colorVal, head_width=.4*flowScalefactor, width=.2*flowScalefactor, length_includes_head=True)
        i += 1

    xr.sort()
    yr.sort()
    xr[-1] += 1
    yr[-1] += 1

    print("Number of Arrows Plotted: " + str(j))

    for i in range(int(xr[0]), int(xr[-1])):
        for j in range(int(yr[0]), int(yr[-1])):
            if (i+j) in elements[j].isotopes:
                axis.text(i, j, str(i+j) + "-" + elements[j].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                    ha='center', va='center')

    minor_locator1 = AutoMinorLocator(2)
    minor_locator2 = AutoMinorLocator(2)
    axis.xaxis.set_minor_locator(minor_locator1)
    axis.yaxis.set_minor_locator(minor_locator2)
    axis.set_xlim(xr[0], xr[-1])
    axis.set_ylim(yr[0], yr[-1])
    axis.set_xticks(range(int(xr[0]), int(xr[-1])))
    axis.set_yticks(range(int(yr[0]), int(yr[-1])))
    axis.grid(which='minor', color='black', linewidth=0.5)
    fig.colorbar(scalarMap, ax=axis)

def readDecays(path, file):
    decays = {}
    with open(path+file, mode='r') as file:
        csvFile = csv.DictReader(file)
        for lines in csvFile:
            if any(i.isdigit() for i in lines['Half-Life']):
                lines['Half-Life'] = float(lines['Half-Life'])
            if lines['Half-Life (Unit)'] == 'y':
                lines['Half-Life (Unit)'] = 'd'
                lines['Half-Life'] = lines['Half-Life'] * 365
            if lines['Half-Life (Unit)'] == 'd':
                lines['Half-Life (Unit)'] = 'h'
                lines['Half-Life'] = lines['Half-Life'] * 24
            if lines['Half-Life (Unit)'] == 'h':
                lines['Half-Life (Unit)'] = 'm'
                lines['Half-Life'] = lines['Half-Life'] * 60
            if lines['Half-Life (Unit)'] == 'm':
                lines['Half-Life (Unit)'] = 's'
                lines['Half-Life'] = lines['Half-Life'] * 60

            if lines['Half-Life (Unit)'] == 'fs':
                lines['Half-Life (Unit)'] = 'ps'
                lines['Half-Life'] = lines['Half-Life'] / 1e3
            if lines['Half-Life (Unit)'] == 'ps':
                lines['Half-Life (Unit)'] = 'ns'
                lines['Half-Life'] = lines['Half-Life'] / 1e3
            if lines['Half-Life (Unit)'] == 'ns':
                lines['Half-Life (Unit)'] = 'us'
                lines['Half-Life'] = lines['Half-Life'] / 1e3
            if lines['Half-Life (Unit)'] == 'us':
                lines['Half-Life (Unit)'] = 'ms'
                lines['Half-Life'] = lines['Half-Life'] / 1e3
            if lines['Half-Life (Unit)'] == 'ms':
                lines['Half-Life (Unit)'] = 's'
                lines['Half-Life'] = lines['Half-Life'] / 1e3

            if lines['Half-Life (Unit)'] == 's':
                decays[lines['Element']+lines['Atomic Mass (A)']] = {'t1/2': lines['Half-Life'],
                                                                'tErr': lines['Half-Life (Error)']
                                                                }
            else:
                if lines['Half-Life (Unit)'] != "":
                    print("Unit Error")
    return decays

def calcBDecay(iso, abund, time):
    times = [f1['flows'][key]['time'][0] for key in list(f1['flows'].keys())]
    times.sort()
    rates = readDecays(walletcardPath, 'walletcards.csv')
    dt = time - times[times.index(find_nearest(times, time))-1]
    if iso in rates.keys():
        rate = rates[iso]
        lam = ln(2) / rate['t1/2']
        return abund * lam * dt
    else:
        return 0

def calcBDecayGen(inp):
    times = [f1['flows'][key]['time'][0] for key in list(f1['flows'].keys())]
    times.sort()
    rates = readDecays(walletcardPath, 'walletcards.csv')
    while True:
        iso = inp[0][:-2]
        abund = inp[1]
        time = inp[2]
        dt = time - times[times.index(find_nearest(times, time))-1]
        if iso in rates.keys():
            rate = rates[iso]
            lam = ln(2) / rate['t1/2']
            inp = yield abund * lam * dt
        else:
            inp = yield 0

def sortFlows(dataset, num, xrange, yrange, minFlow, growOnly=True):
    length = len(dataset['flows'][num]['p_in'][()])
    filt = np.zeros((length))
    for each in yrange:
        filt += np.where(dataset['flows'][num]['p_in'][()] == each, 1, 0)
        filt += np.where(dataset['flows'][num]['p_out'][()] == each, 1, 0)
    for each in xrange:
        filt += np.where(dataset['flows'][num]['n_in'][()] == each, 1, 0)
        filt += np.where(dataset['flows'][num]['n_out'][()] == each, 1, 0)
    filt = filt * np.where(dataset['flows'][num]['flow'][()] >= minFlow, 1, 0)
    if growOnly:
        filt = filt * np.where(dataset['flows'][num]['p_in'][()]+dataset['flows'][num]['n_in'][()] <= dataset['flows'][num]['p_out'][()]+dataset['flows'][num]['n_out'][()], 1, 0)
    filt = np.where(filt > 0)[0]
    return filt

def graphFlow(dataset, temp, xrange, yrange, fig, axis, minFlowshow, flowScalefactor, cbar=[]):
    temp = float(temp)
    keyTemps ={}
    keys = list(dataset['flows'].keys())
    for key in keys:
        keyTemps[key] = dataset['flows'][key]['temp'][0]
    keyTemps = dict(sorted(keyTemps.items(), key=lambda item: item[1]))
    num = list(keyTemps.keys())[list(keyTemps.values()).index(find_nearest(list(keyTemps.values()), temp))]
    arrow = np.array([0 for i in range(6)])

    xr = list(xrange[1:-1])
    yr = list(yrange[1:-1])
    data = procFlow(dataset, num)
    data = np.hstack([data, np.zeros((len(data[:,0]),1))])
    for i in range(len(dataset['flows'][num]['n_out'][:])-1):
        pin = data[i, 0]
        nin = data[i, 1]
        pout = data[i, 2]
        nout = data[i, 3]
        #if pin >= 16:
            #print("Pin >= 16")
        if pin in yr and pout in yr and nin in xr and nout in xr and data[i,4] > minFlowshow:
            #print("Pin = %5.2f, Pout = %5.2f, Nin = %5.2f, Nout = %5.2f" % (pin, pout, nin, nout))
            arrow = np.vstack([arrow, [data[i]]])

    if np.any(arrow):
        arrow = np.delete(arrow, (0), axis=0)
        adjData = arrow[:, 4].copy()
        adjData = np.vstack([adjData, adjData])
        if isinstance(flowScalefactor, np.ndarray):
            for i in range(0, len(arrow)):
                if flowScalefactor[int(arrow[i, 1]), int(arrow[i, 0])] != 0:
                    adjData[0, i] = np.where(adjData[0, i] >= minFlowshow, adjData[0, i]/np.abs(np.log(flowScalefactor[int(arrow[i, 1]), int(arrow[i, 0])])), 1e-10)
                    adjData[1, i] = np.where(adjData[1, i] >= minFlowshow, adjData[1, i], 0)
                else:
                    adjData[0, i] = 1e-10
                    adjData[1, i] = 0
        else:
            adjData = np.where(adjData >= minFlowshow, adjData*flowScalefactor, 0)
        arrow[:, 5] = adjData[0, :]
        arrow[:, 4] = adjData[1, :]

        """    for i in yrange[1:-1]:
            for j in xrange[1:-1]:
                axis.text(j, i, str(i+j) + "-" + elements[i].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                    ha='center', va='center')"""

        cmap = plt.cm.jet
        if cbar == []:
            cNorm = colors.LogNorm(vmin=1e-6, vmax=np.max(arrow[:,5]))
        else:
            cNorm = colors.LogNorm(vmin=cbar[0], vmax=cbar[1])
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

        maxF = np.max(arrow[:, 5])
        flowScale = maxF / 2.5
        i = 0
        while i < len(arrow):
            if arrow[i,4] == 0:
                arrow = np.delete(arrow, i, axis=0)
            else:
                i += 1

        colorVal = scalarMap.to_rgba(arrow[5])
        axis.quiver(arrow[1:0], [arrow[3]-arrow[1]], [arrow[2]-arrow[0]], color = colorVal)

        minor_locator1 = AutoMinorLocator(2)
        minor_locator2 = AutoMinorLocator(2)
        axis.xaxis.set_minor_locator(minor_locator1)
        axis.yaxis.set_minor_locator(minor_locator2)
        axis.set_xlim(xrange[0], xrange[-1])
        axis.set_ylim(yrange[0], yrange[-1])
        axis.grid(which='minor', color='black', linewidth=0.5)
        fig.colorbar(scalarMap, ax=axis)

    else:
        minor_locator1 = AutoMinorLocator(2)
        minor_locator2 = AutoMinorLocator(2)
        axis.xaxis.set_minor_locator(minor_locator1)
        axis.yaxis.set_minor_locator(minor_locator2)
        axis.set_xlim(xrange[0], xrange[-1])
        axis.set_ylim(yrange[0], yrange[-1])
        axis.grid(which='minor', color='black', linewidth=0.5)
        return None, None


    return cNorm.vmin, cNorm.vmax

def graphFlowGen(dataset, temp, xrange, yrange, fig, axis, minFlowshow, flowScalefactor, cbar=[], pRange=[80,100]):
    temp = float(temp)
    keyTemps ={}
    keys = list(dataset['flows'].keys())
    for key in keys:
        keyTemps[key] = dataset['flows'][key]['temp'][0]
    keyTemps = dict(sorted(keyTemps.items(), key=lambda item: item[1]))
    num = list(keyTemps.keys())[list(keyTemps.values()).index(find_nearest(list(keyTemps.values()), temp))]
    xr = list(xrange[1:-1])
    yr = list(yrange[1:-1])
    cnt = 1000
    logrng = np.logspace(-15, -3, cnt)
    cmap = plt.cm.jet

    filter = sortFlows(dataset, num, xrange, yrange, minFlowshow)

    if cbar == []:
        cNorm = colors.LogNorm(vmin=1e-6, vmax=1e-1)
    else:
        if cbar[0] == -1:
            cNorm = colors.LogNorm(vmin=minFlowshow, vmax=cbar[1])
        else:
            cNorm = colors.LogNorm(vmin=cbar[0], vmax=cbar[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

    pmin = yr[0]
    pmax = yr[-1]
    nmin = xr[0]
    nmax = xr[-1]

    for data in procFlowGenMasked(dataset, num, filter):
        """    for i in yrange[1:-1]:
            for j in xrange[1:-1]:
                axis.text(j, i, str(i+j) + "-" + elements[i].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                    ha='center', va='center')"""

        val = find_nearest(logrng, data[4])
        indx = list(logrng).index(val)/cnt
        if indx >= find_nearest_index(logrng, cbar[0])/cnt:
            colorVal = scalarMap.to_rgba(data[4])
            axis.arrow(data[1], data[0], data[3]-data[1], data[2]-data[0], color = colorVal, head_width=indx*.4, width=indx*.2, length_includes_head=True)

    if xrange[0] == 0 or xrange[-1] == 200-1:
        n = np.where(flowScalefactor > 1e-20)[1]
        n = np.max(n)
        axis.set_xlim(xrange[0], n+1)
    else:
        axis.set_xlim(xrange[0], xrange[-1])

    if yrange[0] == 0 or yrange[-1] == 200-1:
        p = np.where(flowScalefactor > 1e-20)[0]
        p = np.max(p)
        axis.set_ylim(xrange[0], p+1)
    else:
        axis.set_ylim(yrange[0], yrange[-1])

    axis.grid(which='minor', color='black', linewidth=0.5)
    fig.colorbar(scalarMap, ax=axis)
    return minFlowshow

def plotWeakRates(dataset, temp, xrange, yrange, fig, axis, minFlowshow, flowScalefactor, cbar=[]):
    temp = float(temp)
    keys = list(dataset['flows'].keys())
    temps = np.sort([dataset['flows'][i]['temp'][0] for i in keys])
    num = keys[np.where(temps == find_nearest(temps, temp))[0][0]]

    arrow = np.array([0 for i in range(6)])
    bdecays = np.array([0 for i in range(6)])

    xr = list(xrange[1:-1])
    yr = list(yrange[1:-1])
    data = procFlow(dataset, num)
    data = np.hstack([data, np.zeros((len(data[:,0]),1))])
    for i in range(len(dataset['flows'][num]['n_out'][:])-1):
        pin = int(data[i, 0])
        nin = int(data[i, 1])
        pout = int(data[i, 2])
        nout = int(data[i, 3])
        #if pin >= 16:
            #print("Pin >= 16")
        if pin in yr and pout in yr and nin in xr and nout in xr:
            #print("Pin = %5.2f, Pout = %5.2f, Nin = %5.2f, Nout = %5.2f" % (pin, pout, nin, nout))
            arrow = np.vstack([arrow, data[i]])
            if pin - pout == 1 and nout - nin == 1:
                iso = str(periodictable.elements[pin].symbol) + str(pin+nin)
                bdecays = np.vstack([bdecays, [pin, nin, pout, nout, calcBDecay(iso, flowScalefactor[nin, pin], dataset['flows'][num]['time'][0]), 0]])
            else:
                bdecays = np.vstack([bdecays, [0,0,0,0,0,0]])


    arrow = np.delete(arrow, (0), axis=0)
    bdecays = np.delete(bdecays, (0), axis=0)

    adjData = arrow[:, 4].copy()
    adjData = np.vstack([adjData, adjData])
    if isinstance(flowScalefactor, np.ndarray):
        for i in range(0, len(arrow)):
            if flowScalefactor[int(arrow[i, 1]), int(arrow[i, 0])] != 0:
                adjData[0, i] = np.where(adjData[0, i] >= minFlowshow, adjData[0, i]/np.abs(np.log(flowScalefactor[int(arrow[i, 1]), int(arrow[i, 0])])), 1e-30)
                adjData[1, i] = np.where(adjData[1, i] >= minFlowshow, adjData[1, i], 0)
            else:
                adjData[0, i] = 1e-30
                adjData[1, i] = 0
    else:
        adjData = np.where(adjData >= minFlowshow, adjData*flowScalefactor, 0)
    arrow[:, 5] = adjData[0, :]
    arrow[:, 4] = adjData[1, :]

    adjData = bdecays[:, 4].copy()
    adjData = np.vstack([adjData, adjData])
    if isinstance(flowScalefactor, np.ndarray):
        for i in range(0, len(bdecays)):
            if flowScalefactor[int(bdecays[i, 1]), int(bdecays[i, 0])] != 0:
                adjData[0, i] = np.where(adjData[0, i] >= minFlowshow, adjData[0, i]/np.abs(np.log(flowScalefactor[int(bdecays[i, 1]), int(bdecays[i, 0])])), 1e-30)
                adjData[1, i] = np.where(adjData[1, i] >= minFlowshow, adjData[1, i], 0)
            else:
                adjData[0, i] = 1e-30
                adjData[1, i] = 0
    else:
        adjData = np.where(adjData >= minFlowshow, adjData*flowScalefactor, 0)
    bdecays[:, 5] = adjData[0, :]
    bdecays[:, 4] = adjData[1, :]

    """    for i in yrange[1:-1]:
        for j in xrange[1:-1]:
            axis.text(j, i, str(i+j) + "-" + elements[i].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                ha='center', va='center')"""

    cmap = plt.cm.jet
    if cbar == []:
        cNorm = colors.Normalize(vmin=0, vmax=100)
    else:
        cNorm = colors.LogNorm(vmin=cbar[0], vmax=cbar[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

    maxF = np.max(bdecays[:, 5])
    flowScale = maxF / 2.5
    for i in range(len(bdecays)):
        if arrow[i, 4] > 1e-30 and bdecays[i, 4] > 1e-30:
            arr = bdecays[i]
            colorVal = scalarMap.to_rgba((arr[5]/(arrow[i, 5]))*100)
            axis.arrow(arr[1], arr[0], arr[3]-arr[1], arr[2]-arr[0], color = colorVal, head_width=arr[5]*.4/flowScale, width=arr[5]*.2/flowScale, length_includes_head=True)

    minor_locator1 = AutoMinorLocator(2)
    minor_locator2 = AutoMinorLocator(2)
    axis.xaxis.set_minor_locator(minor_locator1)
    axis.yaxis.set_minor_locator(minor_locator2)
    axis.set_xlim(xrange[0], xrange[-1])
    axis.set_ylim(yrange[0], yrange[-1])
    axis.grid(which='minor', color='black', linewidth=0.5)
    fig.colorbar(scalarMap, ax=axis)
    return flowScale

def plotWeakRatesGen(dataset, temp, xrange, yrange, fig, axis, minFlowshow, flowScalefactor, cbar=[]):
    temp = float(temp)
    keyTemps ={}
    keys = list(dataset['flows'].keys())
    for key in keys:
        keyTemps[key] = dataset['flows'][key]['temp'][0]
    keyTemps = dict(sorted(keyTemps.items(), key=lambda item: item[1]))
    num = list(keyTemps.keys())[list(keyTemps.values()).index(find_nearest(list(keyTemps.values()), temp))]
    xr = list(xrange[1:-1])
    yr = list(yrange[1:-1])
    cnt = 1000
    logrng = np.logspace(-20, -5, cnt)
    cmap = plt.cm.jet
    if cbar == []:
        cNorm = colors.Normalize(vmin=0, vmax=100)
    else:
        cNorm = colors.LogNorm(vmin=cbar[0], vmax=cbar[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    calcGen = calcBDecayGen(('H1', 0, 0))
    beta = next(calcGen)

    for data in procFlowGen(dataset,num):
        pin = data[0]
        nin = data[1]
        pout = data[2]
        nout = data[3]
        time = data[4]
        #if pin >= 16:
            #print("Pin >= 16")
        if pin in yr and pout in yr and nin in xr and nout in xr:
            #print("Pin = %5.2f, Pout = %5.2f, Nin = %5.2f, Nout = %5.2f" % (pin, pout, nin, nout))
            arrow = [pin, nin, pout, nout, data[4], 0]
            if np.any(arrow):
                if isinstance(flowScalefactor, np.ndarray):
                    if flowScalefactor[int(pin), int(nin)] != 0:
                        arrow[5] = arrow[4]/np.abs(np.log(flowScalefactor[int(pin), int(nin)]))
                        if pin - pout == 1 and nout - nin == 1:
                            iso = str(periodictable.elements[pin].symbol) + str(pin+nin)
                            beta = calcGen.send((iso, flowScalefactor[int(pin), int(nin)], time))
                            if beta >= minFlowshow:
                                bdecays = [pin, nin, pout, nout, beta, 0]
                            else:
                                bdecays = [0, 0, 0, 0, 0, 0]
                            bdecays[5] = bdecays[4] / np.abs(np.log(flowScalefactor[int(pin), int(nin)]))
                        else:
                            bdecays = [0,0,0,0,0,0]
                    else:
                        arrow[5] = 1e-30
                else:
                    arrow[5] = arrow[4] * flowScalefactor
                    if pin - pout == 1 and nout - nin == 1:
                        iso = str(periodictable.elements[pin].symbol) + str(pin+nin)
                        bdecays = [pin, nin, pout, nout, calcBDecay(iso, flowScalefactor, dataset['flows'][num]['time'][0]), 0]
                        bdecays[5] = bdecays[4] * flowScalefactor

                val = find_nearest(logrng, bdecays[4])
                indx = list(logrng).index(val) / cnt
                if bdecays[4] != 0 and bdecays[4]*100/arrow[4] >= .1:
                    colorVal = scalarMap.to_rgba(bdecays[4]*100/arrow[4])
                    axis.arrow(bdecays[1], bdecays[0], bdecays[3] - bdecays[1], bdecays[2] - bdecays[0], color=colorVal, head_width=indx * .4,
                               width=indx * .2, length_includes_head=True)

    minor_locator1 = AutoMinorLocator(2)
    minor_locator2 = AutoMinorLocator(2)
    axis.xaxis.set_minor_locator(minor_locator1)
    axis.yaxis.set_minor_locator(minor_locator2)
    if xrange[0] == 0:
        n = np.where(flowScalefactor > 1e-25)[1]
        n = np.max(n)
        axis.set_xlim(0, n+1)
    else:
        axis.set_xlim(xrange[0], xrange[-1])

    if yrange[0] == 0:
        p = np.where(flowScalefactor > 1e-25)[0]
        p = np.max(p)
        axis.set_ylim(0, p+1)
    else:
        axis.set_ylim(yrange[0], yrange[-1])

    axis.grid(which='minor', color='black', linewidth=0.5)
    fig.colorbar(scalarMap, ax=axis)
    return

def procNuclei(dataset,t):
    maxSize  = dataset['tracked_nuclei']['A'].shape[0]
    data = np.zeros((maxSize, 5))
    data[:, 0] = dataset['tracked_nuclei']['A'][:]
    data[:, 1] = dataset['tracked_nuclei']['N'][:]
    data[:, 2] = dataset['tracked_nuclei']['Z'][:]
    data[:, 3] = dataset['tracked_nuclei']['Y'][t, :]
    data[0, 4] = dataset['tracked_nuclei']['time'][t]
    return data