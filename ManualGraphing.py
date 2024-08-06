import numpy
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import math
from tqdm import tqdm
import time
import h5py
from scipy.optimize import curve_fit


path1 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0017.dat/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0137.dat/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0257.dat/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0377.dat/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/traj_0747.dat/"
filename = "WinNet_data.h5"

trajectorypath = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"
traj1 = "traj_0017.dat"
traj2 = "traj_0137.dat"
traj3 = "traj_0257.dat"
traj4 = "traj_0377.dat"
traj5 = "traj_0747.dat"

f1 = h5py.File(path1 + filename, 'r')
f2 = h5py.File(path2 + filename, 'r')
f3 = h5py.File(path3 + filename, 'r')
f4 = h5py.File(path4 + filename, 'r')
f5 = h5py.File(path5 + filename, 'r')

def procNuclei(dataset,t):
    maxSize  = dataset['tracked_nuclei']['A'].shape[0]
    data = np.zeros((maxSize, 5))
    data[:, 0] = dataset['tracked_nuclei']['A'][:]
    data[:, 1] = dataset['tracked_nuclei']['N'][:]
    data[:, 2] = dataset['tracked_nuclei']['Z'][:]
    data[:, 3] = dataset['tracked_nuclei']['Y'][t, :]
    data[0, 4] = dataset['tracked_nuclei']['time'][t]
    return data

def sepData(data, row):
    out = []
    for each in data:
        out.append(each[row])
    return out

def loadArrayFile(path, filename):
    data = []
    tags = []
    with open(path+filename, encoding="utf-8") as f:
        line = f.readline()
        while line[0] == "#":
            tags.append(line)
            line = f.readline()
        while line != "":
            line = line.split(" ")
            i = len(line) - 1
            while i >= 0:
                if line[i] == '':
                    del line[i]
                i -= 1
            line[-1] = line[-1][:-1]
            line = list(map(float, line))
            data.append(line)
            line = f.readline()
    data = np.array(data)
    return data, tags

def findColumn(data, item):
    if item in data:
        return data.index(item)
    else:
        return -1

def initGraph(dim):
    fig, axis = plt.subplots(dim[0],dim[1], layout="constrained")
    return fig, axis

def find_nearest(array,value):
    array = np.sort(array)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

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

def findReactionFlow(dataset, key, isoIn, isoOut):
    pin = dataset['flows'][key]['p_in'][:]
    nin = dataset['flows'][key]['n_in'][:]
    pout = dataset['flows'][key]['p_out'][:]
    nout = dataset['flows'][key]['n_out'][:]

    maska = np.where(pin == isoIn[0], 1, 0)
    maskb = np.where(nin == isoIn[1], 1, 0)
    maskc = np.where(pout == isoOut[0], 1, 0)
    maskd = np.where(nout == isoOut[1], 1, 0)
    mask = maska*maskb*maskc*maskd

    if np.any(mask):
        result = np.where(mask > 0)[0][0]
        flow = dataset['flows'][key]['flow'][result]
        return flow

    return None


"""f1times = [f1['flows'][i]['temp'][0] for i in list(f1['flows'].keys())]
f2times = [f2['flows'][i]['temp'][0] for i in list(f2['flows'].keys())]
f3times = [f3['flows'][i]['temp'][0] for i in list(f3['flows'].keys())]
f4times = [f4['flows'][i]['temp'][0] for i in list(f4['flows'].keys())]
f5times = [f5['flows'][i]['temp'][0] for i in list(f5['flows'].keys())]
f1times.sort()
f2times.sort()
f3times.sort()
f4times.sort()
f5times.sort()

temp = 3.5
temp1 = f1times[np.where(f1times == find_nearest(f1times, temp))[0][0]]
time1 = f1['mainout']['time'][np.where(f1['mainout']['temp'] == temp1)[0][0]]
time1 = np.where(f1['tracked_nuclei']['time'][:] == find_nearest(f1['tracked_nuclei']['time'][:], time1))[0]

temp2 = f2times[np.where(f2times == find_nearest(f2times, temp))[0][0]]
time2 = f2['mainout']['time'][np.where(f2['mainout']['temp'] == temp2)[0][0]]
time2 = np.where(f2['tracked_nuclei']['time'][:] == find_nearest(f2['tracked_nuclei']['time'][:], time2))[0]

temp3 = f3times[np.where(f3times == find_nearest(f3times, temp))[0][0]]
time3 = f3['mainout']['time'][np.where(f3['mainout']['temp'] == temp3)[0][0]]
time3 = np.where(f3['tracked_nuclei']['time'][:] == find_nearest(f3['tracked_nuclei']['time'][:], time3))[0]

temp4 = f4times[np.where(f4times == find_nearest(f4times, temp))[0][0]]
time4 = f4['mainout']['time'][np.where(f4['mainout']['temp'] == temp4)[0][0]]
time4 = np.where(f4['tracked_nuclei']['time'][:] == find_nearest(f4['tracked_nuclei']['time'][:], time4))[0]

temp5 = f5times[np.where(f5times == find_nearest(f5times, temp))[0][0]]
time5 = f5['mainout']['time'][np.where(f5['mainout']['temp'] == temp5)[0][0]]
time5 = np.where(f5['tracked_nuclei']['time'][:] == find_nearest(f5['tracked_nuclei']['time'][:], time5))[0]

data = procNuclei(f1, time1)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'red')

data = procNuclei(f2, time2)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'black', label = "Y @ 3GK",)

data = procNuclei(f3, time3)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'green')

data = procNuclei(f4, time4)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'blue')

data = procNuclei(f5, time5)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'orange')

temp = 1.0
temp1 = f1times[np.where(f1times == find_nearest(f1times, temp))[0][0]]
time1 = f1['mainout']['time'][np.where(f1['mainout']['temp'] == temp1)[0][0]]
time1 = np.where(f1['tracked_nuclei']['time'][:] == find_nearest(f1['tracked_nuclei']['time'][:], time1))[0]

temp2 = f2times[np.where(f2times == find_nearest(f2times, temp))[0][0]]
time2 = f2['mainout']['time'][np.where(f2['mainout']['temp'] == temp2)[0][0]]
time2 = np.where(f2['tracked_nuclei']['time'][:] == find_nearest(f2['tracked_nuclei']['time'][:], time2))[0]

temp3 = f3times[np.where(f3times == find_nearest(f3times, temp))[0][0]]
time3 = f3['mainout']['time'][np.where(f3['mainout']['temp'] == temp3)[0][0]]
time3 = np.where(f3['tracked_nuclei']['time'][:] == find_nearest(f3['tracked_nuclei']['time'][:], time3))[0]

temp4 = f4times[np.where(f4times == find_nearest(f4times, temp))[0][0]]
time4 = f4['mainout']['time'][np.where(f4['mainout']['temp'] == temp4)[0][0]]
time4 = np.where(f4['tracked_nuclei']['time'][:] == find_nearest(f4['tracked_nuclei']['time'][:], time4))[0]

temp5 = f5times[np.where(f5times == find_nearest(f5times, temp))[0][0]]
time5 = f5['mainout']['time'][np.where(f5['mainout']['temp'] == temp5)[0][0]]
time5 = np.where(f5['tracked_nuclei']['time'][:] == find_nearest(f5['tracked_nuclei']['time'][:], time5))[0]

data = procNuclei(f1, time1)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'red', linestyle = 'dashed')

data = procNuclei(f2, time2)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'black', label = "Y @ 1GK", linestyle = 'dashed')

data = procNuclei(f3, time3)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'green', linestyle = 'dashed')

data = procNuclei(f4, time4)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'blue', linestyle = 'dashed')

data = procNuclei(f5, time5)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 0], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'orange', linestyle = 'dashed')

#axis.plot(finab6['A'][:], finab6['Y'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis[0][0].set_xlabel("Mass Number")
axis[0][0].set_ylabel("Abundance")
axis[0][0].set_xscale("linear")
axis[0][0].set_yscale("log")
axis[0][0].legend()

dYp, dYn, dYe = [], [], []"""

"""for source in [f1,f2,f3,f4,f5]:
    Ypi = source['mainout']['yp'][np.where(source['mainout']['temp'][()] == find_nearest(source['mainout']['temp'][()], 3.5))[0][0]]
    Ypf = source['mainout']['yp'][np.where(source['mainout']['temp'][()] == find_nearest(source['mainout']['temp'][()], 1))[0][0]]
    dYp.append((Ypi-Ypf)/Ypi)

    Yni = source['mainout']['yn'][np.where(source['mainout']['temp'][()] == find_nearest(source['mainout']['temp'][()], 3.5))][0]
    Ynf = source['mainout']['yn'][np.where(source['mainout']['temp'][()] == find_nearest(source['mainout']['temp'][()], 1))][0]
    dYn.append((Yni-Ynf)/Yni)

    Yei = source['mainout']['ye'][np.where(source['mainout']['temp'][()] == find_nearest(source['mainout']['temp'][()], 3.5))][0]
    Yef = source['mainout']['ye'][np.where(source['mainout']['temp'][()] == find_nearest(source['mainout']['temp'][()], 1))][0]
    dYe.append((Yei-Yef)/Yei)

axis[0][1].plot(f1['mainout']['temp'][()], f1['mainout']['yp'][()],  color = 'red', label = "dYp = {:.0%}".format(dYp[0]))
axis[0][1].plot(f2['mainout']['temp'][()], f2['mainout']['yp'][()],  color = 'black', label = "dYp = {:.0%}".format(dYp[1]))
axis[0][1].plot(f3['mainout']['temp'][()], f3['mainout']['yp'][()],  color = 'green', label = "dYp = {:.0%}".format(dYp[2]))
axis[0][1].plot(f4['mainout']['temp'][()], f4['mainout']['yp'][()],  color = 'blue', label = "dYp = {:.0%}".format(dYp[3]))
axis[0][1].plot(f5['mainout']['temp'][()], f5['mainout']['yp'][()],  color = 'orange', label = "dYp = {:.0%}".format(dYp[4]))
axis[0][1].axvline(x=3.5, color = 'black', linestyle= 'solid')
axis[0][1].axvline(x=1, color = 'black', linestyle= 'dashed')
"""

plots =[2, 2]
fig, axis = initGraph(plots)

Ni56 = (28, 28)
Co56 = (27, 29)
Ni57 = (28, 29)
Co57 = (27, 30)
Ni58 = (28, 30)
Cu59 = (29, 30)
Zn60 = (30, 30)

temps = []
Ni56Co56Flows = []
Co56Ni57Flows = []
Ni57Co57Flows = []
Co57Ni58Flows = []
Ni58Cu59Flows = []
Cu59Ni56Flows = []
Cu59Zn60Flows = []

keys = list(f1['flows'].keys())
keyTimes = {}
for val in keys:
    keyTimes[val] = f1['flows'][val]['time'][0]
keyTimes = dict(sorted(keyTimes.items(), key=lambda item: item[1]))
keys = keyTimes.keys()

for key in tqdm(keys):
    temps.append(f1['flows'][key]['temp'][0])
    Ni56Co56Flow = findReactionFlow(f1, key, Ni56, Co56)
    Co56Ni57Flow = findReactionFlow(f1, key, Co56, Ni57)
    Ni57Co57Flow = findReactionFlow(f1, key, Ni57, Co57)
    Co57Ni58Flow = findReactionFlow(f1, key, Co57, Ni58)
    Ni58Cu59Flow = findReactionFlow(f1, key, Ni58, Cu59)
    Cu59Ni56Flow = findReactionFlow(f1, key, Cu59, Ni56)
    Cu59Zn60Flow = findReactionFlow(f1, key, Cu59, Zn60)

    if Ni56Co56Flow != None:
        Ni56Co56Flows.append(Ni56Co56Flow)
    else:
        Ni56Co56Flows.append(0)

    if Co56Ni57Flow != None:
        Co56Ni57Flows.append(Co56Ni57Flow)
    else:
        Co56Ni57Flows.append(0)

    if Ni57Co57Flow != None:
        Ni57Co57Flows.append(Ni57Co57Flow)
    else:
        Ni57Co57Flows.append(0)

    if Co57Ni58Flow != None:
        Co57Ni58Flows.append(Co57Ni58Flow)
    else:
        Co57Ni58Flows.append(0)

    if Ni58Cu59Flow != None:
        Ni58Cu59Flows.append(Ni58Cu59Flow)
    else:
        Ni58Cu59Flows.append(0)

    if Cu59Ni56Flow != None:
        Cu59Ni56Flows.append(Cu59Ni56Flow)
    else:
        Cu59Ni56Flows.append(0)

    if Cu59Zn60Flow != None:
        Cu59Zn60Flows.append(Cu59Zn60Flow)
    else:
        Cu59Zn60Flows.append(0)



axis[0][0].plot(temps, Ni56Co56Flows, color = 'red', label = "Ni56->Co56")
axis[0][0].plot(temps, Co56Ni57Flows, color = 'blue', label = "Co56->Ni57")
axis[0][0].plot(temps, Ni57Co57Flows, color = 'black', label = "Ni57->Co57")
axis[0][0].plot(temps, Co57Ni58Flows, color = 'green', label = "Co57->Ni58")
axis[0][0].plot(temps, Ni58Cu59Flows, color = 'orange', label = "Ni58->Cu59")
axis[0][0].plot(temps, Cu59Ni56Flows, color = 'pink', label = "Cu59->Ni56")
axis[0][0].plot(temps, Cu59Zn60Flows, color = 'lime', label = "Cu59->Zn60")
axis[0][0].set_xscale('linear')
axis[0][0].set_yscale('log')
axis[0][0].set_ylim(1e-20, 1e1)
axis[0][0].set_xlim(0, 7.1)
axis[0][0].grid()
axis[0][0].legend()


"""temps = []
Ni56Co56Flows = []
Co56Ni57Flows = []
Ni57Co57Flows = []
Co57Ni58Flows = []
Ni58Cu59Flows = []
Cu59Ni56Flows = []

for key in f2['flows'].keys():
    temps.append(f2['flows'][key]['temp'][0])
    Ni56Co56Flow = findReactionFlow(f2, key, Ni56, Co56)
    Co56Ni57Flow = findReactionFlow(f2, key, Co56, Ni57)
    Ni57Co57Flow = findReactionFlow(f2, key, Ni57, Co57)
    Co57Ni58Flow = findReactionFlow(f2, key, Co57, Ni58)
    Ni58Cu59Flow = findReactionFlow(f2, key, Ni58, Cu59)
    Cu59Ni56Flow = findReactionFlow(f2, key, Cu59, Ni56)

    if Ni56Co56Flow != None:
        Ni56Co56Flows.append(Ni56Co56Flow)
    else:
        Ni56Co56Flows.append(0)

    if Co56Ni57Flow != None:
        Co56Ni57Flows.append(Co56Ni57Flow)
    else:
        Co56Ni57Flows.append(0)

    if Ni57Co57Flow != None:
        Ni57Co57Flows.append(Ni57Co57Flow)
    else:
        Ni57Co57Flows.append(0)

    if Co57Ni58Flow != None:
        Co57Ni58Flows.append(Co57Ni58Flow)
    else:
        Co57Ni58Flows.append(0)

    if Ni58Cu59Flow != None:
        Ni58Cu59Flows.append(Ni58Cu59Flow)
    else:
        Ni58Cu59Flows.append(0)

    if Cu59Ni56Flow != None:
        Cu59Ni56Flows.append(Cu59Ni56Flow)
    else:
        Cu59Ni56Flows.append(0)



axis[0][1].plot(temps, Ni56Co56Flows, color = 'red', label = "Ni56Co56Flows")
axis[0][1].plot(temps, Co56Ni57Flows, color = 'blue', label = "Co56Ni57Flows")
axis[0][1].plot(temps, Ni57Co57Flows, color = 'black', label = "Ni57Co57Flows")
axis[0][1].plot(temps, Co57Ni58Flows, color = 'green', label = "Co57Ni58Flows")
axis[0][1].plot(temps, Ni58Cu59Flows, color = 'orange', label = "Ni58Cu59Flows")
axis[0][1].plot(temps, Cu59Ni56Flows, color = 'pink', label = "Cu59Ni56Flows")
axis[0][1].plot(temps, Cu59Zn60Flows, color = 'pink', label = "Cu59Ni56Flows")
axis[0][1].set_xscale('linear')
axis[0][1].set_yscale('log')
axis[0][1].set_xlim(0, 7.1)
axis[0][1].grid()
axis[0][1].legend()
"""

axis[1][0].plot(f1['mainout']['temp'][()], f1['mainout']['ye'][()],  color = 'grey', linestyle='dashed', label="Ye")
axis[1][0].plot(f1['mainout']['temp'][()], f1['mainout']['yn'][()],  color = 'grey', linestyle='dotted', label="Yn")
axis[1][0].plot(f1['mainout']['temp'][()], f1['mainout']['yp'][()],  color = 'grey', label="Yp")
axis[1][0].set_xscale('linear')
axis[1][0].set_yscale('log')
axis[1][0].set_xlim(0, 7.1)
axis[1][0].grid()
axis[1][0].legend()






temps = []
Ni56Co56Flows = []
Co56Ni57Flows = []
Ni57Co57Flows = []
Co57Ni58Flows = []
Ni58Cu59Flows = []
Cu59Ni56Flows = []
Cu59Zn60Flows = []

keys = list(f2['flows'].keys())
keyTimes = {}
for val in keys:
    keyTimes[val] = f2['flows'][val]['time'][0]
keyTimes = dict(sorted(keyTimes.items(), key=lambda item: item[1]))
keys = keyTimes.keys()

for key in tqdm(keys):
    temps.append(f2['flows'][key]['temp'][0])
    Ni56Co56Flow = findReactionFlow(f2, key, Ni56, Co56)
    Co56Ni57Flow = findReactionFlow(f2, key, Co56, Ni57)
    Ni57Co57Flow = findReactionFlow(f2, key, Ni57, Co57)
    Co57Ni58Flow = findReactionFlow(f2, key, Co57, Ni58)
    Ni58Cu59Flow = findReactionFlow(f2, key, Ni58, Cu59)
    Cu59Ni56Flow = findReactionFlow(f2, key, Cu59, Ni56)
    Cu59Zn60Flow = findReactionFlow(f2, key, Cu59, Zn60)

    if Ni56Co56Flow != None:
        Ni56Co56Flows.append(Ni56Co56Flow)
    else:
        Ni56Co56Flows.append(0)

    if Co56Ni57Flow != None:
        Co56Ni57Flows.append(Co56Ni57Flow)
    else:
        Co56Ni57Flows.append(0)

    if Ni57Co57Flow != None:
        Ni57Co57Flows.append(Ni57Co57Flow)
    else:
        Ni57Co57Flows.append(0)

    if Co57Ni58Flow != None:
        Co57Ni58Flows.append(Co57Ni58Flow)
    else:
        Co57Ni58Flows.append(0)

    if Ni58Cu59Flow != None:
        Ni58Cu59Flows.append(Ni58Cu59Flow)
    else:
        Ni58Cu59Flows.append(0)

    if Cu59Ni56Flow != None:
        Cu59Ni56Flows.append(Cu59Ni56Flow)
    else:
        Cu59Ni56Flows.append(0)

    if Cu59Zn60Flow != None:
        Cu59Zn60Flows.append(Cu59Zn60Flow)
    else:
        Cu59Zn60Flows.append(0)



axis[0][1].plot(temps, Ni56Co56Flows, color = 'red', label = "Ni56->Co56")
axis[0][1].plot(temps, Co56Ni57Flows, color = 'blue', label = "Co56->Ni57")
axis[0][1].plot(temps, Ni57Co57Flows, color = 'black', label = "Ni57->Co57")
axis[0][1].plot(temps, Co57Ni58Flows, color = 'green', label = "Co57->Ni58")
axis[0][1].plot(temps, Ni58Cu59Flows, color = 'orange', label = "Ni58->Cu59")
axis[0][1].plot(temps, Cu59Ni56Flows, color = 'pink', label = "Cu59->Ni56")
axis[0][1].plot(temps, Cu59Zn60Flows, color = 'lime', label = "Cu59->Zn60")
axis[0][1].set_xscale('linear')
axis[0][1].set_yscale('log')
axis[0][1].set_ylim(1e-20, 1e1)
axis[0][1].set_xlim(0, 7.1)
axis[0][1].grid()
axis[0][1].legend()

axis[1][1].plot(f2['mainout']['temp'][()], f2['mainout']['ye'][()],  color = 'grey', linestyle='dashed', label="Ye")
axis[1][1].plot(f2['mainout']['temp'][()], f2['mainout']['yn'][()],  color = 'grey', linestyle='dotted', label="Yn")
axis[1][1].plot(f2['mainout']['temp'][()], f2['mainout']['yp'][()],  color = 'grey', label="Yp")
axis[1][1].set_xscale('linear')
axis[1][1].set_yscale('log')
axis[1][1].set_xlim(0, 7.1)
axis[1][1].grid()
axis[1][1].legend()

#axis[0][1].plot(f2['mainout']['temp'][()], f2['mainout']['yn'][()],  color = 'black', linestyle='dashed')
#axis[0][1].plot(f3['mainout']['temp'][()], f3['mainout']['yn'][()],  color = 'green', linestyle='dashed')
#axis[0][1].plot(f4['mainout']['temp'][()], f4['mainout']['yn'][()],  color = 'blue', linestyle='dashed')
#axis[0][1].plot(f5['mainout']['temp'][()], f5['mainout']['yn'][()],  color = 'orange', linestyle='dashed')

"""

axis[0][1].plot(f1['mainout']['temp'][()], f1['mainout']['ye'][()],  color = 'red', linestyle='dotted', label = "dYe = {:.0%}".format(dYe[0]))#, label=f1['mainout']['lnuebar'][0]
axis[0][1].plot(f2['mainout']['temp'][()], f2['mainout']['ye'][()],  color = 'black', linestyle='dotted', label = "dYe = {:.0%}".format(dYe[1]))
axis[0][1].plot(f3['mainout']['temp'][()], f3['mainout']['ye'][()],  color = 'green', linestyle='dotted', label = "dYe = {:.0%}".format(dYe[2]))
axis[0][1].plot(f4['mainout']['temp'][()], f4['mainout']['ye'][()],  color = 'blue', linestyle='dotted', label = "dYe = {:.0%}".format(dYe[3]))
axis[0][1].plot(f5['mainout']['temp'][()], f5['mainout']['ye'][()],  color = 'orange', linestyle='dotted', label = "dYe = {:.0%}".format(dYe[4]))

axis[0][1].set_xlabel("Temp")
axis[0][1].set_ylabel("Yp, Ye")
axis[0][1].set_xscale("linear")
axis[0][1].set_yscale("linear")
axis[0][1].legend()


calculationsPath = "/home/thana1dr/Research_Results&Data/Week6/Calculations/"

f1Ncap, tags = loadArrayFile(calculationsPath, "NcapRatef1.txt")
f2Ncap, tags = loadArrayFile(calculationsPath, "NcapRatef2.txt")
f3Ncap, tags = loadArrayFile(calculationsPath, "NcapRatef3.txt")
f4Ncap, tags = loadArrayFile(calculationsPath, "NcapRatef4.txt")
f5Ncap, tags = loadArrayFile(calculationsPath, "NcapRatef5.txt")

axis[1][0].plot(f1['mainout']['temp'][()], f1Ncap[:,0],  color = 'red')
axis[1][0].plot(f2['mainout']['temp'][()], f2Ncap[:,0],  color = 'black')
axis[1][0].plot(f3['mainout']['temp'][()], f3Ncap[:,0],  color = 'green')
axis[1][0].plot(f4['mainout']['temp'][()], f4Ncap[:,0],  color = 'blue')
axis[1][0].plot(f5['mainout']['temp'][()], f5Ncap[:,0],  color = 'orange')
axis[1][0].axvline(x=3.5, color = 'black', linestyle= 'solid')
axis[1][0].axvline(x=1, color = 'black', linestyle= 'dashed')
axis[1][0].set_xlabel("Temp")
axis[1][0].set_ylabel("Average (n,y) Rate")
axis[1][0].set_xscale("linear")
axis[1][0].set_yscale("log")
axis[1][0].legend()


f1Vcap, tags = loadArrayFile(calculationsPath, "VcapRatef1.txt")
f2Vcap, tags = loadArrayFile(calculationsPath, "VcapRatef2.txt")
f3Vcap, tags = loadArrayFile(calculationsPath, "VcapRatef3.txt")
f4Vcap, tags = loadArrayFile(calculationsPath, "VcapRatef4.txt")
f5Vcap, tags = loadArrayFile(calculationsPath, "VcapRatef5.txt")

axis[1][1].plot(f1['mainout']['temp'][()], f1Vcap[:,0],  color = 'red')
axis[1][1].plot(f2['mainout']['temp'][()], f2Vcap[:,0],  color = 'black')
axis[1][1].plot(f3['mainout']['temp'][()], f3Vcap[:,0],  color = 'green')
axis[1][1].plot(f4['mainout']['temp'][()], f4Vcap[:,0],  color = 'blue')
axis[1][1].plot(f5['mainout']['temp'][()], f5Vcap[:,0],  color = 'orange')
axis[1][1].axvline(x=3.5, color = 'black', linestyle= 'solid')
axis[1][1].axvline(x=1, color = 'black', linestyle= 'dashed')
axis[1][1].set_xlabel("Temp")
axis[1][1].set_ylabel("Average (n,p) Rate")
axis[1][1].set_xscale("linear")
axis[1][1].set_yscale("log")
axis[1][1].legend()"""

plt.savefig("/home/thana1dr/Research_Results&Data/Week9/Plots/NiRegionFlows.png", dpi=200)
plt.show()