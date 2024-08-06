import numpy
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import h5py
from scipy.optimize import curve_fit

plots =[3,3]
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

def func(x,a,b):
    return a/(np.power(x, b))
def initGraph(dim):
    if dim[0] != 1 or dim[1] != 1:
        fig, axis = plt.subplots(dim[0],dim[1], layout="constrained", figsize=(20,10))
    else:
        fig, axis = plt.subplots(dim[0], dim[1], layout="constrained")
    return fig, axis

def readTrajFile(path, filename):
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
            line = list(map(float, line[0:-1]))
            data.append(line)
            line = f.readline()
    data = np.array(data)
    return data, tags

def loadFile(path, filename):
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

def splitNLast(item):
    return item[0][-1], item[1][-1]

def find_nearest(array,value):
    array = np.sort(array)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def procNuclei(dataset,t):
    maxSize  = dataset['tracked_nuclei']['A'].shape[0]
    data = np.zeros((maxSize, 5))
    data[:, 0] = dataset['tracked_nuclei']['A'][:]
    data[:, 1] = dataset['tracked_nuclei']['N'][:]
    data[:, 2] = dataset['tracked_nuclei']['Z'][:]
    data[:, 3] = dataset['tracked_nuclei']['Y'][t, :]
    data[0, 4] = dataset['tracked_nuclei']['time'][t]
    return data


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


def track_nuclei(dataset, Z, N):
    t1 = np.where(dataset['tracked_nuclei']['Z'][:] == Z)
    t2 = np.where(dataset['tracked_nuclei']['N'][:] == N)
    t3 = 0
    for each in t1:
        for entry in each:
            if entry in t2[0]:
                t3 = entry
                break

    time = dataset['tracked_nuclei']['time'][:]
    Ynuc = [dataset['tracked_nuclei']['Y'][:, t3]][0]
    return time, Ynuc

def massCalc(r, u, up, P, Pp, ptot):
    top = -r**2*(u*up*P + u*up*ptot - Pp(u**2+1))
    bottom = 2*Pp*r + P + ptot
    return top/bottom


fig, axis = initGraph(plots)


f1times = [f1['flows'][str(i)]['temp'][0] for i in range(1, len(list(f1['flows'].keys())))]
f2times = [f2['flows'][str(i)]['temp'][0] for i in range(1, len(list(f2['flows'].keys())))]
f3times = [f3['flows'][str(i)]['temp'][0] for i in range(1, len(list(f3['flows'].keys())))]
f4times = [f4['flows'][str(i)]['temp'][0] for i in range(1, len(list(f4['flows'].keys())))]
f5times = [f5['flows'][str(i)]['temp'][0] for i in range(1, len(list(f5['flows'].keys())))]

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
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'red')

data = procNuclei(f2, time2)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'black', label = "Y @ 3GK",)

data = procNuclei(f3, time3)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'green')

data = procNuclei(f4, time4)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'blue')

data = procNuclei(f5, time5)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
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
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'red', linestyle = 'dashed')

data = procNuclei(f2, time2)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'black', label = "Y @ 1GK", linestyle = 'dashed')

data = procNuclei(f3, time3)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'green', linestyle = 'dashed')

data = procNuclei(f4, time4)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'blue', linestyle = 'dashed')

data = procNuclei(f5, time5)
t = data[0, 4]  # Time Stamp
combined = np.vstack([data[:, 2], data[:, 3]]) # [mass, Abundance]
combined = combineLikeRows(combined, 0, 1)
axis[0][0].plot(combined[0], combined[1], color = 'orange', linestyle = 'dashed')

#axis.plot(finab6['A'][:], finab6['Y'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis[0][0].set_xlabel("Mass Number")
axis[0][0].set_ylabel("Abundance")
axis[0][0].set_xscale("linear")
axis[0][0].set_yscale("log")
axis[0][0].legend()

"""# reordering the labels
handles, labels = axis[0][0].get_legend_handles_labels()
# specify order
yeValues = [f1['mainout']['ye'][0],f2['mainout']['ye'][0],f3['mainout']['ye'][0],f4['mainout']['ye'][0],f5['mainout']['ye'][0]] #f6['mainout']['ye'][0]
order =[sorted(yeValues).index(i) for i in yeValues]

out1 = [x for x in range(0, len(yeValues))]
out2 = [x for x in range(0, len(yeValues))]
for i in order:
    out1[order[i]] = handles[i]
    out2[order[i]] = labels[i]
axis[0][0].legend(handles = out1, labels = out2)"""

f1out = f1['mainout']
f2out = f2['mainout']
f3out = f3['mainout']
f4out = f4['mainout']
f5out = f5['mainout']
#f6out = f6['mainout']
"""
fig1, axis1 = initGraph(plots)
axis1.plot(f1out['time'][:], f1out['temp'][:], color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis1.plot(f2out['time'][:], f2out['temp'][:], color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis1.plot(f3out['time'][:], f3out['temp'][:], color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis1.plot(f4out['time'][:], f4out['temp'][:], color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis1.plot(f5out['time'][:], f5out['temp'][:], color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
axis1.plot(f6out['time'][:], f6out['temp'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis1.set_xlabel("Time")
axis1.set_ylabel("Temperaure")
axis1.set_xscale("linear")
axis1.set_yscale("log")
axis1.legend(handles = out1, labels = out2)

fig3, axis3 = initGraph(plots)
axis3.plot(f1out['time'][:], f1out['dens'][:], color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis3.plot(f2out['time'][:], f2out['dens'][:], color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis3.plot(f3out['time'][:], f3out['dens'][:], color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis3.plot(f4out['time'][:], f4out['dens'][:], color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis3.plot(f5out['time'][:], f5out['dens'][:], color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
axis3.plot(f6out['time'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis3.set_xlabel("Time")
axis3.set_ylabel("Density")
axis3.set_xscale("linear")
axis3.set_yscale("log")
axis3.legend(handles = out1, labels = out2)

fig2, axis2 = initGraph(plots)
axis2.plot(f1out['time'][:], f1out['rad'][:], color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis2.plot(f2out['time'][:], f2out['rad'][:], color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis2.plot(f3out['time'][:], f3out['rad'][:], color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis2.plot(f4out['time'][:], f4out['rad'][:], color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis2.plot(f5out['time'][:], f5out['rad'][:], color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
axis2.plot(f6out['time'][:], f6out['rad'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis2.set_xlabel("Time")
axis2.set_ylabel("Radius")
axis2.set_xscale("linear")
axis2.set_yscale("log")
axis2.legend(handles = out1, labels = out2)
"""

mass1 = (((f1out['rad'][:]*1000)**3)*((f1out['dens'][:]*1000))*5e-31)
mass2 = (((f2out['rad'][:]*1000)**3)*((f2out['dens'][:]*1000))*5e-31)
mass3 = (((f3out['rad'][:]*1000)**3)*((f3out['dens'][:]*1000))*5e-31)
mass4 = (((f4out['rad'][:]*1000)**3)*((f4out['dens'][:]*1000))*5e-31)
mass5 = (((f5out['rad'][:]*1000)**3)*((f5out['dens'][:]*1000))*5e-31)

vel1 = numpy.gradient(f1out['rad'][:], f1out['time'][:])
vel2 = numpy.gradient(f2out['rad'][:], f2out['time'][:])
vel3 = numpy.gradient(f3out['rad'][:], f3out['time'][:])
vel4 = numpy.gradient(f4out['rad'][:], f4out['time'][:])
vel5 = numpy.gradient(f5out['rad'][:], f5out['time'][:])

mof1 = 4*np.pi*(np.power(f1out['rad'][:]*1000*1000,2))*f1out['dens'][:]*vel1*5e-34
mof2 = 4*np.pi*(np.power(f2out['rad'][:]*1000*1000,2))*f2out['dens'][:]*vel2*5e-34
mof3 = 4*np.pi*(np.power(f3out['rad'][:]*1000*1000,2))*f3out['dens'][:]*vel3*5e-34
mof4 = 4*np.pi*(np.power(f4out['rad'][:]*1000*1000,2))*f4out['dens'][:]*vel4*5e-34
mof5 = 4*np.pi*(np.power(f5out['rad'][:]*1000*1000,2))*f5out['dens'][:]*vel5*5e-34

axis[0][1].plot(f1out['rad'][:], mof1, color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis[0][1].plot(f2out['rad'][:], mof2, color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis[0][1].plot(f3out['rad'][:], mof3, color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis[0][1].plot(f4out['rad'][:], mof4, color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis[0][1].plot(f5out['rad'][:], mof5, color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
#axis[0][1].plot(f6out['rad'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis[0][1].set_xlabel("Radius")
axis[0][1].set_ylabel("Mass Outflow Rate")
axis[0][1].set_xscale("linear")
axis[0][1].set_yscale("linear")
axis[0][1].legend()

axis[0][2].plot(f1['tracked_nuclei']['time'][()],np.sum(f1['tracked_nuclei']['Y'][()]/5e6, axis=1))

axis[0][2].plot(f1out['time'][:], (((f1out['rad'][:]*1000)**3)*((f1out['dens'][:]*1000))*5e-31), color = 'red', label = "Mi = " + str(round((((f1out['rad'][0]*1000)**3)*((f1out['dens'][0]*1000))*5e-31), 10)))
axis[0][2].plot(f2out['time'][:], (((f2out['rad'][:]*1000)**3)*((f2out['dens'][:]*1000))*5e-31), color = 'black', label = "Mi = " + str(round((((f2out['rad'][0]*1000)**3)*((f2out['dens'][0]*1000))*5e-31), 10)))
axis[0][2].plot(f3out['time'][:], (((f3out['rad'][:]*1000)**3)*((f3out['dens'][:]*1000))*5e-31), color = 'green', label = "Mi = " + str(round((((f3out['rad'][0]*1000)**3)*((f3out['dens'][0]*1000))*5e-31), 10)))
axis[0][2].plot(f4out['time'][:], (((f4out['rad'][:]*1000)**3)*((f4out['dens'][:]*1000))*5e-31), color = 'blue', label = "Mi = " + str(round((((f4out['rad'][0]*1000)**3)*((f4out['dens'][0]*1000))*5e-31), 10)))
axis[0][2].plot(f5out['time'][:], (((f5out['rad'][:]*1000)**3)*((f5out['dens'][:]*1000))*5e-31), color = 'orange', label = "Mi = " + str(round((((f5out['rad'][0]*1000)**3)*((f5out['dens'][0]*1000))*5e-31), 10)))
#axis3.plot(f6out['time'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis[0][2].set_xlabel("Time")
axis[0][2].set_ylabel("Mass (Solar Masses)")
axis[0][2].set_xscale("linear")
axis[0][2].set_yscale("log")
axis[0][2].legend()



V1s = 5.67e8*np.sqrt(f1out['temp'][:]*.0862)*np.sqrt(f1out['entr'][:])*1e-5
V2s = 5.67e8*np.sqrt(f2out['temp'][:]*.0862)*np.sqrt(f2out['entr'][:])*1e-5
V3s = 5.67e8*np.sqrt(f3out['temp'][:]*.0862)*np.sqrt(f3out['entr'][:])*1e-5
V4s = 5.67e8*np.sqrt(f4out['temp'][:]*.0862)*np.sqrt(f4out['entr'][:])*1e-5
V5s = 5.67e8*np.sqrt(f5out['temp'][:]*.0862)*np.sqrt(f5out['entr'][:])*1e-5



axis[1][0].plot(f1out['rad'][:], vel1,  color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis[1][0].plot(f2out['rad'][:], vel2,  color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis[1][0].plot(f3out['rad'][:], vel3,  color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis[1][0].plot(f4out['rad'][:], vel4,  color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis[1][0].plot(f5out['rad'][:], vel5,  color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
axis[1][0].plot(f1out['rad'][:], V1s, color = 'red', linestyle= 'dashed')
axis[1][0].plot(f2out['rad'][:], V2s, color = 'black', linestyle= 'dashed')
axis[1][0].plot(f3out['rad'][:], V3s, color = 'green', linestyle= 'dashed')
axis[1][0].plot(f4out['rad'][:], V4s, color = 'blue', linestyle= 'dashed')
axis[1][0].plot(f5out['rad'][:], V5s, color = 'orange', linestyle= 'dashed')
V1s = vel1[np.argwhere(np.diff(np.sign(vel1 - V1s))).flatten()[0]]
axis[1][0].axhline(y=V1s, color = 'red', linestyle= 'dashed')
V2s = vel2[np.argwhere(np.diff(np.sign(vel2 - V2s))).flatten()[0]]
axis[1][0].axhline(y=V2s, color = 'black', linestyle= 'dashed')
V3s = vel3[np.argwhere(np.diff(np.sign(vel3 - V3s))).flatten()[0]]
axis[1][0].axhline(y=V3s, color = 'green', linestyle= 'dashed')
V4s = vel4[np.argwhere(np.diff(np.sign(vel4 - V4s))).flatten()[0]]
axis[1][0].axhline(y=V4s, color = 'blue', linestyle= 'dashed')
V5s = vel5[np.argwhere(np.diff(np.sign(vel5 - V5s))).flatten()[0]]
axis[1][0].axhline(y=V5s, color = 'orange', linestyle= 'dashed')
axis[1][0].set_xlabel("Radius")
axis[1][0].set_ylabel("Velocity")
axis[1][0].set_xscale("log")
axis[1][0].set_yscale("log")
axis[1][0].legend()

axis[1][1].plot(f1out['time'][:], f1out['rad'][:], color = 'red')
axis[1][1].plot(f2out['time'][:], f2out['rad'][:], color = 'black')
axis[1][1].plot(f3out['time'][:], f3out['rad'][:], color = 'green')
axis[1][1].plot(f4out['time'][:], f4out['rad'][:], color = 'blue')
axis[1][1].plot(f5out['time'][:], f5out['rad'][:], color = 'orange')
axis[1][1].set_xlabel("Time")
axis[1][1].set_ylabel("Radius")
axis[1][1].set_xscale("linear")
axis[1][1].set_yscale("linear")

"""axis[1][1].plot(f1out['rad'][:]**3, f1out['temp'][:]/f1out['rad'][:]**3, color = 'red')
axis[1][1].plot(f2out['rad'][:]**3, f2out['temp'][:]/f2out['rad'][:]**3, color = 'black')
axis[1][1].plot(f3out['rad'][:]**3, f3out['temp'][:]/f3out['rad'][:]**3, color = 'green')
axis[1][1].plot(f4out['rad'][:]**3, f4out['temp'][:]/f4out['rad'][:]**3, color = 'blue')
axis[1][1].plot(f5out['rad'][:]**3, f5out['temp'][:]/f5out['rad'][:]**3, color = 'orange')
#axis1.plot(f6out['time'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))

popt, pcov = curve_fit(func, f1out['rad'][:]**3, f1out['temp'][:]/f1out['rad'][:]**3)
# residual sum of squares
ss_res = np.sum((f1out['temp'][:]/f1out['rad'][:]**3 - func(f1out['rad'][:]**3, *popt)) ** 2)
# total sum of squares
ss_tot = np.sum((f1out['temp'][:]/f1out['rad'][:]**3 - np.mean(f1out['temp'][:]/f1out['rad'][:]**3)) ** 2)
# r-squared
r2 = 1 - (ss_res / ss_tot)
axis[1][1].plot(f1out['rad'][:]**3, func(f1out['rad'][:]**3, *popt), color='pink', linestyle='dashed', label='fit1: R=%5.3f, a=%5.3f, b=%5.3f' % (r2,popt[0],popt[1]))
popt, pcov = curve_fit(func, f2out['rad'][:]**3, f2out['temp'][:]/f2out['rad'][:]**3)
# residual sum of squares
ss_res = np.sum((f2out['temp'][:]/f2out['rad'][:]**3 - func(f2out['rad'][:]**3, *popt)) ** 2)
# total sum of squares
ss_tot = np.sum((f2out['temp'][:]/f2out['rad'][:]**3 - np.mean(f2out['temp'][:]/f2out['rad'][:]**3)) ** 2)
# r-squared
r2 = 1 - (ss_res / ss_tot)
axis[1][1].plot(f2out['rad'][:]**3, func(f2out['rad'][:]**3, *popt), color='gray', linestyle='dashed', label='fit2: R=%5.3f, a=%5.3f, b=%5.3f' % (r2,popt[0],popt[1]))
popt, pcov = curve_fit(func, f3out['rad'][:]**3, f3out['temp'][:]/f3out['rad'][:]**3)
# residual sum of squares
ss_res = np.sum((f3out['temp'][:]/f3out['rad'][:]**3 - func(f3out['rad'][:]**3, *popt)) ** 2)
# total sum of squares
ss_tot = np.sum((f3out['temp'][:]/f3out['rad'][:]**3 - np.mean(f3out['temp'][:]/f3out['rad'][:]**3)) ** 2)
# r-squared
r2 = 1 - (ss_res / ss_tot)
axis[1][1].plot(f3out['rad'][:]**3, func(f3out['rad'][:]**3, *popt), color='lime', linestyle='dashed', label='fit3: R=%5.3f, a=%5.3f, b=%5.3f' % (r2,popt[0],popt[1]))
popt, pcov = curve_fit(func, f4out['rad'][:]**3, f4out['temp'][:]/f4out['rad'][:]**3)
# residual sum of squares
ss_res = np.sum((f4out['temp'][:]/f4out['rad'][:]**3 - func(f4out['rad'][:]**3, *popt)) ** 2)
# total sum of squares
ss_tot = np.sum((f4out['temp'][:]/f4out['rad'][:]**3 - np.mean(f4out['temp'][:]/f4out['rad'][:]**3)) ** 2)
# r-squared
r2 = 1 - (ss_res / ss_tot)
axis[1][1].plot(f4out['rad'][:]**3, func(f4out['rad'][:]**3, *popt), color='cyan', linestyle='dashed', label='fit4: R=%5.3f, a=%5.3f, b=%5.3f' % (r2,popt[0],popt[1]))
popt, pcov = curve_fit(func, f5out['rad'][:]**3, f5out['temp'][:]/f5out['rad'][:]**3)
# residual sum of squares
ss_res = np.sum((f5out['temp'][:]/f5out['rad'][:]**3 - func(f5out['rad'][:]**3, *popt)) ** 2)
# total sum of squares
ss_tot = np.sum((f5out['temp'][:]/f5out['rad'][:]**3 - np.mean(f5out['temp'][:]/f5out['rad'][:]**3)) ** 2)
# r-squared
r2 = 1 - (ss_res / ss_tot)
axis[1][1].plot(f5out['rad'][:]**3, func(f5out['rad'][:]**3, *popt), color='yellow', linestyle='dashed', label='fit5: R=%5.3f, a=%5.3f, b=%5.3f' % (r2,popt[0],popt[1]))
axis[1][1].set_xlabel("Volume")
axis[1][1].set_ylabel("Pressure")
axis[1][1].set_xscale("linear")
axis[1][1].set_yscale("linear")
axis[1][1].set_xlim(3500,14000)
axis[1][1].set_ylim(.000,.006)
axis[1][1].legend()"""

axis[1][2].plot(f1out['time'][:], f1out['entr'][:], color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis[1][2].plot(f2out['time'][:], f2out['entr'][:], color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis[1][2].plot(f3out['time'][:], f3out['entr'][:], color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis[1][2].plot(f4out['time'][:], f4out['entr'][:], color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis[1][2].plot(f5out['time'][:], f5out['entr'][:], color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
#axis2.plot(f6out['time'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
E1R = f1out['time'][np.where(f1out['entr'] == max(f1out['entr'][:]))][0]
E1F = f1out['time'][np.where(f1out['entr'] == min(f1out['entr'][:]))][0]
T1F = f1out['temp'][np.where(f1out['entr'] == min(f1out['entr'][:]))][0]
T1R = f1out['temp'][np.where(f1out['entr'] == max(f1out['entr'][:]))][0]
axis[1][2].axvline(x=E1F, color = 'red', linestyle= 'dashed')
axis[1][2].axvline(x=E1R, color = 'red', linestyle= 'dashed')
#axis[2][2].axhline(y=T1R, color = 'red', linestyle= 'dashed', label="T@MaxE = %5.2f" % (T1R))
E2R = f2out['time'][np.where(f2out['entr'] == max(f2out['entr'][:]))][0]
E2F = f2out['time'][np.where(f2out['entr'] == min(f2out['entr'][:]))][0]
T2F = f2out['temp'][np.where(f2out['entr'] == min(f2out['entr'][:]))][0]
T2R = f2out['temp'][np.where(f2out['entr'] == max(f2out['entr'][:]))][0]
axis[1][2].axvline(x=E2F, color = 'black', linestyle= 'dashed')
axis[1][2].axvline(x=E2R, color = 'black', linestyle= 'dashed')
#axis[2][2].axhline(y=T2R, color = 'black', linestyle= 'dashed')
E3R = f3out['time'][np.where(f3out['entr'] == max(f3out['entr'][:]))][0]
E3F = f3out['time'][np.where(f3out['entr'] == min(f3out['entr'][:]))][0]
T3F = f3out['temp'][np.where(f3out['entr'] == min(f3out['entr'][:]))][0]
T3R = f1out['temp'][np.where(f3out['entr'] == max(f3out['entr'][:]))][0]
axis[1][2].axvline(x=E3F, color = 'green', linestyle= 'dashed')
axis[1][2].axvline(x=E3R, color = 'green', linestyle= 'dashed')
#axis[2][2].axhline(y=T3R, color = 'green', linestyle= 'dashed')
E4R = f4out['time'][np.where(f4out['entr'] == max(f4out['entr'][:]))][0]
E4F = f4out['time'][np.where(f4out['entr'] == min(f4out['entr'][:]))][0]
T4F = f4out['temp'][np.where(f4out['entr'] == min(f4out['entr'][:]))][0]
T4R = f4out['temp'][np.where(f4out['entr'] == max(f4out['entr'][:]))][0]
axis[1][2].axvline(x=E4F, color = 'blue', linestyle= 'dashed')
axis[1][2].axvline(x=E4R, color = 'blue', linestyle= 'dashed')
#axis[2][2].axhline(y=T4R, color = 'blue', linestyle= 'dashed')
E5R = f5out['time'][np.where(f5out['entr'] == max(f5out['entr'][:]))][0]
E5F = f5out['time'][np.where(f5out['entr'] == min(f5out['entr'][:]))][0]
T5F = f5out['temp'][np.where(f5out['entr'] == min(f5out['entr'][:]))][0]
T5R = f5out['temp'][np.where(f5out['entr'] == max(f5out['entr'][:]))][0]
axis[1][2].axvline(x=E5F, color = 'orange', linestyle= 'dashed')
axis[1][2].axvline(x=E5R, color = 'orange', linestyle= 'dashed')
#axis[2][2].axhline(y=T5R, color = 'orange', linestyle= 'dashed', label="T@MaxE = %5.2f" % (T5R))
axis[1][2].set_xlabel("Time")
axis[1][2].set_ylabel("Entropy")
axis[1][2].set_xscale("linear")
axis[1][2].set_yscale("linear")
axis[1][2].legend()

acc1 = numpy.gradient(vel1, f1out['rad'][:])
acc2 = numpy.gradient(vel2, f2out['rad'][:])
acc3 = numpy.gradient(vel3, f3out['rad'][:])
acc4 = numpy.gradient(vel4, f4out['rad'][:])
acc5 = numpy.gradient(vel5, f5out['rad'][:])

t1, y1 = track_nuclei(f1, 1, 0)
t2, y2 = track_nuclei(f2, 1, 0)
t3, y3 = track_nuclei(f3, 1, 0)
t4, y4 = track_nuclei(f4, 1, 0)
t5, y5 = track_nuclei(f5, 1, 0)

axis[2][0].plot(f1['mainout']['temp'][()], y1,  color = 'red')
axis[2][0].plot(f2['mainout']['temp'][()], y2,  color = 'black')
axis[2][0].plot(f3['mainout']['temp'][()], y3,  color = 'green')
axis[2][0].plot(f4['mainout']['temp'][()], y4,  color = 'blue')
axis[2][0].plot(f5['mainout']['temp'][()], y5,  color = 'orange')
axis[2][0].axvline(x=3.5, color = 'black', linestyle= 'solid')
axis[2][0].axvline(x=1, color = 'black', linestyle= 'dashed')

t1, y1 = track_nuclei(f1, 0, 1)
t2, y2 = track_nuclei(f2, 0, 1)
t3, y3 = track_nuclei(f3, 0, 1)
t4, y4 = track_nuclei(f4, 0, 1)
t5, y5 = track_nuclei(f5, 0, 1)

lnuebar1 = readTrajFile(trajectorypath, traj1)[0][:,6]
lnuebar2 = readTrajFile(trajectorypath, traj2)[0][:,6]
lnuebar3 = readTrajFile(trajectorypath, traj3)[0][:,6]
lnuebar4 = readTrajFile(trajectorypath, traj4)[0][:,6]
lnuebar5 = readTrajFile(trajectorypath, traj5)[0][:,6]

axis[2][0].plot(f1['mainout']['temp'][()], y1,  color = 'red', linestyle='dashed')#, label=f1['mainout']['lnuebar'][0]
axis[2][0].plot(f2['mainout']['temp'][()], y2,  color = 'black', linestyle='dashed')
axis[2][0].plot(f3['mainout']['temp'][()], y3,  color = 'green', linestyle='dashed')
axis[2][0].plot(f4['mainout']['temp'][()], y4,  color = 'blue', linestyle='dashed')
axis[2][0].plot(f5['mainout']['temp'][()], y5,  color = 'orange', linestyle='dashed')

axis[2][0].set_xlabel("Temp")
axis[2][0].set_ylabel("Yp, Yn")
axis[2][0].set_xscale("linear")
axis[2][0].set_yscale("log")
axis[2][0].legend()

axis[2][1].plot(f1out['entr'][:], f1out['yn'][:], color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis[2][1].plot(f2out['entr'][:], f2out['yn'][:], color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis[2][1].plot(f3out['entr'][:], f3out['yn'][:], color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis[2][1].plot(f4out['entr'][:], f4out['yn'][:], color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis[2][1].plot(f5out['entr'][:], f5out['yn'][:], color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
#axis2.plot(f6out['time'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis[2][1].set_xlabel("Radius")
axis[2][1].set_ylabel("Yn")
axis[2][1].set_xscale("linear")
axis[2][1].set_yscale("log")
axis[2][1].legend()

axis[2][2].plot(f1out['rad'][:], f1out['entr'][:], color = 'red', label = "Ye = " + str(round(f1['mainout']['ye'][0], 10)))
axis[2][2].plot(f2out['rad'][:], f2out['entr'][:], color = 'black', label = "Ye = " + str(round(f2['mainout']['ye'][0], 10)))
axis[2][2].plot(f3out['rad'][:], f3out['entr'][:], color = 'green', label = "Ye = " + str(round(f3['mainout']['ye'][0], 10)))
axis[2][2].plot(f4out['rad'][:], f4out['entr'][:], color = 'blue', label = "Ye = " + str(round(f4['mainout']['ye'][0], 10)))
axis[2][2].plot(f5out['rad'][:], f5out['entr'][:], color = 'orange', label = "Ye = " + str(round(f5['mainout']['ye'][0], 10)))
#axis2.plot(f6out['time'][:], f6out['dens'][:], color = 'cyan', label = "Ye = " + str(round(f6['mainout']['ye'][0], 10)))
axis[2][2].set_xlabel("Radius")
axis[2][2].set_ylabel("Entropy")
axis[2][2].set_xscale("log")
axis[2][2].set_yscale("log")
axis[2][2].legend()


"""
fig2, axis2 = initGraph(plots)
sc2 = axis2.scatter(f1out['time'][:],f1out['temp'][:], c=f1out['ylight'][:], cmap='magma', norm=mpl.colors.Normalize(0, 0.0005))
sc3 = axis2.scatter(f5out['time'][:],f5out['temp'][:], c=f5out['ylight'][:], cmap=mpl.cm.cool, norm=mpl.colors.Normalize(0, 0.0005))
sc2.set_label("Low Entropy")
sc3.set_label("High Entropy")
plt.colorbar(sc2)
plt.colorbar(sc3)
axis2.set_xlabel("Time")
axis2.set_ylabel("Temp")
axis2.set_xscale("linear")
axis2.set_yscale("linear")

fig3, axis3 = initGraph(plots)
sc3 = axis3.scatter(f5out['temp'][:], f5out['time'][:], c=f5out['yheavy'][:])
plt.colorbar(sc3)
axis3.set_xlabel("Temperature")
axis3.set_ylabel("Time")
axis3.set_xscale("log")
axis3.set_yscale("log")"""

for i in axis:
    for j in i:
        j.grid()
plt.show()
