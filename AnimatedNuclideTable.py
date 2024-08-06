import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cmx
import numpy as np
import mendeleev
import time
import tqdm
import h5py
import math

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

xr = (0,150)
yr = (0,90)

sourceFile1 = f4

data = np.zeros((200,200))
abund = procNuclei(sourceFile1, 0)
for iso in abund:
    if iso[1] < data.shape[0] or iso[2] < data.shape[1]:
        data[int(iso[2]),int(iso[1])] = iso[3]
    else:
        break

duration = len(sourceFile1['tracked_nuclei']['time'][:])

# create discrete colormap
cmap = plt.cm.jet
norm = colors.LogNorm()
scalarMap = cmx.ScalarMappable(norm=norm, cmap=cmap)

stable = ([],[])
for elem in mendeleev.get_all_elements():
    for iso in elem.isotopes:
        if not iso.is_radioactive and iso.abundance != None:
            stable[1].append(iso.atomic_number)
            stable[0].append(iso.mass_number-iso.atomic_number)

fig, ax = plt.subplots()
ax.scatter(*stable, color="black", s=5)
im = ax.imshow(data, cmap=cmap, norm=norm)
fig.colorbar(scalarMap, ax=ax)

minT = np.min([round(source['tracked_nuclei']['time'][0], 4) for source in [sourceFile1]])
maxT = np.max([round(source['tracked_nuclei']['time'][-1], 4) for source in [sourceFile1]])

def update(frame):
    time = np.linspace(minT, maxT, duration - 1)[frame]
    t4 = np.where(sourceFile1['tracked_nuclei']['time'][()] == find_nearest(sourceFile1['tracked_nuclei']['time'][()], time))[0][0]
    tempT = np.where(sourceFile1['mainout']['time'][()] == find_nearest(sourceFile1['tracked_nuclei']['time'][()], time))[0][0]
    abund = procNuclei(sourceFile1, t4)
    for iso in abund:
        if iso[1] < data.shape[0] or iso[2] < data.shape[1]:
            data[int(iso[2]), int(iso[1])] = iso[3]
        else:
            break
    im.set(data=data, cmap=cmap, norm=norm)
    # draw gridlines

    ax.set_title("Time = " + str(round(time ,4)) + "   Temp = " + str(round(sourceFile1['mainout']['temp'][tempT] ,4)) + " GK")
    #ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
    return (im)

ani = animation.FuncAnimation(fig=fig, func=update, frames=duration-1, interval=30)
#ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=1)
ax.set_xticks(range(xr[0], xr[1], 5))
ax.set_yticks(range(xr[0], xr[1], 4))
ax.set_xlim(xr[0]-1, xr[1])
ax.set_ylim(yr[0]-1, yr[1])
plt.tight_layout()
plt.show()
if input("Save Animation: (Yes/no)").lower()[0] == "y":
    print("Saving Animation (This might take some time) .  .  .")
    ani.save(filename="TimeNuclidePlot" +str(time.time())+".mp4", writer="ffmpeg")