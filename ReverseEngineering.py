import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.integrate
from periodictable import elements
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from itertools import islice, batched
import numpy as np
import math
from tqdm import tqdm
from tqdm import trange
import h5py

"""
trajectorypath = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"
traj1 = "traj_0017.dat"
traj2 = "traj_0137.dat"
traj3 = "traj_0257.dat"
traj4 = "traj_0377.dat"
traj5 = "traj_0747.dat"

path1 = "/home/thana1dr/WinNet/runs/vp_proc_detailed/" + traj1 + "/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_detailed/" + traj2 + "/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_detailed/" + traj3 + "/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_detailed/" + traj4 + "/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_detailed/" + traj5 + "/"
filename = "WinNet_data.h5"
"""


trajectorypath = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"
traj1 = "traj_0017.dat"
traj2 = "traj_0137.dat"
traj3 = "traj_0257.dat"
traj4 = "traj_0377.dat"
traj5 = "traj_0747.dat"

path1 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/" + traj1 + "/"
path2 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/" + traj2 + "/"
path3 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/" + traj3 + "/"
path4 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/" + traj4 + "/"
path5 = "/home/thana1dr/WinNet/runs/vp_proc_2_flow/" + traj5 + "/"
filename = "WinNet_data.h5"


SummedFlowsPath = "/home/thana1dr/Research_Results&Data/Week8/Automatic/ReactionLists/"

f1 = h5py.File(path1 + filename, 'r')
f2 = h5py.File(path2 + filename, 'r')
f3 = h5py.File(path3 + filename, 'r')
f4 = h5py.File(path4 + filename, 'r')
f5 = h5py.File(path5 + filename, 'r')

path = "/home/thana1dr/Research_Results&Data/Week8/Automatic/ReactionLists/"

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

def procNuclei(dataset,t):
    maxSize  = dataset['tracked_nuclei']['A'].shape[0]
    data = np.zeros((len(t), maxSize, 5))
    data[0, :, 0] = dataset['tracked_nuclei']['A'][:]
    data[0, :, 1] = dataset['tracked_nuclei']['N'][:]
    data[0, :, 2] = dataset['tracked_nuclei']['Z'][:]
    data[:, :, 3] = dataset['tracked_nuclei']['Y'][t, :]
    data[:, 0, 4] = dataset['tracked_nuclei']['time'][t]
    return data

def NucleiDataExists(dataset, Z, N):
    maska = np.where(dataset['tracked_nuclei']['N'][:] == N, 1, 0)
    maskb = np.where(dataset['tracked_nuclei']['Z'][:] == Z, 1, 0)
    mask = np.where(maska*maskb > 0)[0]

    if np.any(mask):
        return True
    else:
        return False

def sepData(data, row):
    out = []
    for each in data:
        out.append(each[row])
    return out

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

def find_nearest_index(array,value):
    array = np.sort(array)
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def findAllReactions(dataset, fileout):
    keys = list(dataset['flows'].keys())
    with open(path + fileout, 'w', encoding="utf-8") as f:
        print("Saving all Reactions: ")
        for key in tqdm(keys):
            pin = dataset['flows'][key]['p_in'][:]
            nin = dataset['flows'][key]['n_in'][:]
            pout = dataset['flows'][key]['p_out'][:]
            nout = dataset['flows'][key]['n_out'][:]
            for i in range(pin.shape[0]):
                f.write(str(pin[i]) + " " + str(nin[i]) + " " + str(pout[i]) + " " + str(nout[i]))
                f.write("\n")
    return

def findUniqueReactions(filein, fileout):
    with open(path+filein, "r") as f:
        with open(path + fileout, "w") as g:
            lines = set()
            t = tqdm()
            print("Processing and Saving Unique Reactions: ")
            for nlines in batched(f, 1000):
                for i, each in enumerate(nlines):
                    if each not in lines:
                        lines.add(each)
                        g.write(each)
                    t.update(1)
            t.close()
            return

def findAllReactionsGen(dataset):
    keys = list(dataset['flows'].keys())
    reacs = []
    for key in tqdm(keys):
        reacCount = dataset['flows'][key]['p_in'].shape[0]
        pin = dataset['flows'][key]['p_in'][:]
        nin = dataset['flows'][key]['n_in'][:]
        pout = dataset['flows'][key]['p_out'][:]
        nout = dataset['flows'][key]['n_out'][:]
        for i in range(reacCount):
            if [int(pin[i]), int(nin[i]), int(pout[i]), int(nout[i])] not in reacs:
                reacs.append([int(pin[i]), int(nin[i]), int(pout[i]), int(nout[i])])
        yield reacs

def MaskReactionGen(dataset, reactions):
    keys = list(dataset['flows'].keys())
    keyTimes = {}
    for val in keys:
        keyTimes[val] = dataset['flows'][val]['time'][0]
    keyTimes = dict(sorted(keyTimes.items(), key=lambda item: item[1]))

    pinAll = [dataset['flows'][key]['p_in'][:] for key in keys]
    ninAll = [dataset['flows'][key]['n_in'][:] for key in keys]
    poutAll = [dataset['flows'][key]['p_out'][:] for key in keys]
    noutAll = [dataset['flows'][key]['n_out'][:] for key in keys]
    for slice in batched(reactions, 10):
        masks = []
        impKeys = []
        for each in tqdm(slice):
            pinw = each[0]
            ninw = each[1]
            poutw = each[2]
            noutw = each[3]
            for i in range(len(keyTimes)):
                key = keys[i]
                pin = pinAll[i]
                nin = ninAll[i]
                pout = poutAll[i]
                nout = noutAll[i]

                maska = np.where(pin == pinw, 1, 0)
                maskb = np.where(nin == ninw, 1, 0)
                maskc = np.where(pout == poutw, 1, 0)
                maskd = np.where(nout == noutw, 1, 0)
                mask = maska*maskb*maskc*maskd
                mask = np.where(mask > 0)[0]
                if len(mask) != 0:
                    impKeys.append(key)
                    masks.append(mask[0])
        yield (masks, impKeys)

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
    for i in range(maxSize):
        data[0] = pin[i]
        data[1] = nin[i]
        data[2] = pout[i]
        data[3] = nout[i]
        data[4] = fout[i]
        data[5] = time
        yield data

def procFlowGenMem(dataset, key):
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
    for i in range(maxSize):
        data[0] = pin[i]
        data[1] = nin[i]
        data[2] = pout[i]
        data[3] = nout[i]
        data[4] = fout[i]
        data[5] = time
        yield data

def procFlowSumMasked(dataset, key, masks):
    keys = list(dataset['flows'].keys())
    keyTimes = {}
    for val in keys:
        keyTimes[val] = dataset['flows'][val]['time'][0]
    keyTimes = dict(sorted(keyTimes.items(), key=lambda item: item[1]))
    if isinstance(key, int):
        key = keys[key]

    pin = []
    nin = []
    pout = []
    nout = []
    fout = 0
    dt = 0
    totFlow = []

    if len(masks) == 0:
        return []
    PINVALUES = [dataset['flows'][key[i]]['p_in'][masks[i]] for i in range(len(masks))]
    NINVALUES = [dataset['flows'][key[i]]['n_in'][masks[i]] for i in range(len(masks))]
    POUTVALUES = [dataset['flows'][key[i]]['p_out'][masks[i]] for i in range(len(masks))]
    NOUTVALUES = [dataset['flows'][key[i]]['n_out'][masks[i]] for i in range(len(masks))]

    PinChange = np.where(np.roll(PINVALUES,1) != PINVALUES)[0]
    NinChange = np.where(np.roll(NINVALUES, 1) != NINVALUES)[0]
    PoutChange = np.where(np.roll(POUTVALUES, 1) != POUTVALUES)[0]
    NoutChange = np.where(np.roll(NOUTVALUES, 1) != NOUTVALUES)[0]

    Changes = np.append(PinChange, NinChange)
    Changes = np.append(Changes, PoutChange)
    Changes = np.append(Changes, NoutChange)
    Changes = np.unique(Changes)[1:-2]

    for i in trange(len(masks)):
        mask = masks[i]
        val = key[i]
        if len(pin) != 0:

            pinNew = dataset['flows'][val]['p_in'][mask]
            ninNew = dataset['flows'][val]['n_in'][mask]
            poutNew = dataset['flows'][val]['p_out'][mask]
            noutNew = dataset['flows'][val]['n_out'][mask]

            if i in Changes:
                pin.append(pinNew)
                nin.append(ninNew)
                pout.append(poutNew)
                nout.append(noutNew)
                totFlow.append(fout)
                fout = 0
                dt = 0
            elif i == len(masks)-1:
                totFlow.append(fout)
                fout = 0
                dt = 0

        else:
            pinNew = dataset['flows'][val]['p_in'][mask]
            ninNew = dataset['flows'][val]['n_in'][mask]
            poutNew = dataset['flows'][val]['p_out'][mask]
            noutNew = dataset['flows'][val]['n_out'][mask]
            pin.append(pinNew)
            nin.append(ninNew)
            pout.append(poutNew)
            nout.append(noutNew)

        flow = (dataset['flows'][val]['flow'][mask])
        time = dataset['flows'][val]['time'][-1]
        if find_nearest_index(list(keyTimes.values()), time) > 1:
            dt = (time - keyTimes[str(find_nearest_index(list(keyTimes.values()), time)-1)])
        else:
            dt = (keyTimes[str(2)] - keyTimes[str(1)])
        fout += (flow)*dt


    """if len(totFlow) <= 0:
        return []"""
    data = np.zeros((5, len(pin)))

    data[0] = pin
    data[1] = nin
    data[2] = pout
    data[3] = nout
    data = np.unique(data, axis=1)
    data[4] = totFlow
    return data

def sumFlowsByTime(dataset, reactions, outPath, outFile):
    keys = list(dataset['flows'].keys())
    keyTimes = {}
    for val in keys:
        keyTimes[val] = dataset['flows'][val]['time'][0]
    keyTimes = dict(sorted(keyTimes.items(), key=lambda item: item[1]))

    flows = np.zeros((len(reactions)))
    flowKeys = [str(reactions[i, :].astype(int).tolist()) for i in range(reactions.shape[0])]
    flowIndexes = [i for i in range(len(flowKeys))]
    flowDict = {}

    for i in range(len(flowKeys)):
        flowDict[flowKeys[i]] = flowIndexes[i]

    for time in tqdm(keyTimes):
        now = keyTimes[time]
        if time == '1' or not np.array_equal(dataset['flows'][time]['p_in'][:], pin):
            pin = dataset['flows'][time]['p_in'][:]
            nin = dataset['flows'][time]['n_in'][:]
            pout = dataset['flows'][time]['p_out'][:]
            nout = dataset['flows'][time]['n_out'][:]
            fout = dataset['flows'][time]['flow'][:]

        if find_nearest_index(list(keyTimes.values()), now) > 1:
            dt = (now - keyTimes[str(find_nearest_index(list(keyTimes.values()), now)-1)])
        else:
            dt = (keyTimes[str(2)] - keyTimes[str(1)])

        for i in range(pin.shape[0]):
            reac = [pin[i], nin[i], pout[i], nout[i]]
            index = flowDict[str(reac)]
            flows[index] += fout[i]*dt

    with open(outPath + outFile, 'w', encoding="utf-8") as f:
        for i in range(len(flows)):
            for each in flowKeys[i].split(','):
                each = ''.join(filter(lambda x: x.isdigit(), each))
                f.write(str(each) + " ")
            f.write(str(flows[i]))
            f.write("\n")
            f.flush()
    return flows

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

def writeFile(path, filename, lines):
    with open(path+filename, 'w', encoding="utf-8") as f:
        for line in lines:
            f.write(str(line))
            f.write("\n")
    return

keys = list(f1['flows'].keys())
keyTimes = {}
for val in tqdm(keys):
    keyTimes[val] = f1['flows'][val]['time'][0]
keyTimes = dict(sorted(keyTimes.items(), key=lambda item: item[1]))
keys = list(keyTimes.keys())

def procIsotopicFlow(dataset, isotope, time):
    pin = dataset['flows'][keys[time]]['p_in'][:]
    nin = dataset['flows'][keys[time]]['n_in'][:]
    pout = dataset['flows'][keys[time]]['p_out'][:]
    nout = dataset['flows'][keys[time]]['n_out'][:]
    fout = dataset['flows'][keys[time]]['flow'][:]
    tempStamp = dataset['flows'][keys[time]]['temp'][0]
    timeStamp = dataset['flows'][keys[time]]['time'][0]

    maska = np.where(pin == isotope[0], 1, 0)
    maskb = np.where(nin == isotope[1], 1, 0)

    maskc = np.where(pout == isotope[0], 1, 0)
    maskd = np.where(nout == isotope[1], 1, 0)

    maske = np.where(pout-pin == isotope[0], 1, 0)
    maskf = np.where(nout-nin == isotope[1], 1, 0)

    mask = np.where(maska*maskb + maskc*maskd + maske*maskf > 0)[0]

    pin = pin[mask]
    nin = nin[mask]
    pout = pout[mask]
    nout = nout[mask]
    fout = fout[mask]

    flow = 0

    for i in range(mask.shape[0]):
        if (pin[i], nin[i]) == isotope:
            flow -= fout[i]
        elif (pout[i], nout[i]) == isotope:
            flow += fout[i]
        elif (pin[i]-pout[i], nin[i]-nout[i]) == isotope:
            flow += fout[i]
        elif(pout[i]-pin[i], nout[i]-nin[i]) == isotope:
            flow -= fout[i]

    if abs(flow) < 1e-25:
        flow = 0

    return flow, timeStamp, tempStamp

def procTotIsoFlows(reactions, iso):
    p = iso[0]
    n = iso[1]
    totFlow = 0

    for each in reactions:
        if each[0] == p and each[1] == n:
            totFlow -= each[4]
        elif each[2] == p and each[3] == n:
            totFlow += each[4]

    return totFlow

def procTotIsoFlowsRegion(reactions, iso, xrange, yrange):
    p = iso[0]
    n = iso[1]
    totFlow = 0

    for each in reactions:
        if each[0] == p and each[1] == n and each[2] in yrange and each[3] in xrange:
            totFlow -= each[4]
        elif each[2] == p and each[3] == n and each[0] in yrange and each[1] in xrange:
            totFlow += each[4]

    return totFlow


def buildReactionMatrix(reactions):
    isos = reactions[:, :2]
    isos = np.vstack([isos, reactions[:, 2:4]])
    isos = np.unique(isos, axis=0)
    isoNum = isos.shape[0]

    isoDict = {}

    for i in range(isos.shape[0]):
        isoDict[str(isos[i])] = i

    ReacMatrix = np.zeros((isoNum, isoNum))
    FlowMatrix = np.zeros((isoNum, isoNum))

    for i in trange(len(reactions)):
        reac = reactions[i]

        if ReacMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] != 0:
            FlowMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] -= reac[4]
            FlowMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] -= reac[4]


            if FlowMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] < 0:
                FlowMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] = FlowMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] * -1
                ReacMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] = ReacMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] * -1

            if FlowMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] < 0:
                FlowMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] = FlowMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] * -1
                ReacMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] = ReacMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] * -1

            #print("Subtracting Inverse Reactions")

        else:
            ReacMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] = -1
            ReacMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] = 1
            FlowMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] = reac[4]
            FlowMatrix[isoDict[str(reac[2:4])], isoDict[str(reac[0:2])]] = reac[4]

    return ReacMatrix, FlowMatrix, isoDict

def buildTestReactionMatrix(reactions):
    xlen = 4
    ylen = 4
    isos = np.zeros((xlen*ylen, 2))

    i = 0
    for p in range(27,27+ylen):
        for n in range(28, 28+xlen):
            isos[i] = [p, n]
            i += 1

    isoNum = isos.shape[0]

    isoDict = {}

    for i in range(isos.shape[0]):
        isoDict[str(isos[i])] = i

    ReacMatrix = np.zeros((isoNum, isoNum))
    FlowMatrix = np.zeros((isoNum, isoNum))

    for i in trange(len(reactions)):
        reac = reactions[i]

        maska = np.zeros((len(isos)))
        maskb = np.zeros((len(isos)))
        for i in range(len(isos)):
            if (reac[0] == isos[i][0] and reac[1] == isos[i][1]):
                maska[i] = 1
            if (reac[2] == isos[i][0] and reac[3] == isos[i][1]):
                maskb[i] = 1

        indxa = np.where(maska > 0)[0]
        indxb = np.where(maskb > 0)[0]
        # Sample:   ReacMatrix[isoDict[str(reac[0:2])], isoDict[str(reac[2:4])]] = -1

        if len(indxa) > 0 and len(indxb) > 0:
            if ReacMatrix[int(indxa[0]), int(indxb[0])] != 0:
                # If the reaction Matrix already contains a value in this position add flow instead of overwriting
                # Reac matrix term insures directionality of flow remains true
                #print("Changing Entry: " + str(FlowMatrix[int(indxa[0]), int(indxb[0])]))
                FlowMatrix[int(indxa[0]), int(indxb[0])] -= ReacMatrix[int(indxa[0]), int(indxb[0])]*reac[4]
                #print("New Value: " + str(FlowMatrix[int(indxa[0]), int(indxb[0])]))

                # After adding/subtracting flow if flow values are now negative multiply Reac and Flow terms by -1 to keep positive
                # Save this processing for the end because it should be simpler then?
                """ if FlowMatrix[int(indxa[0]), int(indxb[0])] < 0:
                    FlowMatrix[int(indxa[0]), int(indxb[0])] *= -1
                    ReacMatrix[int(indxa[0]), int(indxb[0])] *= -1
                    FlowMatrix[int(indxb[0]), int(indxa[0])] *= -1
                    ReacMatrix[int(indxb[0]), int(indxa[0])] *= -1"""

            ReacMatrix[int(indxa[0]), int(indxb[0])] = 1
            FlowMatrix[int(indxa[0]), int(indxb[0])] += reac[4]

            ReacMatrix[int(indxb[0]), int(indxa[0])] = -1
            FlowMatrix[int(indxb[0]), int(indxa[0])] += reac[4]

    for i in range(FlowMatrix.shape[0]):
        for j in range(FlowMatrix.shape[1]):
            if FlowMatrix[i, j] < 0:
                print("Negative Value Found: " + str(FlowMatrix[i, j]))
                FlowMatrix[i, j] *= -1
                ReacMatrix[i, j] *= -1

    return ReacMatrix, FlowMatrix, isoDict

def graphSummedFlow(dataset, xrange, yrange, fig, axis, minFlowshow, flowScalefactor, cbar=[]):
    xr = list(xrange[1:-1])
    yr = list(yrange[1:-1])
    cnt = 1000
    logrng = np.logspace(-15, -3, cnt)
    cmap = plt.cm.jet

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

    for i in yrange[1:-1]:
        for j in xrange[1:-1]:
            axis.text(j, i, str(i+j) + "-" + elements[i].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                ha='center', va='center')

    for data in dataset:


        val = find_nearest(logrng, data[4])
        indx = list(logrng).index(val)/cnt
        if data[4] >= minFlowshow:
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

def plotSolvedRates(dataset, xrange, yrange, fig, axis, minFlowshow, flowScalefactor, cbar=[]):
    xr = list(xrange[1:-1])
    yr = list(yrange[1:-1])
    cnt = 1000
    logrng = np.logspace(-15, -3, cnt)
    cmap = plt.cm.jet

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

    for i in yrange[1:-1]:
        for j in xrange[1:-1]:
            axis.text(j, i, str(i+j) + "-" + elements[i].symbol, fontsize="xx-small", bbox={'facecolor':'white','alpha':1,'edgecolor':'none','pad':1},
                ha='center', va='center')

    for data in dataset:
        val = find_nearest(logrng, data[4])
        indx = list(logrng).index(val)/cnt
        if data[4] >= minFlowshow:
            colorVal = scalarMap.to_rgba(data[4])
            if data[3] == data[1]:
                axis.arrow(data[1], data[2], data[3] - data[1], data[0] - data[2], color=colorVal, head_width=indx * .4,
                           width=indx * .2, length_includes_head=True)
            else:
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



plots = [2,1]
xrange = range(27, 33)
yrange = range(26, 32)
fig, axis = initGraph(plots)
minFlowShow = 1e-3
maxFlowShow = 1e-1


reacs = loadFile(path, "SummedFlowsF1.dat")[0]
reacs = np.unique(reacs, axis=0)
ReacMatrix, FlowMatrix, isos = buildTestReactionMatrix(reacs)
ReacMatrix[-1, -1] = 1
isoList = list(isos.keys())
TotMatrix = ReacMatrix * FlowMatrix

FlowPercMatrix = np.zeros((ReacMatrix.shape[0], ReacMatrix.shape[1]))

for j in range(ReacMatrix.shape[1]):
    rowSum = np.sum(FlowMatrix[j])
    for i in range(ReacMatrix.shape[0]):
        FlowPercMatrix[j][i] = FlowMatrix[j][i] / rowSum

EndYs = np.zeros((ReacMatrix.shape[0]))
f1Sums = loadFile(SummedFlowsPath, "SummedFlowsF1.dat")[0]
for i in range(ReacMatrix.shape[0]):
    iso = isoList[i]
    iso = iso.replace('.', '').strip('[').strip(']').split(" ")
    iso = list(map(int, iso))
    p = iso[0]
    n = iso[1]

    totFlow = procTotIsoFlowsRegion(f1Sums, (p,n), xrange, yrange)
    EndYs[i] = abs(totFlow)

p = 28
n = 28

EndYs[0] -= np.sum(EndYs[1:])
EndYs[isos['['+str(p)+'. '+str(n)+'.]']] *= (1/3)

MatrixTotals = np.sum(TotMatrix, axis=1)
diff = EndYs - MatrixTotals

solution = np.linalg.solve(TotMatrix, EndYs)
#solution = np.ones((solution.shape[0]))
print(solution)

# Updating matrix with Solution
solvedRates = TotMatrix
# solvedRates = np.matmul(TotMatrix, solution)

plotData = np.zeros((solvedRates.shape[0]*solvedRates.shape[1], 5))
reacsConsidered= {}
k = 0
for i in range(solvedRates.shape[0]):
    for j in range(solvedRates.shape[1]):
        inp = isoList[i].replace('.', '').strip('[').strip(']').split(" ")
        out = isoList[j].replace('.', '').strip('[').strip(']').split(" ")
        inp = list(map(int, inp))
        out = list(map(int, out))

        if solvedRates[i, j] < 0: # If flow is negative swap direction
            plotData[k, :] = [out[0], out[1], inp[0], inp[1], -solvedRates[i, j]]
        elif inp != out:
            plotData[k] = [inp[0], inp[1], out[0], out[1], solvedRates[i, j]]

        k += 1

# Remove inverse reactions
reacs = plotData[:, :4]
reacStrings = [str(plotData[i, 0]) + " " + str(plotData[i, 1]) + " , " + str(plotData[i, 2]) + " " + str(plotData[i, 3]) for i in range(plotData.shape[0])]
for i in range(reacs.shape[0]): # Iterate through reactions
    reac = reacs[i]
    search = str(reac[2]) + " " + str(reac[3]) + " , " + str(reac[0]) + " " + str(reac[1])
    if search in reacStrings:
        index = reacStrings.index(search)
        if plotData[i, 4] > plotData[index, 4]:
            plotData[i, 4] -= plotData[index, 4]
            reacStrings[index] = ""
            plotData[index] = [0, 0, 0, 0, 0]
        else:
            plotData[index, 4] -= plotData[i, 4]
            reacStrings[i] = ""
            plotData[i] = [0, 0, 0, 0, 0]

plotData = np.unique(plotData, axis=0)
plotSolvedRates(plotData, xrange, yrange, fig, axis[0], minFlowShow, 1, cbar=(minFlowShow, maxFlowShow))

minFlowshow1 = graphSummedFlow(f1Sums, xrange, yrange, fig, axis[1], minFlowShow, 1, cbar=(minFlowShow, maxFlowShow))





minor_locator1 = AutoMinorLocator(2)
minor_locator2 = AutoMinorLocator(2)
axis[0].xaxis.set_minor_locator(minor_locator1)
axis[0].yaxis.set_minor_locator(minor_locator2)
axis[1].xaxis.set_minor_locator(minor_locator1)
axis[1].yaxis.set_minor_locator(minor_locator2)
plt.show()
print("Done")



"""reacs = loadFile(path, "SummedFlowsF1.dat")[0]
reacs = np.unique(reacs, axis=0)
ReacMatrix, FlowMatrix, isos = buildReactionMatrix(reacs)
TotMatrix = ReacMatrix * FlowMatrix
SupposedYs = np.sum(TotMatrix, axis=1)

BeginningYs = []

abund = procNuclei(f1, 0)
for each in list(isos.keys()):
    each = each.strip().split(" ")

    i = 0
    while i < len(each):
        string = each[i]
        res = any(chr.isdigit() for chr in string)
        if not res:
            each.pop(i)
        else:
            i += 1

    if len(each) == 2:
        p = int(''.join(filter(lambda x: x.isdigit(), each[0].strip())))
        n = int(''.join(filter(lambda x: x.isdigit(), each[1].strip())))

        maska = np.where(abund[:, 2] == p, 1, 0)
        maskb = np.where(abund[:, 1] == n, 1, 0)
        if np.any(maska*maskb):
            indx = np.where(maska*maskb > 0)[0][0]
            BeginningYs.append(abund[indx, 3])
        else:
            BeginningYs.append(0)

BeginningYs = np.array(BeginningYs)
ActualYs = []

abund = procNuclei(f1, -1)
for each in list(isos.keys()):
    each = each.strip().split(" ")

    i = 0
    while i < len(each):
        string = each[i]
        res = any(chr.isdigit() for chr in string)
        if not res:
            each.pop(i)
        else:
            i += 1

    if len(each) == 2:
        p = int(''.join(filter(lambda x: x.isdigit(), each[0].strip())))
        n = int(''.join(filter(lambda x: x.isdigit(), each[1].strip())))

        maska = np.where(abund[:, 2] == p, 1, 0)
        maskb = np.where(abund[:, 1] == n, 1, 0)
        if np.any(maska*maskb):
            indx = np.where(maska*maskb > 0)[0][0]
            ActualYs.append(abund[indx, 3])
        else:
            ActualYs.append(0)

ActualYs = np.array(ActualYs)

PercentDiff = (ActualYs-SupposedYs)/ActualYs*100

PercentDiff = PercentDiff[np.where(np.isfinite(PercentDiff))]

AvgDiff = np.average(PercentDiff)"""

"""reacs = loadFile(path, "SummedFlowsF1.dat")[0]
isos = reacs[:, :2]
isos = np.vstack([isos, reacs[:, 2:4]])
isos = np.unique(isos, axis=0)
"""



"""j = 0

trueTimes = f1['mainout']['time'][:]
trueTemps = f1['mainout']['temp']

AbunLen = len(f1['flows']) + 1
abund = procNuclei(f1, [i for i in range(AbunLen)])
times = abund[:, 0, 4]


for each in isos: # 300:
    p = int(each[0])
    n = int(each[1])
    print("Starting: " + str(p + n) + "-" + periodictable.elements[p].symbol)

    maska = np.where(abund[0, :, 2] == p, 1, 0)
    maskb = np.where(abund[0, :, 1] == n, 1, 0)

    if np.any(maska * maskb):
        indx = np.where(maska * maskb > 0)[0][0]
        YwithTime = abund[:, indx, 3]
        initialT = np.where(YwithTime > 1e-25)[0]
    else:
        YwithTime = np.zeros((AbunLen))

    Yinitial = 0

    FlowWithTime = np.zeros((AbunLen))
    FlowsTime = np.zeros((AbunLen))
    FlowsTemps = np.zeros((AbunLen))
    IntFT = np.zeros((AbunLen))
    IntY = np.zeros((AbunLen))
    dY = np.zeros((AbunLen))

    temps = np.zeros((AbunLen))

    if np.any(initialT):
        for i in trange(AbunLen-1):
            FlowIndex = find_nearest_index(list(keyTimes.values()), times[i])

            flow, timeStamp, tempStamp = procIsotopicFlow(f1, (p, n), FlowIndex)
            
            FlowWithTime[i] = (flow)
            FlowsTime[i] = (timeStamp)
            FlowsTemps[i] = (tempStamp)

            index = find_nearest_index(trueTimes, times[i])
            temps[i] = trueTemps[index]

            if i > 0:
                #IntFT[i] = (IntFT[i-1] + FlowWithTime[i] * (times[i] - times[i-1]))
                IntFT[i] = Yinitial + (IntY[i-1] + (FlowWithTime[i])*(FlowsTime[i]-FlowsTime[i-1]))
            else:
                IntFT[i] = 0

            if i > 0:
                IntY[i] = Yinitial + (IntY[i-1] + (FlowWithTime[i])*(FlowsTime[i]-FlowsTime[i-1])) #
            else:
                IntY[i] = (0)

            if IntY[i] < 1e-25:
                IntY[i] = 0

            if IntFT[i] < 1e-25:
                IntFT[i] = 0


    EvalTime = -10
    percDiff =  round(abs(((YwithTime[EvalTime] - IntFT[EvalTime]) / YwithTime[EvalTime])*100), 0)
    avgDiff = round(np.average(abs((YwithTime-IntFT)*100/(YwithTime))), 0)

    if percDiff > 1000:
        percDiff = ">1000"
    elif np.isfinite(percDiff):
        percDiff = str(int(percDiff))
    else:
        percDiff = str(percDiff)

    if avgDiff > 1000:
        avgDiff = ">1000"
    elif np.isfinite(avgDiff):
        avgDiff = str(int(avgDiff))
    else:
        avgDiff = str(avgDiff)

    YwithTime[-1] = YwithTime[-2]
    FlowWithTime[-1] = FlowWithTime[-2]
    FlowsTime[-1] = FlowsTime[-2]
    IntFT[-1] = IntFT[-2]
    IntY[-1] = IntY[-2]
    dY[-1] = dY[-2]
    times[-1] = times[-2]

    NegFlows = np.where(FlowWithTime < 0, -FlowWithTime, 0)

    plt.plot(temps, YwithTime, label="Abundance")
    plt.plot(FlowsTemps, FlowWithTime, label="Flows")
    plt.plot(FlowsTemps, IntFT, label="Integral Flows")
    plt.plot(FlowsTemps, IntY, label="Experimental Integral")
    plt.plot(FlowsTemps, NegFlows, label="Negative Flows")
    #plt.plot(FlowsTime, -IntFT, label="Negative Integral Flows")
    plt.xscale("linear")
    plt.yscale("log")
    plt.ylim(1e-30, 1)
    plt.grid()
    plt.legend()
    plt.title(str(p+n)+"-"+periodictable.elements[p].symbol + "\n     %Diff @ End = " + percDiff + "%     Avg % Diff = " + avgDiff + "%")
    plt.savefig(("/home/thana1dr/Research_Results&Data/Week9/Automatic/FlowFitting/LargerH5Files/1Traj{:04d}.png".format(j)), dpi=200)
    #plt.show()
    plt.close()
    j += 1
"""
