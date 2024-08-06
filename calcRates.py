from matplotlib import pyplot as plt
from itertools import islice, batched
import numpy as np
import math
from tqdm import tqdm
from tqdm import trange
import subprocess
import threading
import time
import h5py
import os

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

def findColumn(data, item):
    if item in data:
        return data.index(item)
    else:
        return -1

def initGraph(dim):
    fig, axis = plt.subplots(dim[0],dim[1], layout="constrained")
    return fig, axis

def listToString(inputList):
    out = ""
    for each in inputList:
        out += str(each) + " "
    return out

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

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


def FindThreadReactions(dataset, keys, fileout, statusout):
    with open(path + fileout, 'w', encoding="utf-8") as f:
        lines = set()
        for key in keys:
            pin = dataset['flows'][str(key)]['p_in'][:]
            nin = dataset['flows'][str(key)]['n_in'][:]
            pout = dataset['flows'][str(key)]['p_out'][:]
            nout = dataset['flows'][str(key)]['n_out'][:]
            reacs = {str(pin[i]) + " " + str(nin[i]) + " " + str(pout[i]) + " " + str(nout[i]) for i in range(pin.shape[0])}
            reacs = reacs - lines

            for each in reacs:
                lines.add(each)
                f.write(each + "\n")

            with open(path + statusout, 'w', encoding="utf-8") as g:
                perc = keys.index(key)/len(keys)*100
                g.write(str(perc))

            f.flush()
    return

def findAllReactionsThreaded(dataset, fileout):
    keys = list(dataset['flows'].keys())
    tmp1 = "_tmp-1.dat"
    tmp2 = "_tmp-2.dat"
    tmp3 = "_tmp-1.dat"
    tmp4 = "_tmp-2.dat"
    numThreads = 10
    tmps = ["_tmp-" + str(i) + ".dat" for i in range(numThreads)]
    statuses = ["Status_" + str(i) + ".dat" for i in range(numThreads)]

    threadKeys = [key for key in chunks(keys, int(len(keys)/numThreads))]

    print("Finding and Saving all Unique Reactions: ")
    threads = [threading.Thread(target=FindThreadReactions, args=(dataset, threadKeys[i], fileout + tmps[i], statuses[i])) for i in range(numThreads)]

    for thread in threads:
        thread.start()

    print("Threads Started")
    print("Waiting .  .  .")

    for thread in threads:
        thread.join()

    print("All Threads Joined")

    total = np.zeros((4))
    for each in tmps:
        total = np.vstack([total, loadFile(path, fileout + each)[0]])
        total = np.unique(total)

    with open(path + fileout, 'w', encoding="utf-8") as f:
        for each in total:
            f.write(each)

    os.remove(path + fileout + tmp1)
    os.remove(path + fileout + tmp2)

    return


def findAllReactions(dataset, fileout):
    Numkeys = len(list(dataset['flows'].keys()))
    with open(path + fileout, 'w', encoding="utf-8") as f:
        lines = set()

        for key in tqdm(range(1, Numkeys+1)):
            pin = dataset['flows'][str(key)]['p_in'][:]
            nin = dataset['flows'][str(key)]['n_in'][:]
            pout = dataset['flows'][str(key)]['p_out'][:]
            nout = dataset['flows'][str(key)]['n_out'][:]
            reacs = {str(pin[i]) + " " + str(nin[i]) + " " + str(pout[i]) + " " + str(nout[i]) for i in range(pin.shape[0])}
            reacs = reacs - lines

            for each in reacs:
                lines.add(each)
                f.write(each + "\n")

            f.flush()
    return

def findUniqueReactions(filein, fileout):
    with open(path+filein, "r") as f:
        with open(path + fileout, "w") as g:
            lines = set()
            t = tqdm()
            print("Processing and Saving Unique Reactions: ")

            pin = f.readline().strip("\n").strip('[').strip(']').split(" ")
            nin = f.readline().strip("\n").strip('[').strip(']').split(" ")
            pout = f.readline().strip("\n").strip('[').strip(']').split(" ")
            nout = f.readline().strip("\n").strip('[').strip(']').split(" ")
            #blank = f.readline()

            while pin != "":
                #print(pin)
                #print(nin)
                #print(pout)
                #print(nout)
                for i in range(len(pin)):
                    each = pin[i] +" "+ nin[i] +" "+ pout[i] +" "+ nout[i]
                    if each not in lines:
                        lines.add(each)
                        g.write(each + "\n")
                t.update(1)

                pin = f.readline()
                while pin.find(']') == -1:
                    pin += f.readline()
                pin = pin.replace('\n', '').strip('[').strip(']').split(" ")
                pin = list(filter(None, pin))

                nin = f.readline()
                while nin.find(']') == -1:
                    nin += f.readline()
                nin = nin.replace('\n', '').strip('[').strip(']').split(" ")
                nin = list(filter(None, nin))

                pout = f.readline()
                while pout.find(']') == -1:
                    pout += f.readline()
                pout = pout.replace('\n', '').strip('[').strip(']').split(" ")
                pout = list(filter(None, pout))

                nout = f.readline()
                while nout.find(']') == -1:
                    nout += f.readline()
                nout = nout.replace('\n', '').strip('[').strip(']').split(" ")
                nout = list(filter(None, nout))

                #blank = f.readline()

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
    keys = keyTimes.keys()



    dts = []
    vals = list(keyTimes.values())
    for i in range(len(keyTimes)):
        if i != 0:
            dts.append(vals[i]-vals[i-1])

    flows = np.zeros((len(reactions)))
    flowKeys = [str(reactions[i, :].astype(int).tolist()) for i in range(reactions.shape[0])]
    flowIndexes = [i for i in range(len(flowKeys))]

    flowDict = {}

    for i in range(len(flowKeys)):
        flowDict[flowKeys[i]] = flowIndexes[i]

    for i in trange(len(list(keyTimes.keys()))):
        time = str(i+1)
        if i+2 <len(list(keyTimes.keys())):
            dt = keyTimes[str(i + 2)] - keyTimes[str(i + 1)]
        else:
            dt = keyTimes[str(i + 1)] - keyTimes[str(i)]

        if time == '1' or not np.array_equal(dataset['flows'][time]['p_in'][:], pin):
            pin = dataset['flows'][time]['p_in'][:]
            nin = dataset['flows'][time]['n_in'][:]
            pout = dataset['flows'][time]['p_out'][:]
            nout = dataset['flows'][time]['n_out'][:]
            fout = dataset['flows'][time]['flow'][:]

        for i in range(pin.shape[0]):
            reac = [pin[i], nin[i], pout[i], nout[i]]

            """if str(reac) not in flowDict.keys():
                maxI = -1
                minI = 0

                if list(keyTimes.keys()).index(time) == 200:
                    print("300 Keys Done")

                rnge = np.where(reactions[minI:, 0] == pin[i])[0]
                maxI = np.max(rnge) + 1
                minI = np.min(rnge)
                rnge = np.where(reactions[minI:maxI, 1] == nin[i])[0]
                maxI = np.max(rnge+minI) + 1
                minI = np.min(rnge+minI)
                rnge = np.where(reactions[minI:maxI, 2] == pout[i])[0]
                maxI = np.max(rnge+minI) + 1
                minI = np.min(rnge+minI)
                index = np.where(reactions[minI:maxI, 3] == nout[i])[0][0] + minI

                flowDict[str(reac)] = index

            else:"""
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

def writeFilePartial(path, filename, lines, func, args):
    with open(path+filename, 'w', encoding="utf-8") as f:
        for line in func(*args):
            f.write(str(line))
            f.write("\n")
    return

def writeRates(path, filename1, filename2, dataset, temps, Alimit):
    with open(path+filename1, 'w', encoding="utf-8") as f:
        with open(path+filename2, 'w', encoding="utf-8") as g:
            for lines in findNCapVcapGen(dataset, temps, Alimit):
                fline = lines[0]
                f.write(str(fline))
                f.write("\n")
                gline = lines[1]
                g.write(str(gline))
                g.write("\n")
    return

def findNCapVcap(dataset, temps):
    avgNCap = []
    avgVCap = []
    for temp in temps:
        print("Temp: " + str(temp) +" GK")
        nCapsum = 0
        nCapcnt = 0
        vCapsum = 0
        vCapcnt = 0
        temp = float(temp)
        keyTemps ={}
        keys = list(dataset['flows'].keys())
        for key in keys:
            keyTemps[key] = dataset['flows'][key]['temp'][0]
        keyTemps = dict(sorted(keyTemps.items(), key=lambda item: item[1]))
        num = list(keyTemps.keys())[list(keyTemps.values()).index(find_nearest(list(keyTemps.values()), temp))]
        for data in procFlowGen(dataset, num):
            pin = int(data[0])
            nin = int(data[1])
            pout = int(data[2])
            nout = int(data[3])
            # if pin >= 16:
            # print("Pin >= 16")
            if (pin == pout and nout-nin == 1):
                nCapsum += data[4]
                nCapcnt += 1
            elif (pin - nin == 1 and nout - nin == 1):
                vCapsum += data[4]
                vCapcnt += 1

        if nCapcnt != 0:
            avgNCap.append(nCapsum/nCapcnt)
        else:
            avgNCap.append(0)
        if vCapcnt != 0:
            avgVCap.append(vCapsum/vCapcnt)
        else:
            avgVCap.append(0)
    return avgNCap, avgVCap

def findNCapVcapGen(dataset, temps, Alimit):
    keyTemps = {}
    keys = list(dataset['flows'].keys())
    for key in keys:
        keyTemps[key] = dataset['flows'][key]['temp'][0]
    keyTemps = dict(sorted(keyTemps.items(), key=lambda item: item[1]))
    for temp in tqdm(temps):
        nCapsum = 0
        nCapcnt = 0
        vCapsum = 0
        vCapcnt = 0
        temp = float(temp)
        num = list(keyTemps.keys())[list(keyTemps.values()).index(find_nearest(list(keyTemps.values()), temp))]
        for data in procFlowGen(dataset, num):
            pin = int(data[0])
            nin = int(data[1])
            pout = int(data[2])
            nout = int(data[3])
            # if pin >= 16:
            # print("Pin >= 16")
            if (pin == pout and nout-nin == -1) or (pout - pin == 2 and nout-nin == 1):
                nCapsum += data[4]
                nCapcnt += 1
            elif (pout == 0 and nout == 1 and pin+nin >= Alimit):
                vCapsum += data[4]
                vCapcnt += 1

        if nCapcnt != 0 and vCapcnt != 0:
            yield (nCapsum/nCapcnt, vCapsum/vCapcnt)
        elif nCapcnt != 0:
            yield (nCapsum/nCapcnt, 0)
        elif vCapcnt != 0:
            yield (0, vCapsum/vCapcnt)
        else:
            yield (0, 0)

print("Starting Calculations: ")
begin = time.time()

findAllReactions(f1, "Uniq_Reac_f1(deta).dat")
reacs = loadFile(path, "Uniq_Reac_f1(deta).dat")[0]
reacs = np.unique(reacs, axis=0)
sumFlowsByTime(f1, reacs, path, "SummedFlowsF1(deta).dat")


findAllReactions(f2, "Uniq_Reac_f2(deta).dat")
reacs = loadFile(path, "Uniq_Reac_f2(deta).dat")[0]
reacs = np.unique(reacs, axis=0)
sumFlowsByTime(f2, reacs, path, "SummedFlowsF2(deta).dat")


findAllReactions(f3, "Uniq_Reac_f3(deta).dat")
reacs = loadFile(path, "Uniq_Reac_f3(deta).dat")[0]
reacs = np.unique(reacs, axis=0)
sumFlowsByTime(f3, reacs, path, "SummedFlowsF3(deta).dat")


findAllReactions(f4, "Uniq_Reac_f4(deta).dat")
reacs = loadFile(path, "Uniq_Reac_f4(deta).dat")[0]
reacs = np.unique(reacs, axis=0)
sumFlowsByTime(f4, reacs, path, "SummedFlowsF4(deta).dat")


findAllReactions(f5, "Uniq_Reac_f5(deta).dat")
reacs = loadFile(path, "Uniq_Reac_f5(deta).dat")[0]
reacs = np.unique(reacs, axis=0)
sumFlowsByTime(f5, reacs, path, "SummedFlowsF5(deta).dat")


"""with open(path + "SummedFlowsf1_3.dat", "w") as f:
    for (mask, key) in MaskReactionGen(f1, reacs):
        data = procFlowSumMasked(f1, key, mask)
        for i in range(data.shape[1]):
            f.write(str(data[0, i]) + " " + str(data[1, i]) + " " + str(data[2, i]) + " " + str(data[3, i]) + " " + str(data[4, i]))
            f.write("\n")
            f.flush()
        #print(data)"""


""""""

"""for (mask, key) in MaskReactionGen(f1, reacs):
    data = procFlowSumMasked(f1, key, mask)
    if len(data) != 0:
        print(data)"""



"""writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NeutronReleasedf1.txt", "SoleNeutronf1.txt", f1, f1['mainout']['temp'][()], 0)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NeutronReleasedf2.txt", "SoleNeutronf2.txt", f2, f2['mainout']['temp'][()], 0)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NeutronReleasedf3.txt", "SoleNeutronf3.txt", f3, f3['mainout']['temp'][()], 0)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NeutronReleasedf4.txt", "SoleNeutronf4.txt", f4, f4['mainout']['temp'][()], 0)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NeutronReleasedf5.txt", "SoleNeutronf5.txt", f5, f5['mainout']['temp'][()], 0)"""


"""writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NcapRate49f1.txt", "VcapRate49f1.txt", f1, f1['mainout']['temp'][()], 50)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NcapRate49f2.txt", "VcapRate49f2.txt", f2, f2['mainout']['temp'][()], 50)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NcapRate49f3.txt", "VcapRate49f3.txt", f3, f3['mainout']['temp'][()], 50)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NcapRate49f4.txt", "VcapRate49f4.txt", f4, f4['mainout']['temp'][()], 50)
writeRates("/home/thana1dr/Research_Results&Data/Week6/Calculations/", "NcapRate49f5.txt", "VcapRate49f5.txt", f5, f5['mainout']['temp'][()], 50)"""


print("Finished.")
print("Runtime: "+ str((time.time()- begin)*1000) +" ms")