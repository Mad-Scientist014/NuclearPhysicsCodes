from beautifultable import BeautifulTable
import numpy as np
import os
import re
import math
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl

directory = "/home/thana1dr/Downloads/Jacobi_nup/"
pathout = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"

plots =[3,3]

# Define the desired initial values
# To ignore value set = True
min_t = 1e-2
max_t = 2e-1

min_Entropy = 60
max_Entropy = 100

min_T = 1.9e10
max_T = 2.1e10

min_r = 2.1e6
max_r = 3.e6

min_ye = 0.55
max_ye = 0.65

min_p = 3e7
max_p= 1e8

min_Lnue = 4e51
max_Lnue = 6e51

min_Lanue = 2.9e51
max_Lanue = 3.1e51

min_Tnue = 8e0
max_Tnue = 9e0

min_Tanue = 5.2e0
max_Tanue = 5.3e0


# Specify the directory where your .dat files are located
directory = "/home/thana1dr/Downloads/Jacobi_nup/"

def chooseTraj():
    opts = []
    fails = {"time":0,"entr":0,"temp":0,"dens":0,"rad":0,"ye":0,"lnue":0,"lanue":0,"tnue":0,"tanue":0,}
    avgFail = {"time":0,"entr":0,"temp":0,"dens":0,"rad":0,"ye":0,"lnue":0,"lanue":0,"tnue":0,"tanue":0,}
    minFail = {"time":math.inf,"entr":math.inf,"temp":math.inf,"dens":math.inf,"rad":math.inf,"ye":math.inf,"lnue":math.inf,"lanue":math.inf,"tnue":math.inf,"tanue":math.inf,}
    maxFail = {"time":0,"entr":0,"temp":0,"dens":0,"rad":0,"ye":0,"lnue":0,"lanue":0,"tnue":0,"tanue":0,}
    # Iterate through the .dat files in the directory and write matching data to the output file

    for filename in os.listdir(directory):
        if filename.endswith(".dat"):
            file_path = os.path.join(directory, filename)

            entropy = None
            ye = None
            r_value = None

            with (open(file_path, "r") as file):
                for i, line in enumerate(file):
                    if i == 0:
                        # Extract the Entropy value using regular expressions
                        match = re.search(r'S =\s*([0-9.]+[Ee][-+]?[0-9]+)', line)
                        if match:
                            entropy = round(float(match.group(1)), 2)
                    elif i == 4:
                        data_line = line.strip().split()
                        if len(data_line) >= 5:
                            t = round(float(data_line[0]), 4)
                            temp = round(float(data_line[1]), 4)
                            p = round(float(data_line[2]), 4)
                            r_value = round(float(data_line[3]), 2)
                            ye = float(data_line[4])
                            lnue = round(float(data_line[5]), 4)
                            lanue = round(float(data_line[6]), 4)
                            tnue = round(float(data_line[7]), 4)
                            tanue = round(float(data_line[8]), 4)
                            if ye >= min_ye and ye <= max_ye:
                                if r_value >= min_r and r_value <= max_r:
                                    if p>=min_p and p<=max_p:
                                        if lnue >= min_Lnue and lnue <= max_Lnue:
                                            if lanue >= min_Lanue and lanue <= max_Lanue:
                                                if t >= min_t and t <= max_t:
                                                    if tnue >= min_Tnue and tnue <= max_Tnue:
                                                        if tanue >= min_Tanue and tanue <= max_Tanue:
                                                            if temp >= min_T and temp <= max_T:
                                                                if entropy >= min_Entropy and entropy <= max_Entropy:
                                                                    opts.append([filename, entropy, t, temp, p, r_value, ye, lnue, lanue, tnue, tanue])
                                                                else:
                                                                    fails["entr"] += 1
                                                                    avgFail["entr"] = ((fails["entr"] - 1) * avgFail["entr"] + entropy) / fails["entr"]
                                                                    if entropy > maxFail["entr"]: maxFail["entr"] = entropy
                                                                    if entropy < minFail["entr"]: minFail["entr"] = entropy
                                                            else:
                                                                fails["temp"] += 1
                                                                avgFail["temp"] = ((fails["temp"] - 1) * avgFail["temp"] + temp) / fails["temp"]
                                                                if temp > maxFail["temp"]: maxFail["temp"] = temp
                                                                if temp < minFail["temp"]: minFail["temp"] = temp
                                                        else:
                                                            fails["tanue"] += 1
                                                            avgFail["tanue"] = ((fails["tanue"] - 1) * avgFail["tanue"] + tanue) / fails["tanue"]
                                                            if tanue > maxFail["tanue"]: maxFail["tanue"] = tanue
                                                            if tanue < minFail["tanue"]: minFail["tanue"] = tanue

                                                    else:
                                                        fails["tnue"] += 1
                                                        avgFail["tnue"] = ((fails["tnue"] - 1) * avgFail["tnue"] + tnue) / fails["tnue"]
                                                        if tnue > maxFail["tnue"]: maxFail["tnue"] = tnue
                                                        if tnue < minFail["tnue"]: minFail["tnue"] = tnue
                                                else:
                                                    fails["time"] += 1
                                                    avgFail["time"] = ((fails["time"] - 1) * avgFail["time"] + t) / fails["time"]
                                                    if t > maxFail["time"]: maxFail["time"] = t
                                                    if t < minFail["time"]: minFail["time"] = t
                                            else:
                                                fails["lanue"] += 1
                                                avgFail["lanue"] = ((fails["lanue"] - 1) * avgFail["lanue"] + lanue) / fails["lanue"]
                                                if lanue > maxFail["lanue"]: maxFail["lanue"] = lanue
                                                if lanue < minFail["lanue"]: minFail["lanue"] = lanue
                                        else:
                                            fails["lnue"] += 1
                                            avgFail["lnue"] = ((fails["lnue"] - 1) * avgFail["lnue"] + lnue) / fails["lnue"]
                                            if lnue > maxFail["lnue"]: maxFail["lnue"] = lnue
                                            if lnue < minFail["lnue"]: minFail["lnue"] = lnue
                                    else:
                                        fails["dens"] += 1
                                        avgFail["dens"] = ((fails["dens"] - 1) * avgFail["dens"] + p) / fails["dens"]
                                        if p > maxFail["dens"]: maxFail["dens"] = p
                                        if p < minFail["dens"]: minFail["dens"] = p
                                else:
                                    fails["rad"] += 1
                                    avgFail["rad"] = ((fails["rad"] - 1) * avgFail["rad"] + r_value) / fails["rad"]
                                    if r_value > maxFail["rad"]: maxFail["rad"] = r_value
                                    if r_value < minFail["rad"]: minFail["rad"] = r_value
                            else:
                                fails["ye"] += 1
                                avgFail["ye"] = ((fails["ye"] - 1) * avgFail["ye"] + ye) / fails["ye"]
                                if ye > maxFail["ye"]: maxFail["ye"] = ye
                                if ye < minFail["ye"]: minFail["ye"] = ye


    return opts, fails, avgFail, minFail, maxFail

def func(x,a,b):
    return a/(np.power(x, b))

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

def rowConvert(data, row, scale):
    tmp = []
    for each in data[:,row]:
        tmp.append(each * (10**scale))
    data[:, row] = tmp
    return data

def WriteFile(path, filename, data, tags):
    os.makedirs(os.path.dirname(path+filename), exist_ok=True)
    with open(path + filename, "w", encoding="utf-8") as f:
        f.writelines(tags)
        for each in data:
            for item in each:
                f.write(str(item) + " ")
            f.write("\n")

def initGraph(dim):
    if dim[0] != 1 or dim[1] != 1:
        fig, axis = plt.subplots(dim[0],dim[1], layout="constrained", figsize=(20,10))
    else:
        fig, axis = plt.subplots(dim[0], dim[1], layout="constrained")
    return fig, axis

table = BeautifulTable()
trajectories, fails, avg, min, max = chooseTraj()
sucAvg = [0 for i in range(11)]
for each in trajectories:
    table.rows.append(each)
    for i in range(1, len(each)):
        sucAvg[i] += each[i]
if len(trajectories) != 0 : sucAvg = [each/len(trajectories) for each in sucAvg]
table.columns.header = ["Filename", "Entropy", "Time", "Temp", "Dens", "Rad", "Ye", "Lnue", "Lanue", "Tnue", "Tanue"]
sucAvg[0] = "Avg"
table.rows.append(sucAvg)
print(table)
print("# Trajcetories = " + str(len(trajectories)))
failTable = BeautifulTable()
failTable.columns.header = ["Time", "Entropy","Temp", "Dens", "Rad", "Ye", "Lnue", "Lanue", "Tnue", "Tanue"]
failTable.rows.append([fails["time"], fails["entr"], fails["temp"], fails["dens"], fails["rad"], fails["ye"], fails["lnue"], fails["lanue"], fails["tnue"], fails["tanue"]])
failTable.rows.append([avg["time"], avg["entr"], avg["temp"], avg["dens"], avg["rad"], avg["ye"], avg["lnue"], avg["lanue"], avg["tnue"], avg["tanue"]])
failTable.rows.append([min["time"], min["entr"], min["temp"], min["dens"], min["rad"], min["ye"], min["lnue"], min["lanue"], min["tnue"], min["tanue"]])
failTable.rows.append([max["time"], max["entr"], max["temp"], max["dens"], max["rad"], max["ye"], max["lnue"], max["lanue"], max["tnue"], max["tanue"]])
failTable.rows.header = ["Fails", "AvgFail", "MinFail", "MaxFail"]
print(failTable)



if len(trajectories) != 0:
    if input("Graph Trajectories: (Yes/no)").lower()[0] == "y":
        fig, axis = initGraph(plots)
        trajData = []
        for traj in trajectories:
            data, tags = loadFile(directory,traj[0])
            data = rowConvert(data, 1, -9)
            data = rowConvert(data, 3, -5)
            trajData.append([traj[0], data])
        for traj in trajData:
            axis[0][0].plot(traj[1][:, 0][:], traj[1][:, 1][:], label=str(traj[0]))
            axis[0][1].plot(traj[1][:, 1][:], traj[1][:, 2][:], label=str(traj[0]))
            axis[0][2].plot(traj[1][:, 0][:], traj[1][:, 2][:], label=str(traj[0]))
            axis[1][0].plot(traj[1][:, 0][:], traj[1][:, 3][:], label=str(traj[0]))

            axis[1][1].plot(traj[1][:, 3][:] ** 3, traj[1][:, 1][:] / traj[1][:, 3][:] ** 3)
            popt, pcov = curve_fit(func, traj[1][:, 3][:] ** 3, traj[1][:, 1][:] / traj[1][:, 3][:] ** 3)
            ss_res = np.sum((traj[1][:, 1][:] / traj[1][:, 3][:] ** 3 - func(traj[1][:, 3][:] ** 3, *popt)) ** 2)
            ss_tot = np.sum((traj[1][:, 1][:] / traj[1][:, 3][:] ** 3 - np.mean(traj[1][:, 1][:] / traj[1][:, 3][:] ** 3)) ** 2)
            r2 = 1 - (ss_res / ss_tot)
            axis[1][1].plot(traj[1][:, 3][:] ** 3, func(traj[1][:, 3][:] ** 3, *popt), linestyle='dashed',label='fit: R=%5.3f, a=%5.3f, b=%5.3f' % (r2, popt[0], popt[1]))

            axis[1][2].plot(traj[1][:, 3][:], traj[1][:, 2][:], label=str(traj[0]))

            axis[2][0].plot(traj[1][:, 0][:], traj[1][:, 3][:] ** 3*traj[1][:, 2][:], label=str(traj[0]))
            axis[2][1].plot(traj[1][:, 1][:] / traj[1][:, 3][:] ** 3, traj[1][:, 2][:], label=str(traj[0]))

            axis[2][2].plot(traj[1][:, 0][:], traj[1][:, 4][:], label=str(traj[0]))


            #axis[2][0].plot(traj[1][:, 0][:], traj[1][:, 1][:], label=str(traj[0]))
            #axis[2][1].plot(traj[1][:, 0][:], traj[1][:, 1][:], label=str(traj[0]))



        axis[0][0].set_xlabel("Time")
        axis[0][0].set_ylabel("Temperature")
        axis[0][0].set_xscale("log")
        axis[0][0].set_yscale("linear")
        axis[0][0].legend()
        axis[0][0].grid()

        axis[0][1].set_xlabel("Temp")
        axis[0][1].set_ylabel("Density")
        axis[0][1].set_xscale("linear")
        axis[0][1].set_yscale("linear")
        axis[0][1].legend()
        axis[0][1].grid()

        axis[0][2].set_xlabel("Time")
        axis[0][2].set_ylabel("Density")
        axis[0][2].set_xscale("log")
        axis[0][2].set_yscale("log")
        axis[0][2].legend()
        axis[0][2].grid()

        axis[1][0].set_xlabel("Time")
        axis[1][0].set_ylabel("Radius")
        axis[1][0].set_xscale("linear")
        axis[1][0].set_yscale("linear")
        axis[1][0].legend()
        axis[1][0].grid()

        axis[1][1].set_xlabel("Volume")
        axis[1][1].set_ylabel("Pressure")
        axis[1][1].set_xscale("log")
        axis[1][1].set_yscale("log")
        axis[1][1].legend(loc="upper right")
        axis[1][1].grid()

        axis[1][2].set_xlabel("Radius")
        axis[1][2].set_ylabel("Density")
        axis[1][2].set_xscale("log")
        axis[1][2].set_yscale("log")
        axis[1][2].legend()
        axis[1][2].grid()

        axis[2][0].set_xlabel("Time")
        axis[2][0].set_ylabel("Mass")
        axis[2][0].set_xscale("linear")
        axis[2][0].set_yscale("log")
        axis[2][0].legend()
        axis[2][0].grid()

        axis[2][1].set_xlabel("Pressure")
        axis[2][1].set_ylabel("Density")
        axis[2][1].set_xscale("log")
        axis[2][1].set_yscale("log")
        axis[2][1].legend()
        axis[2][1].grid()

        axis[2][2].set_xlabel("Time")
        axis[2][2].set_ylabel("Ye")
        axis[2][2].set_xscale("linear")
        axis[2][2].set_yscale("linear")

        axis[2][2].set_ylim(0.4, 0.8)
        axis[2][2].grid()
        axis[2][2].legend()

        mpl.rcParams['agg.path.chunksize'] = 10000
        plt.show()

    if input("Add Trajectories: (Yes/no)").lower()[0] == "y":
        showChanges = input("Show files added / removed? (Yes/no)").lower()[0] == "y"
        if input("Remove Old Trajectories: (Yes/no)").lower()[0] == "y":
            for filename in os.listdir(pathout):
                if filename.startswith("traj_") and filename.endswith(".dat"):
                    os.remove(pathout+filename)
                    if showChanges:
                        print("Deleting: " + pathout + filename)
        for each in trajectories:
            data, tags = loadFile(directory,each[0])
            data = rowConvert(data, 1, -9)
            data = rowConvert(data, 3, -5)
            WriteFile(pathout,each[0],data, tags)
            if showChanges:
                print("Adding: " + pathout + each[0])