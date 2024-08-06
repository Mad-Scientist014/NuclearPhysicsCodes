from matplotlib import pyplot as plt
import numpy as np
import math
import random

plots =[1,1]
compPlot =[1,1]
path = "/home/thana1dr/WinNet/runs/Example_AGB_Test_cescutti/trajectory_TP/" #TP C13pocket

# 1:iteration 2:time[s], 3:T[GK], 4:rho[g/cm3], 5:Ye
# 6:R[km], 7:Y_n, 8:Y_p, 9:Y_alpha, 10:Y_lights
# 11:Y_heavies 12:<Z> 13:<A> 14:entropy [kB/baryon] (15:Sn [MeV])
def loadFile(path, filename):
    data = []
    tags = []
    with open(path+filename, encoding="utf-8") as f:
        line = f.readline()
        while line[0] == "#":
            for each in line.split(" ")[1:]:
                tags.append(each[2:].replace(",", ""))
            i = 0
            while i < len(tags):
                if tags[i].find("\n") != -1:
                    tags[i].replace("\n", "")
                i += 1
            line = f.readline()
        while line != "":
            line = line.split(" ")
            i = len(line) - 1
            while i >= 0:
                if line[i] == '':
                    del line[i]
                i -= 1
            line[-1] = line[-1][:-1]
            data.append(line)
            line = f.readline()
    return data, tags

def loadNucleiFile(path, filename):
    data = []
    tags = []
    with open(path+filename, encoding="utf-8") as f:
        line = f.readline()
        while line[0] == "#":
            for each in line.split(" ")[1:-1]:
                if each != "":
                    tags.append(each.replace(",", ""))
            i = len(tags) - 1
            while i > 0:
                test = tags[i]
                if test.find('(') != -1 or test.find(')') != -1:
                    del tags[i]

                i -= 1
            line = f.readline()
        while line != "":
            line = line.split(" ")
            i = len(line) - 1
            while i >= 0:
                if line[i] == '':
                    del line[i]
                i -= 1
            line[-1] = line[-1][:-1]
            data.append(line)
            line = f.readline()
    return data, tags

def procData(data):
    out = []
    for each in data:
        row = []
        for num in each:
            row.append(round(float(num),50))
        out.append(row)
    return out

def addData(data, tags, new, name):
    i = 0
    while i < len(data):
        data[i].append(new[i])
    tags.append(name)
    return data, tags

def calcNewRow(data, tags, op, r1, r2, name):
    i = 0
    r1 = sepData(data, r1)
    r2 = sepData(data, r2)
    if op == "*":
        while i < len(data):
            data[i].append(r1[i]*r2[2])
            i += 1
        tags.append(name)
    elif op == "/":
        while i < len(data):
            data[i].append(r1[1]/r2[2])
            i += 1
        tags.append(name)
    return data, tags
def sepData(data, row):
    out = []
    for each in data:
        out.append(each[row])
    return out

def initGraph(dim):
    fig, axis = plt.subplots(dim[0],dim[1])
    return fig, axis

def addGraphManual(fig, axis, n, plots, title, label, data, scale):
    if plots[1] != 1:
        k = n % plots[0]
        j = math.floor(n / plots[1])
        axis[j,k].plot(data[0], data[1])
        axis[j, k].set_title(title)
        axis[j, k].set_xlabel(label[0])
        axis[j, k].set_ylabel(label[1])
        axis[j, k].xscale = scale[0]
        axis[j, k].yscale = scale[1]
    elif(plots[0] != 1):
        j = n % plots[0]
        axis[j].plot(data[0], data[1])
        axis[j].set_title(title)
        axis[j].set_xlabel(label[0])
        axis[j].set_ylabel(label[1])
        axis[j].xscale = scale[0]
        axis[j].yscale = scale[1]
    else:
        axis.plot(data[0], data[1])
        axis.set_title(title)
        axis.set_xlabel(label[0])
        axis.set_ylabel(label[1])
        axis.xscale = scale[0]
        axis.yscale = scale[1]

def addGraphAuto(fig, axis, n, plots, title, source, dataIn, scale, extra = []):
    label = (tags[source[0]],tags[source[1]])
    data = (sepData(dataIn, source[0]),sepData(dataIn, source[1]))
    xticks = 100
    yticks = 100
    i = 0
    hexadecimal_alphabets = '0123456789ABCDEF'

    color = ["#" + ''.join([random.choice(hexadecimal_alphabets) for j in
                            range(6)]) for i in range(len(extra))]

    twins = []
    if plots[1] != 1:
        k = n % plots[0]
        j = math.floor(n / plots[1])
        axis[j,k].plot(data[0], data[1])
        if len(extra) != 0:
            while i < len(extra):
                twins.append(axis[j,k].twinx())
                twins[i].plot(sepData(dataIn, extra[i][0]), sepData(dataIn, extra[i][1]), c=color[i])
                twins[i].set_ylabel(tags[extra[i][1]])
                twins[i].set_xscale(scale[0])
                twins[i].set_yscale("log")
                i += 1
        axis[j, k].set_title(title)
        axis[j, k].set_xlabel(label[0])
        axis[j, k].set_ylabel(label[1])
        axis[j, k].set_xscale(scale[0])
        axis[j, k].set_yscale(scale[1])
    elif(plots[0] != 1):
        j = n % plots[0]
        axis[j].plot(data[0], data[1])
        if len(extra) != 0:
            while i < len(extra):
                twins.append(axis[j].twinx())
                twins[i].plot(sepData(dataIn, extra[i][0]), sepData(dataIn, extra[i][1]), c=color[i])
                twins[i].set_ylabel(tags[extra[i][1]])
                twins[i].set_xscale(scale[0])
                twins[i].set_yscale("log")
                i += 1
        axis[j].set_title(title)
        axis[j].set_xlabel(label[0])
        axis[j].set_ylabel(label[1])
        axis[j].set_xscale(scale[0])
        axis[j].set_yscale(scale[1])
    else:
        axis.plot(data[0], data[1])
        if len(extra) != 0:
            while i < len(extra):
                twins.append(axis.twinx())
                twins[i].plot(sepData(dataIn, extra[i][0]), sepData(dataIn, extra[i][1]), c=color[i])
                twins[i].set_ylabel(tags[extra[i][1]])
                twins[i].set_xscale(scale[0])
                twins[i].set_yscale("linear")
                i += 1
        axis.set_title(title)
        axis.set_xlabel(label[0])
        axis.set_ylabel(label[1])
        axis.set_xscale(scale[0])
        axis.set_yscale(scale[1])

def addCompPlot(fig, axis, title, data, scale):
    i = 1
    while i < len(elem):
        axis.plot(sepData(data, 0), sepData(data, i), label = elem[i])
        i += 1
    axis.set_title(title)
    axis.set_xlabel(elem[0])
    axis.set_ylabel("Abunduncies")
    axis.set_xscale(scale[0])
    axis.set_yscale(scale[1])

data, tags = loadFile(path, "mainout.dat")
data = procData(data)
fig, axis = initGraph(plots)
addGraphAuto(fig,axis,0,plots,"First", (1, 14), data, ("log", "linear"))
addGraphAuto(fig,axis,1,plots,"Second", (1, 3), data, ("log", "linear"))
#addGraphAuto(fig,axis,2,plots,"Third", (2, 3), data, ("linear", "linear"))
#addGraphAuto(fig,axis,3,plots,"Fourth", (1, 6), data, ("log", "linear"))
fig.tight_layout()

abund, elem = loadNucleiFile(path, "tracked_nuclei.dat")
abund = procData(abund)
fig2,axis2 = initGraph(compPlot)
addCompPlot(fig2,axis2,"CompPlot", abund, ("log", "log"))
plt.legend()
plt.show()
