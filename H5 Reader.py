from matplotlib import pyplot as plt
import periodictable
import numpy as np
import math
import random
import h5py

plots =[1,1]
path = "/home/thana1dr/WinNet/runs/vp-proc"
filename = "WinNet_data.h5"
tags = ['abar', 'dens', 'entr', 'iteration', 'rad', 'temp', 'time', 'ya', 'ye', 'yheavy', 'ylight', 'yn', 'yp', 'zbar']
f = h5py.File(path + filename, 'r')

def procNuclei(nuclei):
    i = 0
    A = nuclei['A']
    Z = nuclei['Z']
    elements = [periodictable.elements[x].symbol for x in Z]
    while i < len(A):
        elements[i] = elements[i] + str(A[i])
        i += 1
    print(elements)
    return elements

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
def addGraph(fig, axis, n, plots, title, source, dataIn, color=[]):
    p = []
    i = 1

    if color == []:
        hexadecimal_alphabets = '0123456789ABCDEF'
        color = ["#" + ''.join([random.choice(hexadecimal_alphabets) for j in range(6)]) for i in range(len(source))]

    twins = []
    if plots[1] != 1:
        k = n % plots[0]
        j = math.floor(n / plots[1])
        axis[j, k].yaxis.set_tick_params(labelleft=False)
        axis[j, k].set_yticks([])
        while i < len(source):
            twins.append(axis[j].twinx())
            p.append(
                twins[i-1].plot(dataIn[source[0][0]][:], dataIn[source[i][0]][:], c=color[i-1], label=source[i][0])[
                    0])
            p[-1].set_label(source[i][0])
            if i > math.floor(len(source) / 2):
                twins[i - 1].spines["right"].set_position(("axes", .08 * (i - (math.floor(len(source) / 2) + 1)) + 1))
                twins[i - 1].spines["right"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('right')
                twins[i - 1].yaxis.set_ticks_position('right')
            else:
                twins[i - 1].spines["left"].set_position(("axes", -.08 * (i - 1)))
                twins[i - 1].spines["left"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('left')
                twins[i - 1].yaxis.set_ticks_position('left')
            twins[i - 1].set_zorder(i)
            twins[i - 1].set_ylabel(source[i])
            twins[i - 1].set_xscale(source[0][1])
            twins[i - 1].set_yscale(source[i][1])
            twins[i - 1].yaxis.label.set_color(p[i - 1].get_color())

    elif(plots[0] != 1):
        j = n % plots[0]
        axis[j].yaxis.set_tick_params(labelleft=False)
        axis[j].set_yticks([])
        while i < len(source):
            twins.append(axis[j].twinx())
            p.append(twins[i-1].plot(dataIn[source[0][0]][:], dataIn[source[i][0]][:], c=color[i-1], label=source[i][0])[0])
            p[-1].set_label(source[i][0])
            if i > math.floor(len(source)/2):
                twins[i-1].spines["right"].set_position(("axes", .08*(i-(math.floor(len(source)/2)+ 1)) + 1))
                twins[i-1].spines["right"].set_visible(True)
                twins[i-1].yaxis.set_label_position('right')
                twins[i-1].yaxis.set_ticks_position('right')
            else:
                twins[i-1].spines["left"].set_position(("axes", -.08 * (i - 1)))
                twins[i-1].spines["left"].set_visible(True)
                twins[i-1].yaxis.set_label_position('left')
                twins[i-1].yaxis.set_ticks_position('left')
            twins[i-1].set_zorder(i)
            twins[i-1].set_ylabel(source[i])
            twins[i-1].set_xscale(source[0][1])
            twins[i-1].set_yscale(source[i][1])
            twins[i-1].yaxis.label.set_color(p[i-1].get_color())
    else:
        axis.yaxis.set_tick_params(labelleft=False)
        axis.set_yticks([])
        while i < len(source):
            twins.append(axis.twinx())
            p.append(twins[i-1].plot(dataIn[source[0][0]][:], dataIn[source[i][0]][:], c=color[i-1], label=source[i][0])[0])
            p[-1].set_label(source[i][0])
            if i > math.floor(len(source)/2):
                twins[i-1].spines["right"].set_position(("axes", .08*(i-(math.floor(len(source)/2)+ 1)) + 1))
                twins[i-1].spines["right"].set_visible(True)
                twins[i-1].yaxis.set_label_position('right')
                twins[i-1].yaxis.set_ticks_position('right')
            else:
                twins[i-1].spines["left"].set_position(("axes", -.08 * (i - 1)))
                twins[i-1].spines["left"].set_visible(True)
                twins[i-1].yaxis.set_label_position('left')
                twins[i-1].yaxis.set_ticks_position('left')
            twins[i-1].set_zorder(i)
            twins[i-1].set_ylabel(source[i])
            twins[i-1].set_xscale(source[0][1])
            twins[i-1].set_yscale(source[i][1])
            twins[i-1].yaxis.label.set_color(p[i-1].get_color())


            i += 1
        axis.set_title(title)
        axis.set_xlabel(source[0][0])

    plt.legend()

def addGraphMan(elements, axis, n, plots, title, source, dataIn, color=[]):
    p = []
    i = 1

    if color == []:
        hexadecimal_alphabets = '0123456789ABCDEF'
        color = ["#" + ''.join([random.choice(hexadecimal_alphabets) for j in range(6)]) for i in range(len(source))]

    twins = []
    if plots[1] != 1:
        k = n % plots[0]
        j = math.floor(n / plots[1])
        axis[j, k].yaxis.set_tick_params(labelleft=False)
        axis[j, k].set_yticks([])
        while i < len(source):
            twins.append(axis[j].twinx())
            p.append(
                twins[i - 1].plot(dataIn[source[0][0]][:], dataIn[source[i][0]][:], c=color[i - 1], label=source[i][0])[
                    0])
            p[-1].set_label(source[i][0])
            if i > math.floor(len(source) / 2):
                twins[i - 1].spines["right"].set_position(("axes", .08 * (i - (math.floor(len(source) / 2) + 1)) + 1))
                twins[i - 1].spines["right"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('right')
                twins[i - 1].yaxis.set_ticks_position('right')
            else:
                twins[i - 1].spines["left"].set_position(("axes", -.08 * (i - 1)))
                twins[i - 1].spines["left"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('left')
                twins[i - 1].yaxis.set_ticks_position('left')
            twins[i - 1].set_zorder(i)
            twins[i - 1].set_ylabel(source[i])
            twins[i - 1].set_xscale(source[0][1])
            twins[i - 1].set_yscale(source[i][1])
            twins[i - 1].yaxis.label.set_color(p[i - 1].get_color())

    elif (plots[0] != 1):
        j = n % plots[0]
        axis[j].yaxis.set_tick_params(labelleft=False)
        axis[j].set_yticks([])
        while i < len(source):
            twins.append(axis[j].twinx())
            p.append(
                twins[i - 1].plot(dataIn[source[0][0]][:], dataIn[source[i][0]][:], c=color[i - 1], label=source[i][0])[
                    0])
            p[-1].set_label(source[i][0])
            if i > math.floor(len(source) / 2):
                twins[i - 1].spines["right"].set_position(("axes", .08 * (i - (math.floor(len(source) / 2) + 1)) + 1))
                twins[i - 1].spines["right"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('right')
                twins[i - 1].yaxis.set_ticks_position('right')
            else:
                twins[i - 1].spines["left"].set_position(("axes", -.08 * (i - 1)))
                twins[i - 1].spines["left"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('left')
                twins[i - 1].yaxis.set_ticks_position('left')
            twins[i - 1].set_zorder(i)
            twins[i - 1].set_ylabel(source[i])
            twins[i - 1].set_xscale(source[0][1])
            twins[i - 1].set_yscale(source[i][1])
            twins[i - 1].yaxis.label.set_color(p[i - 1].get_color())
    else:
        axis.yaxis.set_tick_params(labelleft=False)
        axis.set_yticks([])
        while i < len(source):
            twins.append(axis.twinx())
            if source[i][0].isnumeric():
                p.append(
                    twins[i - 1].plot(dataIn[source[0][0]][:], dataIn['Y'][:, int(source[i][0])], c=color[i - 1],
                                      label=elements[int(source[i][0])])[
                        0])
                twins[i - 1].set_ylabel(elements[int(source[i][0])])
            else:
                p.append(
                    twins[i - 1].plot(dataIn[source[0][0]][:], dataIn[source[i][0]][:], c=color[i - 1],
                                      label=source[i][0])[
                        0])
                twins[i - 1].set_ylabel(source[i])
            p[-1].set_label(source[i][0])
            if i > math.floor(len(source) / 2):
                twins[i - 1].spines["right"].set_position(("axes", .08 * (i - (math.floor(len(source) / 2) + 1)) + 1))
                twins[i - 1].spines["right"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('right')
                twins[i - 1].yaxis.set_ticks_position('right')
            else:
                twins[i - 1].spines["left"].set_position(("axes", -.08 * (i - 1)))
                twins[i - 1].spines["left"].set_visible(True)
                twins[i - 1].yaxis.set_label_position('left')
                twins[i - 1].yaxis.set_ticks_position('left')
            twins[i - 1].set_zorder(i)
            twins[i - 1].set_xscale(source[0][1])
            twins[i - 1].set_yscale(source[i][1])
            twins[i - 1].yaxis.label.set_color(p[i - 1].get_color())

            i += 1
        axis.set_title(title)
        axis.set_xlabel(source[0][0])
    return


grp = f['mainout']
nuc = f['tracked_nuclei']
fig, axis = initGraph(plots)
elements = procNuclei(nuc)
addGraphMan(elements,axis,0,plots,"Time Evolution of AGB Star",
             (('time', 'linear'),
              (str(findColumn(elements, 'Th232')), 'linear'),
    (str(findColumn(elements, 'U235')), 'linear'),
              ),nuc)

plt.show()

"""
fig, axis = initGraph(plots)
addGraph(fig,axis,0,plots,"Time Evolution of AGB Star",
             (('time', 'log'), ('temp', 'linear'), ('yp', 'linear'),('yn', 'linear'), ('ye', 'linear')), grp,
             color=('red', 'blue', 'green', 'black'))
fig.get_layout_engine().set(w_pad=4 / 72, h_pad=4 / 72, hspace=0,
                            wspace=0)

plt.show()
"""