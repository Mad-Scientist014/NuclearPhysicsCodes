import numpy as np

def readRate(file):
    Rate = np.zeros((30, 2))
    i = 0
    with open(file, "r") as f:
        for line in f.readlines():
            line = line.strip().replace("\t", "").split(" ")
            Rate[i] = line
            i += 1
    return Rate

def modifyRate(path, rateFile, reaction, change):
    with open(path+rateFile, "r+") as f:
        line = f.readline().strip()
        if line == str(reaction):
            rates = f.readline().strip().split()


Ni56Co56 = "Ni56(n,p)rate.txt"


mod = 2
currentRate = "1.338398e+07 1.521653e+07 1.688445e+07 1.997938e+07 2.292435e+07 2.581363e+07 2.869419e+07 3.159353e+07 3.452954e+07 3.751478e+07 4.055858e+07 5.687031e+07 7.536579e+07 9.638768e+07 1.202003e+08 1.470332e+08 1.770946e+08 2.105752e+08 2.476483e+08 3.331734e+08 4.346602e+08 5.527122e+08 6.874803e+08 8.386064e+08"
inp = currentRate.strip().split()
out = ""
for each in inp:
    out += '%.7E ' % (float(each)*mod)
print(out)

"""Ni56Co56Rate = readRate(Ni56Co56)
print(Ni56Co56Rate)"""