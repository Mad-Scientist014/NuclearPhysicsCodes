from tqdm import trange
from tqdm import tqdm
import numpy as np
import requests

site = "https://reaclib.jinaweb.org/"

outString = "1\n\n2\n\n3\n\n4\n\n5\n\n6\n\n7\n\n8\n\n9\n\n10\n\n11\n\n"
path = "/home/thana1dr/WinNet/data/"

def findRateIndex(reaction):
    if reaction.find(')') - reaction.find('(') <= 2:
        # Reaction is a weak reaction
        reaction = 'weak/'+reaction
    url = site+reaction+"/"
    response = requests.get(url)
    if response.status_code == 200:
        # Proc data to find id
        # Once id is found make request to download table
        print("Reaction Found .  .  .")
        data = response.text
        lookfor = ".php?rateindex="
        lookfor2 = "Rate Details:"
        if lookfor in data:
            index = data.index(lookfor)
            end = index+len(lookfor)+6
            val = data[index+len(lookfor):end]

            index2 = data.index(lookfor2)
            end2 = index2+len(lookfor2) + data[index2+len(lookfor2):].find("</h1>")
            detName = data[index2+len(lookfor2):end2]

            index3 = data.index(url)
            end3 = index3+len(url) + data[index3+len(url):].find("/")
            version = data[index3+len(url):end3]

            index4 = data.index("<th>Q (MeV):</th>")
            end4 = index4+len("<th>Q (MeV):</th>") + data[index4+len("<th>Q (MeV):</th>"):].find("</td>")
            Q = data[index4+len("<th>Q (MeV):</th>")+13:end4-3]

            data = requests.get(site+"/sets.php?rateindex=" + val).text

            if "<td>r</td>" in data:
                res = "r"
            elif "<td>n</td>" in data:
                res = "n"
            elif "<td>w</td>" in data:
                res = "w"

            if "<td>v</td>" in data:
                dir = "v"
            else:
                dir = " "

            print("Rate Index = " + val)
            return detName, version, res, dir, Q, val
    else:
        print("Server Response Error \n Try different Reaction")

def findTabulatedRates(temps, rateIndex):
    url = site + "datapoints.php?rateindex=" + rateIndex + "&ratetemp="
    rates = np.zeros((len(temps)))
    for i in trange(len(temps)):
        each = temps[i]
        tot = url + str(each)
        response = requests.get(tot)
        data = response.text
        lookfor = "<em>Results:</em>"
        lookfor2 = "cm<sup>3</sup> mol<sup>-1</sup> s<sup>-1</sup>"
        if lookfor in data:
            index = data.index(lookfor)
            end = data.index(lookfor2)
            rate = data[index+len(lookfor):end]
            rates[i] = rate
        else:
            print("Rate Not Found")
    return rates

def appendRatestoOutfile(detname, vers, res, dir, Q, rates, outString):
    # Determine Chapter Number
    detname = detname.replace('+', '').split('&rarr;')
    inp = detname[0].split(' ')
    out = detname[1].split(' ')

    i = 0
    while i < len(inp):
        if len(inp[i]) == 0:
            inp.pop(i)
        else:
            i += 1

    i = 0
    while i < len(out):
        if len(out[i]) == 0:
            out.pop(i)
        else:
            i += 1

    if len(inp) == 1 and len(out) == 1:
        chapNum = 1
    elif len(inp) == 1 and len(out) == 2:
        chapNum = 2
    elif len(inp) == 1 and len(out) == 3:
        chapNum = 3
    elif len(inp) == 2 and len(out) == 1:
        chapNum = 4
    elif len(inp) == 2 and len(out) == 2:
        chapNum = 5
    elif len(inp) == 2 and len(out) == 3:
        chapNum = 6
    elif len(inp) == 2 and len(out) == 4:
        chapNum = 7
    elif len(inp) == 3 and len(out) == 1:
        chapNum = 8
    elif len(inp) == 3 and len(out) == 2:
        chapNum = 9
    elif len(inp) == 4 and len(out) == 2:
        chapNum = 10
    else:
        chapNum = 11

    Is = ["     " for i in range(6)]

    for i in range(len(inp)):
        Is[i] = inp[i]

    for j in range(len(out)):
        Is[i+j+1] = out[j]

    for i in range(len(Is)):
        each = Is[i]
        if len(each) < 5:
            diff = 5-len(each)
            Is[i] = " " * diff + each

    addIndex = outString.find("\n"+str(chapNum+1)+"\n")+1
    indRow = " " * 5

    for each in Is:
        indRow += each + " "

    indRow += " " * 8
    indRow += vers+res+dir
    indRow += " " * 3
    if 'E' not in Q:
        Q += 'E+00'
    indRow += Q
    indRow += " " * 10
    indRow += "\n"

    dataRow = ""
    for each in rates:
        if each < 1e-99:
            each = 0
        if each != 0:
            dataRow += "%.3E" % each + " "
        else:
            dataRow += "1.000E-99 "
    dataRow += "\n"
    all = indRow + dataRow
    outString = outString[:addIndex] + all + outString[addIndex:]

    return outString

def cleanOutString(outString):
    save = ""
    for i in range(1, 11):
        if outString.find("\n"+str(i+1)+"\n") - outString.find("\n"+str(i)+"\n") <= 5:
            outString = outString[outString.find("\n"+str(i+1)+"\n")+1:]
        else:
            save += outString[outString.find("\n"+str(i)+"\n")+1:outString.find("\n"+str(i+1)+"\n")+1]
    save = save.strip()
    return save

def writeOutput(filename, data):
    with open(path+filename, 'w') as f:
        for each in data:
            f.write(each)
        f.flush()
    return

tempVals = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 1.5e-1, 2e-1, 2.5e-1,
            3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, 1e0, 1.5e0, 2.0e0,
            2.5e0, 3e0, 3.5e0, 4e0, 5e0, 6e0, 7e0, 8e0, 9e0, 1e1]

detName, vers, res, dir, Q, index = findRateIndex("cu59(p,a)ni56")
rates = findTabulatedRates(tempVals, index)
rates *= 2
outString = appendRatestoOutfile(detName, vers, res, dir, Q, rates, outString)

detName, vers, res, dir, Q, index = findRateIndex("cu59(p,g)zn60")
rates = findTabulatedRates(tempVals, index)
rates *= 2
outString = appendRatestoOutfile(detName, vers, res, dir, Q, rates, outString)

detName, vers, res, dir, Q, index = findRateIndex("ni56(n,p)co56")
rates = findTabulatedRates(tempVals, index)
rates *= 2
outString = appendRatestoOutfile(detName, vers, res, dir, Q, rates, outString)


outString = cleanOutString(outString)
writeOutput("vp_proc_tabulatedRatesVp.dat", outString)

print(outString)
