import periodictable

path = "/home/thana1dr/WinNet/data/Research_Data/vp-proc/"
filename = "track_nuclei"

start = 1
stop = 1000
def a2bE(start, stop):
    totrack = []
    for element in periodictable.elements:
        if element.mass >= start and element.mass <= stop:
            if int(element.mass) == 1 and element.symbol == "n":
                totrack.append(("n" + "\n").lower())
            elif int(element.mass) == 1 and element.symbol == "H":
                totrack.append(("p" + "\n").lower())
            elif int(element.mass) == 2 and element.symbol == "H":
                totrack.append(("d" + "\n").lower())
            elif int(element.mass) == 3 and element.symbol == "H":
                totrack.append(("t" + "\n").lower())
            else:
                totrack.append((str(element.symbol) + str(int(element.mass)) + "\n").lower())
    return totrack

def a2bI(start, stop):
    totrack = []
    for element in periodictable.elements:
        for iso in element:
            if iso.mass >= start and iso.mass <= stop:
                if element.mass >= start and element.mass <= stop:
                    if int(iso.mass) == 1 and iso.symbol == "n":
                        totrack.append((str(iso.symbol) + "\n").lower())
                    elif int(iso.mass) == 1 and iso.symbol == "H":
                        totrack.append(("p" + "\n").lower())
                    elif int(iso.mass) == 2 and iso.symbol == "D":
                        totrack.append((str(iso.symbol) + "\n").lower())
                    elif int(iso.mass) == 3 and iso.symbol == "T":
                        totrack.append((str(iso.symbol) + "\n").lower())
                    else:
                        totrack.append((str(iso.symbol) + str(int(iso.mass)) + "\n").lower())
    return totrack

def fromNumberList(list):
    totrack = []
    for each in list:
        element = periodictable.elements[each]
        for iso in element:
            if iso.mass >= start and iso.mass <= stop:
                if element.mass >= start and element.mass <= stop:
                    if int(iso.mass) == 1 and iso.symbol == "n":
                        totrack.append((str(iso.symbol) + "\n").lower())
                    elif int(iso.mass) == 1 and iso.symbol == "H":
                        totrack.append(("p" + "\n").lower())
                    elif int(iso.mass) == 2 and iso.symbol == "D":
                        totrack.append((str(iso.symbol) + "\n").lower())
                    elif int(iso.mass) == 3 and iso.symbol == "T":
                        totrack.append((str(iso.symbol) + "\n").lower())
                    else:
                        totrack.append((str(iso.symbol) + str(int(iso.mass)) + "\n").lower())
    return totrack

def WriteFile(path, filename, lines):
    with open(path + filename, "w", encoding="utf-8") as f:
        f.writelines(lines)

#totrack = a2bE(start, stop)
totrack = a2bI(start, stop)
#totrack = fromNumberList(sorted([38,39,40]))

print(totrack)
WriteFile(path, filename, totrack)

