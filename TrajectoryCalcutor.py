import periodictable
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import mendeleev
import numpy as np
from numpy import pi
import h5py
import re
from scipy.optimize import curve_fit
from periodictable import elements
from matplotlib.ticker import AutoMinorLocator
from sklearn.preprocessing import MinMaxScaler
from numpy import log as ln
import scipy
from scipy import optimize
from scipy.misc import derivative
from scipy.constants import N_A
#from decimal import *
import math
from bisect import bisect
from datetime import datetime
from tqdm import tqdm
from sys import getsizeof
import time
import csv
import gc

kb = 1.38e-23  # J/K
h = 1.054e-34  # J/s
c = 299792458  # m/s
G = 6.675e-11  # N*m^2/kg^2

lp = np.sqrt(h*G/c**3)
mp = np.sqrt(h*c/G)
tp = np.sqrt(h*G/c**5)
Tp = np.sqrt(h*c**5/(G*kb**2))
Ep = np.sqrt(h*c**5/G)

SolarMassKg = 1.989e30

Ktomev = .086
mevtoK = 1/Ktomev
MeVtoK = 1e9*mevtoK
KtoMeV = 1/MeVtoK
MeVtoGK = MeVtoK/1e9
MeVtoJ = 1/1.60218e-13
mevtoJ = 1/6.2415e21
JtoMeV = 1/MeVtoJ


# ////////////////////////////////////
# Input Paramaters
# ////////////////////////////////////

Mns = 1.4  # Solar Masses
Rns = 10 # * 1e5  # Km
Mout = 5.25e-6*SolarMassKg/(mp/tp)
termRadius = 100

# ////////////////////////////////////

stepNum = 0
Ti = 0
Ni = 0
SimTime = 0
uValues = []
rValues = []
qValues = []
TValues = []
rhoValues = []
q1Val = []
q2Val = []
q3Val = []
q4Val = []
q5Val = []

Rnsm = Rns * 1e3 # * 1e5  # Km
Rnscm = Rnsm * 1e2
Rv6 = Rnscm / 1e6  # 1e6 cm
Ye = 0.5

rhoi = 1e10 # g/cm^3
Lve51 = 1
Lave51 = 1
Lvu51 = 1
Eve = 12
Eave = 22
Evu = 34
Ev = Eav = 34

# USE PLANCK UNITS BELOW
Mp = 1.672e-27 / mp  # Mass in Planck Units
Mn = 1.674e-27 / mp  # Mass in Planck Units
mN = (Mp + Mn)/2

Mnsp = Mns * SolarMassKg / mp
Rnsp = Rns * 1e3 / lp

path = "/home/thana1dr/Research_Results&Data/Week8/Automatic/TrajCalc/First/"

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

def MOutRate(r, rhob, u):
    return 4*np.pi*r**2*(rhob)*u

def uByMout(r, rhob):
    return (Mout/(4*np.pi*(r**2)*(rhob)))

def rhobByMout(r, u):
    return Mout/(4*np.pi*r**2*u)

def dRhodrNum(r, u):
    func = lambda r: rhobByMout(Mout, r, u)
    derv = derivative(func, r, dx=1e-12, order=37)
    """dervs = [derivative(func, r, dx=1e-12, order=2 * i + 1) for i in range(1, 70)]
    dervs = np.array(dervs)"""
    return derv

def dudr(u, rhotot, P, dPdr, r):
    first = (1/u)
    second = (1/(rhotot + P))
    third = (1+u**2-2*Mnsp/r)
    last = (Mnsp/r**2)
    tot = first*(second*third*dPdr-last)
    return tot

def dudrsolved(r, u):
    if u < 0 or u == np.nan:
        return np.inf

    placement = bisect(rValues, r)
    rhob = rhobByMout(r, u)
    if placement != len(rValues):
        while len(rValues) > placement and len(rValues) > 0:
            rValues.pop()
            uValues.pop()
            rhoValues.pop()

    rValues.append(r)
    uValues.append(u[0])
    rhoValues.append(rhob[0])


    T = TGasLa(r, u)
    TMeV = T * Tp / MeVtoK
    E = e(T, rhob)
    Pres = P(T, rhob)
    rhoT = rhotot(rhob, E)
    if r == Rnsp:
        Q = 0
        qValues.append(Q)
        TValues.append(Ti)
        q1Val.append(0)
        q2Val.append(0)
        q3Val.append(0)
        q4Val.append(0)
        q5Val.append(0)
    else:
        qr = r * lp * 1e2
        qrho = rhob[0] * (mp / lp ** 3) / 1e3
        #print(qrho)
        qT = TMeV[0]

        Q = qSum(qr, qT, qrho) * MeVtoK / (Tp / (tp * mp))


        if placement != len(qValues):

            while len(qValues) > len(rValues) and len(qValues) > 0:
                qValues.pop()
                TValues.pop()
                q1Val.pop()
                q2Val.pop()
                q3Val.pop()
                q4Val.pop()
                q5Val.pop()

        qValues.append(Q / MeVtoK * ((Tp / (tp * mp))))
        TValues.append(TMeV[0])
        q1Val.append(q1(qr))
        q2Val.append(q2(qT))
        q3Val.append(q3(qT, qr, qrho))
        q4Val.append(q4(qT, qrho))
        q5Val.append(q5(qr, qrho))



    u2rhot = (1+u**2-2*Mnsp/r)/(rhoT + Pres)
    dpdt = (11*np.pi**2*T**3)/45 + rhob/mN
    dPdrho = T/mN
    dedrhob = -(11*pi**2*T**4)/(60*rhob**2)
    dedT = (11*np.pi**2*T**3)/(15*rhob) + 3/(2*mN)

    first = dpdt*(u2rhot)*(Q/u)
    second = 2*Mout/(4*pi*u*r**3)
    part2 = dPdrho*dedT-dpdt*(Pres/rhob**2-dedrhob)
    third = Mnsp/r**2
    top = first + second*part2 - third
    bot1 = u*dedT
    bot2 = (Mout/(4*pi*u**2*r**2))*u2rhot*part2
    bottom = bot1 + bot2
    tot = top/bottom
    if len(rValues) % 100 == 0:
        print("R = " + str(r / 1e3 * lp) + "  " + "U = " + str(u * (lp/tp)*100) + "  " + "Q = " + str(Q * (lp/tp)*100))

    if tot == np.nan:
        tot = np.inf

    return tot

def dTdrSolved2(q, u, T, rhoTotal, rho, r, P, dRhodr):
    rhobottom = rhoTotal + P
    tfirst = q*(rho**2/P)*u**2*r**3
    tsecond = (T/mN)*dRhodr*(r/rhobottom + r*u**2/rhobottom - 2*Mnsp/rhobottom)
    tthird = Mnsp/r
    tfourth = 2*u**2
    top = tfirst + tsecond - tthird + tfourth

    bottom1 = u**3*r**3*(q*rho**2/(P*u) + dRhodr)
    bottom2 = 11*np.pi**2/45*T**3 + rho/mN
    bottom = bottom1 - bottom2

    total = top / bottom
    return total

def dTdrSolved(q, u, T, rhoTotal, rhob, r, P):
    heating = q*r**4*rhob**2*u**2
    mid = P*Mnsp
    last = -2*P*r*u**2
    top = heating + mid + last

    bot1 = rhob**2*u**3*r**4*((11*np.pi**2/(15*rhob))*T**3+3/(2*mN))
    bot2 = (P/(rhoTotal+P))*((11*np.pi**2/45)*T**3+rhob/mN)*(1+u**2-2*Mnsp/r)
    bottom = bot1+bot2
    tot = top / bottom
    return tot # Currently outpus ~2e-80 then 2e-130 now 1e-39

def dTdr(T, u, r, rhob, P, q):
    last = (-(Mnsp/r**2)+2*u**2/r-(45/(11*np.pi**2))*(u*rhob/T**4)*q)
    prefactor = (1/(1+u**2-(2*Mnsp/r)))
    mid = ((rhob+P)/(4*P))
    tot = T*prefactor*mid*last
    return tot

def TGasLa(r, u):
    rhob = rhobByMout(r, u)

    V = 4 / 3 * pi * r ** 3
    if len(uValues) != 0:
        t = scipy.integrate.simpson(1/np.abs(np.array(uValues)), x=np.array(rValues))
    else:
        t = 0
    beta = ((Mout*t)/mN+Ni)/(Ni)
    if not np.isfinite(beta) or beta < 1:
        beta = 1

    first = (11*pi**2/180*(Ti/KtoMeV/Tp)**3 + (rhoi* 1e3 / (mp / lp**3))/(mN))
    second = 180*4/3*pi*Rnsp**3*beta/(11*pi**2*V)
    third = 180*rhob/(11*pi**2*(mN))
    T3 = (first)*second-third
    T = (T3)**(1/3)
    return T
def dPdr(T, rhob, dTdr, dRhodr):
    dPdt = ((11*np.pi**2/45)*T**3+rhob/mN)
    first = ((dPdt) * dTdr) # *PlanckPtoSI/1.616e-35
    dPdrho = (T/mN)
    second = dPdrho*dRhodr
    tot = first + second
    return tot

def dPdrGrad(pres, rad):
    pres = np.array(pres)
    pres = pres[np.where(pres != 0)[0]]
    tot = np.gradient(pres, rad)
    return tot

def dEdr(T, rhob, dTdr):
    return ((11*np.pi**2/15)*T**3/rhob)*dTdr

def dRhodr(rhob, u, p, q, dEdr):
    val = (q/u - dEdr)*(-rhob**2/p)
    return val

def q(u, dedr, P, rhob, drhobdr):
    return u*(dedr-(P/rhob**2)*(drhobdr))

def P(T, rhob):
    first = (11*np.pi**2/180)*T**4
    second = (rhob/mN)*T
    tot = (first + second)
    return tot

def e(T, rhob):
    first = (11*np.pi**2/60)*T**4/rhob
    second = (3*T/(2*mN))
    tot = (first + second)
    return tot

def rhotot(rhob, e):
    return rhob + rhob*e

def phi(r):
    top = 1 - 2 * Mns * SolarMassKg * 1e3 / (Rnscm)
    bot = 1 - 2 * Mns * SolarMassKg * 1e3 / r

    phi2 = top/bot
    tot = phi2**(1/2)
    return tot

def g1(r):
    first = ((Rnsm * 1e2 / r) ** 2)
    secondtop = 1 - 2 * Mns * SolarMassKg * 1e3 / r
    secondbot = 1 - 2 * Mns * SolarMassKg * 1e3 / (Rnscm)
    second = secondtop/secondbot
    tot2 = 1-first*second
    tot = tot2**(1/2)
    return tot

def g2(r):
    first = (1-g1(r)) ** 4
    second = g1(r) ** 2 + 4*g1(r) + 5
    tot = first*second
    return tot

def q1(r): # r in
    first = (1-Ye)*Lve51*(Eve**2)
    second = Ye*Lave51*(Eave**2)
    firsthalf = 9.65*N_A*(first+second)
    last = (((1-g1(r))/(Rv6**2))*((phi(r))**6))
    return firsthalf*last

def q2(T): # T in Mev output in mev/g*s
    return (2.27*N_A*T**6)

def q3(T, r, rho8): # T in Mev  out in mev/g*s
    first = Lve51*Eve
    second = Lave51*Eave
    third = (6/7)*Lvu51*Evu
    mult = 2.17*N_A*T**4/rho8
    paren = first + second + third
    last = ((1-g1(r))/Rv6**2)*phi(r)**5
    return mult*paren*last

def q4(T, rho8): # T in Mev  out in mev/g*s
    top = T**9
    first = top/rho8
    tot = .144*N_A*first
    return tot

def q5(r, rho8):
    first = Lve51*Lave51*(Eve + Eave)
    second = (6/7)*Lvu51**2*Evu
    paren = first + second
    last = g2(r)/(rho8*Rv6**4)*phi(r)**9
    return 12*N_A*paren*last

def qSum(r, T, rhob):
    print(rhob, rhob/1e8)
    rho8 = rhob / 1e8
    convFactor = 1
    a = 2*q1(r) - 1e27
    b = q2(T) * convFactor
    c = q3(T, r, rho8) * convFactor
    d = q4(T, rho8) * convFactor
    e = q5(r, rho8) * convFactor
    all = (a - b + c - d + e)
    #print(all)
    return all

def initialQ(T):
    r = Rnsm * 1e2
    rho8 = rhoi
    return qSum(r, T, rho8)

def findTi():
    root = optimize.brentq(initialQ, 0, 100)
    return root

def findNi(Ti, Rhobi):
    Vi = 4/3*pi*Rnsp**3
    Pi = P(Ti, Rhobi)
    N = Pi*Vi/Ti
    return N

def solveSystem(r, U, T, P, rhob, E, q, Mout, dPdr, dEdr, dTdr, dRhodr, dUdr):
    eqns = [U*dUdr-(1/(rhob+rhob*e+P))*dPdr(1+U**2-2*Mns/r)+Mns/r**2,
            Mout-4*np.pi*r**2*rhob*U,
            q - U(dEdr- (P/rhob**2)*dRhodr),
            P - (11*np.pi**2/180)*T**4+(rhob/mN)*T,
            dPdr/dTdr - (11*np.pi**2/45)*T**3+(rhob/mN),
            dEdr/dTdr - (11*np.pi**2/15)*T**3+(3/(2*mN))]

def initalConditions():
    R = float(Rnsp)
    t = 0
    Q = 0
    TMev = findTi()
    Ti = TMev
    Temp = TMev * MeVtoK
    Tplank = Temp / Tp
    Ni = findNi(Temp, rhoi)
    rhoPlank = rhoi * 1e3 / (mp / lp**3)
    Pres = P(Tplank, rhoPlank)
    Energy = e(Tplank, rhoPlank)
    rhoTotal = rhotot(rhoPlank, Energy)
    U = uByMout(R, rhoPlank)

    return R, U, Ni, Ti

def solveWithSolver(Y0):
    sol = scipy.integrate.solve_ivp(dudrsolved, (10.00001*1e3/lp, termRadius*1e3/lp), [Y0])
    return sol

def setupCalculations(Mout):
    """RRun = [0 for i in range(1000)]
    PresRun = [0 for i in range(1000)]
    URun = [0 for i in range(1000)]
    TRun = [0 for i in range(1000)]
    RhoRun = [0 for i in range(1000)]
    QRun = [0 for i in range(1000)]
    timeRun = [0 for i in range(1000)]"""
    R = float(Rnsp)
    t = 0
    Q = 0
    TMev = findTi()
    Ni = findNi(Ti, rhoi)
    Temp = TMev * MeVtoK
    Tplank = Temp / Tp
    rhoPlank = rhoi * 1e3 / (mp / lp**3)
    Pres = P(Tplank, rhoPlank, mN)
    Energy = e(Tplank, rhoPlank, mN)
    rhoTotal = rhotot(rhoPlank, Energy)
    U = uByMout(Mout, R, rhoPlank)

    """RRun[0] = (R)
    PresRun[0] = Pres
    URun[0] = (U)
    TRun[0] = (TMev)
    RhoRun[0] = (rhoPlank)
    QRun[0] = (Q)
    timeRun[0] = (t)"""

    RhoPrime = 0
    TempPrime = 0
    PresPrime = 0
    VelPrime = 0

    rStep = 1e-25
    k = 0
    i = 1
    with open(path + "traj_005.dat", 'w', encoding="utf-8") as f:
        while R*lp < 20.5*1e3:
            tStep = rStep/U
            t += tStep
            R += rStep


            Tplank += TempPrime * rStep
            Temp = Tplank * Tp
            TMev = Temp * KtoMeV

            U += VelPrime * rStep

            rhoPlank = rhobByMout(Mout, R, U)

            Energy = e(Tplank, rhoPlank, mN)

            Q = np.absolute(qSum((R * lp * 1e3), TMev, rhoPlank*(mp/lp**3 * 1e3))) * (MeVtoK / (Tp / (tp * mp)))
            Pres = P(Tplank, rhoPlank, mN)
            rhoTotal = rhotot(rhoPlank, Energy)

            if i % 1 == 0:
                f.write(str(R*lp/1e3) + "   " + str(U * (lp/tp)*1e2/1e7) + "    " + str(TMev*10) + "    " + str(rhoPlank*(mp/lp**3)/1e3 / 1e8) + "  " + str(Q / (MeVtoK / (Tp / (tp * mp)))) + "   " + str(t*tp) + "\n")

            """RRun[i % 1000] = (R)
            PresRun[i % 1000] = (Pres)
            URun[i % 1000] = (U)
            TRun[i % 1000] = (Tplank)
            RhoRun[i % 1000] = (rhoPlank)
            QRun[i % 1000] = (Q)
            timeRun[i % 1000] = (t)"""
            RhoPrime = dRhodrNum(Mout, R, U)
            TempPrime = -dTdrSolved2(Q, U, Tplank, rhoTotal, rhoPlank, R, Pres, RhoPrime)
            PresPrime = dPdr(Tplank, rhoPlank, TempPrime, RhoPrime)
            VelPrime = -dudr(U, rhotot(rhoPlank, Energy), Pres, PresPrime, R)
            #print(Rdec)

            i += 1
            if i % 1 == 0:
                print(str(R*lp/1e3) + "   " + str(U * (lp/tp)*1e2/1e7) + "    " + str(TMev*10) + "    " + str(rhoPlank*(mp/lp**3)/1e3 / 1e8) + "  " + str(Q / (MeVtoK / (Tp / (tp * mp)))) + "   " + str(t*tp) + "\n")
                print(str(TempPrime) + "   " + str(PresPrime) + "   " + str(VelPrime) + "\n\n")
                """print("R = " + str(R*lp/1e3))
                print("U = " + str((U * (lp/tp))*1e2/1e7))
                print("Q = " + str(Q / (MeVtoK / (Ep / (tp * mp)))))
                print("Acceleration = " + str(VelPrime))
                print("PresPrime = " + str(PresPrime))
                print("TempPrime = " + str(TempPrime))
                print("TMeV = " + str(TMev))
                print("T = " + str(Temp/1e9))
                print("Time = " + str(t) + "\n\n")
                fig, axs = plt.subplots(2,2)
                axs[0][0].plot(RRun, QRun, label="Q")
                axs[0][0].set_ylim(0, 2e27)
                axs[0][0].set_yscale("linear")
                axs[0][0].set_ylabel("Q")
                axs[0][0].set_xlabel("Rad")
    
                axs[1][0].plot(RRun, URun, label="T")
                axs[1][0].set_ylim(0.001, 1000)
                axs[1][0].set_yscale("log")
                axs[1][0].set_ylabel("U")
                axs[1][0].set_xlabel("Rad")
    
                axs[0][1].plot(RRun, TRun, label="Rho")
                axs[0][1].set_ylim(0.001, 1000)
                axs[0][1].set_yscale("log")
                axs[0][1].set_ylabel("Temp")
                axs[0][1].set_xlabel("Rad")

                axs[1][1].plot(RRun, RhoRun, label="R")
                axs[1][1].set_ylim(0.001, 1000)
                axs[1][1].set_yscale("linear")
                axs[1][1].set_ylabel("Rho")
                axs[1][1].set_xlabel("Rad")
    
                for l in axs:
                    for m in l:
                        m.set_xscale("log")
                plt.savefig(("/home/thana1dr/Research_Results&Data/Week8/Automatic/TrajCalc/First/Traj{:03d}.png".format(k)), dpi=200)
                k += 1
                #plt.show()
                plt.close()"""



    # It is now running semi-stabely. Something is up with the dtdr function may have to make my own :(

    print("Initiated . . .")



R, U, Ni, Ti = initalConditions()
print("TInitial = "  +str(Ti))
solution = solveWithSolver(U)

print("TInitial = "  +str(Ti))

uValues = np.array(uValues) * (lp/tp) * 1e2
rValues = np.array(rValues) * lp / 1e3
qValues = np.array(qValues)
TValues = np.array(TValues)
rhoValues = np.array(rhoValues) * (mp/lp**3) / 1e3

offset=1


plt.plot(rValues, TValues[offset:], label = "TMeV")
plt.ylabel("Temp (MeV)")
plt.xlabel("Radius (Km)")
plt.show()
plt.close()

plt.plot(rValues, uValues, label = "TMeV")
plt.ylabel("Velocity (cm/s)")
plt.xlabel("Radius (Km)")
plt.show()
plt.close()

plt.plot(rValues, rhoValues, label = "TMeV")
plt.ylabel("Density (g/cm^3)")
plt.xlabel("Radius (Km)")
plt.show()
plt.close()


#plt.plot(np.array(rValues) * lp / 1e3, np.array(qValues)[:], label="QTotal")
plt.plot(rValues, q1(rValues*1e5), label = "Q1")
plt.plot(rValues, q2(TValues[offset:]) * 1e-9, label = "Q2", linestyle = "dashed")
plt.plot(rValues, q3(TValues[offset:], rValues*1e5, rhoValues * 1e-8), label = "Q3", linestyle = "dotted")
plt.plot(rValues, q4(TValues[offset:], rhoValues * 1e-8), label = "Q4")
plt.plot(rValues, q5(rValues*1e5, rhoValues * 1e-8), label = "Q5")
plt.ylabel("Neutrino Heating (1e27 MeV/(gs)")
plt.xlabel("Radius (Km)")
plt.yscale("linear")
plt.xscale("log")
plt.legend()
plt.show()
plt.close()

plt.plot(np.array(rValues), np.array(qValues)[offset:], label="QTotal")
plt.ylabel("Total Neutrino Heating (1e27 MeV/(gs)")
plt.xlabel("Radius (Km)")
plt.yscale("linear")
plt.xscale("log")
plt.legend()
plt.show()
plt.close()


print("Done")