import numpy as np
import scipy.integrate as intgr
import matplotlib.pyplot as plot

def IntegrationFunction(x):
    Ez = np.sqrt(OM*(1+x)**3 + Ok*(1+x)**2 + OL)
    func = 1 / Ez
    return func

def CalcDCByDH(z):
    integral = intgr.quad(IntegrationFunction, 0, z)
    return integral[0]

def CalcDAByDH(z):
    dcdh = CalcDCByDH(z)
    dadh = 0
    if (Ok < 0):
        dadh = (1 / (np.sqrt(Ok) * (1+z) ) ) * np.sinh(np.sqrt(Ok) * dcdh)
    elif (Ok > 0):
        dadh = (1 / (np.sqrt(np.abs(Ok)) * (1+z) ) ) * np.sin(np.sqrt(np.abs(Ok)) * dcdh)
    else:
        dadh = dcdh / (1+z)
    return dadh

def CalcDLByDH(z):
    dadh = CalcDAByDH(z)
    dldh = dadh * ((1+z)**2)
    return dldh

def BisectToGetzMax(zl, zu, dlmax): #zl is the lower limit of interval, zu is upper limit, dlmax is the calculated value of DLmax to compare to
    midpoint = (zl+zu)/2
    dlcalc = CalcDLByDH(midpoint) * DH
    diff = dlcalc - dlmax #if +ve, midpoint gives zmax as too big, if -ve, midpoint gives zmax too small
    percDiff = abs(diff) / dlmax
    if (percDiff <= 0.00001):
        return midpoint
    elif (diff > 0):
        return BisectToGetzMax(zl, midpoint, dlmax)
    else:
       return BisectToGetzMax(midpoint, zu, dlmax)
        

OMValues = [1.0, 0.27]
OLValues = [0.0, 0.73]
ind = 1
OM = OMValues[ind]
OL = OLValues[ind]
Ok = 1 - (OM + OL)
h = 0.72
tH = 3.0856 * (10**17) / h
c = 299792458
DH = c * tH
zvalues = []
fbyf0values = []
file = open("quasar.dat", "r") #quasar.dat here has the column titles removed
Varray = []
ratioArray = []
newDataArray = []
for r in file:
    k = r.split()
    zvalues = np.append(zvalues, k[0])
    fbyf0values = np.append(fbyf0values, k[1])
for i in range(len(zvalues)):
    z = float(zvalues[i])
    ff0 = float(fbyf0values[i])
    dcdh = CalcDCByDH(z)
    dc = dcdh * DH
    dldh = CalcDLByDH(z)
    dl = dldh * DH
    dlmax = np.sqrt(ff0) * dl
    zmax = BisectToGetzMax(0, 20, dlmax)
    #print("z = " + str(z) + " zmax = " + str(zmax))
    V = (4 * np.pi / 3) * ((dl/(1+z)) ** 3)
    Vmax = (4 * np.pi / 3) * ((dlmax/(1+zmax)) ** 3)
    VbyVmax = V / Vmax
    propConst = VbyVmax * (ff0 ** 1.5)
    ratioArray = np.append(ratioArray, propConst)
    newDataArray = np.append(newDataArray, VbyVmax)
m = np.mean(newDataArray)
print("<V/Vmax> = " + str(m) + "\n")
#uncomment block below to print data to demonstrate euclidean limit
"""
for i in range(len(ratioArray)):
    print("z = " + str(zvalues[i]) + " (V/Vmax) * (f/f0)^3/2 = " + str(ratioArray[i]))
"""