import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

now = datetime.now()

cs = 340

def EqualTurbulentCells(rMin, rMax, lc, res=1000):

    r = np.linspace(rMin, rMax, res)

    eddyV = SoundSpeed(r)*np.sin((r*2*np.pi)/lc)

    plt.plot(r, eddyV)
    plt.xlabel('Distance from centre of black hole')
    plt.ylabel('Speed')
    plt.savefig(f'{now}.png')
    plt.show()

    return [r, eddyV]

def RandomTurbulentCells(rMin, rMax, met='kolmogorov', lc=1, sd=0.1, res=1000):

    sign = -1

    r = np.linspace(rMin, rMax, res)

    eddyV = np.ones(len(r))

    currentStart = rMin
    eddySize = np.random.normal(lc, sd)
    while eddySize > lc:
        eddySize = np.random.normal(lc, sd)
    currentEnd = r[find_nearest(r, currentStart + eddySize)]
    

    while rMax > currentEnd:

        minIndex = np.where(r==currentStart)[0][0]
        maxIndex = np.where(r==currentEnd)[0][0]

        domain = np.linspace(0, np.pi, (maxIndex-minIndex)+1)
        midPoint = int(np.floor((minIndex+maxIndex)/2))
        eddyV[minIndex:maxIndex+1] = sign * SpeedScale(r[midPoint], method=met) * np.sin(domain)

        sign *= -1
        currentStart = currentEnd
        eddySize = np.random.normal(lc, sd)
        while eddySize > lc:
            eddySize = np.random.normal(lc, sd)

        currentEnd = r[find_nearest(r, currentStart + eddySize)]

    eddySize = rMax-currentStart
    currentEnd = rMax

    minIndex = np.where(r==currentStart)[0][0]
    maxIndex = np.where(r==currentEnd)[0][0]

    domain = np.linspace(0, np.pi, (maxIndex-minIndex)+1)
    midPoint = int(np.floor((minIndex+maxIndex)/2))
    eddyV[minIndex:maxIndex+1] = sign * SpeedScale(r[midPoint], method=met) * np.sin(domain)


    plt.plot(r, eddyV)
    plt.xlabel('Distance from centre of black hole')
    plt.ylabel('Speed')
    plt.savefig(f'{now}.png')
    plt.show()


def VarHeightTurbulentCells(rMin, rMax, met = 'kolmogorov', sd=0.1, res=1000):

    sign = -1

    r = np.linspace(rMin, rMax, res)

    eddyV = np.ones(len(r))

    currentStart = rMin
    sizeLimit = Height(currentStart)
    eddySize = np.random.normal(sizeLimit, sd)
    while eddySize > sizeLimit:
        eddySize = np.random.normal(sizeLimit, sd)
    currentEnd = r[find_nearest(r, currentStart + eddySize)]
    

    while rMax > currentEnd:

        minIndex = np.where(r==currentStart)[0][0]
        maxIndex = np.where(r==currentEnd)[0][0]

        domain = np.linspace(0, np.pi, (maxIndex-minIndex)+1)
        midPoint = int(np.floor((minIndex+maxIndex)/2))
        eddyV[minIndex:maxIndex+1] = sign * SpeedScale(r[midPoint], method=met) * np.sin(domain)

        sign *= -1
        currentStart = currentEnd
        sizeLimit = Height(currentStart)
        eddySize = np.random.normal(sizeLimit, sd)
        while eddySize > sizeLimit:
            eddySize = np.random.normal(sizeLimit, sd)
        currentEnd = r[find_nearest(r, currentStart + eddySize)]

    eddySize = rMax-currentStart
    currentEnd = rMax

    minIndex = np.where(r==currentStart)[0][0]
    maxIndex = np.where(r==currentEnd)[0][0]

    domain = np.linspace(0, np.pi, (maxIndex-minIndex)+1)
    midPoint = int(np.floor((minIndex+maxIndex)/2))
    eddyV[minIndex:maxIndex+1] = sign * SpeedScale(r[midPoint], method=met) * np.sin(domain)


    plt.plot(r, eddyV)
    plt.xlabel('Distance from centre of black hole')
    plt.ylabel('Speed')
    plt.show()



def find_nearest(array, value):
    
    return (np.abs(array - value)).argmin()


def SpeedScale(r, method='kolmogorov', sd=3):

    if method == 'sound':
        return SoundSpeed(r, sd)
    elif method == 'kolmogorov':
        return KolmogorovSpeed(r)


def SoundSpeed(r, sd = 3):

    centre = 100-0.5*r
    return np.random.normal(centre, sd)


def KolmogorovSpeed(r):

    return 100 + 0.5*r


def Height(r):

    return 2+0.1*r


VarHeightTurbulentCells(0, 100)