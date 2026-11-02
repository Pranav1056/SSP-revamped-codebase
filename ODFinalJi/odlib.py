import numpy as np
from astropy.io import fits 
import scipy as scipy
#want v in gaussian
def getMomentum(r, v):
    return np.cross(r, v)

#returns degrees
import math as m
def getAngle(cosine, sine):
    initialAcos = m.acos(cosine)
    initialAsin = m.asin(sine)
    if (abs(cosine**2 + sine**2 - 1) > 0.1):
        return "error"
    else:
        result = 0
        if (sine == 0 and cosine == 1):
            result = 0
        elif (sine == 1 and cosine == 0):
            result = m.pi*0.5
        elif (sine == 0 and cosine == -1):
            result = m.pi
        elif (sine == -1 and cosine == 0):
            result = 1.5*m.pi
        elif (sine > 0 and cosine > 0):
            result = abs(initialAsin)
        elif (sine > 0 and cosine < 0):
            result = initialAcos
        elif (sine < 0 and cosine < 0):
            result = 1*m.pi + abs(initialAsin)
        elif (sine < 0 and cosine > 0):
            result = 2*m.pi-abs(initialAsin)
    return m.degrees(result)

def getPerihelion(rDotVec, hVec):
    return np.cross(-1*hVec, rDotVec) - np.linalg.norm(rDotVec)

def getA(rVec, rDotVec):
    return 1/((2/np.linalg.norm(rVec))-np.linalg.norm(rDotVec)**2)

#RETURNS IN DEGREES!!!!
def getInclination(hVec):
    return np.degrees(np.arctan(np.sqrt(hVec[0]**2 + hVec[1]**2)/hVec[2]))

def getEccentricity(hVec, a):
    return np.sqrt(1-(np.linalg.norm(hVec)**2/a))

# RETURNS IN DEGREES
def getLongitudeOfAscendingNode(hVec, iDegrees):
    cosOmega = -hVec[1]/(np.linalg.norm(hVec)*np.sin(np.radians(iDegrees)))
    sinOmega = hVec[0]/(np.linalg.norm(hVec)*np.sin(np.radians(iDegrees)))
    return getAngle(cosOmega, sinOmega)

def getTrueAnomaly(a, e, hVec, rVec, rDotVec):
    rMag = np.linalg.norm(rVec)
    hMag = np.linalg.norm(hVec)
    sinV = ((a*(1 - e**2))/(e*hMag))*((np.dot(rVec, rDotVec))/rMag)
    # print("sinV:", sinV)
    cosV = (1/e)*(((a*(1 - e**2))/rMag) - 1)
    # print("cosV:", cosV)
    return getAngle(cosV, sinV)

def getArgOfPerihelion(tAnomaly, rVec, LAN, i):
    x = rVec[0]
    y = rVec[1]
    z = rVec[2]
    rMag = np.linalg.norm(rVec)
    omega = np.radians(LAN)

    cosU = (x*np.cos(omega) + y*np.sin(omega))/rMag
    sinU = z/(rMag*np.sin(np.radians(i)))
    # print(sinU)

    u = getAngle(cosU, sinU)
    initialArg = np.radians(u)-np.radians(tAnomaly)
    while initialArg < 0:
        initialArg = initialArg + 2*np.pi
    while initialArg > 2*np.pi:
        initialArg = initialArg - 2*np.pi
    return np.degrees(initialArg)
    
def lspr(fileName, x, y):
    file = open(fileName)
    info = np.loadtxt(file, str)

    #getting the RA column and turning it into a list
    iStarRAs = info[:, 2:3]
    iStarRAs = np.transpose(iStarRAs)[0]

    starRAs = []
    #converting into decimal hours
    for i in range(0, len(iStarRAs)):
        # sum = float(iStarRAs[i][:2])*15 + float(iStarRAs[i][3:5])*0.25 + float(iStarRAs[i][6:])*0.25/60
        # starRAs.append(float(sum))
        starRAs.append(float(iStarRAs[i]))
    # print(len(starRAs))

    iStarDecs = info[:, 3:]
    iStarDecs = np.transpose(iStarDecs)[0]
    
    starDecs = []
    #converting into decimal degrees
    for i in range(0, len(iStarDecs)):
        # coef = 1
        # if (iStarDecs[i][0] == "-"):
        #     coef = -1
        # sum = float(iStarDecs[i][1:3]) + float(iStarDecs[i][4:6])/60 + float(iStarDecs[i][7:])/3600
        starDecs.append(float(iStarDecs[i]))
    # print(starDecs)

    #converting star x's into a list
    iStarXs = info[:, :1]
    iStarXs = np.transpose(iStarXs)[0]
    starXs = []
    for i in range(0, len(iStarXs)):
        starXs.append(float(iStarXs[i]))
    # print(starXs)

    #converting star y's into a list
    iStarYs = info[:, 1:2]
    iStarYs = np.transpose(iStarYs)[0]
    starYs = []
    for i in range(0, len(iStarYs)):
        starYs.append(float(iStarYs[i]))
    # print(starYs)

    #making x square and y square sums
    xSquareSum = np.sum(np.array(starXs)*np.array(starXs))

    ySquareSum = np.sum(np.array(starYs)*np.array(starYs))

    xySum = np.sum(np.array(starXs)*np.array(starYs))

    starYs = np.array(starYs)
    starXs = np.array(starXs)
    starRAs = np.array(starRAs)
    starDecs = np.array(starDecs)
    # print(starRAs)

    #creating first matrix
    xyMatrix = np.array([[len(starRAs), np.sum(starXs), np.sum(starYs)], [np.sum(starXs), xSquareSum, xySum], 
                                   [np.sum(starYs), xySum, ySquareSum]])

    xyMatrix = np.linalg.inv(xyMatrix)
    
    raMatrix = np.array([[np.sum(starRAs)], [np.sum(starRAs*starXs)], 
                        [np.sum(starRAs*starYs)]])
    
    decMatrix = np.array([[np.sum(starDecs)], [np.sum(starDecs*starXs)], 
                        [np.sum(starDecs*starYs)]])
    
    raAnswers = np.dot(xyMatrix, raMatrix)
    decAnswers = np.dot(xyMatrix, decMatrix)

    # print(raAnswers)
    # print(decAnswers)

    RA = raAnswers[0] + raAnswers[1]*x + raAnswers[2]*y
    Dec = decAnswers[0] + decAnswers[1]*x + decAnswers[2]*y

    uSumRA = 0
    for i in range(0, len(starRAs)):
        uSumRA = uSumRA + (starRAs[i] - raAnswers[0] - raAnswers[1]*starXs[i] - raAnswers[2]*starYs[i])**2
    
    uSumDec = 0
    errorList = []
    for i in range(0, len(starDecs)):
        uSumDec = uSumDec + (starDecs[i] - decAnswers[0] - decAnswers[1]*starXs[i] - decAnswers[2]*starYs[i])**2
        errorList.append((starDecs[i] - decAnswers[0] - decAnswers[1]*starXs[i] - decAnswers[2]*starYs[i])*3600)
    
    # print(errorList, "<------ the error list")
    raU = np.sqrt(uSumRA/(len(starRAs) - 3))
    decU = np.sqrt(uSumDec/(len(starDecs) - 3))
    
    # print(decU*3600, "decu")
    # print(raU*3600, "raU")

    # print(RA)
    # print(Dec)

    RA = RA[0]
    Dec = Dec[0]

    RA = RA*24/360
    # print(RA)
    RAhours = int(RA)
    RAminutes = int((RA - RAhours)*60)
    RAseconds = round(3600*(RA - RAhours - RAminutes/60), 3)
    raString = str(RAhours)+":"+str(RAminutes)+":"+str(RAseconds)

    decHours = int(Dec)
    decMinutes = int((Dec - decHours)*60)
    decSeconds = round(3600*(Dec - decHours - decMinutes/60), 3)
    decString = str(decHours)+":"+str(decMinutes)+":"+str(decSeconds)

    print("*"*14)
    print("plate constants")
    print("*"*14)
    print("b1:", raAnswers[0][0], "deg")
    print("b2:", decAnswers[0][0], "deg")
    print("a11:", raAnswers[1][0], "deg/pix")
    print("a12:", raAnswers[2][0], "deg/pix")
    print("a21:", decAnswers[1][0], "deg/pix")
    print("a22:", decAnswers[2][0], "deg/pix")
    print("*"*11)
    print("uncertainty")
    print("*"*11)
    print(" RA:", round(3600*raU[0], 3), "arcsec")
    print("Dec:", round(3600*decU[0], 3), "arcsec")
    print("*"*25)
    print("astrometry for")
    print("(x, y) =", x, y)
    print("*"*25)
    print(" RA =", raString)
    print("Dec =", decString)

# mInitial = 0.42
# eInitial = 0.8
# tolerance = 0.00000000001
def newtonRaphson(m, e, tol):
    oldEtest = m
    ETest = oldEtest - (oldEtest - e*np.sin(oldEtest) - m)/(-e*np.cos(oldEtest) + 1)
    iterationNum = 0
    while (abs(ETest - oldEtest) > tol):
        iterationNum = iterationNum + 1
        oldEtest = ETest
        ETest = oldEtest - (oldEtest - e*np.sin(oldEtest) - m)/(-e*np.cos(oldEtest) + 1)
    # print("E: ", ETest)
    # print("Tolerance: ", tol)
    # print("Iterations: ", iterationNum)
    return ETest

#time is a string- hh:mm:ss
def julian(isLocalTime, yr, month, day, time):
    term1 = (month+9)//12
    term2 = 7*(term1 + yr)//4
    julianAt0 = 367*yr - term2 + 275*month//9 + day + 1721013.5
    UT = float(time[0:2])+ float(time[3:5])/60 + float(time[6:])/3600
    if isLocalTime:
        UT = UT + 4
    return julianAt0  + UT/24
# print(julian(False, 2018, 7, 13, "00:00:00"))

#returns in degrees!!!!!!!!
def getMeanAnomaly(rVec, trueAnomaly, a, e):
    r = np.linalg.norm(rVec)
    E0 = np.arccos((1/e)*(1-(r/a)))
    if(np.radians(trueAnomaly) > np.pi):
        return 360 - np.degrees(E0 - e*np.sin(E0))
    return np.degrees(E0 - e*np.sin(E0))

def getTimeOfLastPerihelion(meanAnomaly, a, julianDate):
    k = 0.01720209895
    return -1*((np.radians(meanAnomaly)/np.sqrt((k**2)/(a**3))) - julianDate)

def convertToDecimalHours(stringToConvert, hasSign, toDegrees):
    coef = 1
    sum = 0
    if(hasSign):
        if (stringToConvert[0] == "-"):
            coef = -1
        sum = coef*(float(stringToConvert[1:3]) + float(stringToConvert[4:6])/60 + float(stringToConvert[7:])/3600)
    else:
        sum = float(stringToConvert[:2]) + float(stringToConvert[3:5])/60 + float(stringToConvert[6:])/3600
    if(toDegrees):
        sum = sum*360/24
    return sum

def convertToDecDegrees(stringToConvert, hasSign, wasHours):
    if (not wasHours):
        coef = 1
        sum = 0
        if(hasSign):
            if (stringToConvert[0] == "-"):
                coef = -1
            sum = coef*(float(stringToConvert[1:3]) + float(stringToConvert[4:6])/60 + float(stringToConvert[7:])/3600)
        else:
            sum = float(stringToConvert[:2]) + float(stringToConvert[3:5])/60 + float(stringToConvert[6:])/3600
        return sum
    else: 
        coef = 1
        sum = 0
        if(hasSign):
            if (stringToConvert[0] == "-"):
                coef = -1
            sum = coef*((float(stringToConvert[1:3])*15) + float(stringToConvert[4:6])*0.25 + float(stringToConvert[7:])*360/86400)
        else:
            sum = (float(stringToConvert[:2])*15) + float(stringToConvert[3:5])*0.25 + float(stringToConvert[6:])*360/86400
        return sum

def getRaAndDec(a, eccentricity, timePerihelion, inclination, LAN, argOfPerihelion, currentTimeJulian, earthSunVec):
    k = 0.01720209895
    n = k*np.sqrt(1/(a**3))
    initialM = n*(currentTimeJulian - timePerihelion)
    E = newtonRaphson(initialM, eccentricity, 0.000000001) #may change tolerance

    inclination = np.radians(inclination)
    LAN = np.radians(LAN)
    argOfPerihelion = np.radians(argOfPerihelion)

    x = a*np.cos(E) - a*eccentricity
    y = a*np.sqrt(1-eccentricity**2)*np.sin(E)
    z = 0

    posVec = np.array([[x], [y], [z]])
    matrix1 = np.array([[np.cos(LAN), -np.sin(LAN), 0], [np.sin(LAN), np.cos(LAN), 0], [0, 0, 1]])
    matrix2 = np.array([[1, 0, 0], [0, np.cos(inclination), -np.sin(inclination)], [0, np.sin(inclination), np.cos(inclination)]])
    matrix3 = np.array([[np.cos(argOfPerihelion), -np.sin(argOfPerihelion), 0], [np.sin(argOfPerihelion), np.cos(argOfPerihelion), 0], [0, 0, 1]])

    newVector = np.matmul(matrix1, np.matmul(matrix2, np.matmul(matrix3, posVec)))
    epsilon = np.radians(23.44)
    rotationMatrix = np.array([[1, 0, 0], [0, np.cos(epsilon), -np.sin(epsilon)], [0, np.sin(epsilon), np.cos(epsilon)]])
    equitorialPosition = np.matmul(rotationMatrix, newVector)
    rangeVec = equitorialPosition + earthSunVec
    # print("rangevec: ", rangeVec)
    pHat = rangeVec/np.linalg.norm(rangeVec)
    # print(pHat)

    Dec = np.arcsin(pHat[2][0])
    RA = getAngle(pHat[0][0]/np.cos(Dec), pHat[1][0]/np.cos(Dec))
    return RA, np.degrees(Dec)

def getCircle(image, x, y, radius, outerRadius):
    square = np.array(image[int(y)-outerRadius:int(y)+outerRadius+1, int(x)-outerRadius:int(x)+outerRadius+1])
    centerCoord = outerRadius
    count = 0
    for i in range(0, square.shape[0]):
        for j in range(0, square.shape[0]):
            if (i-centerCoord)**2 + (j-centerCoord)**2 > radius**2:
                square[i][j] = 0
            else:
                count = count + 1
    return square, count


def aperturePhotometry(x, y, fileName, radius, innerSkyRadius, outerSkyRadius):
    data = fits.getdata(fileName)
    outer, counter1 = getCircle(data, x, y, outerSkyRadius, outerSkyRadius)
    innerI, counter2 = getCircle(data, x, y, innerSkyRadius, outerSkyRadius)
    ringMatrix = outer-innerI
    skyAvg = ringMatrix.sum()/(counter1-counter2)
    inner = np.zeros((innerI.shape[0], innerI.shape[0]))

    for i in range(0, innerI.shape[0]):
        for j in range(0, innerI.shape[0]):
            if (innerI[i][j] > 0):
                inner[i][j] = innerI[i][j]-skyAvg

    xVec = range(0, outerSkyRadius*2+1)
    yVec = range(0, outerSkyRadius*2+1)
    xMatrix, yMatrix = np.meshgrid(xVec, yVec)
    xCentroid = (inner*xMatrix).sum()/inner.sum() + x - outerSkyRadius
    yCentroid = (inner*yMatrix).sum()/inner.sum() + y - outerSkyRadius

    outer, counter1 = getCircle(data, x, y, outerSkyRadius, outerSkyRadius)
    innerI, counter2 = getCircle(data, x, y, innerSkyRadius, outerSkyRadius)
    ringMatrix = outer-innerI
    #shift annulus
    aperture, count3 = getCircle(data, xCentroid, yCentroid, radius, outerSkyRadius)
    # print(count3) #should be 78
    outer1, counter12 = getCircle(data, x, y, outerSkyRadius, outerSkyRadius)
    innerI1, counter22 = getCircle(data, x, y, innerSkyRadius, outerSkyRadius)
    ringMatrix2 = outer1-innerI1
    ringAsList = []

    for i in range(0, ringMatrix2.shape[0]):
        for j in range(0, ringMatrix2.shape[0]):
            if ringMatrix2[i][j] > 0:
                ringAsList.append(ringMatrix2[i][j])
    skyMedian = np.median(np.array(ringAsList))
    signal = aperture.sum() - (skyMedian)*count3
    mInst = -2.5*np.log10(signal)
    SNR = signal/(np.sqrt(signal + count3*(1+(count3/(counter12-counter22)))*(skyMedian+10+(11)**2)))
    mUncertainty = 1.0857/SNR
    signalUncertainty = signal/SNR
    return xCentroid, yCentroid, signal, signalUncertainty, SNR, mInst, mUncertainty
    
# print(aperturePhotometry(490, 293, "D:/Ji_OD/aptest.fit", 5, 8, 13))


def fg(tau, r2, r2dot, order):
    r = np.linalg.norm(r2)
    a = getA(r2, r2dot)
    mu = 1
    n = np.sqrt(mu/(a**3))
    e = getEccentricity(getMomentum(r2, r2dot), a)
    sign = (np.dot(r2, r2dot)/(n*(a**2)))*np.cos(n*tau - (np.dot(r2, r2dot))/(n*(a**2)))
    + (1-r/a)*np.sin(n*tau - (np.dot(r2, r2dot)/(n*(a**2))))
    if (sign > 0):
        x01 = n*tau + 0.85*e - (np.dot(r2, r2dot)/(n*(a**2)))
    else:
        x01 = n*tau - 0.85*e - (np.dot(r2, r2dot)/(n*(a**2)))
    lambdaThing = lambda x: x - (1-(r)/a)*np.sin(x) + (np.dot(r2, r2dot)/(n*(a**2)))*(1-np.cos(x)) - (n*tau)
    E = scipy.optimize.newton(lambdaThing, x0 = x01, tol=1E-12)
    if order == 5:
        return 1-(a/(r))*(1-np.cos(E)), tau + (1/n)*(np.sin(E) - E)
    u = mu/(r**3)
    z = (np.dot(r2, r2dot)/(r**2))
    q = (np.dot(r2dot, r2dot)/(r**2)) - u
    tau = tau
    if(order == 3):
        f1 = 1 - 0.5*u*(tau**2) + 0.5*u*z*(tau**3)
        g1 = tau - (1/6)*mu*(tau**3)/r**3 
    if(order == 4):
        f1 = 1 - 0.5*mu*(tau**2) + 0.5*u*z*(tau**3) + (1/24)*(3*u*q - 15*u*(z**2) + u**2)
        g1 = tau - (1/6)*mu*(tau**3)/r**3  + (1/4)*u*z*(tau**4)
    return f1, g1

def initialFg(r, tau):
    mu = 1.0
    u = mu/r**3
    f1 = 1 - 0.5*u*(tau**2)
    g1 = tau - (1/6)*mu*(tau**3)/r**3 
    return f1, g1

def scaleEqLagrange(tauVec, sunVec, rhoHatVec, dsVec):
    mu = 1
    tau = tauVec[2]
    tau1 = tauVec[0]
    tau3 = tauVec[1]
    A1 = tau3/tau
    B1 = (A1/6)*((tau**2) - (tau3**2))
    A3 = -tau1/tau
    B3 = (A3/6)*((tau**2) - (tau1**2))

    A = (A1*dsVec[1] - dsVec[2] + A3*dsVec[3])/(-dsVec[0])
    B = (B1*dsVec[1] + B3*dsVec[3])/(-dsVec[0])

    E = (np.dot(rhoHatVec, sunVec))*-2
    F = (np.linalg.norm(sunVec))**2
    a = -((A**2) + A*E  + F)
    b = -mu*(2*A*B + B*E)
    c = -(mu**2)*(B**2)

    # print([c, 0, 0, b[0], 0, 0, a[0], 0, 1])
    # print(a, b, c)

    answers = np.polynomial.polynomial.polyroots([c, 0, 0, b, 0, 0, a, 0, 1])
    rootFinal = []
    for root in answers:
        if np.imag(root) == 0 and np.real(root) > 0:
            rootFinal.append(np.real(root))
    rhos = []
    for root2 in rootFinal:
        rhos.append(A + (mu*B)/(root2**3))

    return rootFinal, rhos
# print(scaleEqLagrange([-0.15481889055, 0.15481889055, 0.3096377811], 
#                 [-0.2398478458274071, 0.9065739917845802, 0.3929623749770952], 
#                 [-0.8518563498182248, -0.2484702599212149, 0.4610892421311239], 
#                 [-0.0010461861084885213, -0.17297581974209159, -0.17201260125558127, -0.16712421570714076]))