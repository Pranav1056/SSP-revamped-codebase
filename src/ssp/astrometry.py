"""Astronomy calculations: astrometry."""

import numpy as np
from .coordinates import convertToDecDegrees

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
