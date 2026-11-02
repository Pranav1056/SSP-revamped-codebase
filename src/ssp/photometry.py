"""Astronomy calculations: photometry."""

import numpy as np

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
    from astropy.io import fits

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
