# Cython routines for nemo

#include "python.pxi"
#include "numpy.pxi"

from astLib import *
import numpy as np
cimport numpy as np
import cython
import math
import time
import sys
from cython cimport parallel

#-------------------------------------------------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def makeDegreesDistanceMap(np.ndarray[np.float64_t, ndim=2] degreesMap, wcs, RADeg, decDeg, maxDistDegrees):
    """Fills (in place) the 2d array degreesMap with distance in degrees from the given position, out to some
    user-specified maximum distance.
    
    Args:
        degreesMap (:obj:`np.ndarray`): Map (2d array) that will be filled with angular distance from the 
            given coordinates. Probably you should feed in an array set to some extreme initial value (e.g.,
            1e6 everywhere) to make it easy to filter for pixels near the object coords afterwards.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to degreesMap.
        RADeg (float): RA in decimal degrees of position of interest (e.g., object location).
        decDeg (float): Declination in decimal degrees of position of interest (e.g., object location).
        maxDistDegrees: The maximum radius out to which distance will be calculated.
    
    Returns:
        A map (2d array) of distance in degrees from the given position,
        (min x, max x) pixel coords corresponding to maxDistDegrees box, 
        (min y, max y) pixel coords corresponding to maxDistDegrees box
    
    Note:
        This routine measures the pixel scale local to the given position, then assumes that it does not 
        change. So, this routine may only be accurate close to the given position, depending upon the WCS
        projection used.
    
    """
    
    # Pixel distance grid            
    cdef float x0, y0, ra0, dec0, ra1, dec1, xPixScale, yPixScale
    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY, xDistPix, yDistPix
    #cdef np.ndarray[np.float64_t, ndim=2] degreesMap
    
    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    
    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=degreesMap.shape[0]
    X=degreesMap.shape[1]
        
    # Real space map of angular distance in degrees, but only consider values near x0, y0
    #degreesMap=np.ones([Y, X], dtype=np.float64)*1e6 # Allocating this was the bottleneck
    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX > X:
        maxX=X
    if minY < 0:
        minY=0
    if maxY > Y:
        maxY=Y
    for y in range(minY, maxY):
        for x in range(minX, maxX):
            yDeg=(y-y0)*yPixScale
            xDeg=(x-x0)*xPixScale
            degreesMap[y, x]=math.sqrt(xDeg*xDeg+yDeg*yDeg)

    return degreesMap, [minX, maxX], [minY, maxY]

#-------------------------------------------------------------------------------------------------------------
def makeXYDegreesDistanceMaps(np.ndarray[np.float64_t, ndim=2] data, wcs, RADeg, decDeg, maxDistDegrees):
    """Returns an array of distance along x, y axes in degrees from given position. maxDistDegrees sets the 
    limit around the given RADeg, decDeg position to which the distance is calculated.
    
    """
    
    # Pixel distance grid            
    cdef float x0, y0, ra0, dec0, ra1, dec1, xPixScale, yPixScale
    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY, xDistPix, yDistPix
    cdef np.ndarray[np.float64_t, ndim=2] xDegreesMap
    cdef np.ndarray[np.float64_t, ndim=2] yDegreesMap

    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    
    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=data.shape[0]
    X=data.shape[1]
    
    # Real space map of angular distance in degrees, but only consider values near x0, y0
    xDegreesMap=np.ones([Y, X], dtype=np.float64)*1e6
    yDegreesMap=np.ones([Y, X], dtype=np.float64)*1e6
    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX >= X:
        maxX=X-1
    if minY < 0:
        minY=0
    if maxY >= Y:
        maxY=Y-1
    for y in range(minY, maxY):
        for x in range(minX, maxX):
            yDeg=(y-y0)*yPixScale
            xDeg=(x-x0)*xPixScale
            xDegreesMap[y, x]=xDeg
            yDegreesMap[y, x]=yDeg
    
    return [xDegreesMap, yDegreesMap]

#-------------------------------------------------------------------------------------------------------------
def standardDeviationFilter(np.ndarray[np.float64_t, ndim=2] data, rPix):
    """Standard deviation filter using pixels within rPix.

    """

    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY
    cdef np.ndarray[np.float64_t, ndim=2] stdDevFilteredMap
    
    Y=data.shape[0]
    X=data.shape[1]   

    stdDevFilteredMap=np.zeros([Y, X], dtype=np.float64)
    for y in xrange(Y):
        print y
        for x in xrange(X):
            minX=x-rPix
            maxX=x+rPix
            minY=y-rPix
            maxY=y+rPix
            if minY < 0:
                minY=0
            if maxY > Y-1:
                maxY=Y
            if minX < 0:
                minX=0
            if maxX > X-1:
                maxX=X
            pixels=[]
            for yi in range(minY, maxY):
                for xi in range(minX, maxX):
                    pixels.append(data[yi, xi])
            stdDevFilteredMap[y, x]=np.std(pixels)
    
    return stdDevFilteredMap

#-------------------------------------------------------------------------------------------------------------
# Rank filter

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int partition(double[:] arr, int left, int right, int pivotIndex) nogil:
    """For implementation of quickselect in numba.
    
    """
    
    cdef unsigned int storeIndex
    cdef float pivotValue
    
    pivotValue=arr[pivotIndex]        
    arr[pivotIndex], arr[right]=arr[right], arr[pivotIndex]
    storeIndex=left
    for i in xrange(left, right-1):
        if arr[i] < pivotValue:
            arr[storeIndex], arr[i]=arr[i], arr[storeIndex]
            storeIndex=storeIndex+1
    arr[right], arr[storeIndex]=arr[storeIndex], arr[right]
    return storeIndex

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double quickSelect(double[:] arr, int left, int right, int n) nogil:
    """Quickselect, in numba for speed.
    
    """
    
    cdef unsigned int pivotIndex
    #cdef Py_ssize_t pivotIndex
    
    if left == right:        
        return arr[left]
    pivotIndex=(left+right)/2
    pivotIndex=partition(arr, left, right, pivotIndex)
    if n == pivotIndex:
        return arr[n]
    elif n < pivotIndex:
        return quickSelect(arr, left, pivotIndex - 1, n)
    else:
        return quickSelect(arr, pivotIndex + 1, right, n)
             
@cython.boundscheck(False)
@cython.wraparound(False)
def rankFilter(np.ndarray[np.float64_t, ndim=2] arr, float rank, np.ndarray[np.int64_t, ndim=2] annulus):
    """Rank is relative, with 0 lowest, 1 max.
    
    """
    
    #cdef unsigned int i, y, x, sizeY, sizeX, numValues, rankIndex
    cdef Py_ssize_t i, y, x, sizeY, sizeX, numValues, rankIndex 
    cdef np.ndarray[np.int64_t, ndim=1] annY, annX
    cdef np.ndarray[np.float64_t, ndim=1] vals
    cdef np.ndarray[np.float64_t, ndim=2] filtData
    
    sizeY=arr.shape[0]
    sizeX=arr.shape[1]
    
    numValues=int(annulus.sum())
    rankIndex=int(round(rank*numValues))
    annY, annX=np.where(annulus != 0)
    annY=annY-annulus.shape[0]/2
    annX=annX-annulus.shape[1]/2
    filtData=np.zeros([sizeY, sizeX], dtype=np.float64)
    vals=np.zeros(numValues, dtype=np.float64)
    
    for y in xrange(annulus.shape[0]/2, arr.shape[0]-annulus.shape[0]/2):
        for x in xrange(annulus.shape[1]/2, arr.shape[1]-annulus.shape[1]/2):
            for i in xrange(numValues):
                vals[i]=arr[annY[i]+y, annX[i]+x]
            filtData[y, x]=quickSelect(vals, 0, len(vals)-1, rankIndex)
    
    return filtData

#-------------------------------------------------------------------------------------------------------------
# Local SN map
@cython.boundscheck(False)
@cython.wraparound(False)
def makeLocalSNMap(np.ndarray[np.float64_t, ndim=2] arr, np.ndarray[np.int64_t, ndim=2] annulus):
    """Makes local S/N map - measuring noise in an annulus using the standard deviation. Measures
    background using mean in annulus. Since there is no sorting involved, this is much quicker
    than using rank filters. But it may not give quite the same results.
    
    """
    
    cdef Py_ssize_t i, y, x, sizeY, sizeX, numValues
    cdef float valSum, numValuesFloat, mu
    cdef np.ndarray[np.int64_t, ndim=1] annY, annX
    cdef np.ndarray[np.float64_t, ndim=1] vals
    cdef np.ndarray[np.float64_t, ndim=2] noiseFiltData, bckSubData, SNData
    
    sizeY=arr.shape[0]
    sizeX=arr.shape[1]
    
    numValues=int(annulus.sum())
    numValuesFloat=float(numValues)
    #rankIndex=int(round(rank*numValues))
    annY, annX=np.where(annulus != 0)
    annY=annY-annulus.shape[0]/2
    annX=annX-annulus.shape[1]/2
    bckSubData=np.zeros([sizeY, sizeX], dtype=np.float64)
    noiseFiltData=np.zeros([sizeY, sizeX], dtype=np.float64)
    SNData=np.zeros([sizeY, sizeX], dtype=np.float64)
    vals=np.zeros(numValues, dtype=np.float64)
    
    for y in xrange(annulus.shape[0]/2, arr.shape[0]-annulus.shape[0]/2):
        for x in xrange(annulus.shape[1]/2, arr.shape[1]-annulus.shape[1]/2):
            valSum=0.
            for i in xrange(numValues):
                valSum=valSum+arr[annY[i]+y, annX[i]+x]
            mu=valSum/numValuesFloat
            bckSubData[y, x]=arr[y, x]-mu

    for y in xrange(annulus.shape[0]/2, arr.shape[0]-annulus.shape[0]/2):
        for x in xrange(annulus.shape[1]/2, arr.shape[1]-annulus.shape[1]/2):
            valSum=0.
            for i in xrange(numValues):
                valSum=valSum+bckSubData[annY[i]+y, annX[i]+x]
                vals[i]=bckSubData[annY[i]+y, annX[i]+x]
            mu=valSum/numValuesFloat
            valSum=0.
            for i in xrange(numValues):
                valSum=valSum+(vals[i]-mu)*(vals[i]-mu)
            valSum=valSum/numValuesFloat
            noiseFiltData[y, x]=math.sqrt(valSum)
            if noiseFiltData[y, x] > 0:
                SNData[y, x]=bckSubData[y, x]/noiseFiltData[y, x]
            
    return SNData

#------------------------------------------------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
def convolve(np.ndarray[np.float64_t, ndim=2] arr, np.ndarray[np.float64_t, ndim=2] kern):
    """Convolves image arr with kernel kern. Hopefully faster than ndimage.convolve.
    
    Returns 2d array
    
    """

    cdef Py_ssize_t y, x, ky, kx, ty, tx, sizeY, sizeX, kSizeX, kSizeY
    cdef np.ndarray[np.float64_t, ndim=2] filtData
    
    # Element-by-element convolution (wraps at edges)
    sizeY=arr.shape[0]
    sizeX=arr.shape[1]
    kSizeY=kern.shape[0]
    kSizeX=kern.shape[1]
    #kY, kX=np.where(k != 0)
    #kY=kY-k.shape[0]/2
    #kX=kX-k.shape[1]/2
    filtData=np.zeros([sizeY, sizeX], dtype=np.float64)
    for y in xrange(sizeY):
        for x in xrange(sizeX):
            for ky in xrange(kSizeY):
                for kx in xrange(kSizeX):
                    ty=y+ky-kSizeY/2
                    tx=x+kx-kSizeX/2
                    if ty > kSizeY-1:
                        ty=ty-kSizeY
                    if tx > kSizeX-1:
                        tx=tx-kSizeX
                    filtData[y, x]=filtData[y, x]+arr[ty, tx]*kern[ky, kx]
                
    return filtData

