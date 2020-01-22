"""

Cross-match cluster catalogs made from individual maps in order to check calibration.

"""

import os
import sys
import glob
import numpy as np
import astropy.table as atpy
from nemo import catalogs
from nemo import plotSettings
from collections import OrderedDict as odict
import pylab as plt
import IPython

#------------------------------------------------------------------------------------------------------------
# Main

# Pick one
#relativeTo='median'
#relativeTo='absolute'
#relativeTo='S16'
relativeTo='S18'

# This only makes sense if relativeTo = 'median', 'S16, or 'S18'
#mode='residualDiff'
mode='ratio'

# Applied to reference catalog only (since now doing forced photometry)
refSNRCut=10.0
print("relativeTo = %s; refSNRCut = %.1f" % (relativeTo, refSNRCut))

plotsDir="plots_%s_%s_%.1f" % (mode, relativeTo, refSNRCut)
os.makedirs(plotsDir, exist_ok = True)

# For testing: multiply ycRef by this
calFactor=1.00
if calFactor != 1.0:
    print("WARNING: calFactor set to %.2f - multiplying ycRef by this factor." % (calFactor))

# Reference catalog - used for cross-match positions across other catalogs
refTab=atpy.Table().read("../MFMF_S18_auto/MFMF_S18_auto_M500.fits")
if relativeTo == 'S16':
    refTab=atpy.Table().read("../MFMF_S16_auto/MFMF_S16_auto_M500.fits")
keepCols=['name', 'RADeg', 'decDeg', 'fixed_SNR', 'fixed_y_c', 'fixed_err_y_c']
removeCols=[]
for k in refTab.keys():
    if k not in keepCols:
        removeCols.append(k)
refTab.remove_columns(removeCols)
refTab=refTab[refTab['fixed_SNR'] > refSNRCut]

# Dictionaries in which to store all cross-matched y_c values
# Each entry will be a list
ycRef={}
ycErrRef={}
yc={}
ycErr={}
labels={}

# Cross match against nemo output
print("Collecting fixed_y_c measurements for each cluster across all maps")
xMatchFiles=glob.glob("outputCatalogs/*.fits")
xMatchFiles.sort()
for f in xMatchFiles:    
    label=os.path.split(f)[-1].split("_optimal")[0]
    print("    %s" % (label))
    tab=atpy.Table().read(f)
    #tab=tab[tab['fixed_SNR'] > SNRCut]
    try:
        refMatched, tabMatched, sep=catalogs.crossMatch(refTab, tab, radiusArcmin = 1.4)
    except:
        raise Exception("Matching probably failed because SNRCut is too low so there were no matches")
    
    for ref, match in zip(refMatched, tabMatched):
        name=ref['name']
        if name not in yc.keys():
            yc[name]=[]
            ycErr[name]=[]
            labels[name]=[]
            ycRef[name]=ref['fixed_y_c']*calFactor
            ycErrRef[name]=ref['fixed_err_y_c']
        yc[name].append(match['fixed_y_c'])
        ycErr[name].append(match['fixed_err_y_c'])
        labels[name].append(label)

nameList=list(yc.keys()); nameList.sort()
for name in nameList:
    yc[name]=np.array(yc[name])
    ycErr[name]=np.array(ycErr[name])

# Check scaled residuals relative to median
# Make a plot for each cluster with a fair number of data points
# Gather a set so we know the typical offset for each map
# We can also make a map of typical offset for each cluster
minPoints=10
resByMap={}
nameList=list(yc.keys()); nameList.sort()
print("Making plots of fixed_y_c for each cluster")
for name in nameList:
    
    # Plot
    if len(yc[name]) < minPoints:
        continue
    if relativeTo == 'absolute':
        res=yc[name]
        resSigma=ycErr[name]
        ylim=None
        ylabel="$\\tilde{y_0}\, (10^{-4})$"
    else:
        if mode == 'relativeDiff':
            if relativeTo == 'median':
                res=(yc[name]-np.median(yc[name]))/ycErr[name]
            else:
                # Relative to S16 or S18 act+planck co-add
                res=(np.array(yc[name])-ycRef[name])/np.sqrt(np.array(ycErr[name])**2+ycErrRef[name]**2)
            resSigma=[1]*len(res)
            ylabel="$\Delta \\tilde{y_0} (\sigma)$"
            ylim=(-4, 4)
        elif mode == 'ratio':
            if relativeTo == 'median':
                res=yc[name]/np.median(yc[name])
                resSigma=(ycErr[name]/yc[name])*res
            else:
                res=np.array(yc[name])/np.array(ycRef[name])
                resSigma=np.sqrt(np.power(ycErr[name]/yc[name], 2)+np.power(ycErrRef[name]/ycRef[name], 2))*res
            ylim=(0, 2)
            ylabel="$\\tilde{y_0} / \\tilde{y_0}$ [%s]" % (relativeTo)
        #medAbsRes=np.median(abs(res))
    
    plotSettings.update_rcParams()
    plt.figure(figsize=(18, 10))
    ax=plt.axes([0.1, 0.4, 0.87, 0.5])
    x=np.arange(1, len(labels[name])+1)
    titleStr=name
    if relativeTo != 'absolute':
        if mode == 'relativeDiff':
            plt.plot(np.linspace(0, len(x)+1, 3), [0]*3, 'k--')
        elif mode == 'ratio':
            plt.plot(np.linspace(0, len(x)+1, 3), [1]*3, 'k--') 
    else:
        # plot median absolute value
        ycMedian=np.median(yc[name])
        plt.plot(np.linspace(0, len(x)+1, 3), [ycMedian]*3, 'k--')
        titleStr=titleStr+" - median %s = %.3f" % (ylabel, ycMedian)
    plt.errorbar(x, res, yerr = [resSigma, resSigma], ecolor = '#AAAAFF', elinewidth = 1.5, fmt = 'D', ms = 6, label = name)
    plt.title(titleStr)#+":   median $\left| \Delta \\tilde{y_0} \\right|$ = %.1f $\sigma$" % (medAbsRes))
    plt.xlim(0, len(labels[name])+1)
    plt.xticks(x, list(labels[name]), rotation = 90)
    plt.ylabel(ylabel)
    if ylim is not None:
        plt.ylim(ylim)
    outPDFName=plotsDir+os.path.sep+name.replace(" ", "_")+".pdf"
    plt.savefig(outPDFName)
    plt.savefig(outPDFName.replace(".pdf", ".png"))
    plt.close()
    print("    %s" % (name))#: median | residual (y0 - median_y_0) | = %.1f sigma" % (name, medAbsRes))
    
    # Gather residuals by map
    for resInMap, mapName in zip(res, labels[name]):
        if mapName not in resByMap:
            resByMap[mapName]=[]
        resByMap[mapName].append(resInMap)

# Plot of mean residual (or ratio) per map +/- std error?
if relativeTo == 'absolute':
    sys.exit()

minPoints=5
medianResByMap=[]
meanResByMap=[]
stdResByMap=[]
stdErrByMap=[]
mapNames=[]
allMapNames=list(resByMap.keys())
allMapNames.sort()
for mapName in allMapNames:
    if len(resByMap[mapName]) < minPoints:
        continue
    medianResByMap.append(np.median(resByMap[mapName]))
    meanResByMap.append(np.mean(resByMap[mapName]))
    stdResByMap.append(np.std(resByMap[mapName]))
    stdErrByMap.append(np.std(resByMap[mapName])/np.sqrt(len(resByMap[mapName])))
    mapNames.append(mapName)
meanResByMap=np.array(meanResByMap)
stdResByMap=np.array(stdResByMap)
stdErrByMap=np.array(stdErrByMap)

#errorBar=stdResByMap
#errorBar=stdErrByMap

plotSettings.update_rcParams()
plt.figure(figsize=(18, 10))
ax=plt.axes([0.1, 0.4, 0.87, 0.5])
x=np.arange(1, len(mapNames)+1)
#plt.title("median $\left| \\rm{mean}\,\Delta \\tilde{y_0} (\sigma) \\right|$ = %.1f" % (np.median(abs(np.array(meanResByMap)))))
titleStr="mean %s by map\nmedian-of-means = %.2f" % (ylabel, np.median(meanResByMap))
if calFactor != 1:
    titleStr=titleStr+" (calFactor = %.3f applied)" % (calFactor)
plt.title(titleStr)
if relativeTo != 'absolute':
    if mode == 'relativeDiff':
        plt.plot(np.linspace(0, len(x)+1, 3), [0]*3, 'k--')
    elif mode == 'ratio':
        plt.plot(np.linspace(0, len(x)+1, 3), [1]*3, 'k--') 
for name, xi in zip(mapNames, x):
    plt.plot([xi]*len(resByMap[name]), resByMap[name], '.', color = '#AAAAAA', zorder = 1)
#ecolor = '#AAAAFF'
plt.errorbar(x, meanResByMap, yerr = [stdResByMap, stdResByMap], ecolor = '#AAAAFF', elinewidth = 3, fmt = 'D', ms = 1,
             zorder = 900)
plt.errorbar(x, meanResByMap, yerr = [stdErrByMap, stdErrByMap], elinewidth = 2, fmt = 'D', ms = 6,
             zorder = 1000)
plt.xlim(0, len(mapNames)+1)
plt.xticks(x, list(mapNames), rotation = 90)
plt.ylabel(ylabel)  # from above...
plt.ylim(ylim)      # from above...
outPDFName=plotsDir+os.path.sep+"meansByMap.pdf"
plt.savefig(outPDFName)
plt.savefig(outPDFName.replace(".pdf", ".png"))
plt.close()

# Text output
outFileName=plotsDir+os.path.sep+"meansByMap.txt"
with open(outFileName, "w") as outFile:
    outFile.write("#mapName\tmean\tsigma\tstdErr\tN\n")
    for mapName, mean, sigma, stdErr in zip(mapNames, meanResByMap, stdResByMap, stdErrByMap):
        outFile.write("%s\t%.3f\t%.3f\t%.3f\t%d\n" % (mapName, mean, sigma, stdErr, len(resByMap[mapName])))
# Wiki table
outFileName=plotsDir+os.path.sep+"wikiMeansByMap.txt"
with open(outFileName, "w") as outFile:
    outFile.write("||border=1 width=80%\n||!mapName ||!mean ||!sigma ||!stdErr ||!N ||\n")
    for mapName, mean, sigma, stdErr in zip(mapNames, meanResByMap, stdResByMap, stdErrByMap):
        outFile.write("|| %s || %.3f || %.3f || %.3f || %d ||\n" % (mapName, mean, sigma, stdErr, len(resByMap[mapName])))

# Individual map histograms on one plot
outPDFName=plotsDir+os.path.sep+"distsByMap.pdf"
binEdges=np.linspace(0, 3, 31) 
numCols=9
numRows=np.ceil(len(mapNames)/numCols)
plotCount=1
plt.figure(figsize=(36, numRows*6))
plt.subplots_adjust(left = 0.02, right = 0.99, bottom = 0.02, top = 0.95, 
                    wspace = 0.2, hspace = 0.2)
for name, mean in zip(mapNames, meanResByMap):
    plt.subplot(numRows, numCols, plotCount)
    plt.hist(resByMap[name], bins = binEdges, density = True)
    yMin, yMax=plt.ylim()
    plt.plot([1]*3, np.linspace(0, yMax, 3), '-')
    plt.plot([mean]*3, np.linspace(0, yMax, 3), ':', label = "%.3f" % (mean))
    plt.ylim(0, yMax)
    #if plotCount > plotCount-numRows:
        #plt.xlabel(ylabel)
    plt.xlim(0, 3)
    plt.legend(loc = 'upper right', fontsize = 8)
    plt.xticks([],[])
    plt.yticks([],[])
    plt.title(name, fontdict = {'size': 9})
    plotCount=plotCount+1
plt.savefig(outPDFName)
plt.savefig(outPDFName.replace(".pdf", ".png"))
plt.close()
