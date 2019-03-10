"""

This module contains tools for handling catalogs, which for us are (generally) astropy Table objects.

"""

from astLib import *
import numpy as np
import operator
import os
import sys
import time
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import IPython

# For adding meta data to output
import datetime
import nemo

#------------------------------------------------------------------------------------------------------------
XMATCH_RADIUS_DEG=1.4/60.0  # catalog matching radius, for sim comparisons

#------------------------------------------------------------------------------------------------------------
# Definitions of table column names and formats that can be written in catalogs
COLUMN_NAMES    = ['name', 
                   'RADeg', 
                   'decDeg', 
                   'SNR', 
                   'numSigPix', 
                   'fractionMapsDetected', 
                   'template', 
                   'galacticLatDeg',
                   'deltaT_c',
                   'err_deltaT_c',
                   'y_c',
                   'err_y_c',
                   'Y500_sr',
                   'err_Y500_sr',
                   'fluxJy',
                   'err_fluxJy',
                   'redshift',
                   'redshiftErr',
                   'ellipse_PA',
                   'ellipse_A',
                   'ellipse_B',
                   'ellipse_x0',
                   'ellipse_y0',
                   'ellipse_e'
                   ]
COLUMN_FORMATS  = ['%s',
                   '%.6f',
                   '%.6f',
                   '%.1f',
                   '%d',
                   '%.2f',
                   '%s',
                   '%.6f',
                   '%.3f',
                   '%.3f',
                   '%.3e',
                   '%.3e',
                   '%.3e',
                   '%.3e',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f',
                   '%.3f'
                   ]

columnsToAdd=[]
formatsToAdd=[]
for k in COLUMN_NAMES:
    if k not in ['name', 'RADeg', 'decDeg', 'galacticLatDeg', 'fractionMapsDetected', 'numSigPix']:
        i=COLUMN_NAMES.index(k)
        formatsToAdd.append(COLUMN_FORMATS[i])
        columnsToAdd.append("fixed_"+k)
COLUMN_NAMES=COLUMN_NAMES+columnsToAdd
COLUMN_FORMATS=COLUMN_FORMATS+formatsToAdd
        
if len(COLUMN_NAMES) != len(COLUMN_FORMATS):
    raise exception("COLUMN_NAMES and COLUMN_FORMATS lists should be same length")
        
#------------------------------------------------------------------------------------------------------------
def makeOptimalCatalog(imageDict, constraintsList = []):
    """Identifies common objects between catalogs in the imageDict and creates a master catalog with
    one entry per object, keeping only the highest S/N detection details.
    
    Args:
        imageDict: dictionary containing filtered maps and associated object catalogs
        constraintsList: an optional list of constraints (for format, see selectFromCatalog)
        
    Returns:
        Nothing - an 'optimalCatalog' key is added to imageDict
    
    """
    
    # Get list of templates - assuming here that all keys that are NOT 'mergedCatalog' are template names
    templates=[]
    for key in imageDict['mapKeys']:
        if key != "mergedCatalog" and key != "optimalCatalog":
            templates.append(key)

    allCatalogs=[]
    for temp in templates:
        if len(imageDict[temp]['catalog']) > 0:
            allCatalogs.append(imageDict[temp]['catalog'])
    if len(allCatalogs) > 0:
        allCatalogs=atpy.vstack(allCatalogs)
        mergedCatalog=allCatalogs.copy()
        mergedCatalog['SNR']=-99.
        mergeRow=0
        usedIndices=[]
        for row in allCatalogs:
            rDeg=astCoords.calcAngSepDeg(row['RADeg'], row['decDeg'], allCatalogs['RADeg'], allCatalogs['decDeg']) 
            xIndices=np.where(rDeg < XMATCH_RADIUS_DEG)[0]
            xMatches=allCatalogs[xIndices]
            xMatchIndex=np.argmax(xMatches['SNR'])
            if xIndices[xMatchIndex] not in usedIndices:
                mergedCatalog[mergeRow]=xMatches[xMatchIndex]
                mergeRow=mergeRow+1
                usedIndices=usedIndices+xIndices.tolist()
        mergedCatalog=mergedCatalog[mergedCatalog['SNR'] > 0]
        mergedCatalog.sort(['RADeg', 'decDeg'])
        mergedCatalog=selectFromCatalog(mergedCatalog, constraintsList)
    else:
        mergedCatalog=[]
    imageDict['optimalCatalog']=mergedCatalog   

#------------------------------------------------------------------------------------------------------------
def catalog2DS9(catalog, outFileName, constraintsList = [], addInfo = [], idKeyToUse = 'name', 
                RAKeyToUse = 'RADeg', decKeyToUse = 'decDeg', color = "cyan", writeNemoInfo = True, coordSys = 'fk5'):
    """Converts a catalog containing object dictionaries into a DS9 region file. 
    
    Args:
        catalog: An astropy Table where each row represents an object.
        outFileName: A file name for the output DS9 region file.
        constraintsList: A list of constraints in the same format as used by `selectFromCatalog`.
        addInfo: A list of dictionaries with keys named `key` and `fmt` (e.g., ``{'key': "SNR", 'fmt': "%.3f"}``).
            These will be added to the object label shown in DS9.
        idKeyToUse: The name of the key in each object dictionary that defines the object's name. Used to 
            label objects in the DS9 region file.
        RAKeyToUse: The name of the key in each object dictionary that contains the RA (decimal degrees) of the
            object.
        decKeyToUse: The name of the key in each object dictionary that contains the declination (decimal 
            degrees) of the object.
        color: The color of the plot symbol used by DS9.
        writeNemoInfo: If True, writes a line with the nemo version and date generated at the top of the 
            DS9 .reg file.
        coordSys: A string defining the coordinate system used for RA, dec, as understood by DS9.
        
    Returns:
        None
    
    """
    
    cutCatalog=selectFromCatalog(catalog, constraintsList) 
        
    with open(outFileName, "w") as outFile:
        timeStamp=datetime.datetime.today().date().isoformat()
        comment="# DS9 region file"
        if writeNemoInfo == True:
            comment=comment+" generated by nemo (version: %s on %s)\n" % (nemo.__version__, timeStamp)
        else:
            comment=comment+"\n"
        outFile.write(comment)
        outFile.write('global dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        for obj in cutCatalog:
            if len(addInfo) > 0:
                infoString=""
                for d in addInfo:
                    if infoString != "":
                        infoString=infoString+" "
                    if obj[d['key']] != None:
                        infoString=infoString+d['fmt'] % (obj[d['key']])
                    else:
                        infoString=infoString+"%s" % (str(obj[d['key']]))
                infoString=" ["+infoString+"]"
            else:
                infoString=""
            if color == 'key':
                colorString=obj['color']
            else:
                colorString=color
            outFile.write("%s;point(%.6f,%.6f) # point=boxcircle color={%s} text={%s%s}\n" \
                        % (coordSys, obj[RAKeyToUse], obj[decKeyToUse], colorString, obj[idKeyToUse], infoString))

#------------------------------------------------------------------------------------------------------------
def makeACTName(RADeg, decDeg, prefix = 'ACT-CL'):
    """Makes ACT cluster name from RADeg, decDeg
    
    """
    
    actName=prefix+" J"+makeRA(RADeg)+makeDec(decDeg)
    
    return actName

#------------------------------------------------------------------------------------------------------------
def makeLongName(RADeg, decDeg, prefix = "ACT-CL"):
    """Makes a long format object name from RADeg, decDeg
    
    """
    
    actName=prefix+" J"+makeLongRA(RADeg)+makeLongDec(decDeg)
    
    return actName
    
#------------------------------------------------------------------------------------------------------------
def makeRA(myRADeg):
    """Makes RA part of ACT names.
    
    """
    hours=(myRADeg/360)*24
    if hours<10:
        sHours="0"+str(hours)[0]
    else:
        sHours=str(hours)[:2]
    
    mins=float(str(hours)[str(hours).index("."):])*60
    if mins<10:
        sMins="0"+str(mins)[:3]
    else:
        sMins=str(mins)[:4]
        
    return (sHours+sMins)#[:-2] # Trims off .x as not used in ACT names
        
#------------------------------------------------------------------------------------------------------------
def makeDec(myDecDeg):
    """Makes dec part of ACT names
    
    """
    
    # Positive
    if myDecDeg>0:
        if myDecDeg<10:
            sDeg="0"+str(myDecDeg)[0]
        else:
            sDeg=str(myDecDeg)[:2]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
        
        return "+"+sDeg+sMins
    else:
        if myDecDeg>-10:
            sDeg="-0"+str(myDecDeg)[1]
        else:
            sDeg=str(myDecDeg)[:3]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
        
        return str(sDeg+sMins)

#-------------------------------------------------------------------------------------------------------------
def makeLongRA(myRADeg):
    """Make a long RA string, i.e. in style of long XCS names
    
    """
    
    hours=(myRADeg/360)*24
    if hours<10:
        sHours="0"+str(hours)[0]
    else:
        sHours=str(hours)[:2]
    
    mins=float(str(hours)[str(hours).index("."):])*60
    if mins<10:
        sMins="0"+str(mins)[0]
    else:
        sMins=str(mins)[:2]
        
    secs=float(str(mins)[str(mins).index("."):])*60
    if secs<10:
        sSecs="0"+str(secs)[:3]
    else:
        sSecs=str(secs)[:4]     
        
    return sHours+sMins+sSecs
        
#-------------------------------------------------------------------------------------------------------------
def makeLongDec(myDecDeg):
    """Make a long dec sting i.e. in style of long XCS names
    
    """
    # Positive
    if myDecDeg>0:
        if myDecDeg<10:
            sDeg="0"+str(myDecDeg)[0]
        else:
            sDeg=str(myDecDeg)[:2]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
            
        secs=float(str(mins)[str(mins).index("."):])*60
        if secs<10:
            sSecs="0"+str(secs)[:3]
        else:
            sSecs=str(secs)[:4]         
        
        return "+"+sDeg+sMins+sSecs
    else:
        if myDecDeg>-10:
            sDeg="-0"+str(myDecDeg)[1]
        else:
            sDeg=str(myDecDeg)[:3]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
            
        secs=float(str(mins)[str(mins).index("."):])*60
        if secs<10:
            sSecs="0"+str(secs)[:3]
        else:
            sSecs=str(secs)[:4]         
        
        return sDeg+sMins+sSecs
        
#-------------------------------------------------------------------------------------------------------------
def selectFromCatalog(catalog, constraintsList):
    """Given a catalog (an astropy Table), return a table of objects matching the
    given constraintsList. Each item in constraintsList is a string in the form:
    
    "key < value", "key > value", etc.
    
    Note that the spaces between key, operator (e.g. '<') and value are essential!
    
    """
            
    passedConstraint=catalog
    for constraintString in constraintsList:
        key, op, value=constraintString.split()
        passedConstraint=passedConstraint[eval("passedConstraint['%s'] %s %s" % (key, op, value))]

    return passedConstraint

#------------------------------------------------------------------------------------------------------------
def catalogListToTab(catalogList, keysToWrite = COLUMN_NAMES):
    """Converts catalog (as a list of dictionaries) into an astropy Table.
    
    Returns astropy Table object
    
    """
    
    availKeys=list(catalogList[0].keys())
    
    # A fudge: we don't know what names y_c_weight keys will have in advance, so they aren't already given
    for key in availKeys:
        if key.find("fixed_y_c_weight") != -1 and key not in keysToWrite:
            keysToWrite.append(key)
    
    tab=atpy.Table()
    for key in keysToWrite:
        if key in availKeys:
            arr=[]
            for obj in catalogList:
                if obj[key] != None:
                    arr.append(obj[key])
                else:
                    arr.append(-99)
            tab.add_column(atpy.Column(arr, key))

    return tab

#------------------------------------------------------------------------------------------------------------
def tabToCatalogList(tab):
    """Converts an astropy Table into a list of dictionaries.

    Returns catalog list
    
    """
    
    catalog=[]
    for row in tab:
        objDict={}
        for k in list(tab.keys()):
            objDict[k]=row[k]
        catalog.append(objDict)
    
    return catalog
    
#------------------------------------------------------------------------------------------------------------
def writeCatalog(catalog, outFileName, constraintsList = []):
    """Writes the catalog (astropy Table) to disk. The output format depends upon the extension given by 
    outFileName - any format supported by the astropy Table object can be used.
    
    NOTE: For .csv format, tab is used as the delimiter. So, to read these using astropy you will need to 
    use something like: tab=atpy.Table().read("test.csv", delimiter = '\t').
    
    constraintsList works as in the selectFromCatalog function.
    
    """

    cutCatalog=selectFromCatalog(catalog, constraintsList)
    if outFileName.split(".")[-1] == 'csv':
        cutCatalog.write(outFileName, format = 'ascii.csv', delimiter = '\t', overwrite = True)
    else:
        cutCatalog.write(outFileName, overwrite = True)

#------------------------------------------------------------------------------------------------------------
def removeDuplicates(tab):
    """Given an astropy Table object, remove duplicate objects - keeping the highest SNR detection for each
    duplicate. This routine is used to clean up output of MPI runs (where we have overlapping tiles). This
    method could be applied elsewhere, if we re-did how we handled catalogs everywhere in nemo.
    
    Returns table with duplicates removed, number of duplicates found, names of duplicates
    
    """
    
    # Find all duplicates
    cat=SkyCoord(ra = tab['RADeg'].data, dec = tab['decDeg'].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat, cat, nthneighbor = 2)
    mask=np.less(rDeg.value, XMATCH_RADIUS_DEG)
    noDupMask=np.greater_equal(rDeg.value, XMATCH_RADIUS_DEG)
    dupTab=tab[mask]
    noDupTab=tab[noDupMask]
    
    # All duplicates removed?
    if mask.sum() == 0:
        return tab, 0, []
    
    # Much faster
    keepMask=np.zeros(len(dupTab), dtype = bool)
    for i in range(len(dupTab)):
        # NOTE: astCoords does not like atpy.Columns sometimes...
        rDeg=astCoords.calcAngSepDeg(dupTab['RADeg'][i], dupTab['decDeg'][i], dupTab['RADeg'].data, dupTab['decDeg'].data)
        mask=np.less(rDeg, XMATCH_RADIUS_DEG)
        indices=np.where(mask == True)[0]
        bestIndex=indices[np.equal(dupTab['SNR'][mask], dupTab['SNR'][mask].max())][0]
        keepMask[bestIndex]=True
    keepTab=dupTab[keepMask]
    
    keepTab=atpy.vstack([keepTab, noDupTab])
    keepTab.sort('RADeg')
    
    return keepTab, len(dupTab), dupTab['name']
    
