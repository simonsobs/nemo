"""

This module contains tools for handling catalogs, which for us are lists of dictionaries.

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
                   'err_fluxJy']
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
                   '%.3f']

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
def mergeCatalogs(imageDict):
    """Identifies common objects between catalogs in the imageDict and creates a master catalog with
    one entry per object, but multiple flux measurements where available.
    
    """
    
    # Properties (apart from position, name info) that we keep for each object, if they are present
    # Fix this properly later
    #wantedKeys=['numSigPix', 'flux_arcmin2', 'fluxErr_arcmin2', 'SNR', 'fluxRadius_arcmin', 
                #'template', 'fluxStatus', 'deltaT_c', 'err_deltaT_c', 'y_c', 'err_y_c', 'Y500_sr', 'err_Y500_sr',
                #'fluxJy', 'err_fluxJy']
    #keysToAdd=[]
    #for k in wantedKeys:
        #keysToAdd.append("fixed_"+k)
    #wantedKeys=wantedKeys+keysToAdd
    
    wantedKeys=[]
    for k in COLUMN_NAMES:
        if k not in ['name', 'RADeg', 'decDeg', 'galacticLatDeg']:
            wantedKeys.append(k)
    
    # Get list of templates - assuming here that all keys that are NOT 'mergedCatalog' are template names
    templates=[]
    for key in imageDict['mapKeys']:
        if key != "mergedCatalog" and key != "optimalCatalog":
            templates.append(key)
            
    # Note we always take the name/position for merges from the first entry in the catalog, which may not 
    # always be the best thing to do
    mergedCatalog=[]
    for temp in templates:
        catalog=imageDict[temp]['catalog']        
        for c in catalog:
            cra=c['RADeg']
            cdec=c['decDeg']
            rMin=1e6
            bestMatch=None
            
            # For faster cross matching
            mRAs=[]
            mDecs=[]
            for m in mergedCatalog:
                mRAs.append(m['RADeg'])
                mDecs.append(m['decDeg'])
            mRAs=np.array(mRAs)
            mDecs=np.array(mDecs)
            if mRAs.shape[0] > 0:
                rs=astCoords.calcAngSepDeg(cra, cdec, mRAs, mDecs)
                rMin=rs.min()
                rMinIndex=np.equal(rs, rMin).nonzero()[0][0]
                bestMatch=mergedCatalog[rMinIndex]
            else:
                bestMatch=None
            
            if bestMatch != None and rMin < XMATCH_RADIUS_DEG:
                for key in wantedKeys:
                    if key in list(bestMatch.keys()) and key in list(c.keys()):
                        bestMatch[key].append(c[key])   
            else:
                # Must be an object not already in list
                nonListKeys=['name', 'RADeg', 'decDeg', 'galacticLatDeg']
                newObj={}
                for key in nonListKeys:
                    if key in list(c.keys()):
                        newObj[key]=c[key]
                for key in wantedKeys:
                    if key in list(c.keys()):
                        newObj[key]=[c[key]]
                mergedCatalog.append(newObj)
                        
    imageDict['mergedCatalog']=mergedCatalog
        
#------------------------------------------------------------------------------------------------------------
def makeOptimalCatalog(imageDict, constraintsList):
    """Identifies common objects between catalogs in the imageDict and creates a master catalog with
    one entry per object, keeping only the highest S/N detection details.
    
    """
    
    # Properties (apart from position, name info) that we keep for each object, if they are present
    # Fix this properly later - duplicated for mergeCatalogs also (urgh)
    #wantedKeys=['name', 'RADeg', 'decDeg', 'galacticLatDeg']
    #wantedKeys=wantedKeys+['numSigPix', 'flux_arcmin2', 'fluxErr_arcmin2', 'SNR', 'fluxRadius_arcmin', 
                           #'template', 'fluxStatus', 'deltaT_c', 'err_deltaT_c', 'y_c', 'err_y_c', 'Y500_sr', 'err_Y500_sr',
                           #'fluxJy', 'err_fluxJy']
    #keysToAdd=[]
    #for k in wantedKeys:
        #if k not in ['name', 'RADeg', 'decDeg', 'galacticLatDeg']:
            #keysToAdd.append("fixed_"+k)
    #wantedKeys=wantedKeys+keysToAdd
    
    wantedKeys=COLUMN_NAMES
    
    # Get list of templates - assuming here that all keys that are NOT 'mergedCatalog' are template names
    templates=[]
    for key in imageDict['mapKeys']:
        if key != "mergedCatalog" and key != "optimalCatalog":
            templates.append(key)
            
    # Note we always take the name/position for merges from the first entry in the catalog, which may not 
    # always be the best thing to do
    mergedCatalog=[]
    for temp in templates:
        catalog=imageDict[temp]['catalog']        
        for c in catalog:
            cra=c['RADeg']
            cdec=c['decDeg']
            rMin=1e6
            bestMatch=None
            for m in mergedCatalog:
                mra=m['RADeg']
                mdec=m['decDeg']
                r=astCoords.calcAngSepDeg(mra, mdec, cra, cdec)
                if r < rMin:
                    rMin=r
                    bestMatch=m
            if bestMatch != None and rMin < XMATCH_RADIUS_DEG:
                # Is this better than the current object?
                if c['SNR'] > bestMatch['SNR']:
                    for key in wantedKeys:
                        if key in list(c.keys()):
                            bestMatch[key]=c[key]
            else:
                # Must be an object not already in list
                newObj={}
                for key in wantedKeys:
                    if key in list(c.keys()):
                        newObj[key]=c[key]
                mergedCatalog.append(newObj)
                    
    # Now cut...
    mergedCatalog=selectFromCatalog(mergedCatalog, constraintsList) 
    
    # Sort by dec, RA
    decSorted=sorted(mergedCatalog, key=operator.itemgetter('decDeg'))
    RASorted=sorted(decSorted, key=operator.itemgetter('RADeg'))
    
    imageDict['optimalCatalog']=RASorted

#------------------------------------------------------------------------------------------------------------
def catalog2DS9(catalog, outFileName, constraintsList = [], addInfo = [], idKeyToUse = 'name', 
                RAKeyToUse = 'RADeg', decKeyToUse = 'decDeg', color = "cyan", writeNemoInfo = True, coordSys = 'fk5'):
    """Converts a catalog containing object dictionaries into a ds9 region file. Objects will be labelled
    in the .reg file according to the idKeyToUse.
    
    If color == 'key', then use 'color' key in object dictionary to set color.
    
    constraintsList works the same way as selectFromCatalog function
    
    """

    # Cut catalog according to constraints
    cutCatalog=selectFromCatalog(catalog, constraintsList) 
    
    outFile=open(outFileName, "w")
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
    outFile.close()

#------------------------------------------------------------------------------------------------------------
def matchCatalogs(catalog1, catalog2, matchKey = 'catalog2Match', xMatchRadiusDeg = 1.0/3600.0):
    """Generic catalog cross matching based on position. Objects in catalog2 that are associated with catalog1
    are added as a dictionary to the corresponding catalog1 object at key matchKey. Doesn't check for things
    which may be spuriously associated at angular distances > 90 deg (this is a flaw in the calcAngSepDeg
    routine in astLib), should add that.
    
    """
    
    # Cross match
    c2RAs=[]
    c2Decs=[]
    for c2 in catalog2:
        c2RAs.append(c2['RADeg'])
        c2Decs.append(c2['decDeg'])
    c2RAs=np.array(c2RAs)
    c2Decs=np.array(c2Decs)
    for c1 in catalog1:
        rMin=1e6
        bestMatch=None        
        if c2RAs.shape[0] > 0:
            rs=astCoords.calcAngSepDeg(c1['RADeg'], c1['decDeg'], c2RAs, c2Decs)
            rMin=rs.min()
            rMinIndex=np.equal(rs, rMin).nonzero()[0][0]
            if rMin < xMatchRadiusDeg:
                deltaRA=abs(c1['RADeg']-catalog2[rMinIndex]['RADeg'])
                deltaDec=abs(c2['decDeg']-catalog2[rMinIndex]['decDeg'])
                if deltaRA < 10 and deltaDec < 10:
                    bestMatch=catalog2[rMinIndex]        
        c1[matchKey]=bestMatch

#------------------------------------------------------------------------------------------------------------
def flagCatalogMatches(catalog, flagCatalog, key, matchRadiusDeg = 2.0/60.0):
    """Flags objects in catalog that have a match in flagCatalog, by adding a key with the given name to
    all dictionaries in catalog. If the object is matched, the value is True.
        
    """
    
    ras=[]
    decs=[]
    for sobj in flagCatalog:
        ras.append(sobj['RADeg'])
        decs.append(sobj['decDeg'])
    ras=np.array(ras)
    decs=np.array(decs)
    matchRadiusDeg=2.0/60.0
    for obj in catalog:
        obj[key]=False
        rDeg=astCoords.calcAngSepDeg(obj['RADeg'], obj['decDeg'], ras, decs)
        rMin=rDeg.min()
        if rMin < matchRadiusDeg:
            rMinIndex=rDeg.tolist().index(rMin)
            deltaRA=abs(obj['RADeg']-ras[rMinIndex])
            deltaDec=abs(obj['decDeg']-decs[rMinIndex])
            if deltaRA < 10.0 and deltaDec < 10.0:
                obj[key]=True
                if 'name' in list(flagCatalog[rMinIndex].keys()):
                    obj[key+" name"]=flagCatalog[rMinIndex]['name']
    
    return catalog

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
    """Given a catalog (list of dictionaries representing objects), return a list of objects matching the
    given constraintsList. Each item in constraintsList is a string in the form:
    
    "key < value", "key > value", etc.
    
    Note that the spaces between key, operator (e.g. '<') and value are essential!
    
    """
    
    passedConstraint=catalog
    for constraintString in constraintsList:
        lastCatalog=passedConstraint
        passedConstraint=[]
        for obj in lastCatalog:         
            key, op, value=constraintString.split()
            if eval("obj['%s'] %s %s" % (key, op, value)) == True:
                passedConstraint.append(obj)
    
    return passedConstraint

#------------------------------------------------------------------------------------------------------------
def catalogToTab(catalog, keysToWrite, keyFormats, constraintsList):
    """Converts a nemo catalog (list of dictionaries) into astropy.table format.
    
    constraintsList works as in the selectFromCatalog function.

    Returns astropy.table object
    
    """
    
    cutCatalog=selectFromCatalog(catalog, constraintsList)                                           
    availKeys=list(cutCatalog[0].keys())
    
    # Write a .fits version (easier for topcatting)
    # NOTE: switched to astropy (v1.3) tables interface
    tab=atpy.Table()
    for key in keysToWrite:
        if key in availKeys:
            arr=[]
            for obj in cutCatalog:
                if obj[key] != None:
                    arr.append(obj[key])
                else:
                    arr.append(-99)
            tab.add_column(atpy.Column(arr, key))
    
    return tab

#------------------------------------------------------------------------------------------------------------
def tabToCatalog(tab):
    """Converts an astropy.table to a nemo catalog (list of dictionaries).

    Returns catalog
    
    """
    
    catalog=[]
    for row in tab:
        objDict={}
        for k in list(tab.keys()):
            objDict[k]=row[k]
        catalog.append(objDict)
    
    return catalog

#------------------------------------------------------------------------------------------------------------
def writeTab(tab, outFileName):
    """Writes an astropy.table object to disk. The file format is determined by the extension of outFileName.
    
    """
    
    if os.path.exists(outFileName) == True:
        os.remove(outFileName)
    tab.write(outFileName)

#------------------------------------------------------------------------------------------------------------
def writeCatalogFromTab(tab, outFileName, keysToWrite, keyFormats, constraintsList, headings = True, 
                        writeNemoInfo = True, extraHeaderText = None):
    """Given an astropy.table (rather than a catalog, i.e., list of dictionaries), write a .csv (actually
    tab-delimited), with meta data at the top, in the style of writeCatalog.
    
    This is provided so that the format of the output provided by nemo doesn't change, even though in places
    we have switched to astropy.table for data storage. We may get rid of this eventually...
    
    """
    
    catalog=tabToCatalog(tab)
    writeCatalog(catalog, outFileName, keysToWrite, keyFormats, constraintsList, headings = headings,
                 writeNemoInfo = writeNemoInfo, extraHeaderText = extraHeaderText)
    
#------------------------------------------------------------------------------------------------------------
def writeCatalog(catalog, outFileName, keysToWrite, keyFormats, constraintsList, headings = True, 
                 writeNemoInfo = True, extraHeaderText = None):
    """Dumps the merged catalog to a .csv.
    
    constraintsList works as in the selectFromCatalog function.
    
    NOTE: Now writing a .fits table too.
    
    NOTE: Not using this to write final nemo output anymore - as we need to take care of duplicates if 
    running under MPI with maps that overlap. But this is still used by other routines (e.g., individual
    catalogs output by the routines in the photometry module).
    
    """
        
    # Cut catalog according to constraints
    cutCatalog=selectFromCatalog(catalog, constraintsList)                                           
    availKeys=list(cutCatalog[0].keys())
    
    outFile=open(outFileName, "w")
    
    # Add meta data
    timeStamp=datetime.datetime.today().date().isoformat()
    if writeNemoInfo == True:
        outFile.write("# Output by nemo (version: %s)\n" % (nemo.__version__))
    outFile.write("# Date generated: %s\n" % (timeStamp))
    if extraHeaderText != None:
        outFile.write(extraHeaderText)
    if headings == True:
        heading="# "
        for k in keysToWrite:
            if heading != "# ":
                heading=heading+"\t"
            if k in availKeys:
                heading=heading+k
        outFile.write(heading+"\tnotes\n")
        
    count=0
    for obj in cutCatalog:
        count=count+1
        obj['id']=count
        line=""
        for k, f in zip(keysToWrite, keyFormats):
            if k in availKeys:
                if line != "":
                    line=line+"\t"
                if type(obj[k]) == list or type(obj[k]) == np.ndarray:
                    if obj[k][0] != None:   # merged cat, just take first item for now
                        line=line+f % (obj[k][0])   
                    else:
                        line=line+str(None)
                else:
                    if obj[k] != None:
                        try:
                            line=line+f % (obj[k]) 
                        except:
                            print("Argh!")
                            ipshell()
                            sys.exit()
                    else:
                        line=line+str(None)
        # Add on a 'notes' column - any key which is just a bool gets added to a , delimited list if True
        notes=""
        for key in list(obj.keys()):
            if type(obj[key]) == bool and obj[key] == True:
                if notes != "":
                    notes=notes+","
                notes=notes+key
        outFile.write(line+"\t"+notes+"\n")
    outFile.close()   
    
    # Write a .fits version (easier for topcatting)
    # NOTE: switched to astropy (v1.3) tables interface
    tab=catalogToTab(catalog, keysToWrite, keyFormats, constraintsList)
    writeTab(tab, outFileName.replace(".csv", ".fits"))

#------------------------------------------------------------------------------------------------------------
def readCatalog(fileName):
    """Reads a nemo .csv catalog, returns a list of dictionaries
    
    """
    
    # We'll grab stuff from the heading line, and we'll understand what to do with it if it appears in this
    typesDict={'name': 'str', 
               'RADeg': 'float', 
               'decDeg': 'float', 
               'SNR': 'float', 
               'fractionMapsDetected': 'float',
               'template': 'str',
               'NED_name': 'str',
               'NED_z': 'float',
               'NED_distArcmin': 'float',
               'NED_RADeg': 'float',
               'NED_decDeg': 'float',
               'deltaT': 'float',
               'deltaT_c': 'float',
               'y_c': 'float',
               'scaleArcmin': 'float',
               'YArcmin2': 'float',
               'reference': 'str',
               'notes': 'str',
               'comment': 'str',
               'id': 'int',
               'logMStar': 'float',
               'z': 'float'
              }
    
    inFile=open(fileName, "r")
    lines=inFile.readlines()
    inFile.close()
    
    keysList=None
    catalog=[]
    for line in lines:
        if line[0] == '#' and line.find("\t") != -1:
            # Scan for column headings - really we should just put in a proper header or save as .fits tables 
            # instead, this could easily break
            keysList=line.lstrip("# ").rstrip("\n").split("\t")
        if line[0] != '#' and len(line) > 3:
            if keysList == None:
                raise Exception("couldn't find column headings in catalog file %s" % (fileName))
            objDict={}
            bits=line.rstrip("\n").split("\t")
            for key, bit in zip(keysList, bits):
                if key in typesDict: # automatically skips what we don't understand (add to typesDict if needed)
                    if bit == 'None':
                        objDict[key]=None
                    elif typesDict[key] == 'str':
                        objDict[key]=eval(typesDict[key]+"('"+bit+"')")
                    else:
                        objDict[key]=eval(typesDict[key]+"("+bit+")")
            catalog.append(objDict)
            
    return catalog

#------------------------------------------------------------------------------------------------------------
def removeDuplicatesFromTab(tab):
    """Given an astropy.table object, remove duplicate objects - keeping the highest SNR detection for each
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
    
    # Keep highest SNR of all pairs
    # NOTE: This is correct, but very slow
    #keepMask=np.zeros(len(dupTab), dtype = bool)
    #for i in range(len(dupTab)):
        #iCoord=SkyCoord(ra = dupTab['RADeg'][i], dec = dupTab['decDeg'][i], unit = 'deg')
        #bestSNR=dupTab['SNR'][i]
        #bestIndex=i
        #for j in range(len(dupTab)):
            #if i != j:
                #jCoord=SkyCoord(ra = dupTab['RADeg'][j], dec = dupTab['decDeg'][j], unit = 'deg')
                #if iCoord.separation(jCoord).value < XMATCH_RADIUS_DEG and dupTab['SNR'][j] > bestSNR:
                    #bestSNR=dupTab['SNR'][j]
                    #bestIndex=j
        #keepMask[bestIndex]=True
    #keepTab=dupTab[keepMask]

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
    
