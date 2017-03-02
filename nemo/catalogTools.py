# -*- coding: utf-8 -*-
"""This module contains tools for handling catalogs, which for us are lists of dictionaries.

"""

from astLib import *
import numpy
import operator
import os
import urllib
import urllib2
import sys
import time
import astropy.table as atpy
import IPython

# For adding meta data to output
import datetime
import nemo

#------------------------------------------------------------------------------------------------------------
XMATCH_RADIUS_DEG=1.4/60.0  # catalog matching radius, for sim comparisons

#------------------------------------------------------------------------------------------------------------
def mergeCatalogs(imageDict):
    """Identifies common objects between catalogs in the imageDict and creates a master catalog with
    one entry per object, but multiple flux measurements where available.
    
    """
       
    # Get list of templates - assuming here that all keys that are NOT 'mergedCatalog' are template names
    templates=[]
    for key in imageDict.keys():
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
            mRAs=numpy.array(mRAs)
            mDecs=numpy.array(mDecs)
            if mRAs.shape[0] > 0:
                rs=astCoords.calcAngSepDeg(cra, cdec, mRAs, mDecs)
                rMin=rs.min()
                rMinIndex=numpy.equal(rs, rMin).nonzero()[0][0]
                bestMatch=mergedCatalog[rMinIndex]
            else:
                bestMatch=None
            
            if bestMatch != None and rMin < XMATCH_RADIUS_DEG:
                keysToAppend=['numSigPix', 'flux_arcmin2', 'fluxErr_arcmin2', 'SNR', \
                              'fluxRadius_arcmin', 'template', 'fluxStatus', 'deltaT_c', 'y_c', 
                              'Y500_sr', 'err_Y500Err']
                for key in keysToAppend:
                    if key in bestMatch.keys():
                        bestMatch[key].append(c[key])   
            else:
                # Must be an object not already in list
                nonListKeys=['name', 'RADeg', 'decDeg', 'galacticLatDeg']
                listKeys=['numSigPix', 'flux_arcmin2', 'fluxErr_arcmin2', 'SNR', 'fluxRadius_arcmin', 
                          'template', 'fluxStatus', 'deltaT_c', 'y_c', 'Y500_sr', 'err_Y500_sr']
                newObj={}
                for key in nonListKeys:
                    if key in c.keys():
                        newObj[key]=c[key]
                for key in listKeys:
                    if key in c.keys():
                        newObj[key]=[c[key]]
                mergedCatalog.append(newObj)
                
    # Add a flag to each object indicating how many filtered maps it was identified in
    for obj in mergedCatalog:
        if 'template' in obj.keys():
            obj['fractionMapsDetected']=float(len(obj['template']))/float(len(templates))
        
    imageDict['mergedCatalog']=mergedCatalog

#------------------------------------------------------------------------------------------------------------
def makeOptimalCatalog(imageDict, constraintsList):
    """Identifies common objects between catalogs in the imageDict and creates a master catalog with
    one entry per object, keeping only the highest S/N detection details.
    
    """
    
    # Get list of templates - assuming here that all keys that are NOT 'mergedCatalog' are template names
    templates=[]
    for key in imageDict.keys():
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
                    wantedKeys=['name', 'RADeg', 'decDeg', 'galacticLatDeg', 'numSigPix', 'flux_arcmin2',
                                'fluxErr_arcmin2', 'SNR', 'fluxRadius_arcmin', 'template', 'fluxStatus',
                                'deltaT_c', 'y_c', 'Y500_sr', 'err_Y500_sr']
                    for key in wantedKeys:
                        if key in c.keys():
                            bestMatch[key]=c[key]
                    bestMatch['fractionMapsDetected']=bestMatch['fractionMapsDetected']+1.0
                else:
                    bestMatch['fractionMapsDetected']=bestMatch['fractionMapsDetected']+1.0
            else:
                # Must be an object not already in list
                wantedKeys=['name', 'RADeg', 'decDeg', 'galacticLatDeg', 'numSigPix', 'flux_arcmin2',
                            'fluxErr_arcmin2', 'SNR', 'fluxRadius_arcmin', 'template', 'fluxStatus',
                            'deltaT_c', 'y_c', 'Y500_sr', 'err_Y500_sr']
                newObj={}
                for key in wantedKeys:
                    if key in c.keys():
                        newObj[key]=c[key]
                newObj['fractionMapsDetected']=1.0  # We'll turn this into a fraction after done matching
                mergedCatalog.append(newObj)
                
    # Add a flag to each object indicating how many filtered maps it was identified in
    for obj in mergedCatalog:
        obj['fractionMapsDetected']=obj['fractionMapsDetected']/float(len(templates))
    
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
    
    outFile=file(outFileName, "w")
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
    c2RAs=numpy.array(c2RAs)
    c2Decs=numpy.array(c2Decs)
    for c1 in catalog1:
        rMin=1e6
        bestMatch=None        
        if c2RAs.shape[0] > 0:
            rs=astCoords.calcAngSepDeg(c1['RADeg'], c1['decDeg'], c2RAs, c2Decs)
            rMin=rs.min()
            rMinIndex=numpy.equal(rs, rMin).nonzero()[0][0]
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
    ras=numpy.array(ras)
    decs=numpy.array(decs)
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
                if 'name' in flagCatalog[rMinIndex].keys():
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
def writeCatalog(catalog, outFileName, keysToWrite, keyFormats, constraintsList, headings = True, 
                 writeNemoInfo = True, extraHeaderText = None):
    """Dumps the merged catalog to a .csv, for now this is only names and object positions.
    
    constraintsList works as in the selectFromCatalog function.
    
    NOTE: Now writing a .fits table too.
        
    """
    
    # Cut catalog according to constraints
    cutCatalog=selectFromCatalog(catalog, constraintsList)                                           
    availKeys=cutCatalog[0].keys()
    
    outFile=file(outFileName, "w")
    
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
                if type(obj[k]) == list or type(obj[k]) == numpy.ndarray:
                    if obj[k][0] != None:   # merged cat, just take first item for now
                        line=line+f % (obj[k][0])   
                    else:
                        line=line+str(None)
                else:
                    if obj[k] != None:
                        try:
                            line=line+f % (obj[k]) 
                        except:
                            print "Argh!"
                            ipshell()
                            sys.exit()
                    else:
                        line=line+str(None)
        # Add on a 'notes' column - any key which is just a bool gets added to a , delimited list if True
        notes=""
        for key in obj.keys():
            if type(obj[key]) == bool and obj[key] == True:
                if notes != "":
                    notes=notes+","
                notes=notes+key
        outFile.write(line+"\t"+notes+"\n")
    outFile.close()   
    
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
    #tab.table_name='ACT'
    fitsOutFileName=outFileName.replace(".csv", ".fits")
    if os.path.exists(fitsOutFileName) == True:
        os.remove(fitsOutFileName)
    tab.write(fitsOutFileName)
    

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
    
    inFile=file(fileName, "r")
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
                raise Exception, "couldn't find column headings in catalog file %s" % (fileName)
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
def readFelipeCatalog(fileName):
    """Reads a Felipe-style catalog.
    
    """
    
    inFile=file(fileName, "r")
    lines=inFile.readlines()
    inFile.close()
    
    catalog=[]
    for line in lines:
        if line[0] != '#' and len(line) > 3:
            objDict={}
            bits=line.rstrip("\n").split()
            objDict['name']=bits[0]
            objDict['RADeg']=astCoords.hms2decimal(bits[1], ":")
            objDict['decDeg']=astCoords.dms2decimal(bits[2], ":")
            catalog.append(objDict)
            
    return catalog

#------------------------------------------------------------------------------------------------------------
def readTobyLatexCatalog(fileName):
    """Read LaTeX format catalog like Toby sent for his paper.
    
    """

    inFile=file(fileName, 'r')
    lines=inFile.readlines()
    inFile.close()
    objList=[]
    for line in lines:
        objDict={}
        bits=line.split("&")
        objDict['name']=bits[0].lstrip(" ").rstrip(" ")
        objDict['RADeg']=astCoords.hms2decimal(bits[1].lstrip(" ").rstrip(" ").replace("$", ""), ":")
        objDict['decDeg']=astCoords.dms2decimal(bits[2].lstrip(" ").rstrip(" ").replace("$", ""), ":")
        objDict['SNR']=float(bits[3])
        objDict['thetaCoreArcmin']=float(bits[4])
        deltaTString=bits[5].replace("$", "").lstrip(" ").rstrip(" ")
        deltaTBits=deltaTString.split("\\pm")
        objDict['deltaT0String']="%d &plusmn; %d" % (int(deltaTBits[0]), int(deltaTBits[1]))
        objDict['deltaT0']=float(deltaTBits[0])
        objDict['deltaT0Err']=float(deltaTBits[1])
        deltaTCorrString=bits[6].replace("$", "").lstrip(" ").rstrip(" ")
        deltaTCorrBits=deltaTCorrString.split("\\pm")
        objDict['deltaT0CorrString']="%d &plusmn; %d" % (int(deltaTCorrBits[0]), int(deltaTCorrBits[1]))
        objDict['deltaT0Corr']=float(deltaTCorrBits[0])
        objDict['deltaT0CorrErr']=float(deltaTCorrBits[1])
        y0String=bits[8].replace("$", "").lstrip(" ").rstrip(" ")
        y0Bits=y0String.split("\\pm")
        objDict['y0String']="%.2f &plusmn; %.2f" % (float(y0Bits[0]), float(y0Bits[1]))
        objDict['y0']=float(y0Bits[0])
        objDict['y0Err']=float(y0Bits[1])
        objDict['z']=float(bits[10])
        objDict['altID']=bits[11].replace("\n", "").replace("\\", "").lstrip(" ").rstrip(" ")
        objList.append(objDict)
        
    return objList
            
#------------------------------------------------------------------------------------------------------------
def addNEDInfo(catalog, nedDir = "nedResults"):
    """Queries NED for matches near each object in the merged catalog, adds the nearest cluster match and 
    its distance in arcmin. We search a box 10' on a side near each object.
        
    """
    
    halfMatchBoxLengthDeg=5.0/60.0
    
    if os.path.exists(nedDir) == False:
        os.makedirs(nedDir)
            
    for obj in catalog:
        
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
        RAMin=RADeg-halfMatchBoxLengthDeg
        RAMax=RADeg+halfMatchBoxLengthDeg
        decMin=decDeg-halfMatchBoxLengthDeg
        decMax=decDeg+halfMatchBoxLengthDeg
                
        outFileName=nedDir+os.path.sep+name.replace(" ", "_")+".txt"        
        if os.path.exists(outFileName) == False:
            print "... fetching NED info for %s ..." % (name)
            try:                
                urllib.urlretrieve("http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.6fd&lat=%.6fd&radius=%.2f&dot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&in_objtypes1=QSO&in_objtypes2=Radio&in_objtypes2=SmmS&in_objtypes2=Infrared&in_objtypes2=Xray&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=ascii_tab&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (RADeg, decDeg, halfMatchBoxLengthDeg*60.0), filename = outFileName)
            except:
                print "WARNING: couldn't get NED info"
                #IPython.embed()
                #sys.exit()
                outFileName=None
                
        nedObjs=parseNEDResult(outFileName)
        
        # Flag matches against clusters - choose nearest one
        rMin=10000
        crossMatchRadiusDeg=2.5/60.0
        clusterMatch={}
        if len(nedObjs['RAs']) > 0:
            obj['NED_allClusterMatches']=[]
            for i in range(len(nedObjs['RAs'])):
                ned=nedObjs
                if ned['sourceTypes'][i] == 'GClstr':
                    r=astCoords.calcAngSepDeg(ned['RAs'][i], ned['decs'][i], RADeg, decDeg)
                    if r < crossMatchRadiusDeg:
                        obj['NED_allClusterMatches'].append(ned['names'][i])
                    if r < rMin and r < crossMatchRadiusDeg:
                        keepName=False
                        if 'name' in clusterMatch:
                            if "ABELL" in clusterMatch['name']:
                                keepName=True
                        if keepName == False:
                            rMin=r
                            clusterMatch['name']=ned['names'][i]
                            if ned['redshifts'][i] != 'N/A':
                                clusterMatch['z']=float(ned['redshifts'][i])
                            else:
                                clusterMatch['z']=None
                            clusterMatch['rArcmin']=rMin*60.0
                            clusterMatch['NED_RADeg']=float(ned['RAs'][i])
                            clusterMatch['NED_decDeg']=float(ned['decs'][i])
        
        if clusterMatch != {}:
            obj['NED_name']=clusterMatch['name']  
            obj['NED_z']=clusterMatch['z']
            obj['NED_distArcmin']=clusterMatch['rArcmin']
            obj['NED_RADeg']=clusterMatch['NED_RADeg']
            obj['NED_decDeg']=clusterMatch['NED_decDeg']
        else:
            obj['NED_name']=None 
            obj['NED_z']=None
            obj['NED_distArcmin']=None
            obj['NED_RADeg']=None
            obj['NED_decDeg']=None
             
#----------------------------------------------------------------------------------------------------
def parseNEDResult(inFileName, onlyObjTypes = None):
    """Parses NED tab-delimited text file query result, returns dictionary.
    
    onlyObjTypes can be a string indicating types of objects only to include e.g. GClstr
    
    """
    
    if inFileName != None and os.path.exists(inFileName):
        inFile=file(inFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    else:
        # Fail safe in case we couldn't contact NED
        lines=[]

    dataStarted=False
    labels=[]
    names=[]
    RAs=[]
    decs=[]
    sourceTypes=[]
    redshifts=[]
    for line in lines:
        bits=line.split("\t")
        if bits[0] == "1":
            dataStarted=True
        if dataStarted == True:
            if onlyObjTypes == str(bits[4]) or onlyObjTypes == None:
                labels.append(bits[0])
                names.append(bits[1])
                RAs.append(float(bits[2]))
                decs.append(float(bits[3]))
                sourceTypes.append(str(bits[4]))
                if bits[6] == '':
                    redshifts.append('N/A')
                else:
                    redshifts.append(str(bits[6]))

    return {'labels': labels, 'names': names, 'RAs': RAs, 'decs': decs, 'sourceTypes': sourceTypes, 'redshifts': redshifts}

#----------------------------------------------------------------------------------------------------
def addExtraMatches(catalog, extrasDict):
    """Matches catalog against an external catalog, which is simply a tab-delimited text file with
    columns name, RADeg, decDeg (nothing else just yet ...)
    
    """
    
    # Load data
    inFile=file(extrasDict['fileName'], "r")
    lines=inFile.readlines()
    inFile.close()
    
    names=[]
    RAs=[]
    decs=[]
    #zs=[]
    for line in lines:
        if line[0] != "#" and len(line) > 3:
            bits=line.split("\t")
            names.append(bits[0])
            RAs.append(float(bits[1]))
            decs.append(float(bits[2]))
            #zs.append(float(bits[3]))
    RAs=numpy.array(RAs)
    decs=numpy.array(decs)
    
    matchRadiusDeg=10.0/60.0
    for obj in catalog:
        obj['extra_%s_name' % (extrasDict['label'])]=None
        obj['extra_%s_RADeg' % (extrasDict['label'])]=None
        obj['extra_%s_decDeg' % (extrasDict['label'])]=None
        obj['extra_%s_distArcmin' % (extrasDict['label'])]=None
        rDeg=astCoords.calcAngSepDeg(obj['RADeg'], obj['decDeg'], RAs, decs)
        rMin=rDeg.min()
        if rMin < matchRadiusDeg:
            rMinIndex=rDeg.tolist().index(rMin)
            deltaRA=abs(obj['RADeg']-RAs[rMinIndex])
            deltaDec=abs(obj['decDeg']-decs[rMinIndex])
            if deltaRA < 10.0 and deltaDec < 10.0:
                obj['extra_%s_name' % (extrasDict['label'])]=names[rMinIndex]
                obj['extra_%s_RADeg' % (extrasDict['label'])]=RAs[rMinIndex]
                obj['extra_%s_decDeg' % (extrasDict['label'])]=decs[rMinIndex]
                obj['extra_%s_distArcmin' % (extrasDict['label'])]=rMin*60.0

    
    outKeys=['extra_%s_name'  % (extrasDict['label']), 'extra_%s_distArcmin'  % (extrasDict['label']), 
             'extra_%s_RADeg'  % (extrasDict['label']), 'extra_%s_decDeg'  % (extrasDict['label'])]
    outFormats=["%s", "%.3f", "%.6f", "%.6f"]
    outLabels=["%s name"  % (extrasDict['label']), "Distance from %s object (arcmin)"  % (extrasDict['label']),
               "%s R.A. (degrees)"  % (extrasDict['label']), "%s Dec. (degrees)"  % (extrasDict['label'])]
    
    return [outKeys, outFormats, outLabels]

#-------------------------------------------------------------------------------------------------------------
def addSDSSRedshifts(catalog, cacheDir = "SDSSQueryResults"):
    """Queries SDSS for redshifts. 
    
    """
    
    print ">>> Adding spec zs from SDSS ..."
    #url = 'http://cas.sdss.org/astrodr7/en/tools/search/x_sql.asp'
    #url = 'http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'
    url = 'http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx'

    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    count=0
    consecutiveQueryCount=0
    for obj in catalog:

        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        
        outFileName=cacheDir+os.path.sep+"%s.csv" % (obj['name'].replace(" ", "_"))
        if os.path.exists(outFileName) == False:
        
            sql="""SELECT
            p.objid,p.ra,p.dec,p.r,
            s.specobjid,s.z, 
            dbo.fSpecZWarningN(s.zWarning) as warning,
            s.plate, s.mjd, s.fiberid
            FROM PhotoObj AS p
            JOIN SpecObj AS s ON s.bestobjid = p.objid
            WHERE 
            p.ra < %.6f+0.1 and p.ra > %.6f-0.1
            AND p.dec < %.6f+0.1 and p.dec > %.6f-0.1
            """ % (obj['RADeg'], obj['RADeg'], obj['decDeg'], obj['decDeg'])

            # Filter SQL so that it'll work
            fsql = ''
            for line in sql.split('\n'):
                fsql += line.split('--')[0] + ' ' + os.linesep;
        
            params=urllib.urlencode({'cmd': fsql, 'format': "csv"})
            response=urllib2.urlopen(url+'?%s' % (params))
            lines=response.read()
            lines=lines.split("\n")

            outFile=file(outFileName, "w")
            for line in lines:
                outFile.write(line+"\n")
            outFile.close()
            
            consecutiveQueryCount=consecutiveQueryCount+1
            if consecutiveQueryCount > 50:
                print "... sleeping to give SDSS server a break ..."
                time.sleep(60)
                consecutiveQueryCount=0
        
        else:
            
            inFile=file(outFileName, "r")
            lines=inFile.readlines()
            inFile.close()
        
        # Parse .csv into catalog
        if lines[0] == "No objects have been found\n":
            obj['SDSSRedshifts']=None
        elif len(lines) > 1 and lines[1] == '"ERROR: Maximum 60 queries allowed per minute. Rejected query: SELECT \n':
            os.remove(outFileName)
            raise Exception, "Exceeded 60 queries/min on SDSS server. Take a breather and rerun nemo (previous queries cached)."
        else:
            obj['SDSSRedshifts']=[]
            for line in lines[2:]: # first line (DR7) always heading, first two lines (DR10) always heading
                if len(line) > 3:
                    zDict={}
                    bits=line.replace("\n", "").split(",")
                    zDict['objID']=bits[0]
                    try:
                        zDict['RADeg']=float(bits[1])
                        zDict['decDeg']=float(bits[2])
                    except:
                        if len(lines) > 1 and lines[1].find('"ERROR: Maximum 60 queries allowed per minute. Rejected query: SELECT') != -1:
                            raise Exception, "Exceeded 60 queries/min on SDSS server. Take a breather and rerun nemo (previous queries cached)."
                        else:
                            print "Hmm. Not able to parse SDSS redshifts"
                            IPython.embed()
                            sys.exit()
                    zDict['rMag']=float(bits[3])
                    zDict['specObjID']=bits[4]
                    try:
                        zDict['z']=float(bits[5])
                    except:
                        print "zDict['z'] problem"
                        IPython.embed()
                        sys.exit()
                    zDict['zWarning']=bits[6]
                    zDict['plate']=bits[7]
                    zDict['mjd']=bits[8]
                    zDict['fiberID']=bits[9]
                    obj['SDSSRedshifts'].append(zDict)

    
