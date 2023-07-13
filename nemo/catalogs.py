"""

This module contains tools for handling catalogs, which are usually :obj:`astropy.table.Table` objects.

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
import astropy.io.fits as pyfits
from scipy import ndimage
from . import maps

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
                   'template', 
                   'tileName',
                   'flags',
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
                   '%s',
                   '%s',
                   '%d',
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
                   '%.3f',
                   ]

columnsToAdd=[]
formatsToAdd=[]
for k in COLUMN_NAMES:
    if k in ['y_c', 'err_y_c', 'deltaT_c', 'err_deltaT_c']:
        i=COLUMN_NAMES.index(k)
        formatsToAdd.append(COLUMN_FORMATS[i])
        columnsToAdd.append("fixed_"+k)
COLUMN_NAMES=COLUMN_NAMES+columnsToAdd
COLUMN_FORMATS=COLUMN_FORMATS+formatsToAdd
        
if len(COLUMN_NAMES) != len(COLUMN_FORMATS):
    raise Exception("COLUMN_NAMES and COLUMN_FORMATS lists should be same length")

#------------------------------------------------------------------------------------------------------------
def _posRecFitFunc(snr, snrFold, pedestal, norm):
    """Fitting function used for position recovery offset (') in terms of fixed_SNR - see
    positionRecovery/positionRecoveryTestDriver.py.
    
    NOTE: Don't use this directly - call checkCrossMatch instead.
    
    """
    return norm*np.exp(-snr/snrFold)+pedestal
    
#------------------------------------------------------------------------------------------------------------
def checkCrossMatch(distArcmin, fixedSNR, z = None, addRMpc = 0.5, fitSNRFold = 1.164, fitPedestal = 0.685,
                    fitNorm = 38.097):
    """Checks the cross match offset between a cluster detection and an external catalog using a model derived
    from source injection sims (see :func:`nemo.maps.positionRecoveryAnalysis`). The position recovery test
    itself only accounts for the effect of noise fluctuations in the maps on the recovered SZ cluster
    positions.
    
    Args:
        distArcmin (:obj:`bool`): Distance of the potential cross match from the ACT position in arcmin.
        fixed_SNR (:obj:`float`): Signal-to-noise at reference filter scale (fixed_SNR) in ACT catalog.
        z (:obj:`float`, optional): If given, addRMpc will be converted to arcmin at this redshift, and then added
            in quadrature to the cross matching radius from the position recovery model.
        addRMpc (:obj:`float`, optional): Accounts for additional positional uncertainty (probably unknown) 
            in the external cross match catalog. This will be added in quadrature.
        fitSNRFold (:obj:`float`, optional): Model fit parameter - e-folding 
            (see :func:`nemo.maps.positionRecoveryAnalysis`).
        fitPedestal (:obj:`float`, optional): Model fit parameter - pedestal level
            (see :func:`nemo.maps.positionRecoveryAnalysis`).
        fitNorm (:obj:`float`, optional): Model fit parameter - normalization
            (see :func:`nemo.maps.positionRecoveryAnalysis`).
    
    Returns:
        True if distArcmin < model offset (+ optional addRMpc in arcmin at z), False if not.

    Note:
        The default values for the fit parameters are from a run on the f090, f150 ACT DR5 co-added maps
        (as used in the `ACT DR5 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2020arXiv200911043H/abstract>`_),
        and describe a function that recovers 99.7% of the inserted model clusters in 
        source injection simulations.
        
    """
    
    maxRadiusArcmin=_posRecFitFunc(fixedSNR, fitSNRFold, fitPedestal, fitNorm)
    addArcmin=0.0
    if z is not None and z > 0:
        addArcmin=np.degrees(addRMpc/astCalc.da(z))*60.0
    maxRadiusArcmin=np.sqrt(maxRadiusArcmin**2 + addArcmin**2)
    if distArcmin < maxRadiusArcmin:
        return True
    else:
        return False

#------------------------------------------------------------------------------------------------------------
def makeOptimalCatalog(catalogDict, constraintsList = []):
    """Identifies common objects between every catalog in the input dictionary of catalogs, and creates a 
    master catalog with one entry per object, keeping only the details of the highest signal-to-noise 
    detection.
    
    Args:
        catalogDict (:obj:`dict`): Dictionary where each key points to a catalog of objects.
        constraintsList (:obj:`list`, optional): A list of constraints (for the format, see 
            :func:`selectFromCatalog`).
        
    Returns:
        None - an ``optimalCatalog`` key is added to ``catalogDict`` in place.
    
    """
    
    allCatalogs=[]
    for key in catalogDict.keys():
        if len(catalogDict[key]['catalog']) > 0:
            allCatalogs.append(catalogDict[key]['catalog'])
    if len(allCatalogs) > 0:
        allCatalogs=atpy.vstack(allCatalogs)
        mergedCatalog=allCatalogs.copy()
        mergedCatalog.add_column(atpy.Column(np.zeros(len(mergedCatalog)), 'toRemove'))
        for row in allCatalogs:
            rDeg=astCoords.calcAngSepDeg(row['RADeg'], row['decDeg'], allCatalogs['RADeg'].data, 
                                         allCatalogs['decDeg'].data) 
            xIndices=np.where(rDeg < XMATCH_RADIUS_DEG)[0]
            if len(xIndices) > 1:
                xMatches=allCatalogs[xIndices]
                xMatchIndex=np.argmax(xMatches['SNR'])
                for index in xIndices:
                    if index != xIndices[xMatchIndex]:
                        mergedCatalog['toRemove'][index]=1
        mergedCatalog=mergedCatalog[mergedCatalog['toRemove'] == 0]
        mergedCatalog.remove_column('toRemove')
        mergedCatalog.sort(['RADeg', 'decDeg'])
        mergedCatalog=selectFromCatalog(mergedCatalog, constraintsList)
    else:
        mergedCatalog=atpy.Table(names=('SNR', 'RADeg', 'decDeg'))
    
    return mergedCatalog

#------------------------------------------------------------------------------------------------------------
def catalog2DS9(catalog, outFileName, constraintsList = [], addInfo = [], idKeyToUse = 'name', 
                RAKeyToUse = 'RADeg', decKeyToUse = 'decDeg', color = "cyan", showNames = True,
                writeNemoInfo = True, coordSys = 'fk5', regionShape = 'point', width = 1):
    """Writes a DS9 region file corresponding to the given catalog. 
    
    Args:
        catalog (:obj:`astropy.table.Table`): An astropy Table where each row represents an object.
        outFileName (:obj:`str`): A file name for the output DS9 region file.
        constraintsList (:obj:`list`, optional): A list of constraints in the same format as used by 
            :func:`selectFromCatalog`.
        addInfo (:obj:`list`, optional): A list of dictionaries with keys named `key` and `fmt` (e.g., 
            ``{'key': "SNR", 'fmt': "%.3f"}``). These will be added to the object label shown in DS9.
        idKeyToUse (:obj:`str`, optional): The name of the key in each object dictionary that defines the 
            object's name. Used to label objects in the DS9 region file.
        RAKeyToUse (:obj:`str`, optional): The name of the key in each object dictionary that contains the 
            RA of the object in decimal degrees.
        decKeyToUse (:obj:`str`, optional): The name of the key in each object dictionary that contains the
            declination of the object in decimal degrees.
        color (:obj:`str`, optional): The color of the plot symbol used by DS9.
        writeNemoInfo (:obj:`bool`, optional): If ``True``, writes a line with the `nemo` version and date 
            generated at the top of the DS9 .reg file.
        coordSys (:obj:`str`, optional): A string defining the coordinate system used for RA, dec, as 
            understood by DS9.
        
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
        outFile.write('global dashlist=8 3 width=%d font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' % (width))
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
            if showNames == True:
                infoString=str(obj[idKeyToUse])+infoString
            if regionShape == 'point':
                outFile.write("%s;point(%.6f,%.6f) # point=cross color={%s} text={%s}\n" \
                            % (coordSys, obj[RAKeyToUse], obj[decKeyToUse], colorString, infoString))
            elif regionShape == 'circle':
                outFile.write('%s;circle(%.6f,%.6f,360") # color={%s} text={%s}\n' \
                            % (coordSys, obj[RAKeyToUse], obj[decKeyToUse], colorString, infoString))                

#------------------------------------------------------------------------------------------------------------
def makeName(RADeg, decDeg, prefix = 'ACT-CL'):
    """Makes an object name string from the given object coordinates, following the IAU convention.
    
    Args:
        RADeg (:obj:`float`): Right ascension of the object in J2000 decimal degrees.
        decDeg (:obj:`float`): Declination of the object in J2000 decimal degrees.
        prefix (:obj:`str`, optional): Prefix for the object name.
    
    Returns:
        Object name string in the format `prefix JHHMM.m+/-DDMM`.
    
    """
    
    actName=prefix+" J"+_makeRA(RADeg)+_makeDec(decDeg)
    
    return actName

#------------------------------------------------------------------------------------------------------------
def makeLongName(RADeg, decDeg, prefix = "ACT-CL"):
    """Makes a long format object name string from the given object coordinates, following the IAU convention.
    
    Args:
        RADeg (:obj:`float`): Right ascension of the object in J2000 decimal degrees.
        decDeg (:obj:`float`): Declination of the object in J2000 decimal degrees.
        prefix (:obj:`str`, optional): Prefix for the object name.
    
    Returns:
        Object name string in the format `prefix JHHMMSS.s+/-DDMMSS`.
    
    """
    
    actName=prefix+" J"+_makeLongRA(RADeg)+_makeLongDec(decDeg)
    
    return actName
    
#------------------------------------------------------------------------------------------------------------
def _makeRA(myRADeg):
    """Makes RA part of ACT names.
    
    """
    hours=(myRADeg/360)*24
    strHours=("%.10f" % (hours))
    if hours<10:
        sHours="0"+strHours[0]
    else:
        sHours=strHours[:2]
    
    mins=float(strHours[strHours.index("."):])*60
    strMins=("%.10f" % (mins))
    if mins < 10:
        sMins="0"+strMins[:3]
    else:
        sMins=strMins[:4]

    return (sHours+sMins)#[:-2] # Trims off .x as not used in ACT names
        
#------------------------------------------------------------------------------------------------------------
def _makeDec(myDecDeg):
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
def _makeLongRA(myRADeg):
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
def _makeLongDec(myDecDeg):
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
    """Return a table of objects matching the given constraints from the catalog. 
    
    Args:
        catalog (:obj:`astropy.table.Table`): The catalog from which objects will be selected.
        constraintsList (:obj:`list`): A list of constraints, where each item is a string of the form
            "key < value", "key > value", etc.. Note that the spaces between the key, operator 
            (e.g. '<'), and value are essential.
    
    Returns:
        An astropy Table object.
    
    """
            
    passedConstraint=catalog
    for constraintString in constraintsList:
        key, op, value=constraintString.split()
        passedConstraint=passedConstraint[eval("passedConstraint['%s'] %s %s" % (key, op, value))]

    return passedConstraint

#------------------------------------------------------------------------------------------------------------
def catalogListToTab(catalogList, keysToWrite = COLUMN_NAMES):
    """Converts a catalog in the form of a list of dictionaries (where each dictionary holds the object
    properties) into an :obj:`astropy.table.Table` object.
    
    Args:
        catalogList (:obj:`list`): Catalog in the form of a list of dictionaries.
        keysToWrite (:obj:`list`, optional): Keys to convert into columns in the output table.
    
    Returns:
        An :obj:`astropy.table.Table` object, where each row corresponds to an object in the catalog.
    
    """
    
    availKeys=list(catalogList[0].keys())
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
    """Converts an :obj:`astropy.table.Table` object into a list of dictionaries.
    
    Args:
        tab (:obj:`astropy.table.Table`): Catalog in the form of an astropy Table object.
    
    Returns:
        A list of dictionaries, where each dictionary represents an object in the catalog.
    
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
    """Writes the catalog to disk, including meta data on the Nemo version used to create the catalog in the
    header.
     
    Args:
        catalog (:obj:`astropy.table.Table`): The object catalog.
        outFileName (:obj:`str`): The output filename. The format depends upon the extension given. Note
            that CSV format (resulting from using the ``.csv`` extension) is written tab-delimited. 
        constraintsList (:obj:`list`, optional): A list of constraints, where each item is a string of the
            form "key < value", "key > value", etc.. Note that the spaces between the key, operator 
            (e.g. '<'), and value are essential.
        
    Returns:
        None
        
    """

    # Deal with any blank catalogs (e.g., in blank tiles - not that we would want such things...)
    if type(catalog) == list and len(catalog) == 0:
        return None
    cutCatalog=selectFromCatalog(catalog, constraintsList)
    cutCatalog.meta['NEMOVER']=nemo.__version__
    if outFileName.split(".")[-1] == 'csv':
        cutCatalog.write(outFileName, format = 'ascii.csv', delimiter = '\t', overwrite = True)
    else:
        cutCatalog.write(outFileName, overwrite = True)

#------------------------------------------------------------------------------------------------------------
def removeDuplicates(tab):
    """Removes duplicate objects from the catalog - keeping the highest SNR detection for each duplicate. 
    This routine is used to clean up the output of MPI runs (where we have overlapping tiles).
    
    Args:
        tab (:obj:`astropy.table.Table`): The object catalog to be checked for duplicates.

    Returns:
        Table with duplicates removed (:obj:`astropy.table.Table`), the number of duplicates found, and a
        list of names for the duplicated objects.
    
    """

    if len(tab) == 1:
        return tab, 1, []
    
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
        mask=np.less_equal(rDeg, XMATCH_RADIUS_DEG)
        if mask.sum() == 0:	# This ought not to be possible but catch anyway
            bestIndex=i
        else:
            indices=np.where(mask == True)[0]
            bestIndex=indices[np.equal(dupTab['SNR'][mask], dupTab['SNR'][mask].max())][0]
        keepMask[bestIndex]=True
    keepTab=dupTab[keepMask]
    
    keepTab=atpy.vstack([keepTab, noDupTab])
    keepTab.sort('RADeg')
    
    return keepTab, len(dupTab), dupTab['name']

#------------------------------------------------------------------------------------------------------------
def flagTileBoundarySplits(tab, xMatchRadiusArcmin = 2.5):
    """Flag objects that are closer than some matching radius but which appear in different tiles. These are
    potentially objects that have been de-blended across tile boundaries (in which case one entry in the
    catalog should be kept and the other discarded), but some objects in this category are genuine close
    pairs. At the moment, this is best resolved by visual inspection. This routine adds a flag column named
    `tileBoundarySplit` to the catalog to make it easy to spot these cases.
        
    Args:
        tab (:obj:`astropy.table.Table`): The object catalog to be checked.
        xMatchRadiusArcmin (float, optional): Cross-match radius in arcmin.

    Returns:
        Catalog (:obj:`astropy.table.Table`) with `tileBoundarySplit` column added - this is True for objects
        that may have been deblended across tiles, and require visual inspection.
    
    """
    
    if len(tab) == 1:
        return tab
    
    xMatchRadiusDeg=xMatchRadiusArcmin/60.
    
    # Find all potential duplicates within a given matching radius
    cat=SkyCoord(ra = tab['RADeg'].data, dec = tab['decDeg'].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat, cat, nthneighbor = 2)
    mask=np.less(rDeg.value, xMatchRadiusDeg)
    noDupMask=np.greater_equal(rDeg.value, xMatchRadiusDeg)
    dupTab=tab[mask]
    noDupTab=tab[noDupMask]
        
    # Identify pairs split across tile boundaries
    tileBoundarySplit=np.zeros(len(dupTab), dtype = bool)
    for i in range(len(dupTab)):
        # NOTE: astCoords does not like atpy.Columns sometimes...
        rDeg=astCoords.calcAngSepDeg(dupTab['RADeg'][i], dupTab['decDeg'][i], dupTab['RADeg'].data, dupTab['decDeg'].data)
        mask=np.less_equal(rDeg, xMatchRadiusDeg)
        if mask.sum() == 0:	# This ought not to be possible but catch anyway
            bestIndex=i
        else:
            indices=np.where(mask == True)[0]
            bestIndex=indices[np.equal(dupTab['SNR'][mask], dupTab['SNR'][mask].max())][0]
            if np.unique(dupTab['tileName'][indices]).shape[0] > 1:
                tileBoundarySplit[indices]=True
    dupTab['tileBoundarySplit']=tileBoundarySplit
    dupTab=dupTab[tileBoundarySplit]
    
    # Flag in the main table
    tab['tileBoundarySplit']=np.zeros(len(tab), dtype = bool)
    for row in dupTab:
        mask=(tab['name'] == row['name'])
        tab['tileBoundarySplit'][mask]=True
        
    return tab

#------------------------------------------------------------------------------------------------------------
def generateRandomSourcesCatalog(mapData, wcs, numSources, seed = None):
    """Generate a random source catalog (with amplitudes in deltaT uK), with random positions within the 
    footprint of the given map (areas where pixel values == 0 are ignored). The distribution of source 
    amplitudes is roughly similar to that seen in the 148 GHz ACT maps, but this routine should only be used
    for tests - it is not a serious attempt at simulating the real extragalactic source population.
    
    Args:
        mapData (:obj:`numpy.ndarray`): Map pixel-data, only used for determining valid area in which sources
            may be randomly placed (pixel values == 0 are ignored).
        wcs (:obj:`astWCS.WCS`): WCS corresponding to the map.
        numSources (:obj:`int`): Number of random sources to put into the output catalog.
        seed (optional, :obj:`int`): If given, generate the catalog using this random seed value. This is
            useful for generating the same realization across maps at different frequencies. The seed will
            be reset after this routine exits.
    
    Returns:
        An :obj:`astropy.table.Table` object containing the catalog.
        
    """
    
    if seed is not None:
        np.random.seed(seed)
        
    deltaT=np.random.lognormal(np.log(600), 1.1, numSources)
    ys, xs=np.where(mapData != 0)
    ys=ys+np.random.uniform(0, 1, len(ys))
    xs=xs+np.random.uniform(0, 1, len(xs))
    indices=np.random.randint(0, len(ys), numSources)
    coords=wcs.pix2wcs(xs[indices], ys[indices])
    coords=np.array(coords)
    tab=atpy.Table()
    tab.add_column(atpy.Column(np.arange(0, numSources)+1, "name"))
    tab.add_column(atpy.Column(coords[:, 0], "RADeg"))
    tab.add_column(atpy.Column(coords[:, 1], "decDeg"))
    tab.add_column(atpy.Column(deltaT, "deltaT_c"))
    
    if seed is not None:
        np.random.seed()
        
    return tab

#------------------------------------------------------------------------------------------------------------
def generateTestCatalog(config, numSourcesPerTile, amplitudeColumnName = 'fixed_y_c', 
                        amplitudeRange = [0.001, 1], amplitudeDistribution = 'linear', selFn = None,
                        avoidanceRadiusArcmin = 20.0, maskDilationPix = 0,
                        tileNames = None):
    """Generate a catalog of objects with random positions and amplitudes. This is for testing purposes - 
    see, e.g., :meth:`nemo.maps.sourceInjectionTest`.
    
    Args:
        config (:obj:`nemo.startup.NemoConfig`): Nemo configuration object.
        numSourcesPerTile (:obj:`int`): The maximum number of sources to insert into each tile. The number 
            of sources actually inserted may be less than this depending on the value of 
            `avoidanceRadiusArcmin`.
        amplitudeColumnName (:obj:`str`): Name of the column in the output catalog in which source (or
            cluster) amplitudes will be stored. Typically this should be "deltaT_c" for sources, and
            "fixed_y_c" for clusters.
        amplitudeRange (:obj:`list`): Range for the random amplitudes, in the form [minimum, maximum].
        amplitudeDistribution (:obj:`str`): Either 'linear' or 'log'.
        selFn (:obj:`nemo.completeness.SelFn`, optional): Nemo selection function object, used to access 
            area masks and coordinates info. If not given, a selFn object will be created using the info 
            in the config. Providing this saves time, as the area mask files don't have to be read from 
            disk.
        avoidanceRadiusArcmin (:obj:`float`): Minimum separation between two objects in the output catalog.
            This should be set large enough to avoid crowding and spurious cross-matching problems.
        maskDilationPix (:obj:`int`): Avoid placing objects within this many pixels of the area mask, by
            dilating the holes in the mask by this number of pixels. Avoids problems with e.g., objects
            split across the mask boundary. If not set, gives a small number of objects with larger than
            expected offsets in recovered position in :meth:`nemo.maps.sourceInjectionTest`).

    Returns:
        An :obj:`astropy.table.Table` object containing the catalog.
        
    """
    
    if selFn is None:
        selFn=completeness.SelFn(config.selFnDir, 4.0, configFileName = config.configFileName, 
                                 enableCompletenessCalc = False, setUpAreaMask = True,
                                 tileNames = tileNames)
    if tileNames is None:
        tileNames=selFn.tileNames

    RAs=[]
    decs=[]
    amps=[]
    catTileNames=[]
    for tileName in tileNames:
        mapData=selFn.areaMaskDict[tileName]
        if mapData.sum() == 0:  # Skip any empty/blank tile
            continue
        wcs=selFn.WCSDict[tileName]
        # Dilate mask to avoid putting objects close to borders or holes
        for i in range(maskDilationPix):
            mapData=1-ndimage.binary_dilation(1-mapData)
        ys, xs=np.where(mapData != 0)
        ys=ys+np.random.uniform(0, 1, len(ys))
        xs=xs+np.random.uniform(0, 1, len(xs))
        indices=np.random.randint(0, len(ys), len(ys))
        coords=wcs.pix2wcs(xs[indices], ys[indices])
        coords=np.array(coords)
        tileRAs=[]
        tileDecs=[]
        keepIndices=[]
        for i in indices:
            rai=coords[i, 0]
            deci=coords[i, 1]
            rDeg=astCoords.calcAngSepDeg(rai, deci, tileRAs, tileDecs)
            keepObj=False
            if len(rDeg) > 0:
                if rDeg.min()*60 > avoidanceRadiusArcmin:
                    keepObj=True
            else:
                keepObj=True
            if keepObj == True:
                tileRAs.append(rai)
                tileDecs.append(deci)
                keepIndices.append(i)
            if len(keepIndices) == numSourcesPerTile:
                break
        # Edge case where we have tiny but non-zero area and didn't manage to insert anything
        if len(coords) == 0:
            continue
        RAs=RAs+coords[keepIndices, 0].tolist()
        decs=decs+coords[keepIndices, 1].tolist()
        if amplitudeDistribution == 'linear':
            amp=np.random.uniform(amplitudeRange[0], amplitudeRange[1], len(keepIndices))
        elif amplitudeDistribution == 'log':
            amp=np.power(10, np.random.uniform(np.log10(amplitudeRange)[0], np.log10(amplitudeRange)[1], len(keepIndices)))
        else:
            raise Exception("Must be either 'linear' or 'log'.")
        amps=amps+amp.tolist()
        catTileNames=catTileNames+[tileName]*len(amp.tolist())
    
    tab=atpy.Table()
    tab.add_column(atpy.Column(np.arange(0, len(amps))+1, "name"))
    tab.add_column(atpy.Column(RAs, "RADeg"))
    tab.add_column(atpy.Column(decs, "decDeg"))
    tab.add_column(atpy.Column(amps, amplitudeColumnName))
    tab.add_column(atpy.Column(catTileNames, 'tileName'))

    return tab

#------------------------------------------------------------------------------------------------------------
def crossMatch(refCatalog, matchCatalog, radiusArcmin = 2.5):
    """Cross matches `matchCatalog` onto `refCatalog` for objects found within some angular radius 
    (specified in arcmin).
    
    Args:
        refCatalog (:obj:`astropy.table.Table`): The reference catalog.
        matchCatalog (:obj:`astropy.table.Table`): The catalog to match onto the reference catalog.
        radiusArcmin (:obj:`float`, optional): Cross-match radius in arcmin.
    
    Returns:
        Cross-matched reference catalog, matchCatalog, and array of angular separation in degrees, for 
        objects in common within the matching radius. The cross matched columns are sorted such that rows in
        each correspond to the matched objects.
    
    """
    
    inTab=refCatalog
    outTab=matchCatalog
    RAKey1, decKey1=getTableRADecKeys(inTab)
    RAKey2, decKey2=getTableRADecKeys(outTab)
    cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
    xMatchRadiusDeg=radiusArcmin/60.
    cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
    mask=np.less(rDeg.value, xMatchRadiusDeg)  
    matched_outTab=outTab[xIndices]
    inTab=inTab[mask]
    matched_outTab=matched_outTab[mask]
    rDeg=rDeg.value[mask]
    
    return inTab, matched_outTab, rDeg

#------------------------------------------------------------------------------------------------------------
def removeCrossMatched(refCatalog, matchCatalog, radiusArcmin = 2.5):
    """Cross matches `matchCatalog` onto `refCatalog` for objects found within some angular radius 
    (specified in arcmin), and returns `refCatalog` with the matching entries removed.
    
    Args:
        refCatalog (:obj:`astropy.table.Table`): The reference catalog.
        matchCatalog (:obj:`astropy.table.Table`): The catalog to match onto the reference catalog.
        radiusArcmin (:obj:`float`, optional): Cross-match radius in arcmin.
    
    Returns:
        Cross-matched reference catalog (:obj:`astropy.table.Table`) with matches to `matchCatalog` removed.
        
    """
        
    inTab=refCatalog
    outTab=matchCatalog
    RAKey1, decKey1=getTableRADecKeys(inTab)
    RAKey2, decKey2=getTableRADecKeys(outTab)
    cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
    xMatchRadiusDeg=radiusArcmin/60.
    cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
    mask=np.greater(rDeg.value, xMatchRadiusDeg)  
    inTab=inTab[mask]
    
    return inTab
    
#------------------------------------------------------------------------------------------------------------
def getTableRADecKeys(tab):
    """Returns the column names in the table in which RA, dec coords are stored, after trying a few possible
    name variations.
    
    Args:
        tab (:obj:`astropy.table.Table`): The table to search.
        
    Returns:
        Name of the RA column, name of the dec. column
    
    """
    RAKeysToTry=['ra', 'RA', 'RADeg']
    decKeysToTry=['dec', 'DEC', 'decDeg', 'Dec']
    RAKey, decKey=None, None
    for key in RAKeysToTry:
        if key in tab.keys():
            RAKey=key
            break
    for key in decKeysToTry:
        if key in tab.keys():
            decKey=key
            break
    if RAKey is None or decKey is None:
        raise Exception("Couldn't identify RA, dec columns in the supplied table.")
    
    return RAKey, decKey

#------------------------------------------------------------------------------------------------------------
def getCatalogWithinImage(tab, shape, wcs, mask = None):
    """Returns the subset of the catalog with coordinates within the image defined by the given `shape`, 
    `wcs`. Optionally, a `mask` may also be applied.
    
    Args:
        tab (:obj:`astropy.table.Table`): Catalog, as an astropy Table object. Must have columns called 
            'RADeg', 'decDeg' that contain object coordinates in decimal degrees.
        shape (:obj:`list`): Shape of the array corresponding to the image / map.
        wcs (:obj:`astWCS.WCS`): WCS of the image.
        mask (optional, :obj:`np.ndarray`): Mask with same dimensions and WCS as the image. Pixels with
            value = 1 indicate valid area, and pixels with value = 0 are considered to be outside the mask.
            If this is given, the returned catalog will contain only objects in the valid area defined by
            this image mask.
    
    Returns:
        An astropy Table containing the subset of objects within the image.
    
    """
    
    xyCoords=np.array(wcs.wcs2pix(tab['RADeg'].tolist(), tab['decDeg'].tolist()))
    selected=[]
    for i in range(len(tab)):
        x, y=xyCoords[i][0], xyCoords[i][1]
        if np.isnan(x) == True or np.isnan(y) == True:
            selected.append(False)
            continue
        if x >= 0 and x < shape[1]-1 and y >= 0 and y < shape[0]-1:
            if mask is not None:
                if mask[int(round(y)), int(round(x))] == 1:
                    selected.append(True)
                else:
                    selected.append(False)
            else:
                selected.append(True)
        else:
            selected.append(False)
    
    return tab[selected]

#------------------------------------------------------------------------------------------------------------
def addFootprintColumnToCatalog(tab, label, areaMask, wcs):
    """Add `footprint_label` column to the catalog, flagging objects found within the valid area of the given
    mask.

    Args:
        tab (:obj:`astropy.table.Table`): Catalog, as an astropy Table object. Must have columns called
            'RADeg', 'decDeg' that contain object coordinates in decimal degrees.
        label (:obj:`str`): A column named `footprint_label` will be added to the catalog. Objects in the
            catalog that fall within the valid area of the given area mask will have `footprint_label`
            set to True.
        areaMask (:obj:`np.ndarray`): Mask image defining the footprint corresponding to the given WCS.
            Pixels with value = 1 indicate valid area, and pixels with value = 0 are considered to be
            outside the mask.
        wcs (:obj:`astWCS.WCS`): WCS of the area mask that defines the footprint.

    Returns:
        An astropy Table with `footprint_label` column added.

    """

    inMask=np.zeros(len(tab['RADeg'].data), dtype = bool)
    coords=wcs.wcs2pix(tab['RADeg'].data, tab['decDeg'].data)
    coords=np.array(np.round(coords), dtype = int)
    mask1=np.logical_and(coords[:, 0] >= 0, coords[:, 1] >= 0)
    mask2=np.logical_and(coords[:, 0] < areaMask.shape[1], coords[:, 1] < areaMask.shape[0])
    mask=np.logical_and(mask1, mask2)
    inMask[mask]=inMask[mask]+areaMask[coords[:, 1][mask], coords[:, 0][mask]]
    tab['footprint_%s' % (label)]=inMask

    return tab
