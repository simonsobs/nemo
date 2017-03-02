# -*- coding: utf-8 -*-
"""This module contains tools for making source browser type webpages and associated plots.

"""

import os
import xmlrpclib
import operator
import urllib
import urllib2
import base64
import glob
from astLib import *
import pyfits
import numpy
import pylab
from scipy import ndimage
import catalogTools
import mapTools
import S82Grabber
import nemo
import sys
import time
import string
from PIL import Image
import IPython
import astropy.table as atpy

#-------------------------------------------------------------------------------------------------------------
# Constants
PLOT_HEIGHT_PIX=1400    # adjust this to change height of all .jpg plots generated at once
PLOT_DISPLAY_HEIGHT_PIX=600   # adjust this to change size the plots are displayed at in the .html files

#-------------------------------------------------------------------------------------------------------------
def makeTablePage(catalog, outDir, linksDir, keysToWrite, keyFormats, keyLabels, objectTypeString = "candidates",
                  surveyAreaDeg2 = None, catalogToLink = None, parFileToLink = None, commentsString = None, 
                  addLinkToDiagnosticsPage = False, addLinkToFITSTarFile = False, imageMapOptions = None):
    """Makes table webpage that serves as index. linksDir specifies the imageDir to use for default object
    page links.
    
    """

    templatePage="""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
    <html>
    <head>
        <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
        <title>Source Browser</title>
    </head>
    <body style="font-family: sans-serif; vertical align: top; justify: full;">
    <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
        <tbody>
            <tr>
                <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 125%;">
                    Source Browser
                </td>
            </tr>
        </tbody>
    </table>
    
    $META_DATA
    $COLOR_CODING
    
    <table frame=border cellspacing=0 cols=6 rules=all border=2 width=100% align=center>
    <tbody>
        <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                text-align: center; vertical-align: middle; font-size: 110%;">"""
    for key in keyLabels:
        templatePage=templatePage+"\n           <td><b>%s</b></td>" % (key)
    templatePage=templatePage+"""
        </tr>
        $TABLE_DATA
    </tbody>
    </table>
        
    </tbody>
    </table>
    <br><br>
    </body>
    </html>
    """
    
    # Useful stuff
    READMEComment="""Matches to other catalogs (e.g. NED) listed on this page are within 2.5' radius of the 
    candidate position."""
    #Integrated Y values not to be trusted at present. &Delta;T<sub>c</sub> and 
    #y<sub>c</sub> values are as measured directly from filtered maps - no correction factors applied.
    if commentsString == None:
        commentsString=READMEComment
    else:
        commentsString=commentsString+" "+READMEComment
    
    html=templatePage
    metaData="""<br><b>Total number of %s</b> = %d
    <br><br>""" % (objectTypeString, len(catalog))
    if surveyAreaDeg2 != None:
        metaData=metaData+"""<b>Survey area</b> = %.1f deg<sup>2</sup>
        <br><br>""" % (surveyAreaDeg2)
    if parFileToLink != None:
        shortParFileName=os.path.split(parFileToLink)[-1]
        metaData=metaData+"""<b>Parameters file:</b> <a href=%s>%s</a><br><br>""" % (shortParFileName, shortParFileName)
        os.system("cp %s %s/" % (parFileToLink, outDir)) 
        os.system("chmod a+r %s/%s" % (outDir, shortParFileName)) # For some reason .par files are defaulting to unreadable
    if catalogToLink != None:
        shortCatalogName=os.path.split(catalogToLink)[-1]
        shortFITSName=shortCatalogName.replace(".csv", ".fits")
        shortRegName=shortCatalogName.replace(".csv", ".reg")
        metaData=metaData+"""<b>Catalog:</b>
        <ul><li><a href=%s>%s</a>   (plain text, tab-delimited)</li>
        <li><a href=%s>%s</a>   (FITS table format)</li>
        <li><a href=%s>%s</a>   (DS9 region file)</li></ul>
        """ % (shortCatalogName, shortCatalogName, shortFITSName, shortFITSName, shortRegName, shortRegName)
        os.system("cp %s %s/" % (catalogToLink, outDir)) 
        os.system("cp %s %s/" % (catalogToLink.replace(".csv", ".fits"), outDir))
        os.system("cp %s %s/" % (catalogToLink.replace(".csv", ".reg"), outDir))
    if imageMapOptions != None:
        metaData=metaData+"<a href=imageMap.html>Clickable full size filtered map image</a><br><br>\n"
        makeImageMapPage(catalog, outDir, imageMapOptions)
    if addLinkToDiagnosticsPage == True:
        metaData=metaData+"<a href=diagnostics/plots.html>Catalog properties</a><br><br>\n"
    if addLinkToFITSTarFile == True:
        metaData=metaData+"<a href=FITSImages.tar.gz>Download .tar.gz containing candidate 148 GHz FITS images</a><br>\n"
    if commentsString != None:
        commentsString=commentsString.replace("$EQUALS", "=") # yes, we really do have to do this to read as actDict
        metaData=metaData+"<b>Comments:</b>\n"
        metaData=metaData+'<p style="margin-left: 1em; margin-right: 1em">'+commentsString+'</p>\n'
    metaData=metaData+"<br>\n"
    html=html.replace("$META_DATA", metaData)
            
    tableData=""
    usedBckColors=[]
    usedBckKeys=[]
    for obj in catalog:
        
        # NED and extra matches - which is nearest?
        # This may well need altering to work with Toby's stuff
        bestAltMatchDistArcmin=1000
        minMatchDistArcmin=2.5  # If greater than this, do not show in NED ID column
        if 'NED_name' in obj.keys() and obj['NED_name'] != None:
            bestAltMatchDistArcmin=obj['NED_distArcmin']
            urlName=obj['NED_name'].replace("+", "%2B").replace(" ", "+")
            urlString="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname="+urlName+"&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"
            clusterMatch="<a href=%s>%s</a> (%.1f)" % (urlString, obj['NED_name'], obj['NED_distArcmin'])
        extraDistArcmin=None
        for key in obj.keys():
            if 'extra_' in key and '_name' in key:
                extraName=obj[key]
            if 'extra_' in key and '_distArcmin' in key:
                extraDistArcmin=obj[key]
        if extraDistArcmin != None and extraDistArcmin < bestAltMatchDistArcmin:
            bestAltMatchDistArcmin=extraDistArcmin
            clusterMatch="%s (%.1f)" % (extraName, extraDistArcmin)
        if bestAltMatchDistArcmin > minMatchDistArcmin:
            clusterMatch="-"  
            
        # Highlighting of rows - obviously, order matters here!
        bckColor="white"
        if 'observed 2009B' in obj.keys() and obj['observed 2009B'] == True:
            bckColor="darkgray"
            bckKey='observed 2009B'
        if 'SPT cluster' in obj.keys() and obj['SPT cluster'] == True:
            bckColor="gold"
            bckKey='SPT cluster'
        if 'ACT 2008 cluster' in obj.keys() and obj['ACT 2008 cluster'] == True:
            bckColor="deeppink"
            bckKey='ACT 2008 cluster'
        if bckColor not in usedBckColors and bckColor != "white":
            usedBckColors.append(bckColor)
            usedBckKeys.append(bckKey)
            
        # Row for each cluster in table
        rowString="<tr>\n"
        for key in keysToWrite:
            htmlKey="$"+string.upper(key)+"_KEY"
            rowString=rowString+"   <td style='background-color: "+bckColor+";' align=center width=10%>"+htmlKey+"</td>\n"
        rowString=rowString+"</tr>\n"
        
        # Insert values - note name is special
        # Need to work out how to handle the clusterMatch column with Toby's stuff
        for key, fmt in zip(keysToWrite, keyFormats):
            htmlKey="$"+string.upper(key)+"_KEY"
            if key == "name":
                nameLink="<a href=\"%s\" target='_blank'>%s</a>" % \
                    (os.path.split(linksDir)[-1]+os.path.sep+obj['name'].replace(" ", "_")+".html", obj['name'])
                rowString=rowString.replace(htmlKey, "%s" % (nameLink))
            else:
                if key == "NED_name" and obj[key] != None:
                    nedName=obj[key]
                    nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedName.replace("+", "%2B").replace(" ", "+"))
                    rowString=rowString.replace(htmlKey, "<a href=%s>%s</a>" % (nedLinkURL, nedName))
                else:
                    rowString=rowString.replace(htmlKey, fmt % (obj[key]))
        tableData=tableData+rowString
        
    html=html.replace("$TABLE_DATA", tableData)
    
    # Colou coding table key
    if len(usedBckColors) > 0:
        colorCoding="""<table frame=border cellspacing=0 cols=2 rules=all border=2 width=60% align=center>
                     <tbody>
                        <tr>
                            <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                                text-align: center; vertical-align: middle; font-size: 110%;" colspan=$NUM_COLORS>Color coding</td>
                        </tr>
                        <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; 
                                text-align: center; vertical-align: middle; font-size: 110%;">
                        $COLOR_KEY
                        </tr>
                    </tbody>
                    </table>
                    <br><br>
                    """
        colorCoding=colorCoding.replace("$NUM_COLORS", str(len(usedBckColors)))
        keyString=""
        for bckColor, bckKey in zip(usedBckColors, usedBckKeys):
            keyString=keyString+'<td style="background-color:'+bckColor+';" width='+str(100.0/len(usedBckColors))+'%>'+bckKey+'</td>\n'
        colorCoding=colorCoding.replace("$COLOR_KEY", keyString)
        html=html.replace("$COLOR_CODING", colorCoding)
    else:
        html=html.replace("$COLOR_CODING", "")
        
    outFile=file(outDir+os.path.sep+"index.html", "w")
    outFile.write(html)
    outFile.close()
    
#-------------------------------------------------------------------------------------------------------------
def makeSourcePages(catalog, outDir, keysToWrite, keyFormats, keyLabels, imageDirs, imageLabels, 
                    imageCaptions, sizeArcmin = 12.0, nedDir = "NEDResults", addSDSSRedshifts = True):
    """Makes webpages for each source - one page for each type of imageDir. These get written in the 
    imageDirs (keeps things tidy).
    
    """
            
    templatePage="""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
    <html>
    <head>
        <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
        <title>$ACT_NAME</title>
    </head>
    <body style="font-family: sans-serif; vertical align: top; justify: full;">
    <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
        <tbody>
            <tr>
                $PREV_LINK_CODE
                <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 125%;">
                    $ACT_NAME
                </td>
                $NEXT_LINK_CODE
            </tr>
        </tbody>
    </table>
    
    <table frame=border cellspacing=0 cols=1 rules=None border=0 width=100%>
    <tbody>
    
    <tr>
        <td align=center>$LINK_CODE</td>
    </tr>
    
    <tr>
        <td align=center><a href="$IMAGE_PATH"><img src="$IMAGE_PATH" align="middle" border=2 height="$PLOT_DISPLAY_HEIGHT_PIX"></a>
        </td>
    </tr>
    <tr><td align=center><i>$CAPTION</i></td></tr>

    <tr>
        <td align=center>$NED_MATCHES_TABLE</td>
    </tr>

    <tr>
        <td align=center>$SDSS_MATCHES_TABLE</td>
    </tr>
    
    <tr>
        <td align=center>$PROPERTIES_TABLE</td>
    </tr>
    
    </tbody>
    </table>
    <br><br>
    </body>
    </html>
    """
    templatePage=templatePage.replace("$PLOT_DISPLAY_HEIGHT_PIX", str(PLOT_DISPLAY_HEIGHT_PIX))

    # Taken this out from above the caption line
    #<tr><td align=center><b>$SIZE_ARC_MIN' x $SIZE_ARC_MIN'</b></td></tr>

    for imgDir in imageDirs:
        for k in range(len(catalog)):
            
            obj=catalog[k]
            name=obj['name']
            pagePath=imgDir+os.path.sep+name.replace(" ", "_")+".html"
            imagePath=name.replace(" ", "_")+".jpg"
            
            html=templatePage
            html=html.replace("$ACT_NAME", name)
            html=html.replace("$IMAGE_PATH", imagePath)
            
            # Slight hack here - we may not know dimensions of external images (e.g. Felipe's)
            if imagePath.find("ExternalImages") == -1:
                html=html.replace("$SIZE_ARC_MIN", "%.1f" % (sizeArcmin))
            else:
                html=html.replace("$SIZE_ARC_MIN", "")
            
            linkCode=""
            for i, label, caption in zip(imageDirs, imageLabels, imageCaptions):
                if i == imgDir:
                    linkCode=linkCode+"<b>%s</b>" % (label)
                    html=html.replace("$CAPTION", "%s" % (caption))
                else:
                    linkPath="../"+os.path.split(i)[-1]+os.path.sep+name.replace(" ", "_")+".html"
                    linkCode=linkCode+"<a href=\"%s\">%s</a>" % (linkPath, label)
                linkCode=linkCode+" - "
            linkCode=linkCode[:-3]
            html=html.replace("$LINK_CODE", linkCode)
            
            # Previous and next object page links
            if k > 0:
                prev=catalog[k-1]
            else:
                prev=None
            if k < len(catalog)-1:
                next=catalog[k+1]
            else:
                next=None
            if prev != None:
                prevLinkCode="""<td style="background-color: rgb(0, 0, 0); font-family: sans-serif; 
                color: rgb(255, 255, 255); text-align: center; vertical-align: middle; font-size: 125%;">
                <a href=$PREV_LINK><b><<</b></a></td>
                """ 
                prevLinkCode=prevLinkCode.replace("$PREV_LINK", prev['name'].replace(" ", "_")+".html")
            else:
                prevLinkCode=""
            if next != None:
                nextLinkCode="""<td style="background-color: rgb(0, 0, 0); font-family: sans-serif; 
                color: rgb(255, 255, 255); text-align: center; vertical-align: middle; font-size: 125%;">
                <a href=$NEXT_LINK><b>>></b></a></td>
                """
                nextLinkCode=nextLinkCode.replace("$NEXT_LINK", next['name'].replace(" ", "_")+".html")
            else:
                nextLinkCode=""
            html=html.replace("$PREV_LINK_CODE", prevLinkCode)
            html=html.replace("$NEXT_LINK_CODE", nextLinkCode)
                
            # NED matches table
            # We should already have the files for this from doing addNEDInfo earlier
            nedFileName=nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
            nedObjs=catalogTools.parseNEDResult(nedFileName)
            if len(nedObjs['RAs']) > 0:
                nedTable="""<br><table frame=border cellspacing=0 cols=6 rules=all border=2 width=80% align=center>
                <tbody>
                <tr>
                    <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 110%;" colspan=6>
                        <b>NED Matches</b>
                    </th>
                </tr>
                <tr>
                    <td><b>ID</b></td>
                    <td><b>Name</b></td>
                    <td><b>R.A.</b></td>
                    <td><b>Dec.</b></td>
                    <td><b>Object Type</b></td>
                    <td><b>Redshift</b></td>
                </tr>
                """                
                for i in range(len(nedObjs['RAs'])):
                    rowString="""<tr>
                        <td align=center width=10%>$ID</td>
                        <td align=center width=10%>$NAME</td>
                        <td align=center width=10%>$RA</td>
                        <td align=center width=10%>$DEC</td>
                        <td align=center width=10%>$TYPE</td>
                        <td align=center width=10%>$REDSHIFT</td>
                    </tr>
                    """
                    nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedObjs['names'][i].replace("+", "%2B").replace(" ", "+"))
                    rowString=rowString.replace("$ID", "%s" % (nedObjs['labels'][i]))
                    rowString=rowString.replace("$NAME", "<a href =%s>%s</a>" % (nedLinkURL, nedObjs['names'][i]))
                    rowString=rowString.replace("$RA", "%.5f" % (nedObjs['RAs'][i]))
                    rowString=rowString.replace("$DEC", "%.5f" % (nedObjs['decs'][i]))
                    rowString=rowString.replace("$TYPE", "%s" % (nedObjs['sourceTypes'][i]))
                    rowString=rowString.replace("$REDSHIFT", "%s" % (nedObjs['redshifts'][i]))
                    nedTable=nedTable+rowString
                nedTable=nedTable+"</tbody></table>"
            else:
                nedTable=""
            html=html.replace("$NED_MATCHES_TABLE", nedTable)

            # SDSS matches table
            if 'SDSSRedshifts' in obj.keys() and addSDSSRedshifts == True and obj['SDSSRedshifts'] != None:
                sdssTable="""<br><table frame=border cellspacing=0 cols=7 rules=all border=2 width=80% align=center>
                <tbody>
                <tr>
                    <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 110%;" colspan=5>
                        <b>SDSS Redshifts</b>
                    </th>
                </tr>
                <tr>
                    <td><b>ID</b></td>
                    <td><b>R.A.</b></td>
                    <td><b>Dec.</b></td>
                    <td><b>z</b></td>
                    <td><b>zWarning</b></td>                    
                </tr>
                """              
                sdssCount=0
                for sdssObj in obj['SDSSRedshifts']:
                    sdssCount=sdssCount+1
                    rowString="""<tr>
                        <td align=center width=10%>$ID</td>
                        <td align=center width=10%>$RA</td>
                        <td align=center width=10%>$DEC</td>
                        <td align=center width=10%>$REDSHIFT</td>
                        <td align=center width=10%>$Z_WARNING</td>
                    </tr>
                    """
                    rowString=rowString.replace("$ID", "%d" % (sdssCount))
                    rowString=rowString.replace("$RA", "%.5f" % (sdssObj['RADeg']))
                    rowString=rowString.replace("$DEC", "%.5f" % (sdssObj['decDeg']))
                    rowString=rowString.replace("$REDSHIFT", "%.3f" % (sdssObj['z']))
                    rowString=rowString.replace("$Z_WARNING", "%s" % (sdssObj['zWarning']))
                    sdssTable=sdssTable+rowString
                sdssTable=sdssTable+"</tbody></table>"
            else:
                sdssTable=""
            html=html.replace("$SDSS_MATCHES_TABLE", sdssTable)

            # Source properties table
            propTable="""<br><table frame=border cellspacing=0 cols=2 rules=all border=2 width=80% align=center>
            <tbody>
            <tr>
                <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;" colspan=2>
                    <b>Source Properties</b>
                </th>
            </tr>
            """
            for pkey, fmt, label in zip(keysToWrite, keyFormats, keyLabels):
                if pkey in obj.keys():
                    rowString="""<tr><td align=left width=50%><b>$KEY_LABEL</b></td>
                    <td align=center width=50%>$KEY_VALUE</td></tr>
                    """
                    rowString=rowString.replace("$KEY_LABEL", label)
                    if obj[pkey] == None:
                        rowString=rowString.replace("$KEY_VALUE", "-")
                    else:
                        if pkey == "NED_name":
                            nedName=obj[pkey]
                            nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedName.replace("+", "%2B").replace(" ", "+"))
                            rowString=rowString.replace("$KEY_VALUE", "<a href=%s>%s</a>" % (nedLinkURL, nedName))
                        else:
                            rowString=rowString.replace("$KEY_VALUE", fmt % (obj[pkey]))
                    propTable=propTable+rowString
            propTable=propTable+"</td></tr></tbody></table>"
            html=html.replace("$PROPERTIES_TABLE", propTable)
            
            outFile=file(pagePath, "w")
            outFile.write(html)
            outFile.close() 

#-------------------------------------------------------------------------------------------------------------
def makeImageMapPage(catalog, outDir, imageMapOptions):
    """Makes a page with a clickable image map of candidates.
    
    """
    
    # We'll allow this map to be downloaded
    imgPath=imageMapOptions['imageDict'][imageMapOptions['filteredMapLabel']]['ycFilteredMap']
    os.system("cp %s %s/" % (imgPath, outDir))
              
    # The plot - make a plain version, one with labelled circles, one without
    img=pyfits.open(imgPath)
    wcs=astWCS.WCS(imgPath)
    mapData=img[0].data
    outImagePath_plain=outDir+os.path.sep+"imageMap_plain.jpg"
    outImagePath_labelledCircles=outDir+os.path.sep+"imageMap.jpg"
    outImagePath_unlabelledCircles=outDir+os.path.sep+"imageMap_noLabels.jpg"
    sigma=numpy.std(mapData)
    numSigma=4.0
    radiusArcmin=4.0

    pylab.figure(figsize = (mapData.shape[1]/100.0, mapData.shape[0]/100.0), dpi = 100)
    p=astPlots.ImagePlot(mapData, wcs, axes=[0., 0., 1., 1.], cutLevels = [-numSigma*sigma, numSigma*sigma],
                         axesLabels = None, colorBar = False)
    p.save(outImagePath_plain)
        
    objRAs=[]
    objDecs=[]
    objLabels=[]
    for obj in catalog:
        objRAs.append(obj['RADeg'])
        objDecs.append(obj['decDeg'])
        objLabels.append(obj['name'])
    p.addPlotObjects(objRAs, objDecs, 'candidates', symbol='circle', size = radiusArcmin*2*60.0, width = 1.0, 
                     color = 'yellow', objLabels = objLabels, objLabelSize = 12.0)
    p.save(outImagePath_labelledCircles)
    p.addPlotObjects(objRAs, objDecs, 'candidates', symbol='circle', size = radiusArcmin*2*60.0, width = 1.0, 
                     color = 'yellow', objLabels = None, objLabelSize = 12.0)
    p.save(outImagePath_unlabelledCircles)
    pylab.close()
    
    # Plot of the point source mask
    psMaskPath=os.path.split(os.path.split(imgPath)[0])[0]+os.path.sep+"diagnostics"+os.path.sep+"psMask_148.fits"
    makePSImageMap=False
    if os.path.exists(psMaskPath) == True:
        makePSImageMap=True
        maskImg=pyfits.open(psMaskPath)
        psMaskData=maskImg[0].data
        outImagePath_psMask=outDir+os.path.sep+"imageMap_psMask.jpg"
        pylab.figure(figsize = (psMaskData.shape[1]/100.0, psMaskData.shape[0]/100.0), dpi = 100)
        p=astPlots.ImagePlot(psMaskData, wcs, axes=[0., 0., 1., 1.], cutLevels = [0, 1],
                            axesLabels = None, colorBar = False)
        p.addPlotObjects(objRAs, objDecs, 'candidates', symbol='circle', size = radiusArcmin*2*60.0, width = 1.0, 
                        color = 'yellow', objLabels = objLabels, objLabelSize = 12.0)
        p.save(outImagePath_psMask)
        pylab.close()
    
    # The pages
    templatePage="""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
    <html>
    <head>
        <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
        <title>Filtered Map ($FILTER_LABEL)</title>
    </head>
    <body style="font-family: sans-serif; vertical align: top; justify: full;">
    <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: $PIX_SIZEpx;">
        <tbody>
            <tr>
                <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 125%;">
                    Filtered Map ($FILTER_LABEL)
                </td>
            </tr>
        </tbody>
    </table>
    
    <table frame=border cellspacing=0 cols=1 rules=None border=0 width=$PIX_SIZEpx>
    <tbody>    
    
    $LINK_CODE
    $PLOT_CODE
    $MAP_CODE
    $DOWNLOAD_CODE
        
    </tbody>
    </table>
    <br><br>
    </body>
    </html>
    """
    
    # Title replacement
    templatePage=templatePage.replace("$FILTER_LABEL", imageMapOptions['filteredMapLabel'])
    templatePage=templatePage.replace("$PIX_SIZE", str(mapData.shape[1]))
    
    # Add links to image maps with/without circles, labels (for folks to download) - and also a .fits map
    jpgPaths=[outImagePath_labelledCircles, outImagePath_unlabelledCircles, outImagePath_plain]
    linkLabels=["Show Name Labels", "Show Candidates", "Hide Candidates"]
    pagePaths=["imageMap.html", "imageMap_unlabelled.html", "imageMap_plain.html"]
    if makePSImageMap == True:
        jpgPaths.append(outImagePath_psMask)
        linkLabels.append("Point Source Mask")
        pagePaths.append("imageMap_psMask.html")
        
    for pageJpgPath, pageLabel, pagePath in zip(jpgPaths, linkLabels, pagePaths):
        
        mapPage=templatePage
        
        linkCode=""
        for jpgPath, label, linkPath in zip(jpgPaths, linkLabels, pagePaths):
            if label == pageLabel:
                linkCode=linkCode+"<b>%s</b>" % (label)
            else:
                linkCode=linkCode+"<a href=\"%s\">%s</a>" % (linkPath, label)
            linkCode=linkCode+" - "
        linkCode=linkCode[:-3]
        mapPage=mapPage.replace("$LINK_CODE", linkCode)
            
        # Add the plot
        plotCode='<img src="%s" width="%d" height="%d" alt="FilteredMap" usemap="#filteredMap">' % (os.path.split(pageJpgPath)[-1], mapData.shape[1], mapData.shape[0])
        mapPage=mapPage.replace("$PLOT_CODE", plotCode)

        mapCode='<map name="filteredMap">\n'
        radiusPix=int(round(radiusArcmin/(wcs.getPixelSizeDeg()*60.0)))
        for obj in catalog:
            x, y=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
            x=int(round(x))
            y=mapData.shape[0]-int(round(y))
            candidateUrl="ACTImages"+os.path.sep+obj['name'].replace(" ", "_")+".html"
            mapCode=mapCode+'<area shape="circle" coords="%d,%d,%d" href="%s" alt="%s" target="_blank">\n' % (x, y, radiusPix, candidateUrl, obj['name'])
        mapCode=mapCode+"</map>\n"
        mapPage=mapPage.replace("$MAP_CODE", mapCode)

        # Add map download link
        downLink='<a href="%s">Download FITS image</a>' % (os.path.split(imgPath)[-1])
        mapPage=mapPage.replace("$DOWNLOAD_CODE", downLink)
        
        outFile=file(outDir+os.path.sep+pagePath, "w")
        outFile.write(mapPage)
        outFile.close()  
    
#-------------------------------------------------------------------------------------------------------------
def makeDiagnosticsPage(outDir, imageDict, diagnosticsDir):
    """Make a page showing plots of filters used and stats on noise, area of filtered maps etc..
    
    Could potentially dump linked sim results here also (e.g. completeness, flux recovery etc.)
    
    """
    
    plotDisplayWidthPix=800
    
    # The filtered maps diagnostics dir
    if diagnosticsDir == None:
        raise Exception, "don't know where the diagnosticsDir is!"
    
    # Make where we're going to store the sourcebrowser diagnostics page and stuff
    browserDiagnosticsDir=outDir+os.path.sep+"diagnostics"
    if os.path.exists(browserDiagnosticsDir) == False:
        os.makedirs(browserDiagnosticsDir)

    # Sort filtered map labels into sensible order
    labels=[]
    for key in imageDict:
        if key != 'optimalCatalog' and key != 'mergedCatalog':
            labels.append(key)
    labels.sort()

    candidatePagesDir=outDir+os.path.sep+"candidates"
    if os.path.exists(candidatePagesDir) == False:
        os.makedirs(candidatePagesDir)
        
    diagnosticsPage="""<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
    <html>
    <head>
        <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
        <title>Catalog Properties</title>
    </head>
    <body style="font-family: sans-serif; vertical align: top; justify: full;">
    <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
        <tbody>
            <tr>
                <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 125%;">
                    Catalog Properties
                </td>
            </tr>
        </tbody>
    </table>
    
    <table frame=border cellspacing=0 cols=1 rules=None border=0 width=100%>
    <tbody>    
    
    $PLOTS_CODE
        
    </tbody>
    </table>
    <br><br>
    </body>
    </html>
    """
    
    # Add each plot in turn, if the plot .png exists
    plotsCode=""
    
    # Cumulative SNR
    plotPath=diagnosticsDir+os.path.sep+"cumulativeSNR.png"
    plotCaption="""Cumulative number of cluster candidates above given SNR (candidates), compared to the<br>
                number found when running the source detection over inverted filtered maps (inverted maps)."""
    plotsCode=plotsCode+makeDiagnosticPlotHTMLCode(plotPath, plotCaption, plotDisplayWidthPix, 
                                                   browserDiagnosticsDir)

    # Contamination estimate vs SNR
    plotPath=diagnosticsDir+os.path.sep+"contaminationEstimateSNR.png"
    plotCaption="""Estimate of contamination above given SNR threshold. This is defined as the ratio of<br>
                the number of sources found when running the source detection over inverted maps compared<br>
                to the number of cluster candidates in the catalog (see top plot)."""
    plotsCode=plotsCode+makeDiagnosticPlotHTMLCode(plotPath, plotCaption, plotDisplayWidthPix, 
                                                   browserDiagnosticsDir)

        
    # Beta model - recovered fraction as fn. delta T
    plotPath=diagnosticsDir+os.path.sep+"fakeSourceSims"+os.path.sep+"recoveredFraction_betaModel_deltaT.png"
    plotCaption="""Recovered fraction of fake clusters (&beta; model) inserted into map as a function of<br>
                decrement size (&Delta;T)."""
    plotsCode=plotsCode+makeDiagnosticPlotHTMLCode(plotPath, plotCaption, plotDisplayWidthPix, 
                                                   browserDiagnosticsDir)    

    # Beta model - recovered fraction as fn. core radius
    plotPath=diagnosticsDir+os.path.sep+"fakeSourceSims"+os.path.sep+"recoveredFraction_betaModel_scaleArcmin.png"
    plotCaption="""Recovered fraction of fake clusters (&beta; model) inserted into map as a function of<br>
                core radius (&theta;<sub>c</sub>)."""
    plotsCode=plotsCode+makeDiagnosticPlotHTMLCode(plotPath, plotCaption, plotDisplayWidthPix, 
                                                   browserDiagnosticsDir)  
                                                   
    # Put it all together
    diagnosticsPage=diagnosticsPage.replace("$PLOTS_CODE", plotsCode)
    outFile=file(browserDiagnosticsDir+os.path.sep+"plots.html", "w")
    outFile.write(diagnosticsPage)
    outFile.close()  

#------------------------------------------------------------------------------------------------------------
def makeDiagnosticPlotHTMLCode(plotPath, plotCaption, plotDisplayWidthPix, browserDiagnosticsDir):
    """Makes HTML code for a plot on the diagnostics page.
    
    """
    
    if os.path.exists(plotPath) == True:
        os.system("cp %s %s/" % (plotPath, browserDiagnosticsDir))
        plotCode="""<tr><td align=center><br></td></tr>
        <tr>
            <td align=center><img src="$PLOT_PATH" align="middle" border=2 width="$PLOT_DISPLAY_WIDTH_PIX">
            </td>
        </tr>
        <tr><td align=center><i>$PLOT_CAPTION</i></td></tr>
        <tr><td align=center><br></td></tr>
        """
        plotCode=plotCode.replace("$PLOT_PATH", os.path.split(plotPath)[-1])
        plotCode=plotCode.replace("$PLOT_DISPLAY_WIDTH", str(plotDisplayWidthPix))
        plotCode=plotCode.replace("$PLOT_CAPTION", plotCaption)
    else:
        plotCode=""
        
    return plotCode
    
#------------------------------------------------------------------------------------------------------------
def clipSmoothedTanResampledACTImage(objDict, mapData, mapWCS, sizeDeg, gaussSmoothArcSecRadius, 
                                     outFileName = None, sizePix = 200):
    """Clips a tan resampled, (optionally smoothed) section around an object in an ACT map, writes it out
    to outFileName, and returns a dictionary containing the clipped map data and WCS. 
        
    """
    
    RADeg=objDict['RADeg']
    decDeg=objDict['decDeg']
    
    # This solves for the RA, dec coords we need to clip to get exactly the right dimensions and dodge
    # the cea projection shenanigans
    tolerance=1e-8  # in degrees on sky
    targetHalfSizeSkyDeg=(sizeDeg*1.1)/2.0  # slightly bigger, trim down afterwards
    funcCalls=["astCoords.calcAngSepDeg(RADeg, decDeg, guess, decDeg)",
                "astCoords.calcAngSepDeg(RADeg, decDeg, guess, decDeg)",
                "astCoords.calcAngSepDeg(RADeg, decDeg, RADeg, guess)",
                "astCoords.calcAngSepDeg(RADeg, decDeg, RADeg, guess)"]
    coords=[RADeg, RADeg, decDeg, decDeg]
    signs=[1.0, -1.0, 1.0, -1.0]
    results=[]
    for f, c, sign in zip(funcCalls, coords, signs):
        # Initial guess range
        maxGuess=sign*targetHalfSizeSkyDeg*2.0
        minGuess=sign*targetHalfSizeSkyDeg/10.0
        guessStep=(maxGuess-minGuess)/10.0
        guesses=numpy.arange(minGuess+c, maxGuess+c, guessStep)
        for i in range(60):
            minSizeDiff=1e6
            bestGuess=None
            for guess in guesses:
                sizeDiff=abs(eval(f)-targetHalfSizeSkyDeg)
                if sizeDiff < minSizeDiff:
                    minSizeDiff=sizeDiff
                    bestGuess=guess
            if minSizeDiff < tolerance:
                break
            else:
                guessRange=abs((maxGuess-minGuess))
                maxGuess=bestGuess+guessRange/4.0
                minGuess=bestGuess-guessRange/4.0
                guessStep=(maxGuess-minGuess)/10.0
                guesses=numpy.arange(minGuess, maxGuess, guessStep)
        results.append(bestGuess)
    RAMax=results[0]
    RAMin=results[1]
    decMax=results[2]
    decMin=results[3]

    # To tan proj, scale, smooth, save
    # We'll make these 1% bigger than 12 arcmin actually to see if we can get rid of annoying edge effect
    # in contour plots
    tanClip=astImages.clipUsingRADecCoords(mapData, mapWCS, RAMin, RAMax, decMin, decMax)
    try:
        tanClip=astImages.resampleToTanProjection(tanClip['data'], tanClip['wcs'], outputPixDimensions = [sizePix, sizePix])
    except:
        print "Hmm - tan reprojection failed"
        ipshell()
        sys.exit()
    #tanClip=astImages.clipImageSectionWCS(tanClip['data'], tanClip['wcs'], RADeg, decDeg, sizeDeg*1.01)
    dataClip=tanClip['data']
    scaleFactor=float(sizePix)/float(tanClip['data'].shape[1])
    tanClip=astImages.scaleImage(tanClip['data'], tanClip['wcs'], scaleFactor)
    if gaussSmoothArcSecRadius != None:
        radPix=(gaussSmoothArcSecRadius/3600.0)/tanClip['wcs'].getPixelSizeDeg()
        tanClip['data']=ndimage.gaussian_filter(tanClip['data'], radPix)                        
    
    if outFileName != None:
        astImages.saveFITS(outFileName, tanClip['data'], tanClip['wcs'])      
    
    return tanClip

#------------------------------------------------------------------------------------------------------------
def makeACTPlots(catalog, imageDict, plotsDir, mapOptionsDict, colorMap = "spectral", verbose = True, 
                 noAxes = False, showLabels = False, figSize = (10, 7), cutLevels = None, sizePix = 200):
    """Makes .jpg images of sizeArcmin on a side centered on sources in catalog. cutLevels is in uK (units
    of maps). Saves the tan-projected, smoothed images we make here as .fits to speed things up later.
    
    Rewritten this to plot AR1 and AR2 side by side, if asked for.
    
    Optionally smoothes with gaussian kernel of radius gaussSmoothArcSecRadius.
    
    If noAxes == True, do not plot coord axes and instead make the plot fill the plot area (for nemomontage).
    
    If showLabels == 'SN' or 'NED', put on cluster name label inside top of plot, and either SNR or NED name
    at bottom. If 'NEDOnly', only show name labels if a known cluster.
    
    """
    
    # Stuff that we've taken out of the options
    gaussSmoothArcSecRadius=mapOptionsDict['gaussSmoothArcsec']
    remakePlots=mapOptionsDict['remakePlots']
        
    # Turn off interactive plotting if it is already on, as drawing the plots on screen takes ages
    plotInteractiveStatus=pylab.isinteractive()
    if plotInteractiveStatus == True:
        pylab.matplotlib.interactive(False)
        
    # Get list of templates - assuming here that all keys that are NOT 'mergedCatalog' are template names
    templatesList=[]
    for key in imageDict.keys():
        if key != "mergedCatalog" and key != "optimalCatalog":
            templatesList.append(key)
       
    # Plot size stuff - e.g. target size to zoom plots up to
    sizeDeg=mapOptionsDict['plotSizeArcmin']/60.0
    
    # Make the plots ...
    if verbose == True: print ">>> Making ACT plots: "
    count=0
    for obj in catalog:
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        name=obj['name']
        plotFileName=plotsDir+os.path.sep+name.replace(" ", "_")+".png"
        
        if remakePlots == True or os.path.exists(plotFileName.replace(".png", ".jpg")) == False:
            # Load the appropriate map, store in these dictionaries with keys corresponding to bands
            wcsDict={}
            dataDict={}
            if mapOptionsDict['mapType'] == 'unfiltered':
                for temp in templatesList:
                    if obj['template'] == temp:
                        unfilteredMapsDictList=imageDict[temp]['unfilteredMapsDictList']
                        for mapDict in unfilteredMapsDictList:
                            if mapDict['obsFreqGHz'] in mapOptionsDict['bands']:
                                img=pyfits.open(mapDict['mapFileName'])
                                wcs=astWCS.WCS(mapDict['mapFileName'])
                                dataDict['%d' % (mapDict['obsFreqGHz'])]=img[0].data
                                wcsDict['%d' % (mapDict['obsFreqGHz'])]=wcs
            elif mapOptionsDict['mapType'] == 'filtered':
                # Load in the appropriate yc filtered template map, convert to T assuming 148 GHz
                img=pyfits.open(imageDict[obj['template']]['ycFilteredMap'])
                wcs=astWCS.WCS(imageDict[obj['template']]['ycFilteredMap'])
                dataDict['148']=mapTools.convertToDeltaT(img[0].data, obsFrequencyGHz = 148.0)
                wcsDict['148']=wcs
                #if os.path.exists("test_%s.fits" % (obj['template'])) == False:
                    #astImages.saveFITS("test_%s.fits" % (obj['template']), dataDict['148'], wcs)
            
            # For stretch of plots
            plotSigma=numpy.std(dataDict['148'])
            
            # Stuff to set up plot
            numBands=len(dataDict.keys())
            cutLevelsDict={}    # If plotting AR1 and AR2 side by side, need to have same scale
            fig=pylab.figure(num = 200, figsize = (10*numBands, 8))
            fig.canvas.set_window_title('ACT plot: %s' % obj['name'])
            name=obj['name']
            RADeg=obj['RADeg']
            decDeg=obj['decDeg']
            bandCount=0    
            if numBands == 2:
                plotWidth=0.35  
                plotStep=0.4
                separateColorBarAxis=True
                cbShrink=0.8
                cbFraction=0.2
                cbLabelOffset=0.05
            elif numBands == 1:
                plotWidth=0.8
                plotStep=0.0
                separateColorBarAxis=False
                cbShrink=1.0
                cbFraction=0.1
                cbLabelOffset=0.10
            
            # Now make plot for each band
            for band in dataDict.keys():
                if 'saveFITS' in mapOptionsDict.keys() and mapOptionsDict['saveFITS'] == True:
                    outFITSDir=os.path.split(plotsDir)[0]+os.path.sep+"FITSImages"
                    if os.path.exists(outFITSDir) == False:
                        os.makedirs(outFITSDir)                  
                    outFITSFileName=outFITSDir+os.path.sep+str(band)+"_"+name.replace(" ", "_").replace("\n", "_")+".fits"
                else:
                    outFITSFileName=None
                tanClip=clipSmoothedTanResampledACTImage(obj, dataDict[band], wcsDict[band], sizeDeg, \
                                                         gaussSmoothArcSecRadius, outFileName = outFITSFileName,
                                                         sizePix = sizePix)   
                
                # Unfiltered maps - relative scaling, with median removed so can see if local decrement
                if mapOptionsDict['mapType'] == 'unfiltered':
                    tanClip['data']=tanClip['data']-numpy.median(tanClip['data'])
                
                # Set cut levels
                if cutLevels == None:
                    if mapOptionsDict['mapType'] == 'unfiltered':
                        zeroLevel=0
                        cutLevels=[-8.0*plotSigma, plotSigma]
                    else:
                        zeroLevel=0
                        cutLevels=[-8.0*plotSigma, plotSigma]
                    if name not in cutLevelsDict.keys():
                        cutLevelsDict[name]=cutLevels
                    else:
                        cutLevels=cutLevelsDict[name]
                                
                scaleTitle="T ($\mu$K)"
                if noAxes == False:
                    axes=[0.1+bandCount*plotStep, 0.1, plotWidth, 0.8]
                    axesLabels="sexagesimal"
                else:
                    axes=[0, 0, 1, 1]
                    axesLabels=None
                try:
                    if 'figSize' in mapOptionsDict.keys(): # overrides if present what we set when calling this thing
                        figSize=mapOptionsDict['figSize']
                    fig=pylab.figure(figsize = figSize)
                    if noAxes == True: # otherwise will be backwards if no coords are plotted
                        tanClip['data']=numpy.fliplr(tanClip['data']) 
                    p=astPlots.ImagePlot(tanClip['data'], tanClip['wcs'], axes = axes, 
                                         axesLabels = axesLabels,cutLevels = cutLevels,
                                         title = "%d GHz" % (int(band)), colorMapName = colorMap)
                    if 'waterMark' in mapOptionsDict.keys() and mapOptionsDict['waterMark'] != None:
                        pylab.text(0.05, 0.05, mapOptionsDict['waterMark'], horizontalalignment = 'left', 
                                   verticalalignment = 'center', fontsize = 14, color = 'black', 
                                   transform = p.axes.transAxes)
                except:
                    print "makeACTPlots - astPlots.ImagePlot failed - eh?"
                    ipshell()
                    sys.exit()
                if 'plotClusterPos' in mapOptionsDict.keys() and mapOptionsDict['plotClusterPos'] == True:
                    p.addPlotObjects([obj['RADeg']], [obj['decDeg']], 'clusterPos', symbol='cross', \
                                      size=sizeDeg/20.0*3600.0, color='white')
                
                
                if 'plotFIRSTObjects' in mapOptionsDict.keys() and mapOptionsDict['plotFIRSTObjects'] == True:
                    firstTab=getFIRSTCatalog(obj['RADeg'], obj['decDeg'], searchRadiusArcmin = sizeDeg*60, 
                                             source = "VIII/71", cacheDir = "VizierCache")
                    if colorMap == 'spectral':
                        firstColor='black'
                    else:
                        firstColor='red'
                    p.addPlotObjects(firstTab['RADeg'], firstTab['decDeg'], 'first', size = 20, width = 2, 
                                     color = firstColor, objLabels = firstTab['name'], objLabelSize = 12)
                    
                if showLabels != False:
                    if showLabels == 'NEDOnly':
                        if 'NED_distArcmin' in obj.keys() and obj['NED_distArcmin'] < 2.0 and obj['NED_name'] != None:
                            pylab.figtext(0.5, 0.075, obj['NED_name'], ha="center", va="center", size = 40)
                            if obj['name'].find("\n") != -1:
                                yTextPos=0.87
                            else:
                                yTextPos=0.9
                        pylab.figtext(0.5, 0.9, obj['name'], ha="center", va="center", size = 60)
                    elif showLabels == 'custom':
                        pylab.figtext(0.5, 0.075, obj['bottomPlotLabel'], ha="center", va="center", size = 40)
                        pylab.figtext(0.5, 0.9, obj['topPlotLabel'], ha="center", va="center", size = 60)
                    else:
                        if obj['name'].find("\n") != -1:
                            yTextPos=0.87
                        else:
                            yTextPos=0.9
                        pylab.figtext(0.5, yTextPos, obj['name'], ha="center", va="center", size = 60)
                        if showLabels == 'NED':
                            if 'NED_distArcmin' in obj.keys() and obj['NED_distArcmin'] < 2.0 and obj['NED_name'] != None:
                                pylab.figtext(0.5, 0.075, obj['NED_name'], ha="center", va="center", size = 40)
                        elif showLabels == 'SN':
                            pylab.figtext(0.5, 0.075, "SNR=%.1f" % (obj['SN']), ha="center", va="center", size = 60)

                bandCount=bandCount+1
            
            if separateColorBarAxis == True and noAxes == False:
                cbAxes=pylab.axes([0.8, 0.0, 0.1, 1.0], frameon = False)
                pylab.xticks("", "")
                pylab.yticks("", "")
            if noAxes == False:
                cb=pylab.colorbar(shrink=cbShrink, fraction=cbFraction)
                a=cb.outline.get_axes()
                pos=a.get_position()
                pylab.figtext(pos.p0[0]+cbLabelOffset, 0.5, scaleTitle, va="center", rotation="vertical")
            pylab.savefig(plotsDir+os.path.sep+name.replace(" ", "_").replace("\n", "_")+".png")
            pylab.close()
                            
    # Convert all plots to .jpg and resize for web friendly display
    pngFiles=glob.glob(plotsDir+os.path.sep+"*.png")
    for p in pngFiles:
        j=p.replace(".png", ".jpg")
        os.system("convert -resize x%d %s %s" % (PLOT_HEIGHT_PIX, p, j))
        os.remove(p)

    pylab.matplotlib.interactive(plotInteractiveStatus)

#------------------------------------------------------------------------------------------------------------
def getFIRSTCatalog(RADeg, decDeg, searchRadiusArcmin = 10.0, source = "VIII/71", cacheDir = "VizierCache"):
    """This queries vizier for a list of object positions. The default source is the FIRST. As this
    uses the vizquery script, it won't work behind a firewall where port 80 is not open or goes through a
    proxy.
    
    """
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    outFileName=cacheDir+os.path.sep+"%.6f_%.6f_%s.cat" % (RADeg, decDeg, source.replace("/", ""))
    
    if os.path.exists(outFileName) == False:
        os.system("vizquery -mime=tsv -site=vizier.u-strasbg.fr -source=%s -c=%.6f%+.6f,rm=%.1f -out=FIRST,_RAJ2000,_DEJ2000,Fpeak >> %s" \
                  % (source, RADeg, decDeg, searchRadiusArcmin, outFileName))
    
    try:
        vizTab=atpy.Table(outFileName, type = 'ascii')
        names=[]
        for n, f in zip(vizTab['FIRST'][2:], vizTab['Fpeak'][2:]):
            names.append(n+" (%.1f mJy)" % (float(f)))
        firstTab=atpy.Table()
        firstTab.add_column("name", names)
        #firstTab.add_column("name", vizTab['FIRST'][2:])
        firstTab.add_column("RADeg", numpy.array(vizTab['_RAJ2000'][2:], dtype = float))
        firstTab.add_column("decDeg", numpy.array(vizTab['_DEJ2000'][2:], dtype = float))
    except:
        firstTab=atpy.Table()
        firstTab.add_column("name", [])
        firstTab.add_column("RADeg", numpy.array([], dtype = float))
        firstTab.add_column("decDeg", numpy.array([], dtype = float))
        
    return firstTab
    
#------------------------------------------------------------------------------------------------------------
def makeSkyviewPlots(catalog, plotsDir, sizeArcmin = 12.0, JPEGFolder = "DSSJPEGs", plotNEDObjects = True, 
                    verbose = True, contourOverlayImageDir = None, contourOverlayBand = 148, 
                    contourSmoothArcSecRadius = 15.0, figSize = (10, 7), 
                    contourLevels = None, remakePlots = False, nedDir = "NEDResults"):
    """Makes RGB plots from DSS images. Can optionally overlay contours from ACT maps.
    
    """

    plotInteractiveStatus=pylab.isinteractive()
    if plotInteractiveStatus == True:
        pylab.matplotlib.interactive(False)
        
    sizeDeg=sizeArcmin/60.0

    if contourOverlayImageDir != None:
        contourImageFiles=glob.glob(contourOverlayImageDir+os.path.sep+"%s_*.fits" % (str(contourOverlayBand)))
    else:
        contourImageFiles=[]
        contourImg=None
        contourWCS=None
    
    if verbose == True and contourOverlayImageDir == None: print ">>> Making DSS plots: "
    if verbose == True and contourOverlayImageDir != None: print ">>> Making DSS plots with contour overlays: "
    count=0
    for obj in catalog:
        count=count+1
        if verbose == True: print "... %s (%d/%d) ... " % (obj['name'], count, len(catalog))
                    
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']

        contourImg=None
        contourWCS=None
        for fileName in contourImageFiles:
            if fileName.find(name.replace(" ", "_")) != -1:
                contourImg=pyfits.open(fileName)
                contourWCS=astWCS.WCS(fileName)
                break
            
        if os.path.exists(plotsDir+os.path.sep+name.replace(" ", "_")+".jpg") == False or remakePlots == True:
                        
            name=obj['name']
        
            # Load data
            channels=["R", "G", "B"]
            loadedAll=True # if one of the images fail, abandon making plot
            channelsData=[]
            for channel in channels:
                inJPGPath=JPEGFolder+os.path.sep+obj['name'].replace(" ", "_")+"_%s.jpg" % (channel)
                if os.path.exists(inJPGPath) == True and checkImageDownloadSuccess(inJPGPath) == True:
                    im=Image.open(inJPGPath)
                    im=im.convert(mode="L")
                    imData=numpy.array(im)
                    imData=numpy.fliplr(imData)
                    imData=numpy.flipud(imData)
                    channelsData.append(imData)
                else:
                    loadedAll=False
            
            # Abort if not got all data - we should maybe do something nicer than this!
            if loadedAll == False:
                print "... couldn't find all images to make plot - aborting ..."
                continue
            
            R=channelsData[0]
            G=channelsData[1]
            B=channelsData[2]
                        
            # Make a WCS
            RADeg, decDeg=obj['RADeg'], obj['decDeg']
            xSizeDeg, ySizeDeg=sizeArcmin/60.0, sizeArcmin/60.0
            xSizePix=R.shape[1]
            ySizePix=R.shape[0]
            xRefPix=xSizePix/2.0
            yRefPix=ySizePix/2.0
            xOutPixScale=xSizeDeg/xSizePix
            yOutPixScale=ySizeDeg/ySizePix
            cardList=pyfits.CardList()
            cardList.append(pyfits.Card('NAXIS', 2))
            cardList.append(pyfits.Card('NAXIS1', xSizePix))
            cardList.append(pyfits.Card('NAXIS2', ySizePix))
            cardList.append(pyfits.Card('CTYPE1', 'RA---TAN'))
            cardList.append(pyfits.Card('CTYPE2', 'DEC--TAN'))
            cardList.append(pyfits.Card('CRVAL1', RADeg))
            cardList.append(pyfits.Card('CRVAL2', decDeg))
            cardList.append(pyfits.Card('CRPIX1', xRefPix+1))
            cardList.append(pyfits.Card('CRPIX2', yRefPix+1))
            cardList.append(pyfits.Card('CDELT1', xOutPixScale))
            cardList.append(pyfits.Card('CDELT2', xOutPixScale))    # Makes more sense to use same pix scale
            cardList.append(pyfits.Card('CUNIT1', 'DEG'))
            cardList.append(pyfits.Card('CUNIT2', 'DEG'))
            newHead=pyfits.Header(cards=cardList)
            wcs=astWCS.WCS(newHead, mode='pyfits')
        
            # new
            #R=R/R.std()
            #G=G/G.std()
            #B=B/B.std()
            #RCut=[R.min(), 10.0]#[R.min(), R.max()]
            #GCut=[G.min(), 10.0] #[G.min(), G.max()]
            #BCut=[B.min(), 10.0]#[B.min(), B.max()]
            
            # old
            minCutR=R.mean()-R.std()
            maxCutR=R.mean()+2.0*R.std()
            minCutG=G.mean()-G.std()
            maxCutG=G.mean()+2.0*G.std()
            minCutB=B.mean()-B.std()
            maxCutB=B.mean()+2.0*B.std()
            RCut=[minCutR, maxCutR]
            GCut=[minCutG, maxCutG]
            BCut=[minCutB, maxCutB]
            
            fig=pylab.figure(figsize = figSize)
            p=astPlots.ImagePlot([R, G, B], wcs, cutLevels=[RCut, GCut, BCut], title = name)
            p.addPlotObjects([RADeg], [decDeg], 'clusterPos', symbol='cross', size=sizeDeg/20.0*3600.0, color='white')
            
            if plotNEDObjects == True:
                # We should already have the files for this from doing addNEDInfo earlier
                nedFileName=nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
                nedObjs=catalogTools.parseNEDResult(nedFileName)
                if len(nedObjs['RAs']) > 0:
                    p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                        size = sizeDeg/40.0*3600.0, color = "#7cfc00")
        
            if contourImg != None:
                p.addContourOverlay(contourImg[0].data, contourWCS, 'actContour', levels = contourLevels, width = 2,     
                                        color = 'yellow', highAccuracy = False)

            pylab.savefig(plotsDir+os.path.sep+name.replace(" ", "_")+".png")
            pylab.close()
    
    # Convert all plots to .jpg
    pngFiles=glob.glob(plotsDir+os.path.sep+"*.png")
    for p in pngFiles:
        j=p.replace(".png", ".jpg")
        os.system("convert -resize x%d %s %s" % (PLOT_HEIGHT_PIX, p, j))
        os.remove(p)

    pylab.matplotlib.interactive(plotInteractiveStatus)
    
#------------------------------------------------------------------------------------------------------------
def fetchSkyviewImages(catalog, JPEGFolder = "DSSJPEGs", surveys = ["DSS2IR", "DSS2R", "DSS2B"], 
                       sizeArcmin = 6.0, refetch = False, verbose = True):
    """Fetches images from Skyview from the given surveys (in order R, G, B) for each object in the catalog 
    and stores them under JPEGFolder. If refetch == False, does not fetch images for which .jpg images already
    exist under JPEGFolder
    
    """
    
    if os.path.exists(JPEGFolder) == False:
        os.makedirs(JPEGFolder)
    sizeDeg=sizeArcmin/60.0
    
    if verbose == True: print ">>> Fetching Skyview images %s: " % (surveys)
    count=0
    channels=["R", "G", "B"]
    for obj in catalog:
        count=count+1
        if verbose == True: 
            print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
        sleepiness=10
        retries=1   # don't flog a dead horse
        for channel, survey in zip(channels, surveys):
            outFileName=JPEGFolder+os.path.sep+name.replace(" ", "_")+"_%s.jpg" % (channel)
            if refetch == True or checkImageDownloadSuccess(outFileName) == False:
                success=False
                attempts=0
                while success == False and attempts < retries:
                    urlString="http://skyview.gsfc.nasa.gov/cgi-bin/images?Position="+str(RADeg)+","+str(decDeg)
                    urlString=urlString+"&Size="+str(sizeDeg)+"&Pixels=1000&Projection=Tan&Survey="
                    urlString=urlString+survey+"&Coordinates=J2000&Return=jpeg"
                    urllib.urlretrieve(urlString, filename = outFileName)
                    success=checkImageDownloadSuccess(outFileName)
                    if success == False:
                        print "... failed to fetch image - trying again in %d sec ..." % (sleepiness)
                        time.sleep(sleepiness)  # be kind to Skyview
                        attempts=attempts+1

#------------------------------------------------------------------------------------------------------------
def checkImageDownloadSuccess(fileName):
    """Checks if the image fileName was successfully downloaded or not. If it wasn't, we'll have an html 
    file masquerading as a .jpg (this may need tweaking to work with e.g. SDSS, works for Skyview),
    
    """
    
    if os.path.exists(fileName) == True:
        inFile=file(fileName, "r")
        line=inFile.readline()
        inFile.close()
        if string.lower(line).find("html") == -1:
            success=True
        else:
            success=False
    else:
        success=False
    
    return success  
    
#------------------------------------------------------------------------------------------------------------
def makeBCSPlots(catalog, plotsDir, opticalOptionsDict, plotSizeArcmin = 12.0, nedDir = "NEDResults", 
                 contourOverlayImageDir = None, contourOverlayBand = None, remakePlots = False):
    """Makes RGB plots from BCS images. Will eventually optionally overlay contours from ACT maps.
    
    """

    # Force turning off interactive plotting, as slows things down lots if draw on screen
    plotInteractiveStatus=pylab.isinteractive()
    if plotInteractiveStatus == True:
        pylab.matplotlib.interactive(False)
        
    if contourOverlayImageDir != None:
        contourImageFiles=glob.glob(contourOverlayImageDir+os.path.sep+"%s_*.fits" % (str(contourOverlayBand)))
    else:
        contourImageFiles=[]
        contourImg=None
        contourWCS=None
    
    rBCSImages=glob.glob(opticalOptionsDict['BCSImagesDir']+os.path.sep+"*"+os.path.sep+"*r.fits")

    count=0
    for obj in catalog:
        count=count+1
        print "... %s (%d/%d) ... " % (obj['name'], count, len(catalog))
                    
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']

        contourImg=None
        contourWCS=None
        for fileName in contourImageFiles:
            if fileName.find(name.replace(" ", "_")) != -1:
                contourImg=pyfits.open(fileName)
                contourWCS=astWCS.WCS(fileName)
                break
        
        if os.path.exists(plotsDir+os.path.sep+name.replace(" ", "_")+".jpg") == False or remakePlots == True:
                        
            # Is this image in the BCS? If so, which field?
            inBCSFile=None
            for bcsFile in rBCSImages:
                wcs=astWCS.WCS(bcsFile)
                if wcs.coordsAreInImage(obj['RADeg'], obj['decDeg']) == True:
                    inBCSFile=bcsFile
            
            if inBCSFile != None:
                print "... making BCS RGB plot (%s) ..." % (obj['name'])
                RGBDict=opticalOptionsDict['BCSRGBDict']
                for key in RGBDict.keys():
                    RGBDict[key]['fileName']=inBCSFile.replace("r.fits", RGBDict[key]['band']+".fits")            
                outFileName=plotsDir+os.path.sep+name.replace(" ", "_")+".png"
                makeRGBPlot(obj, outFileName, RGBDict, plotSizeArcmin, nedDir = nedDir)
                
    # Convert all plots to .jpg
    pngFiles=glob.glob(plotsDir+os.path.sep+"*.png")
    for p in pngFiles:
        j=p.replace(".png", ".jpg")
        os.system("convert -resize x%d %s %s" % (PLOT_HEIGHT_PIX, p, j))
        os.remove(p)

    pylab.matplotlib.interactive(plotInteractiveStatus)

#------------------------------------------------------------------------------------------------------------
def makeRGBPlot(objDict, outFileName, RGBDict, plotSizeArcmin, plotNEDObjects = True, 
                nedDir = "NEDResults", contourImg = None, figSize = (10, 7)):
    """Makes an RGB plot from the info in the RGBDict, which must contain dictionaries in R, G, B keys.
    Each of these must have keys 'cutLevels' and 'fileName'. 'cutLevels' can be set to "auto", in which
    case we use a simple algorithm to come up with a set of cuts that look okay-ish. We clip the image to
    the given dimensions around the object position.
    
    We assume here that the images are already registered and have N at top, E at left.
    
    """
    
    sizeDeg=plotSizeArcmin/60.0
    
    channels=['R', 'G', 'B']
    cutLevels=[]
    dataList=[]
    clippedSection=None
    wcs=None
    for channel in channels:
        img=pyfits.open(RGBDict[channel]['fileName'])
        channelWCS=astWCS.WCS(RGBDict[channel]['fileName'])
        data=img[0].data
        if clippedSection == None:
            clip=astImages.clipImageSectionWCS(data, channelWCS, objDict['RADeg'], objDict['decDeg'], \
                                               plotSizeArcmin/60.0)
            clippedSection=clip['clippedSection']
            dataList.append(clip['data'])
            wcs=clip['wcs']
        else:
            dataList.append(data[clippedSection[2]:clippedSection[3], clippedSection[0]:clippedSection[1]])
        if RGBDict[channel]['cutLevels'] == 'auto':
            minCut=data[-1].mean()-data[-1].std()
            maxCut=data[-1].mean()+2.0*data[-1].std()
            cutLevels.append([minCut, maxCut])
        else:
            cutLevels.append(RGBDict[channel]['cutLevels'])

    fig=pylab.figure(figsize = figSize)
    p=astPlots.ImagePlot(dataList, wcs, cutLevels = cutLevels, title = objDict['name'])
    p.addPlotObjects([objDict['RADeg']], [objDict['decDeg']], 'clusterPos', symbol='cross', \
                      size=sizeDeg/20.0*3600.0, color='white')
    
    if plotNEDObjects == True:
        # We should already have the files for this from doing addNEDInfo earlier
        nedFileName=nedDir+os.path.sep+objDict['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName)
        if len(nedObjs['RAs']) > 0:
            p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                size = sizeDeg/40.0*3600.0, color = "#7cfc00")

    if contourImg != None:
        p.addContourOverlay(contourImg[0].data, contourWCS, 'actContour', levels = contourLevels, width = 2,     
                                color = 'yellow', highAccuracy = False)

    pylab.savefig(outFileName)
    pylab.close()

#-------------------------------------------------------------------------------------------------------------
def fetchStripe82ColorImages(catalog, sizeArcmin = 12.0, jpegFolder = "jpegImages", refetch = False):
    """This uses the S82Grabber class to fetch the .fits images directly from the DAS, then mosaic them, 
    before trimming them down to size and saving them as .jpgs like we would if we just grabbed .jpgs in the
    first place.
    
    """
   
    grabber=S82Grabber.S82Grabber()

    if os.path.exists(jpegFolder) == False:
        os.makedirs(jpegFolder)
    
    # in r, g, b order
    bands=['i', 'r', 'g']
    
    print ">>> Fetching and making Stripe 82 color images ..."
    count=0
    for obj in catalog:
 
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
 
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
                
        outFileName=jpegFolder+os.path.sep+name.replace(" ", "_")+".jpg"      
        
        if os.path.exists(outFileName) == False:
            
            fieldsList=grabber.getNearestFieldsList(RADeg, decDeg, number = 5)
            # chuck in to a 'no image' image if not in stripe 82
            if len(fieldsList) == 0:
                noDataImagePath=nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
                os.system("cp %s %s" % (noDataImagePath, outFileName))
            else:               
                for f in fieldsList:
                    grabber.fetchImages(f)
                
                # Store everything (including mosaiced .fits images) under output jpegs dir for now
                # Fail gracefully (i.e. with no data image) if mosaicing goes wrong. It certainly does when
                # crossing 0h RA
                rgb=[]
                cutLevelsList=[[-5, 50], [-5, 40], [-5, 30]]
                fullWCS=None
                for b, cutLevels in zip(bands, cutLevelsList):
                    bandFITSFileName=outFileName.replace(".jpg", "_%s.fits" % (b))
                    if os.path.exists(bandFITSFileName) == False:
                        grabber.mosaicFields(fieldsList, bandFITSFileName, band = b)
                    if os.path.exists(bandFITSFileName) == True:
                        if fullWCS == None:
                            fullWCS=astWCS.WCS(bandFITSFileName)
                        img=pyfits.open(bandFITSFileName)
                        clip=astImages.clipImageSectionWCS(img[0].data, fullWCS, RADeg, decDeg, sizeArcmin/60.0)
                        rotated=ndimage.rotate(clip['data'], -90)
                        rotated=astImages.normalise(rotated, cutLevels)
                        rotated=rotated.transpose()
                        rgb.append(rotated)
                        #os.remove(bandFITSFileName) # save on disk space
                
                rgb=numpy.array(rgb)
                rgb=rgb.transpose()
                if len(rgb) > 0:             
                    fig=pylab.figure(figsize=(12, 12))   
                    pylab.axes([0,0,1,1])
                    try:
                        pylab.imshow(rgb, interpolation="bilinear", origin='lower')
                    except:
                        print "What happened here?"
                        ipshell()
                        sys.exit()
                    pylab.savefig(outFileName.replace(".jpg", ".png"))
                    pylab.close()
                
                    # We make sure the .jpg is at the native SDSS resolution
                    os.system("convert -resize %dx %s %s" % (rgb.shape[1], outFileName.replace(".jpg", ".png"), outFileName))

                else:
                    noDataImagePath=nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
                    os.system("cp %s %s" % (noDataImagePath, outFileName))
                    
#-------------------------------------------------------------------------------------------------------------
def fetchStripe82Images(catalog, sizeArcmin = 12.0, jpegFolder = "jpegImages", refetch = False):
    """Fetches the Stripe 82 .jpgs for the given image size using the casjobs webservice, stores them under
    jpegFolder. makeSDSSPlots loads these jpegs in, and use matplotlib to make them into plots with
    coord axes etc.
    
    """

    if os.path.exists(jpegFolder) == False:
        os.makedirs(jpegFolder)
    
    print ">>> Fetching Stripe 82 images ..."
    count=0
    for obj in catalog:
 
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
 
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
                
        outFileName=jpegFolder+os.path.sep+name.replace(" ", "_")+".jpg"
        
        SDSSWidth=int(round((sizeArcmin*60.0)/0.396127))
        if os.path.exists(outFileName) == False or refetch == True:
            try:
                urllib.urlretrieve("http://casjobs.sdss.org/ImgCutoutStripe82/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)+"&scale=0.396127&width="+str(SDSSWidth)+"&height="+str(SDSSWidth), filename = outFileName)
            except:
                print "... WARNING: couldn't get Stripe 82 image ..."
                outFileName=None

#------------------------------------------------------------------------------------------------------------
def makeStripe82Plots(catalog, plotsDir, opticalOptionsDict, plotSizeArcmin = 12.0, remakePlots = False, 
                  plotNEDObjects = True, nedDir = "NEDResults", plotSDSSObjects = True, 
                  SDSSQueryFolder = "SDSSQueryResults", figSize = (10, 7)):
    """Makes astPlot plots out of mono Stripe 82 JPEG images, overplots NED objects + SDSS spec objects.
    
    """
    
    print ">>> Making Stripe 82 plots ..."
    plotInteractiveStatus=pylab.isinteractive()
    if plotInteractiveStatus == True:
        pylab.matplotlib.interactive(False)
    
    sizeDeg=plotSizeArcmin/60.0
        
    count=0
    for obj in catalog:
        
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        
        name=obj['name']
        
        # Load data
        inJPGPath=opticalOptionsDict['Stripe82ImagesDir']+os.path.sep+obj['name'].replace(" ", "_")+".jpg"
        if os.path.exists(inJPGPath) == False:
            print "... couldn't find .jpg image - aborting ..."
            continue
        im=Image.open(inJPGPath)
        im=im.convert(mode="L")
        data=numpy.array(im)
        data=numpy.flipud(data)
        data=numpy.fliplr(data)
        
        # Make a WCS
        RADeg, decDeg=obj['RADeg'], obj['decDeg']
        xSizeDeg, ySizeDeg=plotSizeArcmin/60.0, plotSizeArcmin/60.0
        xSizePix=data.shape[1]
        ySizePix=data.shape[0]
        xRefPix=xSizePix/2.0
        yRefPix=ySizePix/2.0
        xOutPixScale=xSizeDeg/xSizePix
        yOutPixScale=ySizeDeg/ySizePix
        cardList=pyfits.CardList()
        cardList.append(pyfits.Card('NAXIS', 2))
        cardList.append(pyfits.Card('NAXIS1', xSizePix))
        cardList.append(pyfits.Card('NAXIS2', ySizePix))
        cardList.append(pyfits.Card('CTYPE1', 'RA---TAN'))
        cardList.append(pyfits.Card('CTYPE2', 'DEC--TAN'))
        cardList.append(pyfits.Card('CRVAL1', RADeg))
        cardList.append(pyfits.Card('CRVAL2', decDeg))
        cardList.append(pyfits.Card('CRPIX1', xRefPix+1))
        cardList.append(pyfits.Card('CRPIX2', yRefPix+1))
        cardList.append(pyfits.Card('CDELT1', xOutPixScale))
        cardList.append(pyfits.Card('CDELT2', xOutPixScale))    # Makes more sense to use same pix scale
        cardList.append(pyfits.Card('CUNIT1', 'DEG'))
        cardList.append(pyfits.Card('CUNIT2', 'DEG'))
        newHead=pyfits.Header(cards=cardList)
        wcs=astWCS.WCS(newHead, mode='pyfits')
        
        # Make plot
        fig=pylab.figure(figsize = figSize)
        p=astPlots.ImagePlot(data, wcs, cutLevels=[data.min(), data.max()], title = name)
        p.addPlotObjects([RADeg], [decDeg], 'clusterPos', symbol='cross', size=sizeDeg/20.0*3600.0, color='white')
                
        if plotNEDObjects == True:
            # We should already have the files for this from doing addNEDInfo earlier
            nedFileName=nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
            nedObjs=catalogTools.parseNEDResult(nedFileName)
            if len(nedObjs['RAs']) > 0:
                p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                    size = sizeDeg/40.0*3600.0, color = "#7cfc00")
    
        if plotSDSSObjects == True:
            if 'SDSSRedshifts' in obj.keys() and obj['SDSSRedshifts'] != None:
                sdssRAs=[]
                sdssDecs=[]
                sdssLabels=[]
                sdssCount=0
                for sdssObj in obj['SDSSRedshifts']:
                    sdssCount=sdssCount+1
                    sdssRAs.append(sdssObj['RADeg'])
                    sdssDecs.append(sdssObj['decDeg'])
                    sdssLabels.append(str(sdssCount))
                if len(sdssRAs) > 0:
                    p.addPlotObjects(sdssRAs, sdssDecs, 'sdssObjects', objLabels = sdssLabels,
                                     size = sizeDeg/40.0*3600.0, symbol = 'box', color = "red")
                                     
        #if contourImg != None:
            #p.addContourOverlay(contourImg[0].data, contourWCS, 'actContour', levels = contourLevels, width = 2,     
                                    #color = 'yellow', highAccuracy = False)

        pylab.savefig(plotsDir+os.path.sep+name.replace(" ", "_")+".png")
        pylab.close()
    
    # Convert all plots to .jpg
    pngFiles=glob.glob(plotsDir+os.path.sep+"*.png")
    for p in pngFiles:
        j=p.replace(".png", ".jpg")
        os.system("convert -resize x%d %s %s" % (PLOT_HEIGHT_PIX, p, j))
        os.remove(p)

    pylab.matplotlib.interactive(plotInteractiveStatus)

#------------------------------------------------------------------------------------------------------------
def fetchSDSSDR7Images(catalog, sizeArcmin = 12.0, JPEGFolder = "SDSSDR7Images", refetch = False):
    """Fetches the SDSS 82 .jpgs for the given image size using the casjobs webservice, stores them under
    jpegFolder. makeSDSSPlots loads these jpegs in, and use matplotlib to make them into plots with
    coord axes etc.
    
    """

    if os.path.exists(JPEGFolder) == False:
        os.makedirs(JPEGFolder)
    
    print ">>> Fetching SDSS DR7 images ..."
    count=0
    for obj in catalog:
 
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
 
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
                
        outFileName=JPEGFolder+os.path.sep+name.replace(" ", "_")+".jpg"
        
        # old
        #SDSSWidth=int(round((sizeArcmin*60.0)/0.396127))
        #SDSSScale=0.396127

        # new - as reducing size of image, no point in getting them at ~6 Mb a pop
        SDSSWidth=1200.0
        SDSSScale=(sizeArcmin*60.0)/SDSSWidth # 0.396127
        
        if os.path.exists(outFileName) == False or refetch == True:
            try:
                urlString="http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)
                urlString=urlString+"&scale="+str(SDSSScale)+"&width="+str(int(SDSSWidth))+"&height="+str(int(SDSSWidth))
                urllib.urlretrieve(urlString, filename = outFileName)
            except:
                print "... WARNING: couldn't get SDSS DR7 image ..."
                outFileName=None

#------------------------------------------------------------------------------------------------------------
def fetchSDSSDR8Images(catalog, sizeArcmin = 12.0, JPEGFolder = "SDSSDR8Images", refetch = False):
    """Fetches the SDSS 82 .jpgs for the given image size using the casjobs webservice, stores them under
    jpegFolder. makeSDSSPlots loads these jpegs in, and use matplotlib to make them into plots with
    coord axes etc.
    
    """

    if os.path.exists(JPEGFolder) == False:
        os.makedirs(JPEGFolder)
    
    print ">>> Fetching SDSS DR8 images ..."
    count=0
    for obj in catalog:
 
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
 
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
                
        outFileName=JPEGFolder+os.path.sep+name.replace(" ", "_")+".jpg"
        
        # old
        #SDSSWidth=int(round((sizeArcmin*60.0)/0.396127))
        #SDSSScale=0.396127

        # new - as reducing size of image, no point in getting them at ~6 Mb a pop
        SDSSWidth=1200.0
        SDSSScale=(sizeArcmin*60.0)/SDSSWidth # 0.396127
        
        if os.path.exists(outFileName) == False or refetch == True:
            try:
                urlString="http://skyservice.pha.jhu.edu/DR8/ImgCutout/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)
                urlString=urlString+"&scale="+str(SDSSScale)+"&width="+str(int(SDSSWidth))+"&height="+str(int(SDSSWidth))
                urllib.urlretrieve(urlString, filename = outFileName)
            except:
                print "... WARNING: couldn't get SDSS DR8 image ..."
                outFileName=None
                
#-------------------------------------------------------------------------------------------------------------
def fetchCFHTImages(catalog, sizeArcmin = 12.0, JPEGFolder = "CFHTLSImages", refetch = False):
    """Retrieves .jpgs from CFHT legacy survey.
    
    """

    if os.path.exists(JPEGFolder) == False:
        os.makedirs(JPEGFolder)
    
    print ">>> Fetching CFHTLS images ..."
    count=0
    for obj in catalog:
 
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
 
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
                
        outFileName=JPEGFolder+os.path.sep+name.replace(" ", "_")+".jpg"
        
        url="http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/cutcfhtls.pl?ra=$RA&dec=$DEC&size=$SIZE_ARCMIN&units=arcminutes&wide=true&deep=true&preview=colour"
        url=url.replace("$RA", str(RADeg))
        url=url.replace("$DEC", str(decDeg))
        url=url.replace("$SIZE_ARCMIN", str(sizeArcmin))
        print url
        
        if os.path.exists(outFileName) == False:
            response=urllib2.urlopen(url)
            lines=response.read()
            lines=lines.split("\n")
            for line in lines:
                if line.find("cutout preview") != -1:
                    break
            imageURL="http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/"+line[line.find("src=")+4:].split('"')[1]
            try:
                urllib.urlretrieve(imageURL, outFileName)
            except:
                noDataImagePath=nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
                os.system("cp %s %s" % (noDataImagePath, outFileName))
            
#------------------------------------------------------------------------------------------------------------
def makeSDSSDR7Plots(catalog, plotsDir, JPEGFolder = "SDSSDR7JPEGs", sizeArcmin = 12.0, remakePlots = False, 
                  plotNEDObjects = True, nedDir = "NEDResults", plotSDSSObjects = True, 
                  SDSSQueryFolder = "SDSSQueryResults", noAxes = False, plotClusterPos = True, 
                  figSize = (10, 7)):
    """Makes astPlot plots out of SDSS JPEG images.
    
    If noAxes == True, do not plot coord axes and instead make the plot fill the plot area (for nemomontage).

    """
    
    print ">>> Making SDSS DR7, DR8 or S82 plots ..."
    plotInteractiveStatus=pylab.isinteractive()
    if plotInteractiveStatus == True:
        pylab.matplotlib.interactive(False)
    
    sizeDeg=sizeArcmin/60.0
        
    count=0
    for obj in catalog:
        
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        
        name=obj['name']
        
        # Load data
        inJPGPath=JPEGFolder+os.path.sep+obj['name'].replace(" ", "_")+".jpg"
        if os.path.exists(inJPGPath) == False:
            print "... couldn't find .jpg image - aborting ..."
            continue
        
        im=Image.open(inJPGPath)
        data=numpy.array(im)
        try:
            data=numpy.flipud(data)
            if noAxes == False:             # If we're not plotting coords, don't flip or it will be backwards
                data=numpy.fliplr(data)
        except:
            "... something odd about image (1d?) - aborting ..."
            continue
        
        R=data[:, :, 0]
        G=data[:, :, 1]
        B=data[:, :, 2]
        cutLevels=[[R.min(), R.max()], [G.min(), G.max()], [B.min(), B.max()]]
        
        # Make a WCS
        RADeg, decDeg=obj['RADeg'], obj['decDeg']
        xSizeDeg, ySizeDeg=sizeArcmin/60.0, sizeArcmin/60.0
        xSizePix=R.shape[1]
        ySizePix=R.shape[0]
        xRefPix=xSizePix/2.0
        yRefPix=ySizePix/2.0
        xOutPixScale=xSizeDeg/xSizePix
        yOutPixScale=ySizeDeg/ySizePix
        cardList=pyfits.CardList()
        cardList.append(pyfits.Card('NAXIS', 2))
        cardList.append(pyfits.Card('NAXIS1', xSizePix))
        cardList.append(pyfits.Card('NAXIS2', ySizePix))
        cardList.append(pyfits.Card('CTYPE1', 'RA---TAN'))
        cardList.append(pyfits.Card('CTYPE2', 'DEC--TAN'))
        cardList.append(pyfits.Card('CRVAL1', RADeg))
        cardList.append(pyfits.Card('CRVAL2', decDeg))
        cardList.append(pyfits.Card('CRPIX1', xRefPix+1))
        cardList.append(pyfits.Card('CRPIX2', yRefPix+1))
        cardList.append(pyfits.Card('CDELT1', xOutPixScale))
        cardList.append(pyfits.Card('CDELT2', xOutPixScale))    # Makes more sense to use same pix scale
        cardList.append(pyfits.Card('CUNIT1', 'DEG'))
        cardList.append(pyfits.Card('CUNIT2', 'DEG'))
        newHead=pyfits.Header(cards=cardList)
        wcs=astWCS.WCS(newHead, mode='pyfits')
        
        # Make plot
        fig=pylab.figure(figsize = figSize)
        if noAxes == True:
            axes=[0, 0, 1, 1]
            axesLabels=None
        else:
            axes=[0.1,0.1,0.8,0.8]
            axesLabels="sexagesimal"
        p=astPlots.ImagePlot([R, G, B], wcs, cutLevels = cutLevels, title = name, axes = axes, 
                             axesLabels = axesLabels)
        if plotClusterPos == True:
            p.addPlotObjects([RADeg], [decDeg], 'clusterPos', symbol='cross', size=sizeDeg/20.0*3600.0, color='white')
                
        if plotNEDObjects == True:
            # We should already have the files for this from doing addNEDInfo earlier
            nedFileName=nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
            nedObjs=catalogTools.parseNEDResult(nedFileName)
            if len(nedObjs['RAs']) > 0:
                p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                    size = sizeDeg/40.0*3600.0, color = "#7cfc00")
    
        if plotSDSSObjects == True:
            if 'SDSSRedshifts' in obj.keys() and obj['SDSSRedshifts'] != None:
                sdssRAs=[]
                sdssDecs=[]
                sdssLabels=[]
                sdssCount=0
                for sdssObj in obj['SDSSRedshifts']:
                    sdssCount=sdssCount+1
                    sdssRAs.append(sdssObj['RADeg'])
                    sdssDecs.append(sdssObj['decDeg'])
                    sdssLabels.append(str(sdssCount))
                if len(sdssRAs) > 0:
                    p.addPlotObjects(sdssRAs, sdssDecs, 'sdssObjects', objLabels = sdssLabels,
                                     size = sizeDeg/40.0*3600.0, symbol = 'box', color = "red")
                                     
        #if contourImg != None:
            #p.addContourOverlay(contourImg[0].data, contourWCS, 'actContour', levels = contourLevels, width = 2,     
                                    #color = 'yellow', highAccuracy = False)

        pylab.savefig(plotsDir+os.path.sep+name.replace(" ", "_")+".png")
        pylab.close()
    
    # Convert all plots to .jpg
    pngFiles=glob.glob(plotsDir+os.path.sep+"*.png")
    for p in pngFiles:
        j=p.replace(".png", ".jpg")
        os.system("convert -resize x%d %s %s" % (PLOT_HEIGHT_PIX, p, j))
        os.remove(p)

    pylab.matplotlib.interactive(plotInteractiveStatus)

#------------------------------------------------------------------------------------------------------------
def makeCFHTPlots(catalog, plotsDir, JPEGFolder = "CFHTLSJPEGs", sizeArcmin = 12.0, remakePlots = False, 
                  plotNEDObjects = True, nedDir = "NEDResults", noAxes = False, plotClusterPos = True, 
                  figSize = (10, 7)):
    """Makes astPlot plots out of CFHT JPEG images.
    
    NOTE: For now, actually this is just copying files from one dir to another (until we fiddle and work
    out pixel scale of CFHTLS .jpgs)
    
    If noAxes == True, do not plot coord axes and instead make the plot fill the plot area (for nemomontage).

    """

    print ">>> Making CFHTLS plots ..."
    count=0
    for obj in catalog:
        
        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        
        name=obj['name']
        
        # Load data
        inJPGPath=JPEGFolder+os.path.sep+obj['name'].replace(" ", "_")+".jpg"
        outPath=plotsDir+os.path.sep+name.replace(" ", "_")+".jpg"
        os.system("cp %s %s" % (inJPGPath, outPath))
    
#------------------------------------------------------------------------------------------------------------
def montage(catalog, numCols, numRows, outFileName, imageDir):
    """Makes a single figure montage using images found in imageDir.
    
    """
    
    if numCols*numRows < len(catalog):
        raise Exception, "numRows, numCols isn't large enough to fit all objects into montage plot"

    # So that fill columns first
    #plotNumbers=[]
    #for c in range(numCols):
        #plotNumbers.append(numpy.arange(numRows)*numCols+(c+1))
    #plotNumbers=numpy.array(plotNumbers).flatten()

    fig=pylab.figure(figsize = (4*numCols, 4*numRows))
    plotCount=1
    pylab.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.02, hspace=0.02)
    for obj in catalog:
        #print obj['name'], plotCount
        inJPGPath=imageDir+os.path.sep+obj['name'].replace(" ", "_").replace("\n", "_")+".jpg"
        im=Image.open(inJPGPath)
        data=numpy.array(im)

        # Flips may vary?
        data=numpy.flipud(data)
        
        R=data[:, :, 0]
        G=data[:, :, 1]
        B=data[:, :, 2]
        
        cutLevels=[[R.min(), R.max()], [G.min(), G.max()], [B.min(), B.max()]]
    
        r=astImages.normalise(R, cutLevels[0])
        g=astImages.normalise(G, cutLevels[1])
        b=astImages.normalise(B, cutLevels[2])
        rgb=numpy.array([r.transpose(), g.transpose(), b.transpose()])
        rgb=rgb.transpose()
        
        pylab.subplot(numRows, numCols, plotCount)
        pylab.imshow(rgb, interpolation="bilinear", origin='lower')
        pylab.xticks([], [])
        pylab.yticks([], [])
        plotCount=plotCount+1

    pylab.savefig(outFileName)
    pylab.close()
    
#------------------------------------------------------------------------------------------------------------
def makeSourceBrowser(imageDict, tableKeys, tableFormats, tableLabels, keysToWrite, keyFormats, keyLabels, 
                      cacheDir, browserOptionsDict, catalogToLink = None, parFileToLink = None,
                      diagnosticsDir = None):
    """Makes source browser webpages with cut out images around each source (both ACT, DSS and contour 
    plots), and info on fluxes, best fitting template, and NED matches etc.
    
    """
    
    
    # Will make the below understand ssh etc. at some point
    outDir=browserOptionsDict['outputPath']
    if os.path.exists(outDir) == False:
        os.makedirs(outDir)
                            
    # Where have we stashed various things?
    nedDir=cacheDir+os.path.sep+"NEDResults"
    
    # Make dirs for images
    outputDirs=[outDir, outDir+os.path.sep+"ACTImages", outDir+os.path.sep+"DSSImages",
                outDir+os.path.sep+"DSSContourImages", outDir+os.path.sep+"BCSImages",
                outDir+os.path.sep+"Stripe82Images", outDir+os.path.sep+"SDSSDR7Images", 
                outDir+os.path.sep+"SDSSDR8Images", outDir+os.path.sep+"Stripe82ColorImages", 
                outDir+os.path.sep+"CFHTLSImages", outDir+os.path.sep+"CFHTLSImages", 
                outDir+os.path.sep+"ExternalImages"]
    for o in outputDirs:
        if os.path.exists(o) == False:
            os.makedirs(o)
    
    # Can place cuts on e.g. S/N on optimal catalog
    catalog=catalogTools.selectFromCatalog(imageDict['optimalCatalog'], browserOptionsDict['constraints'])
    
    # Lists of plot stuff to go into the browser pages - entries in these all correspond to each other
    imageDirs=[]        # where .jpg format images stored
    imageLabels=[]      # labels at top of each source page that allow us to select image to view
    imageCaptions=[]    # caption that goes under image shown on the source pages
    
    # ACT plots
    # We need to store clipped out .fits images around each source for making contour plot - that's not updated yet
    mapOptionsDict=browserOptionsDict['ACTPlotOptions']
    makeACTPlots(catalog, imageDict, outDir+os.path.sep+"ACTImages", mapOptionsDict, 
                 colorMap = "jet")  # was spectral  
    if len(mapOptionsDict['bands']) == 2:
        positions=[" (left)", " (right)"]
    elif len(mapOptionsDict['bands']) == 1:
        positions=[""]
    ACTCaption=""
    for band, pos in zip(mapOptionsDict['bands'], positions):
        if ACTCaption != "":
            ACTCaption=ACTCaption+" and "
        ACTCaption=ACTCaption+"%d GHz%s" % (int(band), pos)
    if mapOptionsDict['mapType'] == 'unfiltered':
        ACTCaption=ACTCaption+" unfiltered ACT map, with median value subtracted and smoothed with a Gaussian"
        ACTCaption=ACTCaption+" kernel with &sigma;=%d\".<br> The source position is marked with the white cross." \
                              % (mapOptionsDict['gaussSmoothArcsec'])
    if mapOptionsDict['mapType'] == 'filtered':
        ACTCaption=ACTCaption+" filtered ACT map, smoothed with a Gaussian"
        ACTCaption=ACTCaption+" kernel with &sigma;=%d\".<br> The source position is marked with the white cross." \
                              % (mapOptionsDict['gaussSmoothArcsec'])
    if len(mapOptionsDict['bands']) > 1:
        ACTCaption=ACTCaption.replace("image", "images")
    imageCaptions.append(ACTCaption)
    imageDirs.append(outDir+os.path.sep+"ACTImages")
    imageLabels.append("ACT")
    
    # Optical plots
    # DSS - from Skyview (could get other surveys via this interface in similar way)
    if browserOptionsDict['addDSS'] == True:
        fetchSkyviewImages(catalog, sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                           surveys = ["DSS2IR", "DSS2R", "DSS2B"], \
                           JPEGFolder = cacheDir+os.path.sep+"DSSJPEGImages", \
                           refetch = browserOptionsDict['refetchDSS'])
        makeSkyviewPlots(catalog, outDir+os.path.sep+"DSSImages", sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                         JPEGFolder = cacheDir+os.path.sep+"DSSJPEGImages", \
                         remakePlots = False, nedDir = nedDir)
        imageLabels.append("DSS2")
        imageDirs.append(outDir+os.path.sep+"DSSImages")
        imageCaptions.append("%.1f' x %.1f' false color DSS2 image. The source position is marked with the white cross." \
                             % (mapOptionsDict['plotSizeArcmin'], mapOptionsDict['plotSizeArcmin']))
    
    # SDSS DR7 - color .jpgs
    if 'addSDSSDR7' in browserOptionsDict.keys() and browserOptionsDict['addSDSSDR7'] == True:
        fetchSDSSDR7Images(catalog, sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                           JPEGFolder = cacheDir+os.path.sep+"SDSSDR7JPEGImages", \
                           refetch = browserOptionsDict['refetchSDSSDR7'])
        makeSDSSDR7Plots(catalog, outDir+os.path.sep+"SDSSDR7Images", sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                         JPEGFolder = cacheDir+os.path.sep+"SDSSDR7JPEGImages", \
                         remakePlots = False, nedDir = nedDir)
        imageLabels.append("SDSS DR7")
        imageDirs.append(outDir+os.path.sep+"SDSSDR7Images")
        imageCaptions.append("%.1f' x %.1f' false color (g,r,i) SDSS DR7 image. The source position is marked with the white cross. <br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR10 spectroscopic redshifts." \
                             % (mapOptionsDict['plotSizeArcmin'], mapOptionsDict['plotSizeArcmin']))

    # SDSS DR8 - color .jpgs
    if 'addSDSSDR8' in browserOptionsDict.keys() and browserOptionsDict['addSDSSDR8'] == True:
        fetchSDSSDR8Images(catalog, sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                           JPEGFolder = cacheDir+os.path.sep+"SDSSDR8JPEGImages", \
                           refetch = browserOptionsDict['refetchSDSSDR8'])
        makeSDSSDR7Plots(catalog, outDir+os.path.sep+"SDSSDR8Images", sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                         JPEGFolder = cacheDir+os.path.sep+"SDSSDR8JPEGImages", \
                         remakePlots = False, nedDir = nedDir)
        imageLabels.append("SDSS DR8")
        imageDirs.append(outDir+os.path.sep+"SDSSDR8Images")
        imageCaptions.append("%.1f' x %.1f' false color (g,r,i) SDSS DR8 image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR10 spectroscopic redshifts." \
                             % (mapOptionsDict['plotSizeArcmin'], mapOptionsDict['plotSizeArcmin']))

    # CFHTLS - color .jpgs
    if 'addCFHTLS' in browserOptionsDict.keys() and browserOptionsDict['addCFHTLS'] == True:
        fetchCFHTImages(catalog, sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                        JPEGFolder = cacheDir+os.path.sep+"CFHTLSJPEGImages", \
                        refetch = browserOptionsDict['refetchCFHTLS'])
        makeCFHTPlots(catalog, outDir+os.path.sep+"CFHTLSImages", sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                      JPEGFolder = cacheDir+os.path.sep+"CFHTLSJPEGImages", \
                      remakePlots = False, nedDir = nedDir)
        imageLabels.append("CFHTLS")
        imageDirs.append(outDir+os.path.sep+"CFHTLSImages")
        imageCaptions.append("%.1f' x %.1f' false color (g,r,i) CFHTLS image centred on the candidate position." \
                             % (mapOptionsDict['plotSizeArcmin'], mapOptionsDict['plotSizeArcmin']))
        
    # BCS - on goliath, painfully slow
    #if opticalOptionsDict['addBCS'] == True:
        #makeBCSPlots(catalog, outDir+os.path.sep+"BCSImages", opticalOptionsDict, \
                     #plotSizeArcmin = mapOptionsDict['plotSizeArcmin'], nedDir = nedDir)
    
    # SDSS Stripe 82 images - grab .fits (and cache), make color .jpgs
    if browserOptionsDict['addStripe82Color'] == True:
        fetchStripe82ColorImages(catalog, sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                                 jpegFolder = cacheDir+os.path.sep+"Stripe82ColorJPEGImages", \
                                 refetch = browserOptionsDict['refetchStripe82Color'])
        makeSDSSDR7Plots(catalog, outDir+os.path.sep+"Stripe82ColorImages", \
                      JPEGFolder = cacheDir+os.path.sep+"Stripe82ColorJPEGImages", \
                      sizeArcmin = mapOptionsDict['plotSizeArcmin'], remakePlots = False, \
                      nedDir = nedDir, plotNEDObjects = False, plotSDSSObjects = False)
        imageLabels.append("SDSS Stripe 82")
        imageDirs.append(outDir+os.path.sep+"Stripe82ColorImages")
        imageCaptions.append("%.1f' x %.1f' (g, r, i) color composite SDSS Stripe 82 image. The source position is marked with the white cross." \
                             % (mapOptionsDict['plotSizeArcmin'], mapOptionsDict['plotSizeArcmin']))

    # SDSS Stripe 82 - mono only .jpgs
    #if browserOptionsDict['addStripe82'] == True:
        #fetchStripe82Images(catalog, sizeArcmin = mapOptionsDict['plotSizeArcmin'], \
                        #jpegFolder = cacheDir+os.path.sep+"Stripe82JPEGImages", \
                        #refetch = browserOptionsDict['refetchStripe82'])
        #makeStripe82Plots(catalog, outDir+os.path.sep+"Stripe82Images", opticalOptionsDict, \
                      #plotSizeArcmin = mapOptionsDict['plotSizeArcmin'], remakePlots = False, \
                      #nedDir = nedDir)
        #imageLabels.append("SDSS Stripe 82")
        #imageDirs.append(outDir+os.path.sep+"Stripe82Images")
        #imageCaptions.append("%.1f' x %.1f' r-band SDSS Stripe 82 image. The source position is marked with the white cross." \
                             #% (mapOptionsDict['plotSizeArcmin'], mapOptionsDict['plotSizeArcmin']))
                            
    # External directory of pre-made images, e.g. from Felipe's paper
    if 'addExternalImages' in browserOptionsDict.keys() and browserOptionsDict['addExternalImages'] == True:
        os.system("cp %s/*.jpg %s" % (browserOptionsDict['externalImagesDir'], outDir+os.path.sep+"ExternalImages"))
        imageLabels.append(browserOptionsDict['externalImagesLabel'])
        imageDirs.append(outDir+os.path.sep+"ExternalImages")
        imageCaptions.append(browserOptionsDict['externalImagesCaption'])
    
    
    # Compress FITS images if we saved them
    if 'saveFITS' in mapOptionsDict.keys() and mapOptionsDict['saveFITS'] == True:
        outFITSDir=outDir+os.path.sep+"FITSImages"
        os.system("tar -cf %s.tar %s" % (outFITSDir, outFITSDir))
        os.system("gzip -f %s.tar" % (outFITSDir))
        os.system("rm -r %s" % (outFITSDir))
        linkToFITSTarFile=True
    else:
        linkToFITSTarFile=False
    
    # Calculate map area, if we can find an area mask
    if diagnosticsDir != None and os.path.exists(diagnosticsDir+os.path.sep+"areaMask.fits") == True:
        img=pyfits.open(diagnosticsDir+os.path.sep+"areaMask.fits")
        mapData=img[0].data
        wcs=astWCS.WCS(diagnosticsDir+os.path.sep+"areaMask.fits")
        mapAreaDeg2=mapTools.getPixelAreaArcmin2Map(mapData, wcs)/3600.0
        surveyAreaDeg2=mapAreaDeg2[numpy.equal(mapData, 1)].sum()
    else:
        surveyAreaDeg2=None
    
    # Make webpages for each source
    makeSourcePages(catalog, outDir, keysToWrite, keyFormats, keyLabels, imageDirs, imageLabels, imageCaptions,
                    nedDir = nedDir, sizeArcmin = mapOptionsDict['plotSizeArcmin'])    
    
    # Make table page (acts as index) - add links on this to the catalog in .csv, .reg, .fits
    if 'imageMapOptions' not in browserOptionsDict.keys():
        browserOptionsDict['imageMapOptions']=None
    else:
        browserOptionsDict['imageMapOptions']['imageDict']=imageDict
    makeTablePage(catalog, outDir, imageDirs[0], tableKeys, tableFormats, tableLabels, 
                  catalogToLink = catalogToLink, parFileToLink = parFileToLink, 
                  surveyAreaDeg2 = surveyAreaDeg2,
                  objectTypeString = browserOptionsDict['objectTypeString'], 
                  commentsString = browserOptionsDict['comments'], 
                  addLinkToDiagnosticsPage = browserOptionsDict['makeDiagnosticsPage'], 
                  addLinkToFITSTarFile = linkToFITSTarFile,
                  imageMapOptions = browserOptionsDict['imageMapOptions'])
    
    # Make diagnostics page if asked for
    if browserOptionsDict['makeDiagnosticsPage'] == True:
        makeDiagnosticsPage(outDir, imageDict, diagnosticsDir)
