"""

Generate .yml configs for various maps and beams.

It's too complicated to understand MFMF maps (because so many pa and set combos)

So, let's start with single frequency only

This runs quick enough that we can stick everything in one batch job

"""

import os
import sys
import numpy as np
import IPython

# Map calibration factors (Steve, via Sigurd ACTPol Slack) --------------------------------------------------
# NOTE: These are just map * factor (pre-S17 needed sqrt, which has already been applied here)
mapCalibFactors={"s13_pa1_f150_deep1":  np.sqrt(0.974),
                 "s13_pa1_f150_deep5":  np.sqrt(1.113),
                 "s13_pa1_f150_deep6":  np.sqrt(0.942),
                 "s14_pa1_f150_deep56": np.sqrt(0.997),
                 "s14_pa2_f150_deep56": np.sqrt(0.978),
                 "s15_pa1_f150_deep56": np.sqrt(0.987),
                 "s15_pa1_f150_boss":   np.sqrt(1.070),
                 "s15_pa1_f150_deep8":  np.sqrt(1.023),
                 "s15_pa2_f150_deep56": np.sqrt(0.958),
                 "s15_pa2_f150_boss":   np.sqrt(1.063),
                 "s15_pa2_f150_deep8":  np.sqrt(1.099),
                 "s15_pa3_f090_deep56": np.sqrt(1.059),
                 "s15_pa3_f090_boss":   np.sqrt(1.125),
                 "s15_pa3_f090_deep8":  np.sqrt(1.110),
                 "s15_pa3_f150_deep56": np.sqrt(1.093),
                 "s15_pa3_f150_boss":   np.sqrt(1.219),
                 "s15_pa3_f150_deep8":  np.sqrt(1.221),
                 "s16_pa2_f150_cmb":    np.sqrt(0.951),
                 "s16_pa3_f090_cmb":    np.sqrt(1.136),
                 "s16_pa3_f150_cmb":    np.sqrt(1.960),
                 "s17_pa4_f150":        1.20995868,
                 "s17_pa4_f220":        0.94868330,
                 "s17_pa5_f090":        0.98843310,
                 "s17_pa5_f150":        1.03778611,
                 "s17_pa6_f090":        1.00000000,
                 "s17_pa6_f150":        1.02176318,
                 "s18_pa4_f150":        1.200,
                 "s18_pa4_f220":        1.081,
                 "s18_pa5_f090":        1.174,
                 "s18_pa5_f150":        1.049,
                 "s18_pa6_f090":        1.214,
                 "s18_pa6_f150":        1.186}

# Configs: maps/beams/labels --------------------------------------------------------------------------------
configs=[{'label':  "s14_deep56",
          'map':    "maps/s14_deep56_pa$PA_$BAND_nohwp_night_3pass_4way_set$SET_map.fits",
          'beam':   "beams/190809/b20190809_s14_pa$PA_$BAND_nohwp_night_beam_profile_jitter_deep56.txt",
          'calKey': "s14_pa$PA_$BAND_deep56",
          'surveyMask': "null"},
         {'label':  "s15_deep56",
          'map':    "maps/s15_deep56_pa$PA_$BAND_nohwp_night_3pass_4way_set$SET_map.fits",
          'beam':   "beams/190809/b20190809_s15_pa$PA_$BAND_nohwp_night_beam_profile_jitter_deep56.txt",
          'calKey': "s15_pa$PA_$BAND_deep56",
          'surveyMask': "null"},
         {'label':  "s15_boss",
          'map':    "maps/s15_boss_pa$PA_$BAND_nohwp_night_3pass_4way_set$SET_map.fits",
          'beam':   "beams/190809/b20190809_s15_pa$PA_$BAND_nohwp_night_beam_profile_jitter_boss.txt",
          'calKey': "s15_pa$PA_$BAND_boss",
          'surveyMask': "null"},
         {'label':  "s16_cmb",
          'map':    "maps/s16_cmb_pa$PA_$BAND_nohwp_night_3pass_2way_set$SET_map.fits",
          'beam':   "beams/190809/b20190809_s16_pa$PA_$BAND_nohwp_night_beam_profile_jitter_cmb.txt",
          'calKey': "s16_pa$PA_$BAND_cmb",
          'surveyMask': "surveyMask_v7_S18.fits"},
         {'label':  "s17_cmb",
          'map':    "maps/s17_cmb_pa$PA_$BAND_nohwp_night_1pass_2way_set$SET_map.fits",
          'beam':   "beams/s17_pa$PA_$BAND_nohwp_night_beam_profile_jitter.txt",
          'calKey': "s17_pa$PA_$BAND",
          'surveyMask': "surveyMask_v7_S18.fits"},
         {'label':  "s18_cmb",
          'map':    "maps/s18_cmb_pa$PA_$BAND_nohwp_night_1pass_2way_set$SET_map.fits",
          'beam':   "beams/s17_pa$PA_$BAND_nohwp_night_beam_profile_jitter.txt",
          'calKey': "s18_pa$PA_$BAND",
          'surveyMask': "surveyMask_v7_S18.fits"},]

# Make Nemo .yml configs and slurm batch script -------------------------------------------------------------         
batchLines=[]
for c in configs:
    for band, freq in zip(['f150', 'f090'], [145.3, 94.1]):
        for paNum in range(1, 7):
            for setNum in range(0, 2):
                with open("template_singlefreq.yml") as inFile:
                    lines=inFile.readlines()
                mapName=c['map'].replace("$PA", str(paNum)).replace("$SET", str(setNum)).replace("$BAND", band)
                #weightName=mapName.replace("_map.fits", "_div.fits")
                beamName=c['beam'].replace("$PA", str(paNum)).replace("$BAND", band)
                if os.path.exists(mapName) == True and os.path.exists(beamName) == True:
                    calibFactor=mapCalibFactors[c['calKey'].replace("$PA", str(paNum)).replace("$BAND", band)]
                    surveyMask=c['surveyMask']
                    outFileName=c['label']+"_%s_pa%d_set%d.yml" % (band, paNum, setNum)
                    with open(outFileName, "w") as outFile:
                        for line in lines:
                            if line.find("$MAP") != -1:
                                outFile.write(line.replace("$MAP", mapName))
                            elif line.find("$FREQ") != -1:
                                outFile.write(line.replace("$FREQ", str(freq)))
                            elif line.find("$BEAM") != -1:
                                outFile.write(line.replace("$BEAM", beamName))
                            elif line.find("$CALIBFACTOR") != -1:
                                outFile.write(line.replace("$CALIBFACTOR", str(calibFactor)))
                            elif line.find("$SURVEYMASK") != -1:
                                outFile.write(line.replace("$SURVEYMASK", surveyMask))
                            else:
                                outFile.write(line)
                    batchLines.append("time mpiexec nemo %s -M" % (outFileName))

with open("clustercal.sh", "w") as batchFile:
    batch="""#!/bin/sh
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64000
#SBATCH --time=23:59:00

source ~/.bashrc
"""
    batchFile.write(batch)
    for line in batchLines:
        batchFile.write(line+"\n")
