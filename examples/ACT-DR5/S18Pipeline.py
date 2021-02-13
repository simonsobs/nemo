#!/usr/bin/env python

import subprocess, os

# Point sources run (two passes)
process=subprocess.run(['sbatch', 'PS_S18_f150_auto.sh'], stdout=subprocess.PIPE)
if (process.returncode == 0):
    jobnum=process.stdout.split(b"job ")[-1].rstrip().decode()
    print("   jobnum = %s" % (jobnum))
else:
    raise Exception("Error submitting point sources run.")

# Clusters run if point sources completes
process=subprocess.run(['sbatch', '--depend=afterok:%s' % (jobnum), 'MFMF_S18_auto.sh'], 
                       stdout=subprocess.PIPE)
if (process.returncode == 0):
    jobnum=process.stdout.split(b"job ")[-1].rstrip().decode()
    print("   jobnum = %s" % (jobnum))
else:
    raise Exception("Error submitting clusters run.")
