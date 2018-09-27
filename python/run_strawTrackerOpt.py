__author__ = 'Mikhail Hushchyn'

import os, sys, shutil
import multiprocessing
import numpy as np
import pandas as pd

output_dir = "output"

if os.path.exists(output_dir):
    print('The directiry '+output_dir+' is already exists. Remove it.')
    sys.exit(2)
else:
    os.mkdir(output_dir)


# Number of random point to generate.
n_points = 1000

# Define random values on the tracking system parameters
min_dist = 3.6

StrawPitch_grid    = np.random.RandomState(12).uniform(min_dist, min_dist, n_points)
YLayerOffset_grid  = np.random.RandomState(13).uniform(min_dist/2., min_dist, n_points)
YPlaneOffset_grid  = np.random.RandomState(14).uniform(min_dist*0.25, min_dist*1.25, n_points)
DeltazLayer_grid   = np.random.RandomState(15).uniform(1, 12, n_points)
DeltazPlane_grid   = np.random.RandomState(16).uniform(1, 12, n_points)
DeltazView_grid    = np.random.RandomState(17).uniform(10, 12, n_points)
ViewAngle_grid     = np.random.RandomState(18).uniform(5, 15, n_points)

points = []

# Baseline point
apoint = [0, 3.60, 1.9, 1.3, 1.6, 4.2, 10., 5]
points.append(apoint)

for i in range(n_points):

    StrawPitch    = StrawPitch_grid[i]
    YLayerOffset  = YLayerOffset_grid[i]
    YPlaneOffset  = YPlaneOffset_grid[i]
    DeltazLayer   = DeltazLayer_grid[i]
    DeltazPlane   = DeltazPlane_grid[i]
    DeltazView    = DeltazView_grid[i]
    ViewAngle     = int(ViewAngle_grid[i])

    apoint = [i+1, StrawPitch, YLayerOffset, YPlaneOffset, DeltazLayer, DeltazPlane, DeltazView, ViewAngle]
    points.append(apoint)



def run_opt(point):

    [i, StrawPitch, YLayerOffset, YPlaneOffset, DeltazLayer, DeltazPlane, DeltazView, ViewAngle] = point

    print ("Run for point: ", i)

    cols = ['N', 'StrawPitch', 'YLayerOffset', 'YPlaneOffset', 'DeltazLayer', 'DeltazPlane', 'DeltazView', 'ViewAngle']
    df = pd.DataFrame(data=[point], columns=cols)

    i_dir = output_dir+"/"+str(i)

    if os.path.exists(i_dir):
        shutil.rmtree(i_dir)
        os.mkdir(i_dir)
    else:
        os.mkdir(i_dir)

    df.to_csv(i_dir + '/point.csv')

    cmd = 'sudo docker run -i -t --rm'
    cmd += ' -v $(pwd)/' + i_dir + ':/working_dir'
    cmd += ' ship_tracker_opt /bin/bash -l -c'
    cmd += ' "source /SHiPBuild/config.sh && python /SHiPBuild/FairShip/python/strawTrackerOpt.py'
    cmd += ' --StrawPitch=' + str(StrawPitch)
    cmd += ' --YLayerOffset=' + str(YLayerOffset)
    cmd += ' --YPlaneOffset=' + str(YPlaneOffset)
    cmd += ' --DeltazLayer=' + str(DeltazLayer)
    cmd += ' --DeltazPlane=' + str(DeltazPlane)
    cmd += ' --DeltazView=' + str(DeltazView)
    cmd += ' --ViewAngle=' + str(ViewAngle)
    cmd += ' --nEvents=10000 --MCFile=/SHiPBuild/Cascade-parp16-MSTP82-1-MSEL4-978Bpot.root >> log.txt"'
    os.system(cmd)

    os.system("rm -f " + i_dir + "/ship.conical.Pythia8-TGeant4.root")
    os.system("rm -f " + i_dir + "/ship.conical.Pythia8-TGeant4_rec.root")

pool = multiprocessing.Pool(20)
pool.map(run_opt, points)
