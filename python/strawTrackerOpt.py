__author__ = 'Mikhail Hushchyn'

import sys, os
import getopt
import json
import ROOT


def objective_function(StrawPitch = 3.60, YLayerOffset = 1.9, YPlaneOffset = 1.3, DeltazLayer = 1.6,
                       DeltazPlane = 4.2, DeltazView = 10., ViewAngle = int(5), nEvents = int(100), MCFile=""):

    FAIRSHIP = os.environ.get('FAIRSHIP')

    # Change geometry config
    with open(FAIRSHIP+'/geometry/geometry_config.py') as f:
        lines = f.read().splitlines()

    lines[152-1] = "     c.strawtubes.StrawPitch         = " + str(StrawPitch) + "*u.cm  "
    lines[153-1] = "     c.strawtubes.DeltazLayer        = " + str(DeltazLayer) + "*u.cm  "
    lines[154-1] = "     c.strawtubes.DeltazPlane        = " + str(DeltazPlane) + "*u.cm  "
    lines[155-1] = "     c.strawtubes.YLayerOffset       = " + str(YLayerOffset) + "*u.cm  "
    lines[156-1] = "     c.strawtubes.YPlaneOffset       = " + str(YPlaneOffset) + "*u.cm  "
    lines[165-1] = "    c.strawtubes.ViewAngle          = " + str(ViewAngle) + "  "
    lines[167-1] = "    c.strawtubes.DeltazView         = " + str(DeltazView) + "*u.cm  "

    with open(FAIRSHIP+'/geometry/geometry_config.py','w') as f:
        f.write('\n'.join(lines))


    # Run MC generation
    cmd = "python $FAIRSHIP/macro/run_simScript.py --caloDesign=3 --tankDesign=6 --muShieldDesign=9 --nuTauTargetDesign=3"
    cmd += " --nEvents "+str(nEvents)
    if MCFile != "":
        cmd += " -f "+str(MCFile)
    os.system(cmd)


    # Run PatRec
    os.system("python $FAIRSHIP/macro/ShipReco.py -f ship.conical.Pythia8-TGeant4.root -g geofile_full.conical.Pythia8-TGeant4.root --realPR=FH")


    # Read quality metrics
    recohists = ROOT.TFile("recohists.root")
    ahist = recohists.Get("Reco_T14")
    metric = ahist.GetBinContent(2)

    return metric



if __name__=='__main__':

    argv = sys.argv[1:]

    #default values for parameters
    StrawPitch = 3.60
    YLayerOffset = 1.9
    YPlaneOffset = 1.3
    DeltazLayer = 1.6
    DeltazPlane = 4.2
    DeltazView = 10.
    ViewAngle = int(5)
    nEvents = 100
    MCFile = ""

    output_file = "output.txt"

    try:
        opts, args = getopt.getopt(argv, "", ["StrawPitch=", "YLayerOffset=", "YPlaneOffset=", "DeltazLayer=", "DeltazPlane=", "DeltazView=", "ViewAngle=", "output=", "nEvents=", "MCFile="])
    except getopt.GetoptError:
        print("Wrong parameters. Available params: StrawPitch, YLayerOffset, YPlaneOffset, DeltazLayer, DeltazPlane, DeltazView, ViewAngle, output, nEvents, MCFile.\n")
        sys.exit(2)

    for opt, arg in opts:
        if opt == "--StrawPitch":
            StrawPitch = arg
        elif opt == "--YLayerOffset":
            YLayerOffset = arg
        elif opt == "--YPlaneOffset":
            YPlaneOffset = arg
        elif opt == "--DeltazLayer":
            DeltazLayer = arg
        elif opt == "--DeltazPlane":
            DeltazPlane = arg
        elif opt == "--DeltazView":
            DeltazView = arg
        elif opt == "--ViewAngle":
            ViewAngle = arg
        elif opt == "--output":
            output_file = arg
        elif opt == "--nEvents":
            nEvents = arg
        elif opt == "--MCFile":
            MCFile = arg


    metric = objective_function(StrawPitch, YLayerOffset, YPlaneOffset, DeltazLayer,
                                DeltazPlane, DeltazView, ViewAngle, nEvents, MCFile)

    with open(output_file, 'w') as tf:
        json.dump(metric, tf)

    print "Metric: ", metric

