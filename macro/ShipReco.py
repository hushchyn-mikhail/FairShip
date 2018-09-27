#!/usr/bin/env python
inputFile = 'ship.conical.Pythia8-TGeant4.root'
geoFile   = None
debug = False
EcalDebugDraw = False
withNoStrawSmearing = None # True   for debugging purposes
nEvents    = 999999
firstEvent = 0
withHists = True
vertexing = True
dy  = None
saveDisk  = False # remove input file
pidProton = False # if true, take truth, if False fake with pion mass
realPR = ''
realPROptions=["Prev", "FH", "AR", "Baseline"]
withT0 = False

import resource
def mem_monitor():
 # Getting virtual memory size 
    pid = os.getpid()
    with open(os.path.join("/proc", str(pid), "status")) as f:
        lines = f.readlines()
    _vmsize = [l for l in lines if l.startswith("VmSize")][0]
    vmsize = int(_vmsize.split()[1])
    #Getting physical memory size  
    pmsize = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print "memory: virtuell = %5.2F MB  physical = %5.2F MB"%(vmsize/1.0E3,pmsize/1.0E3)

import ROOT,os,sys,getopt
import __builtin__ as builtin
import rootUtils as ut
import shipunit as u
import shipRoot_conf

shipRoot_conf.configure()

try:
        opts, args = getopt.getopt(sys.argv[1:], "o:D:FHPu:n:f:g:c:hqv:sl:A:Y:i:",\
           ["ecalDebugDraw","inputFile=","geoFile=","nEvents=","noStrawSmearing","noVertexing","saveDisk","realPR=","withT0"])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter --inputFile=  --geoFile= --nEvents=  --firstEvent=,'
        print ' noStrawSmearing: no smearing of distance to wire, default on'
        print ' outputfile will have same name with _rec added'  
        print ' --realPR= defines track pattern recognition. Possible options: ',realPROptions, "if no option given, fake PR is used."
        sys.exit()
for o, a in opts:
        if o in ("noVertexing",):
            vertexing = False
        if o in ("noStrawSmearing",):
            withNoStrawSmearing = True
        if o in ("--withT0",):
            withT0 = True
        if o in ("-f", "--inputFile",):
            inputFile = a
        if o in ("-g", "--geoFile",):
            geoFile = a
        if o in ("-n", "--nEvents=",):
            nEvents = int(a)
        if o in ("-Y",):
            dy = float(a)
        if o in ("--ecalDebugDraw",):
            EcalDebugDraw = True
        if o in ("--saveDisk",):
            saveDisk = True
        if o in ("--realPR",):
            realPR = a
            if not realPR in realPROptions:
              print "wrong option given for realPR,",a," should be one of ",realPROptions
              exit(1)
if EcalDebugDraw: ROOT.gSystem.Load("libASImage")

# need to figure out which geometry was used
if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')
  try:
    dy = float( tmp[1]+'.'+tmp[2] )
  except:
    dy = None
realPRoption = realPR
if realPR=='':realPRoption='No' 
print 'configured to process ',nEvents,' events from ' ,inputFile, \
      ' starting with event ',firstEvent, ' with option Yheight = ',dy,' with vertexing',vertexing,' and real pattern reco ',realPRoption
if not inputFile.find('_rec.root') < 0: 
  outFile   = inputFile
  inputFile = outFile.replace('_rec.root','.root') 
else:
  outFile = inputFile.replace('.root','_rec.root') 
# outfile should be in local directory
  tmp = outFile.split('/')
  outFile = tmp[len(tmp)-1]
  if inputFile[:7]=="root://" : os.system('xrdcp '+inputFile+' '+outFile)
  elif saveDisk: os.system('mv '+inputFile+' '+outFile)
  else :       os.system('cp '+inputFile+' '+outFile)

if not geoFile:
 tmp = inputFile.replace('ship.','geofile_full.')
 geoFile = tmp.replace('_rec','')

fgeo = ROOT.TFile.Open(geoFile)

from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
#load Shipgeo dictionary
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')
ecalGeoFile = ShipGeo.ecal.File

h={}
log={}
if withHists:
 ut.bookHist(h,'distu','distance to wire',100,0.,5.)
 ut.bookHist(h,'distv','distance to wire',100,0.,5.)
 ut.bookHist(h,'disty','distance to wire',100,0.,5.)
 ut.bookHist(h,'nmeas','nr measuerements',100,0.,50.)
 ut.bookHist(h,'chi2','Chi2/DOF',100,0.,20.)

 ut.bookHist(h,'hits_true_y12','Number of hits per MC track. Stations 1&2, Y views', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_stereo12','Number of hits per MC track. Stations 1&2, Stereo views', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_12','Number of hits per MC track. Stations 1&2', 33, -0.5, 32.5)
 ut.bookHist(h,'hits_true_y34','Number of hits per MC track. Stations 3&4, Y views', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_stereo34','Number of hits per MC track. Stations 3&3, Stereo views', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_34','Number of hits per MC track. Stations 3&4', 33, -0.5, 32.5)
 ut.bookHist(h,'hits_true_1234','Number of hits per MC track. Stations 1-4', 65, -0.5, 64.5)
 ut.bookHist(h,'hits_true_1','Number of hits per MC track. Station 1', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_2','Number of hits per MC track. Station 2', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_3','Number of hits per MC track. Station 3', 17, -0.5, 16.5)
 ut.bookHist(h,'hits_true_4','Number of hits per MC track. Station 4', 17, -0.5, 16.5)

 ut.bookHist(h,'pdg','Number of particles', 6, -0.5, 6.5)
 h['pdg'].GetXaxis().SetBinLabel(1,"Electrons")
 h['pdg'].GetXaxis().SetBinLabel(2,"Muons")
 h['pdg'].GetXaxis().SetBinLabel(3,"Pions")
 h['pdg'].GetXaxis().SetBinLabel(4,"Kaons")
 h['pdg'].GetXaxis().SetBinLabel(5,"Protons")
 h['pdg'].GetXaxis().SetBinLabel(6,"Others")

 ut.bookHist(h,'p/pt_truth','P vs Pt (GeV), MC Truth',100,0.,200.,100,0.,5.)
 ut.bookHist(h,'p_truth','P (GeV), MC Truth',100,0.,200.)
 ut.bookHist(h,'pt_truth','Pt (GeV), MC Truth',100,0.,5.)
 ut.bookHist(h,'pz_truth','Pz (GeV), MC Truth',100,0.,200.)
 ut.bookHist(h,'px_truth','Px (GeV), MC Truth',100,0.,5.)
 ut.bookHist(h,'py_truth','Py (GeV), MC Truth',100,0.,5.)

 ut.bookHist(h,'p/pt_reco','P vs Pt (GeV), Reconstructed',100,0.,200.,100,0.,5.)
 ut.bookHist(h,'p_reco','P (GeV), Reconstructed',100,0.,200.)
 ut.bookHist(h,'pt_reco','Pt (GeV), Reconstructed',100,0.,5.)
 ut.bookHist(h,'pz_reco','Pz (GeV), Reconstructed',100,0.,200.)
 ut.bookHist(h,'px_reco','Px (GeV), Reconstructed',100,0.,5.)
 ut.bookHist(h,'py_reco','Py (GeV), Reconstructed',100,0.,5.)

 ut.bookHist(h,'p_rel_error','(P_reco - P_true) / P_true',100,-0.5,0.5)
 ut.bookHist(h,'pt_rel_error','(Pt_reco - Pt_true) / Pt_true',100,-0.5,0.5)
 ut.bookHist(h,'px_rel_error','(Px_reco - Px_true) / Px_true',100,-0.5,0.5)
 ut.bookHist(h,'py_rel_error','(Py_reco - Py_true) / Py_true',100,-0.5,0.5)
 ut.bookHist(h,'pz_rel_error','(Pz_reco - Pz_true) / Pz_true',100,-0.5,0.5)

 ut.bookHist(h,'tracks_per_event_true','Number of tracks per MC event.', 11, -0.5, 10.5)


 ut.bookHist(h,'Reco_y12','Number of recognized tracks, clones and ghosts in y12 views.', 5, -0.5, 4.5)
 h['Reco_y12'].GetXaxis().SetBinLabel(1,"N total")
 h['Reco_y12'].GetXaxis().SetBinLabel(2,"N recognized tracks")
 h['Reco_y12'].GetXaxis().SetBinLabel(3,"N clones")
 h['Reco_y12'].GetXaxis().SetBinLabel(4,"N ghosts")
 h['Reco_y12'].GetXaxis().SetBinLabel(5,"N others")

 ut.bookHist(h,'Reco_stereo12','Number of recognized tracks, clones and ghosts in stereo12 views.', 5, -0.5, 4.5)
 h['Reco_stereo12'].GetXaxis().SetBinLabel(1,"N total")
 h['Reco_stereo12'].GetXaxis().SetBinLabel(2,"N recognized tracks")
 h['Reco_stereo12'].GetXaxis().SetBinLabel(3,"N clones")
 h['Reco_stereo12'].GetXaxis().SetBinLabel(4,"N ghosts")
 h['Reco_stereo12'].GetXaxis().SetBinLabel(5,"N others")

 ut.bookHist(h,'Reco_y34','Number of recognized tracks, clones and ghosts in y34 views.', 5, -0.5, 4.5)
 h['Reco_y34'].GetXaxis().SetBinLabel(1,"N total")
 h['Reco_y34'].GetXaxis().SetBinLabel(2,"N recognized tracks")
 h['Reco_y34'].GetXaxis().SetBinLabel(3,"N clones")
 h['Reco_y34'].GetXaxis().SetBinLabel(4,"N ghosts")
 h['Reco_y34'].GetXaxis().SetBinLabel(5,"N others")

 ut.bookHist(h,'Reco_stereo34','Number of recognized tracks, clones and ghosts in stereo34 views.', 5, -0.5, 4.5)
 h['Reco_stereo34'].GetXaxis().SetBinLabel(1,"N total")
 h['Reco_stereo34'].GetXaxis().SetBinLabel(2,"N recognized tracks")
 h['Reco_stereo34'].GetXaxis().SetBinLabel(3,"N clones")
 h['Reco_stereo34'].GetXaxis().SetBinLabel(4,"N ghosts")
 h['Reco_stereo34'].GetXaxis().SetBinLabel(5,"N others")

 ut.bookHist(h,'Reco_T14','Number of recognized tracks, clones and ghosts for tracks with at least one hit in T1-4.', 5, -0.5, 4.5)
 h['Reco_T14'].GetXaxis().SetBinLabel(1,"N total")
 h['Reco_T14'].GetXaxis().SetBinLabel(2,"N recognized tracks")
 h['Reco_T14'].GetXaxis().SetBinLabel(3,"N clones")
 h['Reco_T14'].GetXaxis().SetBinLabel(4,"N ghosts")
 h['Reco_T14'].GetXaxis().SetBinLabel(5,"N others")

 ut.bookHist(h,'Reco_3hits_per_views','Number of recognized tracks, clones and ghosts for tracks with at least 3 hits in y12, y34, stereo12 and stereo34.', 5, -0.5, 4.5)
 h['Reco_3hits_per_views'].GetXaxis().SetBinLabel(1,"N total")
 h['Reco_3hits_per_views'].GetXaxis().SetBinLabel(2,"N recognized tracks")
 h['Reco_3hits_per_views'].GetXaxis().SetBinLabel(3,"N clones")
 h['Reco_3hits_per_views'].GetXaxis().SetBinLabel(4,"N ghosts")
 h['Reco_3hits_per_views'].GetXaxis().SetBinLabel(5,"N others")

import shipDet_conf
run = ROOT.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetOutputFile("dummy")  # Output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for creating VMC field
rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
modules = shipDet_conf.configure(run,ShipGeo)
run.Init()
import geomGeant4

if hasattr(ShipGeo.Bfield,"fieldMap"):
  fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True)

# make global variables
builtin.debug    = debug
builtin.pidProton = pidProton
builtin.withT0 = withT0
builtin.realPR = realPR
builtin.vertexing = vertexing
builtin.ecalGeoFile = ecalGeoFile
builtin.ShipGeo = ShipGeo
builtin.modules = modules
builtin.EcalDebugDraw  = EcalDebugDraw
builtin.withNoStrawSmearing = withNoStrawSmearing
builtin.h    = h
builtin.log  = log
iEvent = 0
builtin.iEvent  = iEvent

# import reco tasks
import shipDigiReco
geoMat =  ROOT.genfit.TGeoMaterialInterface()  # if only called in ShipDigiReco -> crash, reason unknown

SHiP = shipDigiReco.ShipDigiReco(outFile,fgeo)
nEvents   = min(SHiP.sTree.GetEntries(),nEvents)
# main loop
for iEvent in range(firstEvent, nEvents):
 if iEvent%100 == 0 or debug: print 'event ',iEvent
 rc    = SHiP.sTree.GetEvent(iEvent) 
 SHiP.digitize()
 SHiP.reconstruct()
 # memory monitoring
 # mem_monitor() 
# end loop over events
SHiP.finish()
