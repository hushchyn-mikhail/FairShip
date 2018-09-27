import os,ROOT,shipVertex,shipDet_conf
if realPR == "Prev": import shipPatRec_prev as shipPatRec # The previous version of the pattern recognition
else: import shipPatRec
import shipunit as u
import rootUtils as ut
from array import array
import sys, math
import TrackExtrapolateTool
stop  = ROOT.TVector3()
start = ROOT.TVector3()

class ShipDigiReco:
    " convert FairSHiP MC hits / digitized hits to measurements"
    def __init__(self,fout,fgeo):
        self.fn = ROOT.TFile.Open(fout,'update')
        self.sTree     = self.fn.cbmsim
        if self.sTree.GetBranch("FitTracks"):
            print "remove RECO branches and rerun reconstruction"
            self.fn.Close()
            # make a new file without reco branches
            f = ROOT.TFile(fout)
            sTree = f.cbmsim
            if sTree.GetBranch("FitTracks"): sTree.SetBranchStatus("FitTracks",0)
            if sTree.GetBranch("goodTracks"): sTree.SetBranchStatus("goodTracks",0)
            if sTree.GetBranch("VetoHitOnTrack"): sTree.SetBranchStatus("VetoHitOnTrack",0)
            if sTree.GetBranch("Particles"): sTree.SetBranchStatus("Particles",0)
            if sTree.GetBranch("fitTrack2MC"): sTree.SetBranchStatus("fitTrack2MC",0)
            if sTree.GetBranch("EcalClusters"): sTree.SetBranchStatus("EcalClusters",0)
            if sTree.GetBranch("EcalReconstructed"): sTree.SetBranchStatus("EcalReconstructed",0)
            if sTree.GetBranch("Pid"): sTree.SetBranchStatus("Pid",0)
            if sTree.GetBranch("Digi_StrawtubesHits"): sTree.SetBranchStatus("Digi_StrawtubesHits",0)
            if sTree.GetBranch("Digi_SBTHits"): sTree.SetBranchStatus("Digi_SBTHits",0)
            if sTree.GetBranch("digiSBT2MC"):   sTree.SetBranchStatus("digiSBT2MC",0)
            if sTree.GetBranch("Digi_TimeDetHits"): sTree.SetBranchStatus("Digi_TimeDetHits",0)
            if sTree.GetBranch("Digi_MuonHits"): sTree.SetBranchStatus("Digi_MuonHits",0)

            rawFile = fout.replace("_rec.root","_raw.root")
            recf = ROOT.TFile(rawFile,"recreate")
            newTree = sTree.CloneTree(0)
            for n in range(sTree.GetEntries()):
                sTree.GetEntry(n)
                rc = newTree.Fill()
            sTree.Clear()
            newTree.AutoSave()
            f.Close()
            recf.Close()
            os.system('cp '+rawFile +' '+fout)
            self.fn = ROOT.TFile(fout,'update')
            self.sTree     = self.fn.cbmsim
        #  check that all containers are present, otherwise create dummy version
        self.dummyContainers={}
        branch_class = {"vetoPoint":"vetoPoint","ShipRpcPoint":"ShipRpcPoint","TargetPoint":"TargetPoint", \
                        "strawtubesPoint":"strawtubesPoint","EcalPointLite":"ecalPoint","HcalPointLite":"hcalPoint", \
                        "splitcalPoint":"splitcalPoint","TimeDetPoint":"TimeDetPoint","muonPoint":"muonPoint"}
        for x in branch_class:
            if not self.sTree.GetBranch(x):
                self.dummyContainers[x+"_array"] = ROOT.TClonesArray(branch_class[x])
                self.dummyContainers[x] = self.sTree.Branch(x,self.dummyContainers[x+"_array"],32000,-1)
                setattr(self.sTree,x,self.dummyContainers[x+"_array"])
                self.dummyContainers[x].Fill()
            #
        if self.sTree.GetBranch("GeoTracks"): self.sTree.SetBranchStatus("GeoTracks",0)
        # prepare for output
        # event header
        self.header  = ROOT.FairEventHeader()
        self.eventHeader  = self.sTree.Branch("ShipEventHeader",self.header,32000,-1)
        # fitted tracks
        self.fGenFitArray = ROOT.TClonesArray("genfit::Track")
        self.fGenFitArray.BypassStreamer(ROOT.kFALSE)
        self.fitTrack2MC  = ROOT.std.vector('int')()
        self.goodTracksVect  = ROOT.std.vector('int')()
        self.mcLink      = self.sTree.Branch("fitTrack2MC",self.fitTrack2MC,32000,-1)
        self.fitTracks   = self.sTree.Branch("FitTracks",  self.fGenFitArray,32000,-1)
        self.goodTracksBranch      = self.sTree.Branch("goodTracks",self.goodTracksVect,32000,-1)
        self.fTrackletsArray = ROOT.TClonesArray("Tracklet")
        self.Tracklets   = self.sTree.Branch("Tracklets",  self.fTrackletsArray,32000,-1)
        #
        self.digiStraw    = ROOT.TClonesArray("strawtubesHit")
        self.digiStrawBranch   = self.sTree.Branch("Digi_StrawtubesHits",self.digiStraw,32000,-1)
        self.digiSBT    = ROOT.TClonesArray("vetoHit")
        self.digiSBTBranch=self.sTree.Branch("Digi_SBTHits",self.digiSBT,32000,-1)
        self.vetoHitOnTrackArray    = ROOT.TClonesArray("vetoHitOnTrack")
        self.vetoHitOnTrackBranch=self.sTree.Branch("VetoHitOnTrack",self.vetoHitOnTrackArray,32000,-1)
        self.digiSBT2MC  = ROOT.std.vector('std::vector< int >')()
        self.mcLinkSBT   = self.sTree.Branch("digiSBT2MC",self.digiSBT2MC,32000,-1)
        self.digiTimeDet    = ROOT.TClonesArray("TimeDetHit")
        self.digiTimeDetBranch=self.sTree.Branch("Digi_TimeDetHits",self.digiTimeDet,32000,-1)
        self.digiMuon    = ROOT.TClonesArray("muonHit")
        self.digiMuonBranch=self.sTree.Branch("Digi_muonHits",self.digiMuon,32000,-1)
        # for the digitizing step
        self.v_drift = modules["Strawtubes"].StrawVdrift()
        self.sigma_spatial = modules["Strawtubes"].StrawSigmaSpatial()

        # setup ecal reconstruction
        self.caloTasks = []
        if self.sTree.GetBranch("EcalPoint"):
            # Creates. exports and fills calorimeter structure
            dflag = 0
            if debug: dflag = 10
            ecalGeo = ecalGeoFile+'z'+str(ShipGeo.ecal.z)+".geo"
            if not ecalGeo in os.listdir(os.environ["FAIRSHIP"]+"/geometry"): shipDet_conf.makeEcalGeoFile(ShipGeo.ecal.z,ShipGeo.ecal.File)
            ecalFiller=ROOT.ecalStructureFiller("ecalFiller", dflag,ecalGeo)
            ecalFiller.SetUseMCPoints(ROOT.kTRUE)
            ecalFiller.StoreTrackInformation()
            self.caloTasks.append(ecalFiller)
            #GeV -> ADC conversion
            ecalDigi=ROOT.ecalDigi("ecalDigi",0)
            self.caloTasks.append(ecalDigi)
            #ADC -> GeV conversion
            ecalPrepare=ROOT.ecalPrepare("ecalPrepare",0)
            self.caloTasks.append(ecalPrepare)
            # Maximums locator
            ecalMaximumFind=ROOT.ecalMaximumLocator("maximumFinder",dflag)
            self.caloTasks.append(ecalMaximumFind)
            # Cluster calibration
            ecalClusterCalib=ROOT.ecalClusterCalibration("ecalClusterCalibration", 0)
            #4x4 cm cells
            ecalCl3PhS=ROOT.TFormula("ecalCl3PhS", "[0]+x*([1]+x*([2]+x*[3]))")
            ecalCl3PhS.SetParameters(6.77797e-04, 5.75385e+00, 3.42690e-03, -1.16383e-04)
            ecalClusterCalib.SetStraightCalibration(3, ecalCl3PhS)
            ecalCl3Ph=ROOT.TFormula("ecalCl3Ph", "[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y")
            ecalCl3Ph.SetParameters(0.000750975, 5.7552, 0.00282783, -8.0025e-05, -0.000823651, 0.000111561)
            ecalClusterCalib.SetCalibration(3, ecalCl3Ph)
            #6x6 cm cells
            ecalCl2PhS=ROOT.TFormula("ecalCl2PhS", "[0]+x*([1]+x*([2]+x*[3]))")
            ecalCl2PhS.SetParameters(8.14724e-04, 5.67428e+00, 3.39030e-03, -1.28388e-04)
            ecalClusterCalib.SetStraightCalibration(2, ecalCl2PhS)
            ecalCl2Ph=ROOT.TFormula("ecalCl2Ph", "[0]+x*([1]+x*([2]+x*[3]))+[4]*x*y+[5]*x*y*y")
            ecalCl2Ph.SetParameters(0.000948095, 5.67471, 0.00339177, -0.000122629, -0.000169109, 8.33448e-06)
            ecalClusterCalib.SetCalibration(2, ecalCl2Ph)
            self.caloTasks.append(ecalClusterCalib)
            # Cluster finder
            ecalClusterFind=ROOT.ecalClusterFinder("clusterFinder",dflag)
            self.caloTasks.append(ecalClusterFind)
            # Calorimeter reconstruction
            ecalReco=ROOT.ecalReco('ecalReco',0)
            self.caloTasks.append(ecalReco)
            # Match reco to MC
            ecalMatch=ROOT.ecalMatch('ecalMatch',0)
            self.caloTasks.append(ecalMatch)
            if EcalDebugDraw:
                # ecal drawer: Draws calorimeter structure, incoming particles, clusters, maximums
                ecalDrawer=ROOT.ecalDrawer("clusterFinder",10)
                self.caloTasks.append(ecalDrawer)
                # add pid reco
            import shipPid
            self.caloTasks.append(shipPid.Task(self))
        # prepare vertexing
        self.Vertexing = shipVertex.Task(h,self.sTree)
        # setup random number generator
        self.random = ROOT.TRandom()
        ROOT.gRandom.SetSeed(13)
        self.PDG = ROOT.TDatabasePDG.Instance()
        # access ShipTree
        self.sTree.GetEvent(0)
        if len(self.caloTasks)>0:
            print "** initialize Calo reconstruction **"
            self.ecalStructure     = ecalFiller.InitPython(self.sTree.EcalPointLite)
            ecalDigi.InitPython(self.ecalStructure)
            ecalPrepare.InitPython(self.ecalStructure)
            self.ecalMaximums      = ecalMaximumFind.InitPython(self.ecalStructure)
            self.ecalCalib         = ecalClusterCalib.InitPython()
            self.ecalClusters      = ecalClusterFind.InitPython(self.ecalStructure, self.ecalMaximums, self.ecalCalib)
            self.EcalClusters = self.sTree.Branch("EcalClusters",self.ecalClusters,32000,-1)
            self.ecalReconstructed = ecalReco.InitPython(self.sTree.EcalClusters, self.ecalStructure, self.ecalCalib)
            self.EcalReconstructed = self.sTree.Branch("EcalReconstructed",self.ecalReconstructed,32000,-1)
            ecalMatch.InitPython(self.ecalStructure, self.ecalReconstructed, self.sTree.MCTrack)
            if EcalDebugDraw: ecalDrawer.InitPython(self.sTree.MCTrack, self.sTree.EcalPoint, self.ecalStructure, self.ecalClusters)
        else:
            ecalClusters      = ROOT.TClonesArray("ecalCluster")
            ecalReconstructed = ROOT.TClonesArray("ecalReconstructed")
            self.EcalClusters = self.sTree.Branch("EcalClusters",ecalClusters,32000,-1)
            self.EcalReconstructed = self.sTree.Branch("EcalReconstructed",ecalReconstructed,32000,-1)
        #
        # init geometry and mag. field
        gMan  = ROOT.gGeoManager
        self.geoMat =  ROOT.genfit.TGeoMaterialInterface()
        #
        self.bfield = ROOT.genfit.FairShipFields()
        self.fM = ROOT.genfit.FieldManager.getInstance()
        self.fM.init(self.bfield)
        ROOT.genfit.MaterialEffects.getInstance().init(self.geoMat)

        # init fitter, to be done before importing shipPatRec
        #fitter          = ROOT.genfit.KalmanFitter()
        #fitter          = ROOT.genfit.KalmanFitterRefTrack()
        self.fitter      = ROOT.genfit.DAF()
        self.fitter.setMaxIterations(50)
        if debug: self.fitter.setDebugLvl(1) # produces lot of printout
        #set to True if "real" pattern recognition is required also
        if debug == True: shipPatRec.debug = 1

        # for 'real' PatRec
        shipPatRec.initialize(fgeo)

    def reconstruct(self):
        ntracks = self.findTracks()
        nGoodTracks = self.findGoodTracks()
        self.linkVetoOnTracks()
        for x in self.caloTasks:
            if hasattr(x,'execute'): x.execute()
            elif x.GetName() == 'ecalFiller': x.Exec('start',self.sTree.EcalPointLite)
            elif x.GetName() == 'ecalMatch':  x.Exec('start',self.ecalReconstructed, self.sTree.MCTrack)
            else : x.Exec('start')
        if len(self.caloTasks)>0:
            self.EcalClusters.Fill()
            self.EcalReconstructed.Fill()
        if vertexing:
            # now go for 2-track combinations
            self.Vertexing.execute()

    def digitize(self):
        self.sTree.t0 = self.random.Rndm()*1*u.microsecond
        self.header.SetEventTime( self.sTree.t0 )
        self.header.SetRunId( self.sTree.MCEventHeader.GetRunID() )
        self.header.SetMCEntryNumber( self.sTree.MCEventHeader.GetEventID() )  # counts from 1
        self.eventHeader.Fill()
        self.digiSBT.Delete()
        self.digiSBT2MC.clear()
        self.digitizeSBT()
        self.digiSBTBranch.Fill()
        self.mcLinkSBT.Fill()
        self.digiStraw.Delete()
        self.digitizeStrawTubes()
        self.digiStrawBranch.Fill()
        self.digiTimeDet.Delete()
        self.digitizeTimeDet()
        self.digiTimeDetBranch.Fill()
        self.digiMuon.Delete()
        self.digitizeMuon()
        self.digiMuonBranch.Fill()

    def digitizeTimeDet(self):
        index = 0
        hitsPerDetId = {}
        for aMCPoint in self.sTree.TimeDetPoint:
            aHit = ROOT.TimeDetHit(aMCPoint,self.sTree.t0)
            if self.digiTimeDet.GetSize() == index: self.digiTimeDet.Expand(index+1000)
            self.digiTimeDet[index]=aHit
            detID = aHit.GetDetectorID()
            if aHit.isValid():
                if hitsPerDetId.has_key(detID):
                    t = aHit.GetMeasurements()
                    ct = aHit.GetMeasurements()
                    # this is not really correct, only first attempt
                    # case that one measurement only is earlier not taken into account
                    # SetTDC(Float_t val1, Float_t val2)
                    if  t[0]>ct[0] or t[1]>ct[1]:
                        # second hit with smaller tdc
                        self.digiTimeDet[hitsPerDetId[detID]].setInvalid()
                        hitsPerDetId[detID] = index
            index+=1

    def digitizeMuon(self):
        index = 0
        hitsPerDetId = {}
        for aMCPoint in self.sTree.muonPoint:
            aHit = ROOT.muonHit(aMCPoint,self.sTree.t0)
            if self.digiMuon.GetSize() == index: self.digiMuon.Expand(index+1000)
            self.digiMuon[index]=aHit
            detID = aHit.GetDetectorID()
            if aHit.isValid():
                if hitsPerDetId.has_key(detID):
                    if self.digiMuon[hitsPerDetId[detID]].GetDigi() > aHit.GetDigi():
                        # second hit with smaller tdc
                        self.digiMuon[hitsPerDetId[detID]].setValidity(0)
                        hitsPerDetId[detID] = index
            index+=1

    def digitizeSBT(self):
        ElossPerDetId    = {}
        tOfFlight        = {}
        listOfVetoPoints = {}
        key=-1
        for aMCPoint in self.sTree.vetoPoint:
            key+=1
            detID=aMCPoint.GetDetectorID()
            if not detID>100000: continue  # not a LiSc or plastic detector
            Eloss=aMCPoint.GetEnergyLoss()
            if not ElossPerDetId.has_key(detID):
                ElossPerDetId[detID]=0
                listOfVetoPoints[detID]=[]
                tOfFlight[detID]=[]
            ElossPerDetId[detID] += Eloss
            listOfVetoPoints[detID].append(key)
            tOfFlight[detID].append(aMCPoint.GetTime())
        index=0
        for seg in ElossPerDetId:
            aHit = ROOT.vetoHit(seg,ElossPerDetId[seg])
            aHit.SetTDC(min( tOfFlight[seg] ) + self.sTree.t0 )
            if self.digiSBT.GetSize() == index:
                self.digiSBT.Expand(index+1000)
            if seg<999999 and ElossPerDetId[seg]<0.045:    aHit.setInvalid()  # threshold for liquid scintillator, source Berlin group
            if seg>999999 and ElossPerDetId[seg]<0.001:    aHit.setInvalid()  # precise threshold for plastic to be determined
            self.digiSBT[index] = aHit
            v = ROOT.std.vector('int')()
            for x in listOfVetoPoints[seg]:
                v.push_back(x)
            self.digiSBT2MC.push_back(v)
            index=index+1
    def digitizeStrawTubes(self):
        # digitize FairSHiP MC hits
        index = 0
        hitsPerDetId = {}
        for aMCPoint in self.sTree.strawtubesPoint:
            aHit = ROOT.strawtubesHit(aMCPoint,self.sTree.t0)
            if self.digiStraw.GetSize() == index: self.digiStraw.Expand(index+1000)
            self.digiStraw[index]=aHit
            if aHit.isValid():
                detID = aHit.GetDetectorID()
                if hitsPerDetId.has_key(detID):
                    if self.digiStraw[hitsPerDetId[detID]].GetTDC() > aHit.GetTDC():
                        # second hit with smaller tdc
                        self.digiStraw[hitsPerDetId[detID]].setInvalid()
                        hitsPerDetId[detID] = index
            index+=1

    def withT0Estimate(self):
        # loop over all straw tdcs and make average, correct for ToF
        n = 0
        t0 = 0.
        key = -1
        SmearedHits = []
        v_drift = modules["Strawtubes"].StrawVdrift()
        modules["Strawtubes"].StrawEndPoints(10002001,start,stop)
        z1 = stop.z()
        for aDigi in self.digiStraw:
            key+=1
            if not aDigi.isValid: continue
            detID = aDigi.GetDetectorID()
            # don't use hits from straw veto
            station = int(detID/10000000)
            if station > 4 : continue
            modules["Strawtubes"].StrawEndPoints(detID,start,stop)
            delt1 = (start[2]-z1)/u.speedOfLight
            t0+=aDigi.GetDigi()-delt1
            SmearedHits.append( {'digiHit':key,'xtop':stop.x(),'ytop':stop.y(),'z':stop.z(),'xbot':start.x(),'ybot':start.y(),'dist':aDigi.GetDigi()} )
            n+=1
        if n>0: t0 = t0/n - 73.2*u.ns
        for s in SmearedHits:
            delt1 = (s['z']-z1)/u.speedOfLight
            s['dist'] = (s['dist'] -delt1 -t0)*v_drift
        return SmearedHits

    def smearHits(self,no_amb=None):
        # smear strawtube points
        SmearedHits = []
        key = -1
        v_drift = modules["Strawtubes"].StrawVdrift()
        modules["Strawtubes"].StrawEndPoints(10002001,start,stop)
        z1 = stop.z()
        for aDigi in self.digiStraw:
            key+=1
            if not aDigi.isValid: continue
            detID = aDigi.GetDetectorID()
            # don't use hits from straw veto
            station = int(detID/10000000)
            if station > 4 : continue
            modules["Strawtubes"].StrawEndPoints(detID,start,stop)
            #distance to wire
            delt1 = (start[2]-z1)/u.speedOfLight
            p=self.sTree.strawtubesPoint[key]
            # use true t0  construction:
            #     fdigi = t0 + p->GetTime() + t_drift + ( stop[0]-p->GetX() )/ speedOfLight;
            smear = (aDigi.GetDigi() - self.sTree.t0  - p.GetTime() - ( stop[0]-p.GetX() )/ u.speedOfLight) * v_drift
            if no_amb: smear = p.dist2Wire()
            SmearedHits.append( {'digiHit':key,'xtop':stop.x(),'ytop':stop.y(),'z':stop.z(),'xbot':start.x(),'ybot':start.y(),'dist':smear, 'detID':detID} )
            # Note: top.z()==bot.z() unless misaligned, so only add key 'z' to smearedHit
            if abs(stop.y())==abs(start.y()): h['disty'].Fill(smear)
            if abs(stop.y())>abs(start.y()): h['distu'].Fill(smear)
            if abs(stop.y())<abs(start.y()): h['distv'].Fill(smear)

        return SmearedHits

    def findTracks(self):

        hitPosLists    = {}
        stationCrossed = {}
        fittedtrackids=[]
        listOfIndices  = {}
        self.fGenFitArray.Delete()
        self.fTrackletsArray.Delete()
        self.fitTrack2MC.clear()

        #
        if withT0:  self.SmearedHits = self.withT0Estimate()
        # old procedure, not including estimation of t0
        else:       self.SmearedHits = self.smearHits(withNoStrawSmearing)

        nTrack = -1
        trackCandidates = []

        if realPR:

            # Do real PatRec
            track_hits = shipPatRec.execute(self.SmearedHits)

            # Create hitPosLists for track fit
            for i_track in track_hits.keys():

                atrack = track_hits[i_track]
                atrack_y12 = atrack['y12']
                atrack_stereo12 = atrack['stereo12']
                atrack_y34 = atrack['y34']
                atrack_stereo34 = atrack['stereo34']
                atrack_smeared_hits = list(atrack_y12) + list(atrack_stereo12) + list(atrack_y34) + list(atrack_stereo34)

                for sm in atrack_smeared_hits:

                    # detID = self.digiMufluxSpectrometer[sm['digiHit']].GetDetectorID()
                    detID = sm['detID']
                    station = int(detID/10000000)
                    trID = i_track

                    # For track fit
                    if not hitPosLists.has_key(trID):
                        hitPosLists[trID] = ROOT.std.vector('TVectorD')()
                        listOfIndices[trID] = []
                        stationCrossed[trID]  = {}
                    m = array('d',[sm['xtop'],sm['ytop'],sm['z'],sm['xbot'],sm['ybot'],sm['z'],sm['dist']])
                    hitPosLists[trID].push_back(ROOT.TVectorD(7,m))
                    listOfIndices[trID].append(sm['digiHit'])
                    if not stationCrossed[trID].has_key(station):
                        stationCrossed[trID][station] = 0
                    stationCrossed[trID][station] += 1

        else:
            # Fake PatRec
            track_hits = {}
            for sm in self.SmearedHits:

                detID = self.digiStraw[sm['digiHit']].GetDetectorID()
                station = int(detID/10000000)
                vnb = (detID - station * 10000000) // 1000000
                is_y12 = ((station == 1) + (station == 2)) * ((vnb == 0) + (vnb == 3))
                is_stereo12 = ((station == 1) + (station == 2)) * ((vnb == 1) + (vnb == 2))
                is_y34 = ((station == 3) + (station == 4)) * ((vnb == 0) + (vnb == 3))
                is_stereo34 = ((station == 3) + (station == 4)) * ((vnb == 1) + (vnb == 2))
                trID = self.sTree.strawtubesPoint[sm['digiHit']].GetTrackID()

                # For track fit
                if not hitPosLists.has_key(trID):
                    hitPosLists[trID] = ROOT.std.vector('TVectorD')()
                    listOfIndices[trID] = []
                    stationCrossed[trID]  = {}
                m = array('d',[sm['xtop'],sm['ytop'],sm['z'],sm['xbot'],sm['ybot'],sm['z'],sm['dist']])
                hitPosLists[trID].push_back(ROOT.TVectorD(7,m))
                listOfIndices[trID].append(sm['digiHit'])
                if not stationCrossed[trID].has_key(station):
                    stationCrossed[trID][station] = 0
                stationCrossed[trID][station] += 1

                # For PatRec metrics
                if not track_hits.has_key(trID):
                    atrack = {'y12': [], 'stereo12': [], 'y34': [], 'stereo34': []}
                    track_hits[trID] = atrack
                if is_y12:
                    track_hits[trID]['y12'].append(sm)
                if is_stereo12:
                    track_hits[trID]['stereo12'].append(sm)
                if is_y34:
                    track_hits[trID]['y34'].append(sm)
                if is_stereo34:
                    track_hits[trID]['stereo34'].append(sm)

        # Save track hit indeces
        for atrack in listOfIndices:
            # make tracklets out of trackCandidates, just for testing, should be output of proper pattern recognition
            nTracks   = self.fTrackletsArray.GetEntries()
            aTracklet  = self.fTrackletsArray.ConstructedAt(nTracks)
            listOfHits = aTracklet.getList()
            aTracklet.setType(3)
            for index in listOfIndices[atrack]:
                listOfHits.push_back(index)

        # Track fit
        for atrack in hitPosLists:
            if atrack < 0:
                continue # these are hits not assigned to MC track because low E cut
            pdg = self.sTree.MCTrack[atrack].GetPdgCode()
            if not self.PDG.GetParticle(pdg):
                continue # unknown particle
            meas = hitPosLists[atrack]
            nM = meas.size()
            if nM < 25 :
                continue # not enough hits to make a good trackfit
            if len(stationCrossed[atrack]) < 3 :
                continue  # not enough stations crossed to make a good trackfit
            if debug:
                mctrack = self.sTree.MCTrack[atrack]
            charge = self.PDG.GetParticle(pdg).Charge()/(3.)
            posM = ROOT.TVector3(0, 0, 0)
            momM = ROOT.TVector3(0,0,3.*u.GeV)
            # approximate covariance
            covM = ROOT.TMatrixDSym(6)
            resolution = self.sigma_spatial
            if withT0: resolution = resolution*1.4 # worse resolution due to t0 estimate
            for  i in range(3):   covM[i][i] = resolution*resolution
            covM[0][0]=resolution*resolution*100.
            for  i in range(3,6): covM[i][i] = ROOT.TMath.Power(resolution / nM / ROOT.TMath.Sqrt(3), 2)
            # trackrep
            rep = ROOT.genfit.RKTrackRep(pdg)
            # smeared start state
            stateSmeared = ROOT.genfit.MeasuredStateOnPlane(rep)
            rep.setPosMomCov(stateSmeared, posM, momM, covM)
            # create track
            seedState = ROOT.TVectorD(6)
            seedCov   = ROOT.TMatrixDSym(6)
            rep.get6DStateCov(stateSmeared, seedState, seedCov)
            theTrack = ROOT.genfit.Track(rep, seedState, seedCov)
            hitCov = ROOT.TMatrixDSym(7)
            hitCov[6][6] = resolution*resolution
            for m in meas:
                tp = ROOT.genfit.TrackPoint(theTrack) # note how the point is told which track it belongs to
                measurement = ROOT.genfit.WireMeasurement(m,hitCov,1,6,tp) # the measurement is told which trackpoint it belongs to
                # print measurement.getMaxDistance()
                measurement.setMaxDistance(ShipGeo.strawtubes.InnerStrawDiameter/2.)
                # measurement.setLeftRightResolution(-1)
                tp.addRawMeasurement(measurement) # package measurement in the TrackPoint
                theTrack.insertPoint(tp)  # add point to Track
                # print "debug meas",atrack,nM,stationCrossed[atrack],self.sTree.MCTrack[atrack],pdg
            trackCandidates.append([theTrack,atrack])

        # Check fitted tracks and save them
        reco_mom = {}
        for entry in trackCandidates:
            #check
            atrack = entry[1]
            theTrack = entry[0]
            if not theTrack.checkConsistency():
                print 'Problem with track before fit, not consistent',atrack,theTrack
                continue
            # do the fit
            try:
                self.fitter.processTrack(theTrack) # processTrackWithRep(theTrack,rep,True)
            except:
                if debug:
                    print "genfit failed to fit track"
                error = "genfit failed to fit track"
                ut.reportError(error)
                continue
            #check
            if not theTrack.checkConsistency():
                if debug:
                    print 'Problem with track after fit, not consistent',atrack,theTrack
                error = "Problem with track after fit, not consistent"
                ut.reportError(error)
                continue
            try:
                fittedState = theTrack.getFittedState()
                fittedMom = fittedState.getMomMag()
                Px,Py,Pz = fittedState.getMom().x(),fittedState.getMom().y(),fittedState.getMom().z()
                P = fittedMom
                Pt = math.sqrt(Px**2 + Py**2)
                reco_mom[atrack] = {'Px': Px, 'Py': Py, 'Pz': Pz, 'P': P, 'Pt': Pt}

                #rc,pos,mom = TrackExtrapolateTool.extrapolateToPlane(theTrack, ShipGeo.Bfield.z)
                #StartX, StartY, StartZ = self.getStartAtOrigin(atrack)
                #print "StartX, StartY, StartZ: ", StartX, StartY, StartZ
            except:
                error = "problem with fittedstate"
                print error
                ut.reportError(error)
                continue
            fitStatus   = theTrack.getFitStatus()
            nmeas = fitStatus.getNdf()
            chi2        = fitStatus.getChi2()/nmeas
            h['chi2'].Fill(chi2)
            h['nmeas'].Fill(nmeas)
            # make track persistent
            nTrack   = self.fGenFitArray.GetEntries()
            if not debug:
                theTrack.prune("CFL")  #  http://sourceforge.net/p/genfit/code/HEAD/tree/trunk/core/include/Track.h#l280
            self.fGenFitArray[nTrack] = theTrack
            self.fitTrack2MC.push_back(atrack)
            if debug:
                print 'save track',theTrack,chi2,nmeas,fitStatus.isFitConverged()
        self.Tracklets.Fill()
        self.fitTracks.Fill()
        self.mcLink.Fill()
        # debug
        if debug:
            print 'save tracklets:'
            for x in self.sTree.Tracklets:
                print x.getType(),x.getList().size()


        #### PatRec metrics
        track_hit_ids = {}
        for trID in track_hits.keys():
            a_track_hits = track_hits[trID]
            a_track_ids = {'y12': [], 'stereo12': [], 'y34': [], 'stereo34': []}
            a_track_ids['y12'] = [self.sTree.strawtubesPoint[hit['digiHit']].GetTrackID() for hit in a_track_hits['y12']]
            a_track_ids['stereo12'] = [self.sTree.strawtubesPoint[hit['digiHit']].GetTrackID() for hit in a_track_hits['stereo12']]
            a_track_ids['y34'] = [self.sTree.strawtubesPoint[hit['digiHit']].GetTrackID() for hit in a_track_hits['y34']]
            a_track_ids['stereo34'] = [self.sTree.strawtubesPoint[hit['digiHit']].GetTrackID() for hit in a_track_hits['stereo34']]
            track_hit_ids[trID] = a_track_ids

        # print "track_hit_ids: ", track_hit_ids
        # print "reco_mom: ", reco_mom

        # Number of hits per track
        n_hits_true_y12 = {}
        n_hits_true_stereo12 = {}
        n_hits_true_12 = {}
        n_hits_true_y34 = {}
        n_hits_true_stereo34 = {}
        n_hits_true_34 = {}
        n_hits_true_1234 = {}
        n_hits_true_1 = {}
        n_hits_true_2 = {}
        n_hits_true_3 = {}
        n_hits_true_4 = {}
        n_tracks_per_event = 0

        for ahit in self.sTree.strawtubesPoint:
            track_id = ahit.GetTrackID()
            pdg = abs(ahit.PdgCode())
            detID = ahit.GetDetectorID()
            statnb = detID // 10000000
            vnb = (detID - statnb * 10000000) // 1000000

            is_y12 = ((statnb == 1) + (statnb == 2)) * ((vnb == 0) + (vnb == 3))
            is_stereo12 = ((statnb == 1) + (statnb == 2)) * ((vnb == 1) + (vnb == 2))
            is_12 = ((statnb == 1) + (statnb == 2))
            is_y34 = ((statnb == 3) + (statnb == 4)) * ((vnb == 0) + (vnb == 3))
            is_stereo34 = ((statnb == 3) + (statnb == 4)) * ((vnb == 1) + (vnb == 2))
            is_34 = ((statnb == 3) + (statnb == 4))
            is_1 = statnb == 1
            is_2 = statnb == 2
            is_3 = statnb == 3
            is_4 = statnb == 4

            if not n_hits_true_1234.has_key(track_id):
                n_hits_true_y12[track_id] = 0
                n_hits_true_stereo12[track_id] = 0
                n_hits_true_12[track_id] = 0
                n_hits_true_y34[track_id] = 0
                n_hits_true_stereo34[track_id] = 0
                n_hits_true_34[track_id] = 0
                n_hits_true_1234[track_id] = 0
                n_hits_true_1[track_id] = 0
                n_hits_true_2[track_id] = 0
                n_hits_true_3[track_id] = 0
                n_hits_true_4[track_id] = 0

                if pdg == 11: h['pdg'].Fill("Electrons", 1)
                elif pdg == 13: h['pdg'].Fill("Muons", 1)
                elif pdg == 211: h['pdg'].Fill("Pions", 1)
                elif pdg == 2212: h['pdg'].Fill("Protons", 1)
                elif pdg == 321: h['pdg'].Fill("Kaons", 1)
                else : h['pdg'].Fill("Others", 1)

                Ptruth, Ptruthx, Ptruthy, Ptruthz = self.getPtruthFirst(track_id)
                Pttruth = math.sqrt(Ptruthx**2 + Ptruthy**2)
                h['p/pt_truth'].Fill(Ptruth,Pttruth)
                h['p_truth'].Fill(Ptruth)
                h['pt_truth'].Fill(Pttruth)
                h['pz_truth'].Fill(Ptruthz)
                h['px_truth'].Fill(Ptruthx)
                h['py_truth'].Fill(Ptruthy)

                n_tracks_per_event += 1



            if is_y12: n_hits_true_y12[track_id] += 1
            if is_stereo12: n_hits_true_stereo12[track_id] += 1
            if is_12: n_hits_true_12[track_id] += 1
            if is_y34: n_hits_true_y34[track_id] += 1
            if is_stereo34: n_hits_true_stereo34[track_id] += 1
            if is_34: n_hits_true_34[track_id] += 1
            if is_1: n_hits_true_1[track_id] += 1
            if is_2: n_hits_true_2[track_id] += 1
            if is_3: n_hits_true_3[track_id] += 1
            if is_4: n_hits_true_4[track_id] += 1
            n_hits_true_1234[track_id] += 1

        for track_id in n_hits_true_1234.keys():
            h['hits_true_y12'].Fill(n_hits_true_y12[track_id])
            h['hits_true_stereo12'].Fill(n_hits_true_stereo12[track_id])
            h['hits_true_12'].Fill(n_hits_true_12[track_id])
            h['hits_true_y34'].Fill(n_hits_true_y34[track_id])
            h['hits_true_stereo34'].Fill(n_hits_true_stereo34[track_id])
            h['hits_true_34'].Fill(n_hits_true_34[track_id])
            h['hits_true_1234'].Fill(n_hits_true_1234[track_id])
            h['hits_true_1'].Fill(n_hits_true_1[track_id])
            h['hits_true_2'].Fill(n_hits_true_2[track_id])
            h['hits_true_3'].Fill(n_hits_true_3[track_id])
            h['hits_true_4'].Fill(n_hits_true_4[track_id])

        h['tracks_per_event_true'].Fill(n_tracks_per_event)

        for i_track in track_hit_ids.keys():
            atrack = track_hit_ids[i_track]
            atrack_y12 = atrack["y12"]
            frac_y12, tmax_y12 = self.fracMCsame(atrack_y12)

            Ptruth, Ptruthx, Ptruthy, Ptruthz = self.getPtruthFirst(tmax_y12)
            Pttruth = math.sqrt(Ptruthx**2 + Ptruthy**2)

            # StartX, StartY, StartZ = self.getStartAtOrigin(tmax_y12)
            # print "StartX, StartY, StartZ: ", StartX, StartY, StartZ

            if reco_mom.has_key(i_track):
                mom = reco_mom[i_track]
                Px, Py, Pz, P, Pt = mom['Px'], mom['Py'], mom['Pz'], mom['P'], mom['Pt']

                h['p/pt_reco'].Fill(P,Pt)
                h['p_reco'].Fill(P)
                h['pt_reco'].Fill(Pt)
                h['pz_reco'].Fill(Pz)
                h['px_reco'].Fill(Px)
                h['py_reco'].Fill(Py)
                h['p_rel_error'].Fill((P - Ptruth) / Ptruth)
                h['pt_rel_error'].Fill((Pt - Pttruth) / Pttruth)
                h['px_rel_error'].Fill((Px - Ptruthx) / Ptruthx)
                h['py_rel_error'].Fill((Py - Ptruthy) / Ptruthy)
                h['pz_rel_error'].Fill((Pz - Ptruthz) / Ptruthz)

        (n_tracks, n_recognized, n_clones, n_ghosts, n_others) = self.track_metrics_per_view(track_hit_ids, view_name='y12')
        h['Reco_y12'].Fill("N total", n_tracks)
        h['Reco_y12'].Fill("N recognized tracks", n_recognized)
        h['Reco_y12'].Fill("N clones", n_clones)
        h['Reco_y12'].Fill("N ghosts", n_ghosts)
        h['Reco_y12'].Fill("N others", n_others)

        (n_tracks, n_recognized, n_clones, n_ghosts, n_others) = self.track_metrics_per_view(track_hit_ids, view_name='stereo12')
        h['Reco_stereo12'].Fill("N total", n_tracks)
        h['Reco_stereo12'].Fill("N recognized tracks", n_recognized)
        h['Reco_stereo12'].Fill("N clones", n_clones)
        h['Reco_stereo12'].Fill("N ghosts", n_ghosts)
        h['Reco_stereo12'].Fill("N others", n_others)

        (n_tracks, n_recognized, n_clones, n_ghosts, n_others) = self.track_metrics_per_view(track_hit_ids, view_name='y34')
        h['Reco_y34'].Fill("N total", n_tracks)
        h['Reco_y34'].Fill("N recognized tracks", n_recognized)
        h['Reco_y34'].Fill("N clones", n_clones)
        h['Reco_y34'].Fill("N ghosts", n_ghosts)
        h['Reco_y34'].Fill("N others", n_others)

        (n_tracks, n_recognized, n_clones, n_ghosts, n_others) = self.track_metrics_per_view(track_hit_ids, view_name='stereo34')
        h['Reco_stereo34'].Fill("N total", n_tracks)
        h['Reco_stereo34'].Fill("N recognized tracks", n_recognized)
        h['Reco_stereo34'].Fill("N clones", n_clones)
        h['Reco_stereo34'].Fill("N ghosts", n_ghosts)
        h['Reco_stereo34'].Fill("N others", n_others)



        (n_tracks, n_recognized, n_clones, n_ghosts, n_others) = self.full_track_metrics(track_hit_ids, mode='one_hit_per_T14')
        h['Reco_T14'].Fill("N total", n_tracks)
        h['Reco_T14'].Fill("N recognized tracks", n_recognized)
        h['Reco_T14'].Fill("N clones", n_clones)
        h['Reco_T14'].Fill("N ghosts", n_ghosts)
        h['Reco_T14'].Fill("N others", n_others)


        (n_tracks, n_recognized, n_clones, n_ghosts, n_others) = self.full_track_metrics(track_hit_ids, mode='3_hits_per_views')
        h['Reco_3hits_per_views'].Fill("N total", n_tracks)
        h['Reco_3hits_per_views'].Fill("N recognized tracks", n_recognized)
        h['Reco_3hits_per_views'].Fill("N clones", n_clones)
        h['Reco_3hits_per_views'].Fill("N ghosts", n_ghosts)
        h['Reco_3hits_per_views'].Fill("N others", n_others)






        return nTrack+1

    def findGoodTracks(self):
        self.goodTracksVect.clear()
        nGoodTracks = 0
        for i,track in enumerate(self.fGenFitArray):
            fitStatus = track.getFitStatus()
            if not fitStatus.isFitConverged(): continue
            nmeas = fitStatus.getNdf()
            chi2  = fitStatus.getChi2()/nmeas
            if chi2<50 and not chi2<0:
                self.goodTracksVect.push_back(i)
                nGoodTracks+=1
        self.goodTracksBranch.Fill()
        return nGoodTracks

    def findVetoHitOnTrack(self,track):
        distMin = 99999.
        vetoHitOnTrack = ROOT.vetoHitOnTrack()
        xx  = track.getFittedState()
        rep   = ROOT.genfit.RKTrackRep(xx.getPDG())
        state = ROOT.genfit.StateOnPlane(rep)
        rep.setPosMom(state,xx.getPos(),xx.getMom())
        for i,vetoHit in enumerate(self.digiSBT):
            vetoHitPos = vetoHit.GetXYZ()
            try:
                rep.extrapolateToPoint(state,vetoHitPos,False)
            except:
                error =  "shipDigiReco::findVetoHitOnTrack extrapolation did not worked"
                ut.reportError(error)
                if debug: print error
                continue
            dist = (rep.getPos(state) - vetoHitPos).Mag()
            if dist < distMin:
                distMin = dist
                vetoHitOnTrack.SetDist(distMin)
                vetoHitOnTrack.SetHitID(i)
        return vetoHitOnTrack

    def linkVetoOnTracks(self):
        self.vetoHitOnTrackArray.Delete()
        index = 0
        for goodTrak in self.goodTracksVect:
            track = self.fGenFitArray[goodTrak]
            if self.vetoHitOnTrackArray.GetSize() == index: self.vetoHitOnTrackArray.Expand(index+1000)
            self.vetoHitOnTrackArray[index] = self.findVetoHitOnTrack(track)
            index+=1
        self.vetoHitOnTrackBranch.Fill()

    def finish(self):
        del self.fitter
        print 'finished writing tree'
        self.sTree.Write()
        ut.errorSummary()
        ut.writeHists(h,"recohists.root")
        if realPR: shipPatRec.finalize()


    #### PatRec metric functions

    def getPtruthFirst(self,mcPartKey):
        Ptruth,Ptruthx,Ptruthy,Ptruthz = -1.,-1.,-1.,-1.
        for ahit in self.sTree.strawtubesPoint:
            if ahit.GetTrackID() == mcPartKey:
                Ptruthx,Ptruthy,Ptruthz = ahit.GetPx(),ahit.GetPy(),ahit.GetPz()
                Ptruth  = ROOT.TMath.Sqrt(Ptruthx**2+Ptruthy**2+Ptruthz**2)
                break
        return Ptruth,Ptruthx,Ptruthy,Ptruthz

    def getStartAtOrigin(self,mcPartKey):
        Ptruth,Ptruthx,Ptruthy,Ptruthz = -1.,-1.,-1.,-1.
        atrack=self.sTree.MCTrack[mcPartKey]
        StartX= atrack.GetStartX()
        StartY= atrack.GetStartY()
        StartZ= atrack.GetStartZ()
        return StartX, StartY, StartZ

    def fracMCsame(self, trackids):

        track = {}
        nh = len(trackids)

        for tid in trackids:
            if track.has_key(tid):
                track[tid] += 1
            else:
                track[tid] = 1

        # now get track with largest number of hits
        if track != {}:
            tmax = max(track, key=track.get)
        else:
            track = {-999:0}
            tmax = -999

        frac = 0.
        if nh > 0:
            frac = float(track[tmax]) / float(nh)

        return frac,tmax

    def track_metrics_per_view(self, track_hit_ids, view_name='y12'):

        n_tracks = 0
        n_recognized = 0
        n_clones = 0
        n_ghosts = 0
        n_others = 0

        true_track_ids = []
        for ahit in self.sTree.strawtubesPoint:
            track_id = ahit.GetTrackID()
            detID = ahit.GetDetectorID()
            station = int(detID/10000000)
            vnb = (detID - station * 10000000) // 1000000
            is_y12 = ((station == 1) + (station == 2)) * ((vnb == 0) + (vnb == 3))
            is_stereo12 = ((station == 1) + (station == 2)) * ((vnb == 1) + (vnb == 2))
            is_y34 = ((station == 3) + (station == 4)) * ((vnb == 0) + (vnb == 3))
            is_stereo34 = ((station == 3) + (station == 4)) * ((vnb == 1) + (vnb == 2))
            is_sel = (is_y12 * (view_name == 'y12')) + \
                     (is_stereo12 * (view_name == 'stereo12'))  + \
                     (is_y34 * (view_name == 'y34')) + \
                     (is_stereo34 * (view_name == 'stereo34'))
            if is_sel and track_id not in true_track_ids:
                true_track_ids.append(track_id)

        n_tracks = len(true_track_ids)
        found_track_ids = []

        min_hits = 3
        min_eff = 0.7

        for i_track in track_hit_ids.keys():

            atrack = track_hit_ids[i_track]
            atrack_view = atrack[view_name]

            if len(atrack_view) == 0:
                continue

            if len(atrack_view) >= min_hits:

                frac_view, tmax_view = self.fracMCsame(atrack_view)

                if frac_view >= min_eff:

                    if tmax_view in true_track_ids and tmax_view not in found_track_ids:
                        n_recognized += 1
                        found_track_ids.append(tmax_view)
                    elif tmax_view in true_track_ids and tmax_view in found_track_ids:
                        n_clones += 1
                    elif tmax_view not in true_track_ids:
                        n_others += 1

                else:
                    n_ghosts += 1

            else:
                n_ghosts += 1

        output = (n_tracks, n_recognized, n_clones, n_ghosts, n_others)

        return output

    def full_track_metrics(self, track_hit_ids, mode='one_hit_per_T14'):

        n_tracks = 0
        n_recognized = 0
        n_clones = 0
        n_ghosts = 0
        n_others = 0

        true_track_ids = []

        if mode == 'one_hit_per_T14':

            stationCrossed = {}
            for ahit in self.sTree.strawtubesPoint:
                track_id = ahit.GetTrackID()
                detID = ahit.GetDetectorID()
                station = int(detID/10000000)
                if not stationCrossed.has_key(track_id):
                    stationCrossed[track_id]  = {}
                if not stationCrossed[track_id].has_key(station):
                    stationCrossed[track_id][station] = 0
                stationCrossed[track_id][station] += 1

            for track_id in stationCrossed.keys():
                crossed = stationCrossed[track_id]
                if crossed.has_key(1) and crossed.has_key(2) and crossed.has_key(3) and crossed.has_key(4):
                    if track_id not in true_track_ids:
                        true_track_ids.append(track_id)

        elif mode == '3_hits_per_views':

            n_true_track_hits_y12 = {}
            n_true_track_hits_stereo12 = {}
            n_true_track_hits_y34 = {}
            n_true_track_hits_stereo34 = {}

            for ahit in self.sTree.strawtubesPoint:

                track_id = ahit.GetTrackID()

                detID = ahit.GetDetectorID()
                statnb = detID // 10000000
                vnb = (detID - statnb * 10000000) // 1000000

                is_y12 = ((statnb == 1) + (statnb == 2)) * ((vnb == 0) + (vnb == 3))
                is_stereo12 = ((statnb == 1) + (statnb == 2)) * ((vnb == 1) + (vnb == 2))
                is_y34 = ((statnb == 3) + (statnb == 4)) * ((vnb == 0) + (vnb == 3))
                is_stereo34 = ((statnb == 3) + (statnb == 4)) * ((vnb == 1) + (vnb == 2))

                if is_y12:
                    if n_true_track_hits_y12.has_key(track_id):
                        n_true_track_hits_y12[track_id] += 1
                    else:
                        n_true_track_hits_y12[track_id] = 1

                if is_stereo12:
                    if n_true_track_hits_stereo12.has_key(track_id):
                        n_true_track_hits_stereo12[track_id] += 1
                    else:
                        n_true_track_hits_stereo12[track_id] = 1

                if is_y34:
                    if n_true_track_hits_y34.has_key(track_id):
                        n_true_track_hits_y34[track_id] += 1
                    else:
                        n_true_track_hits_y34[track_id] = 1

                if is_stereo34:
                    if n_true_track_hits_stereo34.has_key(track_id):
                        n_true_track_hits_stereo34[track_id] += 1
                    else:
                        n_true_track_hits_stereo34[track_id] = 1

            min_hits = 3
            for key in n_true_track_hits_y12.keys():
                if n_true_track_hits_y12[key] >= min_hits:
                    if n_true_track_hits_stereo12.has_key(key):
                        if n_true_track_hits_stereo12[key] >= min_hits:
                            if n_true_track_hits_y34.has_key(key):
                                if n_true_track_hits_y34[key] >= min_hits:
                                    if n_true_track_hits_stereo34.has_key(key):
                                        if n_true_track_hits_stereo34[key] >= min_hits:
                                            true_track_ids.append(key)



        n_tracks = len(true_track_ids)
        found_track_ids = []

        min_hits = 1
        min_eff = 0.7

        for i_track in track_hit_ids.keys():

            atrack = track_hit_ids[i_track]
            atrack_y12 = atrack['y12']
            atrack_stereo12 = atrack['stereo12']
            atrack_y34 = atrack['y34']
            atrack_stereo34 = atrack['stereo34']

            if len(atrack_y12) >= min_hits and len(atrack_stereo12) >= min_hits and len(atrack_y34) >= min_hits and len(atrack_stereo34) >= min_hits:

                frac_y12, tmax_y12 = self.fracMCsame(atrack_y12)
                frac_stereo12, tmax_stereo12 = self.fracMCsame(atrack_stereo12)
                frac_y34, tmax_y34 = self.fracMCsame(atrack_y34)
                frac_stereo34, tmax_stereo34 = self.fracMCsame(atrack_stereo34)

                if tmax_y12 == tmax_stereo12 and tmax_y12 == tmax_y34 and tmax_y12 == tmax_stereo34:
                    if frac_y12 >= min_eff and frac_stereo12 >= min_eff and frac_y34 >= min_eff and frac_stereo34 >= min_eff:

                        if tmax_y12 in true_track_ids and tmax_y12 not in found_track_ids:
                            n_recognized += 1
                            found_track_ids.append(tmax_y12)
                        elif tmax_y12 in true_track_ids and tmax_y12 in found_track_ids:
                            n_clones += 1
                        elif tmax_y12 not in true_track_ids:
                            n_others += 1

                    else:
                        n_ghosts += 1
                else:
                    n_ghosts += 1

            else:
                n_ghosts += 1

        output = (n_tracks, n_recognized, n_clones, n_ghosts, n_others)

        return output