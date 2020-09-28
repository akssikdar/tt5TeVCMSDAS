from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
import os.path

import bamboo.treedecorators as btd
class ReadVariableVarWithSuffix(btd.ReadVariableVarWithSuffix):
    def getVarName(self, branchName, collgrpname=None):
        variNm = btd.normVarName(branchName[len(self.prefix):].lstrip(self.sep))
        return self.prefix, f"{self.systName}{variNm}" if variNm else self.nomName
class ReadMuScaleCorrection(btd.NanoSystematicVarSpec):
    def appliesTo(self, name):
        return name == "Muon"
    def nomName(self, name):
        return ""
    def getVarName(self, name, collgrpname=None):
        if name.startswith("corrected"):
            vari, v = name[len("corrected"):].split("_")
            return v, f"{self.systName}{btd.normVarName(vari)}" if vari else self.nomName(name)
class ReadElScaleCorrection(btd.NanoSystematicVarSpec):
    def appliesTo(self, name):
        return name == "Electron"
    def getVarName(self, name, collgrpname=None):
        if name.split("_")[0] in ("pt", "mass") and len(name.split("_")) >= 2:
            return name.split("_")[0], "_".join(name.split("_")[1:])
    def nomName(self, name):
        return "nom"
_descr_5TeV_removedGroups = ["CaloMET_", "ChsMET_", "RawMET_", "TkMET_"]
_descr_5TeV_removedGroups_MC = _descr_5TeV_removedGroups + ["Generator_", "HTXS_"]
_descr_5TeV_removedCollections = ["nFatJet", "nIsoTrack", "nOtherPV", "nSV", "nSoftActivityJet", "nSubJet", "nTau", "nTrigObj"]
_descr_5TeV_removedCollections_MC = _descr_5TeV_removedCollections + ["nGenJetAK8", "nGenVisTau", "nSubGenJetAK8"]
description_nanov5_5TeV_data = btd.NanoAODDescription.get("v5", year="2018", isMC=False,
        removeGroups=_descr_5TeV_removedGroups, removeCollections=_descr_5TeV_removedCollections,
        systVariations=[
            ReadMuScaleCorrection("muScale"),
            ReadElScaleCorrection("elScale"), ## FIXME implement on the fly
            btd.ReadJetMETVar("Jet", "MET", jetsNomName="nom", metNomName="nom")
            ])
description_nanov5_5TeV_MC = btd.NanoAODDescription.get("v5", year="2018", isMC=True,
        removeGroups=_descr_5TeV_removedGroups_MC, removeCollections=_descr_5TeV_removedCollections_MC,
        addGroups=["Pileup_"], addCollections=["nLHEPart"],
        systVariations=[
            ReadVariableVarWithSuffix("puWeight"),
            ReadVariableVarWithSuffix("PrefireWeight"),
            ReadMuScaleCorrection("muScale"),
            ReadElScaleCorrection("elScale"), ## FIXME implement on the fly
            btd.ReadJetMETVar("Jet", "MET", jetsExclVars=["raw"], metNomName="nom", metExclVars=["raw", "nom"])
            ])

## no need to propagate the scale uncertainty into the scalefactor, always use nominal
import bamboo.scalefactors
binningVariables_nano_noScaleSyst = dict(bamboo.scalefactors.binningVariables_nano)
binningVariables_nano_noScaleSyst.update({ "Pt" : lambda obj : obj.pt.op.arg.wrapped.result[obj.idx] })
## for use with bamboo.scalefactors.get_scalefactor
## e.g. get_scalefactor("lepton", "Muon_RecoToLoose", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="muLoose")
scalefactors_lepMVA = {
    f"{lFlav}_{ratioSel}" : os.path.join(os.path.dirname(__file__), "data", f"{lFlav}_{ratioSel}SF.json")
        for ratioSel in ("RecoToLoose", "LooseToTight")
        for lFlav in ("Electron", "Muon")
    }

class Nano5TeVBase(NanoAODModule):
    """ Base module for postprocessed 5TeV NanoAODv5 samples """
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        ## Decorate the tree
        tree,noSel,be,lumiArgs = super(Nano5TeVBase, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg,
                description=(description_nanov5_5TeV_MC if self.isMC(sample) else description_nanov5_5TeV_data))
        return tree,noSel,be,lumiArgs

class Nano5TeVHistoModule(Nano5TeVBase, NanoAODHistoModule):
    pass

class HelloWorld(Nano5TeVHistoModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)
        if isMC:
            noSel = noSel.refine("mcWeight", weight=[ t.genWeight ])
        noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        plots = []

        muons = op.select(t.Muon, lambda mu : mu.pt > 20.)
        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])
        plots.append(Plot.make1D("dimu_M",
            op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel, EqB(100, 20., 120.),
            title="Dimuon invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        return plots
ttnomname = 'TT'

### Channel to ints
class ch():
  ElMu = 0
  MuMu = 1
  ElEl = 2
  Muon = 3
  Elec = 4
chan = {ch.ElMu:'ElMu', ch.MuMu:'MuMu', ch.ElEl:'ElEl'}

### Levels to ints
class lev():
  dilepton = 0
  ZVeto    = 1
  MET      = 2
  jets2    = 3
  jets2nomet = 4
  btag1    = 5
  ww       = 6
level = {lev.dilepton:'dilepton', lev.ZVeto:'ZVeto', lev.MET:'MET', lev.jets2:'2jets', lev.jets2nomet:'2jetsnomet', lev.btag1:'1btag', lev.ww:'ww'}
invlevel = {'dilepton':lev.dilepton, 'ZVeto':lev.ZVeto, 'MET':lev.MET, '2jets':lev.jets2, '2jetsnomet':lev.jets2nomet, '1btag':lev.btag1, 'ww':lev.ww}

### Systematic uncertainties
class systematic():
  nom       = -1 # Nominal
  MuonEffUp = 0
  MuonEffDo = 1
  ElecEffUp = 2
  ElecEffDo = 3
  JESUp = 4
  JESDo = 5
  JERUp = 6
  JERDo = 7
  PUUp = 8
  PUDo = 9
  TrigEffUp = 10
  TrigEffDo = 11
  PrefireUp = 12
  PrefireDo = 13
  BtagUp = 14
  BtagDo = 15
  MisTagUp = 16
  MisTagDo = 17
  ISRUp = 18
  ISRDo = 19
  FSRUp = 20
  FSRDo = 21
#systlabel = {systematic.nom:'', systematic.MuonEffUp:'MuonEffUp', systematic.MuonEffDo:'MuonEffDown', systematic.ElecEffUp:'ElecEffUp', systematic.ElecEffDo:'ElecEffDown', systematic.JESUp:'JESUp', systematic.JESDo:'JESDown', systematic.JERUp:'JERUp', systematic.JERDo:'JERDown', systematic.PUUp:'PUUp', systematic.PUDo:'PUDown', systematic.TrigEffUp:'TrigEffUp', systematic.TrigEffDo:'TrigEffDown', systematic.PrefireUp:'PrefireUp', systematic.PrefireDo:'PrefireDown', systematic.BtagUp:'BtagUp', systematic.BtagDo:'BtagDown', systematic.MisTagUp:'MisTagUp', systematic.MisTagDo:'MisTagDown', systematic.ISRUp:'ISRUp', systematic.ISRDo:'ISRDown', systematic.FSRUp:'FSRUp', systematic.FSRDo:'FSRDown'}
systlabel = {systematic.nom:'', systematic.MuonEffUp:'MuonEffUp', systematic.MuonEffDo:'MuonEffDown', systematic.ElecEffUp:'ElecEffUp', systematic.ElecEffDo:'ElecEffDown', systematic.TrigEffUp:'TrigEffUp', systematic.TrigEffDo:'TrigEffDown', systematic.PrefireUp:'PrefireUp', systematic.PrefireDo:'PrefireDown', systematic.BtagUp:'BtagUp', systematic.BtagDo:'BtagDown', systematic.MisTagUp:'MisTagUp', systematic.MisTagDo:'MisTagDown'}
class datasets():
  SingleMuon = 0
  SingleElec = 1
  DoubleMuon = 2
  DoubleElec = 3
  MuonEG     = 4
dataset = {datasets.SingleElec:'HighEGJet', datasets.SingleMuon:'SingleMuon', datasets.DoubleMuon:'DoubleMuon'}

class Dilepton(Nano5TeVHistoModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)
        if isMC:
            noSel = noSel.refine("mcWeight", weight=[ t.genWeight ])
        noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        plots = []

    def GetName(self, var, ichan, ilevel = '', isyst = ''):
       ''' Crafts the name for a histo '''
       if isinstance(ichan,  int): ichan  = chan[ichan]
       if isinstance(ilevel, int): ilevel = level[ilevel]
       if isinstance(isyst,  int): isyst  = systlabel[isyst]
       name = var + ('_' + ichan if ichan != '' else '') + ('_' + ilevel if ilevel != '' else '') + ('_'+isyst if isyst != '' else '')
       return name

    def NewHisto(self, var, chan, level, syst, nbins, bin0, binN):
        ''' Used to create the histos following a structure '''
        self.CreateTH1F(self.GetName(var, chan, level, syst), "", nbins, bin0, binN)

    def GetHisto(self, var, chan, level = '', syst = ''):
        ''' Get a given histo using the tthisto structure '''
        return self.obj[self.GetName(var, chan, level, syst)]

    def createHistos(self):
        ''' Create all the histos for the analysis '''
        self.isTTnom = True if self.outname      == 'TT' else False
        self.isTT    = True if self.outname[0:2] == 'TT' else False
        self.isDY    = True if 'DY' in self.sampleName      else False

        self.NewHisto('PUWeights', '', '', '', 100,0,5)
        self.NewHisto('MET_ElEl_2jetsnomet', '', '', '', 30,0,150)
        self.NewHisto('MET_MuMu_2jetsnomet', '', '', '', 30,0,150)

        ### Yields histos
        for key_chan in chan:
          ichan = chan[key_chan]
        if self.isTT: self.NewHisto('FiduEvents', ichan, '', '', 5,0,5)
        for key_syst in systlabel.keys():
          if not self.doSyst and key_syst != systematic.nom: continue
          isyst = systlabel[key_syst]
          self.NewHisto('Yields',   ichan, '', isyst, 5, 0, 5)
          self.NewHisto('YieldsSS', ichan, '', isyst, 5, 0, 5)


        ### Analysis histos
        for key_chan in chan:
          ichan = chan[key_chan]
          for key_level in level:
            ilevel = level[key_level]
            if ichan == 'ElMu' and (ilevel == 'ZVeto' or ilevel == 'MET'): continue
            # Create histos for PDF and scale systematics
            if self.isTTnom:
              self.NewHisto('PDFweights',ichan,ilevel,'',33,0.5,33.5)
              self.NewHisto('ScaleWeights',ichan,ilevel,'',9,0.5,9.5)
            for key_syst in systlabel.keys():
              if key_syst != systematic.nom and self.isData: continue
              if not self.doSyst and key_syst != systematic.nom: continue
              #if key_syst in [systematic.PUUp, systematic.PUDo, systematic.PrefireUp, systematic.PrefireDo, systematic.ISRUp, systematic.ISRDo, systematic.FSRUp, systematic.FSRDo, systematic.BtagDo, systematic.BtagUp]: continue
            
              isyst = systlabel[key_syst]
              # Event
              self.NewHisto('HT',   ichan,ilevel,isyst, 80, 0, 400)
              self.NewHisto('MET',  ichan,ilevel,isyst, 30, 0, 150)
              self.NewHisto('NJets',ichan,ilevel,isyst, 8 ,-0.5, 7.5)
              self.NewHisto('Btags',ichan,ilevel,isyst, 4 ,-0.5, 3.5)
              self.NewHisto('Vtx',  ichan,ilevel,isyst, 10, -0.5, 9.5)
              self.NewHisto('NBtagNJets', ichan,ilevel,isyst, 7, -0.5, 6.5)
              self.NewHisto('NBtagNJets_3bins', ichan,ilevel,isyst, 3, -0.5, 2.5)
              # Leptons
              self.NewHisto('Lep0Pt', ichan,ilevel,isyst, 20, 20, 120)
              self.NewHisto('Lep1Pt', ichan,ilevel,isyst, 16, 10, 90)
              self.NewHisto('Lep0Eta', ichan,ilevel,isyst, 50, -2.5, 2.5)
              self.NewHisto('Lep1Eta', ichan,ilevel,isyst, 50, -2.5, 2.5)
              self.NewHisto('Lep0Phi', ichan,ilevel,isyst, 20, -1, 1)
              self.NewHisto('Lep1Phi', ichan,ilevel,isyst, 20, -1, 1)
              self.NewHisto('DilepPt', ichan,ilevel,isyst, 40, 0, 200)
              self.NewHisto('InvMass', ichan,ilevel,isyst, 60, 0, 300)
              self.NewHisto('DYMass',  ichan,ilevel,isyst, 200, 70, 110)
              self.NewHisto('DYMassBB',  ichan,ilevel,isyst, 200, 70, 110)
              self.NewHisto('DYMassBE',  ichan,ilevel,isyst, 200, 70, 110)
              self.NewHisto('DYMassEB',  ichan,ilevel,isyst, 200, 70, 110)
              self.NewHisto('DYMassEE',  ichan,ilevel,isyst, 200, 70, 110)
              self.NewHisto('DeltaPhi',  ichan,ilevel,isyst, 20, 0, 1)
              if ichan == chan[ch.ElMu]:
                self.NewHisto('ElecEta', 'ElMu',ilevel,isyst, 50, -2.5, 2.5)
                self.NewHisto('MuonEta', 'ElMu',ilevel,isyst, 50, -2.5, 2.5)
                self.NewHisto('ElecPt', 'ElMu',ilevel,isyst, 20, 10, 110)
                self.NewHisto('MuonPt', 'ElMu',ilevel,isyst, 20, 10, 110)
                self.NewHisto('ElecPhi', 'ElMu',ilevel,isyst, 20, -1, 1)
                self.NewHisto('MuonPhi', 'ElMu',ilevel,isyst, 20, -1, 1)

    def FillHistograms(self, leptons, jets, pmet, ich, ilev, isys):
        ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
        if self.SS: return               # Do not fill histograms for same-sign events
        if not len(leptons) >= 2: return # Just in case
        self.SetWeight(isys)

        # Re-calculate the observables
        lep0  = leptons[0]; lep1 = leptons[1]
        l0pt  = lep0.Pt();  l1pt  = lep1.Pt()
        l0eta = lep0.Eta(); l1eta = lep1.Eta()
        l0phi = lep0.Phi(); l1phi = lep1.Phi()
        dphi  = DeltaPhi(lep0, lep1)
        mll   = InvMass(lep0, lep1)
        dipt  = DiPt(lep0, lep1)
        mupt  = 0; elpt  = 0
        mueta = 0; eleta = 0
        muphi = 0; elphi = 0
        if ich == ch.ElMu:
           if lep0.IsMuon():
              mu = lep0
              el = lep1
        else:
            mu = lep1
            el = lep0
        elpt  = el.Pt();  mupt  = mu.Pt()
        eleta = el.Eta(); mueta = mu.Eta()
        elphi = el.Phi(); muphi = mu.Phi()
        tWmvaVal1j1b=self.tWmvaVal1j1b      
        tWmvaVal2j1b=self.tWmvaVal2j1b               
        met = pmet.Pt()
        ht = 0; 

    
        ### Fill the histograms
        #if ich == ch.ElMu and ilev == lev.dilepton: print 'Syst = ', isys, ', weight = ', self.weight
        self.GetHisto('HT',   ich,ilev,isys).Fill(ht, self.weight)
        self.GetHisto('MET',  ich,ilev,isys).Fill(met, self.weight)
        self.GetHisto('NJets',ich,ilev,isys).Fill(njet, self.weight)
        self.GetHisto('Btags',ich,ilev,isys).Fill(nbtag, self.weight)
        self.GetHisto('Vtx',  ich,ilev,isys).Fill(self.nvtx, self.weight)
        self.GetHisto("InvMass", ich, ilev, isys).Fill(mll, self.weight)

    
        # Leptons
        self.GetHisto('Lep0Pt', ich,ilev,isys).Fill(l0pt, self.weight)
        self.GetHisto('Lep1Pt', ich,ilev,isys).Fill(l1pt, self.weight)
        self.GetHisto('Lep0Eta', ich,ilev,isys).Fill(l0eta, self.weight)
        self.GetHisto('Lep1Eta', ich,ilev,isys).Fill(l1eta, self.weight)
        self.GetHisto('Lep0Phi', ich,ilev,isys).Fill(l0phi/3.141592, self.weight)
        self.GetHisto('Lep1Phi', ich,ilev,isys).Fill(l1phi/3.141592, self.weight)
        self.GetHisto('DilepPt', ich,ilev,isys).Fill(dipt, self.weight)
        self.GetHisto('DeltaPhi',  ich,ilev,isys).Fill(dphi/3.141592, self.weight)
        self.GetHisto('InvMass', ich,ilev,isys).Fill(mll, self.weight)
        if mll > 70 and mll < 110: 
           self.GetHisto('DYMass',  ich,ilev,isys).Fill(mll, self.weight)
           l0eta = abs(l0eta); l1eta = abs(l1eta)
           if ich == ch.ElEl:
              if   l0eta <= 1.479 and l1eta <= 1.479: self.GetHisto('DYMassBB',  ich,ilev,isys).Fill(mll, self.weight)
              elif l0eta <= 1.479 and l1eta  > 1.479: self.GetHisto('DYMassBE',  ich,ilev,isys).Fill(mll, self.weight)
              elif l0eta  > 1.479 and l1eta <= 1.479: self.GetHisto('DYMassEB',  ich,ilev,isys).Fill(mll, self.weight)
              elif l0eta  > 1.479 and l1eta  > 1.479: self.GetHisto('DYMassEE',  ich,ilev,isys).Fill(mll, self.weight)
        if ich == ch.ElMu:
            self.GetHisto('ElecEta', ich,ilev,isys).Fill(eleta, self.weight)
            self.GetHisto('ElecPt',  ich,ilev,isys).Fill(elpt, self.weight)
            self.GetHisto('ElecPhi', ich,ilev,isys).Fill(elphi, self.weight)
            self.GetHisto('MuonEta', ich,ilev,isys).Fill(mueta, self.weight)
            self.GetHisto('MuonPt',  ich,ilev,isys).Fill(mupt, self.weight)
            self.GetHisto('MuonPhi', ich,ilev,isys).Fill(muphi, self.weight)
     
    def FillAll(self, ich, ilev, isyst, leps, jets, pmet):
        ''' Fill histograms for a given variation, channel and level '''
        #self.FillYieldsHistos(ich, ilev, isyst)
        #if isyst not in [systematic.PUUp, systematic.PUDo, systematic.PrefireUp, systematic.PrefireDo, systematic.ISRUp, systematic.ISRDo, systematic.FSRUp, systematic.FSRDo, systematic.BtagDo, systematic.BtagUp]:
        self.FillHistograms(leps, jets, pmet, ich, ilev, isyst)



    def SetVariables(self, isyst):
        leps = self.selLeptons
        jets = self.selJets
        pmet = self.pmet
        if   isyst == systematic.nom: return leps, jets, pmet
        elif isyst == systematic.JESUp:
           jets = self.selJetsJESUp
           pmet = self.pmetJESUp
        elif isyst == systematic.JESDo:
           jets = self.selJetsJESDo
           pmet = self.pmetJESDo
        elif isyst == systematic.JERUp:
           jets = self.selJetsJERUp
           pmet = self.pmetJERUp
        elif isyst == systematic.JERDo:
           jets = self.selJetsJERDo
           pmet = self.pmetJERDo
        return leps, jets, pmet



    #===================
    #   Muon Selection
    #===================
    def insideLoop(self, t):
        for i in range(t.nMuon):
          p = TLorentzVector()
          muonpt = t.Muon_pt[i] if not self.doMuonScaleCorrections else t.Muon_corrected_pt[i]
          p.SetPtEtaPhiM(muonpt, t.Muon_eta[i], t.Muon_phi[i], t.Muon_mass[i])
          charge = t.Muon_charge[i]
          dxy = abs(t.Muon_dxy[i])
          dz  = abs(t.Muon_dz[i] )
          passLepMVAID = True
          genPartFlav  = t.Muon_genPartFlav[i] if not self.isData else 1

          if not self.doLepMVA:
            # Tight ID
            if not t.Muon_tightId[i]              : continue #if not t.Muon_mediumId[i]: continue
            if not t.Muon_pfRelIso04_all[i] < 0.15: continue # Tight ISO, RelIso04 < 0.15
            if dxy > 0.2 or dz > 0.5              : continue # IP

          else:
            if not int(t.Muon_mediumPromptId[i])    : continue
            if t.Muon_sip3d[i] > 8                  : continue
            if dxy > 0.05 or dz > 0.1               : continue
            if t.Muon_jetIdx[i] >= 0:
              if     t.Jet_btagDeepB[t.Muon_jetIdx[i]] > 0.1522: passLepMVAID = False
              if not t.Jet_btagDeepB[t.Muon_jetIdx[i]] > -999.  : passLepMVAID = False
            if t.Muon_miniPFRelIso_all[i] > 0.325 : passLepMVAID = False
            if t.Muon_mvaTTH[i] < 0.55            : passLepMVAID = False
            #if t.Muon_miniIsoId[i] < 4: passLepMVAID = False

          # pT < 12 GeV, |eta| < 2.4
          if p.Pt() < self.LepPtCut or abs(p.Eta()) > 2.4: continue
          if (self.doLepMVA and passLepMVAID) or not self.doLepMVA:
            lep = lepton(p, charge, 13, genID = genPartFlav)
            setattr(lep, 'jetId', t.Muon_jetIdx[i])
            self.selLeptons.append(lep)
       
        ##### Electron
        for i in range(t.nElectron):
          p      = TLorentzVector()
          pt     = t.Electron_pt[i]
          eta    = t.Electron_eta[i]
          phi    = t.Electron_phi[i]
          m      = t.Electron_mass[i]
          ecorr  = t.Electron_eCorr[i]
          r9     = t.Electron_r9[i]
          if self.doElecScaleCorrections:
             run = 306936
             self.elecES.SetPrevCor(ecorr)
             if self.isData:
                pt, m = self.elecES.GetScaleCorr(run, pt, eta, phi, m, r9)
             else: 
                pt, m = self.elecES.GetSmearCorr(run, pt, eta, phi, m, r9)

          p.SetPtEtaPhiM(pt, eta, phi, m)
          charge   = t.Electron_charge[i]
          etaSC    = abs(p.Eta());
          dEtaSC   = t.Electron_deltaEtaSC[i]
          convVeto = t.Electron_convVeto[i]
          R9       = t.Electron_r9[i]
          dxy      = abs(t.Electron_dxy[i])
          dz       = abs(t.Electron_dz[i] )
          genPartFlav  = t.Electron_genPartFlav[i] if not self.isData else 1
          passLepMVAID =  True

          if not self.doLepMVA:
            if not t.Electron_cutBased[i] >= 3: continue # 4 Tightcut-based Id
            if ord(t.Electron_lostHits[i]) > 1: continue
            relIso03 = t.Electron_pfRelIso03_all[i]
            #if   etaSC <= 1.479 and relIso03 > 0.0361: continue
            #elif etaSC >  1.479 and relIso03 > 0.094:  continue
            if   etaSC <= 1.479 and (dxy > 0.05 or dz > 0.1): continue 
            elif etaSC >  1.479 and (dxy > 0.10 or dz > 0.2): continue 

          else: # LepMVA
            if ord(t.Electron_lostHits[i]) > 0: continue
            #if not(t.Electron_mvaFall17V2Iso_WPL[i]): continue
            if (t.Electron_sip3d[i] > 8):                   passLepMVAID = False
            if dxy > 0.05 or dz > 0.1:                   passLepMVAID = False
            if (t.Electron_mvaTTH[i] < 0.125):           passLepMVAID = False
            if (t.Electron_miniPFRelIso_all[i] > 0.085): passLepMVAID = False
            if t.Electron_jetIdx[i] >= 0:
               if t.Jet_btagDeepB[t.Electron_jetIdx[i]] > 0.1522: passLepMVAID = False
               if not t.Jet_btagDeepB[t.Electron_jetIdx[i]] > -999.: passLepMVAID = False

      # pT > 12 GeV, |eta| < 2.4
          if not convVeto: continue
          if p.Pt() < self.LepPtCut or abs(p.Eta()) > 2.5: continue
          if (self.doLepMVA and passLepMVAID) or (not self.doLepMVA):
            lep = lepton(p, charge, 11, genID = genPartFlav)
            setattr(lep, 'jetId', t.Electron_jetIdx[i])
            self.selLeptons.append(lep)

          leps = self.selLeptons
          pts  = [lep.Pt() for lep in leps]
          self.selLeptons = [lep for _,lep in sorted(zip(pts,leps))]
          self.selLeptons.reverse()
          leps = self.selLeptons

        ### Set dilepton channel
          nLep = len(self.selLeptons)
          if nLep < 2: return
          l0 = self.selLeptons[0]
          l1 = self.selLeptons[1]
          totId = l0.GetPDGid() + l1.GetPDGid()
          ich = -1
          if   totId == 24: ich = ch.ElMu
          elif totId == 22: ich = ch.ElEl
          elif totId == 26: ich = ch.MuMu
          ### Trigger
          ###########################################
          trigger = {
          ch.Elec:t.HLT_HIEle17_WPLoose_Gsf,
          ch.Muon:t.HLT_HIMu17, #HLT_HIL3Mu20
          ch.ElMu:t.HLT_HIMu17 or t.HLT_HIEle17_WPLoose_Gsf, #t.HLT_HIL3Mu20 or t.HLT_HIEle20_WPLoose_Gsf,
          ch.ElEl:t.HLT_HIEle17_WPLoose_Gsf,
          ch.MuMu:t.HLT_HIMu17# or t.HLT_HIL3DoubleMu0
          }

          passTrig = trigger[ich]

        ### Remove overlap events in datasets
        # In tt @ 5.02 TeV,
 
        if self.isData:
          if   self.sampleDataset == datasets.SingleElec:
            if   ich == ch.ElEl: passTrig = trigger[ch.ElEl]
            elif ich == ch.ElMu: passTrig = trigger[ch.ElEl] and not trigger[ch.MuMu]
            else:                passTrig = False
          elif self.sampleDataset == datasets.SingleMuon:
            if   ich == ch.MuMu: passTrig = trigger[ch.MuMu]
            elif ich == ch.ElMu: passTrig = trigger[ch.MuMu]
            else:                passTrig = False
          elif self.sampleDataset == datasets.DoubleMuon:
            if   ich == ch.MuMu: passTrig = trigger[ich]
            else:                passTrig = False

        ### Event selection
        ###########################################
        self.SetWeight(systematic.nom)
        weight = self.weight
    
        if not len(leps) >= 2: return
        l0 = leps[0]; l1 = leps[1]
        self.SS = l0.charge*l1.charge > 0
   
        if not passTrig: return
        for isyst in systlabel.keys():
          if not self.doSyst and isyst != systematic.nom: continue
          if self.isData and isyst != systematic.nom: continue
        leps, jets, pmet = self.SetVariables(isyst)
        nJets = len(jets)
        nBtag = self.GetNBtagJets(jets, isyst)

        ### Dilepton pair
        #if l0.Pt() < self.Lep0PtCut and l1.Pt() < self.Lep0PtCut: continue
        #if GetMinInvMass(leps) < 20: continue
        self.FillAll(ich, lev.dilepton, isyst, leps, jets, pmet)
        #if self.isTTnom and isyst == systematic.nom: self.FillLHEweights(t, ich, lev.dilepton)




