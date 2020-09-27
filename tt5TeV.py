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

class chan():
    ElEl =1
    MuMu =2
    ElMu =3
channel = {chan.ElEl:'ElEl', chan.MuMu:'MuMu',chan.ElMu:'ElMu'}

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
    
    def GetName(self, var, ichan):
    ''' Crafts the name for a histo '''
        if isinstance(ichan,  int): ichan  = channel[ichan]
        name = var + ('_' + ichan if ichan != '' else '')
        return name

    def NewHisto(self, var, channnel, nbins, bin0, binN):
    ''' Used to create the histos following a structure '''
        self.CreateTH1F(self.GetName(var, channel), "", nbins, bin0, binN)

    def GetHisto(self, var, channel):
    ''' Get a given histo using the tthisto structure '''
        return self.obj[self.GetName(var, channel)]

    def createHistos(self):
    ''' Create all the histos for the analysis '''

        ### Analysis histos
        for key_chan in channel:
        ichan = channnel[key_chan]

        # Event
        self.NewHisto('HT',   ichan, 80, 0, 400)
        self.NewHisto('MET',  ichan, 30, 0, 150)
        self.NewHisto('NJets',ichan, 8 ,-0.5, 7.5)
        self.NewHisto('Btags',ichan, 4 ,-0.5, 3.5)
        self.NewHisto('Vtx',  ichan, 10, -0.5, 9.5)
        self.NewHisto('NBtagNJets', ichan, 7, -0.5, 6.5)
        self.NewHisto('NBtagNJets_3bins', ichan, 3, -0.5, 2.5)
        # Leptons
        self.NewHisto('Lep0Pt', ichan, 20, 20, 120)
        self.NewHisto('Lep1Pt', ichan, 16, 10, 90)
        self.NewHisto('Lep0Eta', ichan, 50, -2.5, 2.5)
        self.NewHisto('Lep1Eta', ichan, 50, -2.5, 2.5)
        self.NewHisto('Lep0Phi', ichan, 20, -1, 1)
        self.NewHisto('Lep1Phi', ichan, 20, -1, 1)
        self.NewHisto('DilepPt', ichan, 40, 0, 200)
        self.NewHisto('InvMass', ichan, 60, 0, 300)
        self.NewHisto('DYMass',  ichan, 200, 70, 110)
        self.NewHisto('DYMassBB',  ichan, 200, 70, 110)
        self.NewHisto('DYMassBE',  ichan, 200, 70, 110)
        self.NewHisto('DYMassEB',  ichan, 200, 70, 110)
        self.NewHisto('DYMassEE',  ichan, 200, 70, 110)
        self.NewHisto('DeltaPhi',  ichan, 20, 0, 1)
        if ichan == chan[ch.ElMu]:
           self.NewHisto('ElecEta', 'ElMu', 50, -2.5, 2.5)
           self.NewHisto('MuonEta', 'ElMu', 50, -2.5, 2.5)
           self.NewHisto('ElecPt', 'ElMu', 20, 10, 110)
           self.NewHisto('MuonPt', 'ElMu', 20, 10, 110)
           self.NewHisto('ElecPhi', 'ElMu', 20, -1, 1)
           self.NewHisto('MuonPhi', 'ElMu', 20, -1, 1)

    def FillHistograms(self, leptons, jets, pmet, ich, ilev, isys):
    ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
        if self.SS: return               # Do not fill histograms for same-sign events
        if not len(leptons) >= 2: return # Just in case
    
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
        if ich == channel.ElMu:
           if lep0.IsMuon():
           mu = lep0
           el = lep1
        else:
           mu = lep1
           el = lep0
        elpt  = el.Pt();  mupt  = mu.Pt()
        eleta = el.Eta(); mueta = mu.Eta()
        elphi = el.Phi(); muphi = mu.Phi()
          
        met = pmet.Pt()
        ht = 0; 
   
    
        ### Fill the histograms
        #if ich == ch.ElMu and ilev == lev.dilepton: print 'Syst = ', isys, ', weight = ', self.weight
        self.GetHisto('HT',   ich).Fill(ht)
        self.GetHisto('MET',  ich).Fill(met)
        self.GetHisto('NJets',ich).Fill(njet)
        self.GetHisto('Btags',ich).Fill(nbtag)
        self.GetHisto('Vtx',  ich).Fill(self.nvtx)
        self.GetHisto("InvMass", ich, ).Fill(mll)

    
       # Leptons
       self.GetHisto('Lep0Pt', ich).Fill(l0pt)
       self.GetHisto('Lep1Pt', ich).Fill(l1pt)
       self.GetHisto('Lep0Eta', ich).Fill(l0eta)
       self.GetHisto('Lep1Eta', ich).Fill(l1eta)
       self.GetHisto('Lep0Phi', ich).Fill(l0phi/3.141592)
       self.GetHisto('Lep1Phi', ich).Fill(l1phi/3.141592)
       self.GetHisto('DilepPt', ich).Fill(dipt, self.weight)
       self.GetHisto('DeltaPhi',  ich).Fill(dphi/3.141592)
       self.GetHisto('InvMass', ich).Fill(mll, self.weight)
       if mll > 70 and mll < 110: 
          self.GetHisto('DYMass',  ich).Fill(mll)
       l0eta = abs(l0eta); l1eta = abs(l1eta)
       if ich == chan.ElEl:
          if   l0eta <= 1.479 and l1eta <= 1.479: self.GetHisto('DYMassBB',  ich).Fill(mll)
          elif l0eta <= 1.479 and l1eta  > 1.479: self.GetHisto('DYMassBE',  ich).Fill(mll)
          elif l0eta  > 1.479 and l1eta <= 1.479: self.GetHisto('DYMassEB',  ich).Fill(mll)
          elif l0eta  > 1.479 and l1eta  > 1.479: self.GetHisto('DYMassEE',  ich).Fill(mll)
      if ich == chan.ElMu:
         self.GetHisto('ElecEta', ich).Fill(eleta)
         self.GetHisto('ElecPt',  ich).Fill(elpt)
         self.GetHisto('ElecPhi', ich).Fill(elphi)
         self.GetHisto('MuonEta', ich).Fill(mueta)
         self.GetHisto('MuonPt',  ich).Fill(mupt)
         self.GetHisto('MuonPhi', ich).Fill(muphi)


        muons = op.select(t.Muon, lambda mu : mu.pt > 20.)
        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])
        plots.append(Plot.make1D("dimu_M",
            op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel, EqB(100, 20., 120.),
            title="Dimuon invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        electrons = op.select(t.Electron, lambda el : el.pt > 22.)
        twoElSel = noSel.refine("twoElectrons", cut=[ op.rng_len(electrons) > 1 ])
        plots.append(Plot.make1D("diel_M",
            op.invariant_mass(electrons[0].p4, electrons[1].p4), twoElSel, EqB(100, 20., 120.),
            title="Dielectron invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))



        return plots

