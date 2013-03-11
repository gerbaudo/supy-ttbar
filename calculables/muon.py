from supy import wrappedChain,utils,calculables
import ROOT as r
##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, isoMax = None, isoVar = "ptcone30", absEtaMax = 1000) :
        self.fixes = collection
        self.stash(["charge","author","ptcone30","loose","P4","cb_d0_exPV","cb_z0_exPV"])
        self.ptMin = ptMin
        self.absEtaMax = absEtaMax
        self.isoMax = isoMax
        self.moreName = "muon pt>%.1f GeV; %s<%.2f"%(ptMin, ISO, isoMax )

    def update(self,ignored) :
        self.value = []        
        p4s    = self.source[self.P4]
        iso = self.source[self.ISO]
        isGlobal = self.source[self.IsGlobalMuon]
        for i in range(p4s.size()) :
            p4 = p4s.at(i)
            if p4.pt() < self.ptMin : continue
            if self.absEtaMax < abs(p4.eta()) : continue
            self.value.append(i)
            # todo: implement isolation etc.
##############################
