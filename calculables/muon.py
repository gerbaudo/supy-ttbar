from supy import wrappedChain,utils,calculables
import ROOT as r
##############################
class Indices(wrappedChain.calculable) :
    def __init__(self, collection = None, ptMin = None, isoMax = None, isoVar='ptcone30', isoType='relative', absEtaMax = 1000) :
        self.fixes = collection
        self.stash(["charge","author",isoVar,"loose","P4","cb_d0_exPV","cb_z0_exPV"])
        self.ptMin = ptMin
        self.absEtaMax = absEtaMax
        self.isoVar = isoVar
        self.isoMax = isoMax
        self.isoType = isoType
        self.isoRel = isoType=='relative'
        iso = isoVar + ('/pt' if self.isoRel else '')
        self.moreName = "muon pt>%.1f GeV; %s<%.2f"%(ptMin, iso, isoMax )

    def update(self,ignored) :
        self.value = []        
        p4s    = self.source[self.P4]
        isoVars = self.source[self.isoVar]
        iso = self.isoMax
        isoRel = self.isoRel
        isoMax = self.isoMax
        for i in range(p4s.size()) :
            p4 = p4s.at(i)
            isoVar = isoVars[i]
            if p4.pt() < self.ptMin : continue
            if self.absEtaMax < abs(p4.eta()) : continue
            if isoMax and ((isoVar/p4.pt()) > isoMax if isoRel else isoVar > isoMax) : continue
            self.value.append(i)
            # todo: implement dist pv etc.
##############################
