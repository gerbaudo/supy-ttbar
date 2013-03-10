from supy import wrappedChain,utils,calculables
import ROOT as r

#___________________________________________________________
class genP4(wrappedChain.calculable) :
    @property
    def name(self) : return "genP4"
    def __init__(self, collection=('mc_','')):
        self.p4=utils.LorentzV
        self.fixes = collection
        self.stash(["E", "pt", "eta", "phi", "m"])
    def update(self, _) :
        self.value = [self.p4(pt, eta, phi, m) for pt,eta,phi,m in zip(self.source[self.pt],
                                                                       self.source[self.eta],
                                                                       self.source[self.phi],
                                                                       self.source[self.m])]
#___________________________________________________________
class genIndices(wrappedChain.calculable) :
    @property
    def name(self) : return "genIndices" + self.label

    def __init__(self, pdgs = [], label = None, status = [], parentPdgs = [], parentIndexLabel='',maxLen=None) :
        self.label = label
        self.PDGs = frozenset(pdgs)
        self.status = frozenset(status)
        self.parentPdgs = frozenset(parentPdgs)
        self.parentIndexLabel = parentIndexLabel
        self.maxLen = maxLen
        self.moreName = "; ".join(["pdgId in %s" %str(list(self.PDGs)),
                                   "status in %s"%str(list(self.status)),
                                   "parentPdgs in %s"%str(list(self.parentPdgs)),
                                   "parentIndexLabel in %s"%self.parentIndexLabel,
                                   "maxLen %s"%str(self.maxLen),
                                   ])
    def update(self,_) :
        pdg = self.source['mc_pdgId']
        status = self.source['mc_status']
        parents = self.source['mc_parent_index']
        parentsPdgs = [[pdg[p] for p in par] for par in parents]
        parentsIndices = self.source[self.parentIndexLabel] if self.parentIndexLabel else []

        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.status) or (status.at(i) in self.status) ) and \
                                 ( (not self.parentPdgs) or any(p in self.parentPdgs for p in parentsPdgs[i])) and \
                                 ( (not self.parentIndexLabel) or any(p in parents[i] for p in parentsIndices))
                                 ,
                             range(pdg.size()) )
        if self.maxLen : self.value = self.value[:self.maxLen] # trim if needed
#___________________________________________________________
class smTopIndex(wrappedChain.calculable) :
    def __init__(self) :
        self.PDGs = frozenset([-6,+6])
        self.childrenPdgs = frozenset([-5,5,-24,24])
        self.moreName = "; ".join(["pdgId in %s" %str(list(self.PDGs)),
                                   "childrenPdgs in %s"%str(list(self.childrenPdgs))])
    def update(self,_) :
        pdg = self.source['mc_pdgId']
        childrenIndices = self.source['mc_child_index']
        childrenPdgs = [[pdg[c] for c in chi] for chi in childrenIndices]
        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.childrenPdgs) or frozenset(childrenPdgs[i]).issubset(self.childrenPdgs)),
                             range(pdg.size()) )
#___________________________________________________________
class genIndicesWqq(wrappedChain.calculable) :
    def update(self,_) :
        ids = self.source['genPdgId']
        mom = self.source['genMotherPdgId']
        self.value = filter(lambda i: abs(mom[i]) is 24 and abs(ids[i]) < 5, range(len(ids)))
