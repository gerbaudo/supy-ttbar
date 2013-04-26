from supy import wrappedChain,utils,calculables
import ROOT as r

# list of leaf names
# ['mc_channel_number', 'mc_event_number', 'mc_event_weight', 'mc_n',
# 'mc_pt', 'mc_m', 'mc_eta', 'mc_phi', 'mc_status', 'mc_barcode',
# 'mc_pdgId', 'mc_charge', 'mc_parents', 'mc_children', 'mc_vx_x',
# 'mc_vx_y', 'mc_vx_z', 'mc_vx_barcode', 'mc_child_index',
# 'mc_parent_index']
#

defaultColl=('mc_','')
#___________________________________________________________
class genP4(wrappedChain.calculable) :
    @property
    def name(self) : return "genP4"
    def __init__(self, collection=defaultColl):
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

    def __init__(self, collection=defaultColl,
                 pdgs = [], label = None, status = [], parentPdgs = [],
                 parentIndexLabel='',maxLen=None) :
        self.fixes = collection
        self.stash(['pdgId','status','parent_index'])
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
        pdg = self.source[self.pdgId]
        status = self.source[self.status]
        parents = self.source[self.parent_index]
        parentsPdgs = [[pdg[p] for p in par] for par in parents]
        parentsIndices = self.source[self.parentIndexLabel] if self.parentIndexLabel else []
        print 'trying to update genIndices ',self.label
        print 'input : %d elem'%len(pdg)
        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.status) or (status.at(i) in self.status) ) and \
                                 ( (not self.parentPdgs) or any(p in self.parentPdgs for p in parentsPdgs[i])) and \
                                 ( (not self.parentIndexLabel) or any(p in parents[i] for p in parentsIndices))
                                 ,
                             range(pdg.size()) )
        print "found %d"%len(self.value)
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

class sherpaTtbarProductsIndices(wrappedChain.calculable) :
    """ The MC record for a ttbar generated with sherpa looks like
    this, with all the decay products coming from two gluons.All the
    interesting particles have status 3, all the other ones have 11
    (mostly).
                 t -> lpos v b, t_ -> lneg b_ v_

    --------------------------------------
    i    st   par         id          name
    --------------------------------------
    21   3 [...]            21         g
    22   3 [...]            21         g
    23   3 [21, 22]         12        ve
    24   3 [21, 22]        -11        e+
    25   3 [21, 22]          5         b
    26   3 [21, 22]         11        e-
    27   3 [21, 22]        -12        ve
    28   3 [21, 22]         -5        /b
    """
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId','status','parent_index'])
    def update(self,_) :
        pdg = self.source[self.pdgId]
        status = self.source[self.status]
        parents = self.source[self.parent_index]
        self.value = filter( lambda i: status[i]==3 and len(parents[i])==2, range(pdg.size()) )
        #pdg = self.source['mc_pdgId']
        #parents = self.source['mc_parent_index']
        #parentsPdgs = [[pdg[p] for p in par] for par in parents]
        #parentsIndices = self.source[self.parentIndexLabel] if self.parentIndexLabel else []

class genIndicesSherpa(wrappedChain.calculable) :
    def __init__(self, pdgs = []) :
        self.pdgs = pdgs
    def update(self,_) :
        pdg = self.source['mc_pdgId']
        sherpaIndices = self.source['sherpaTtbarProductsIndices']
        self.value = filter( lambda i: pdg[i] in self.pdgs, sherpaIndices )

class genIndiceslpos(genIndicesSherpa) :
    def __init__(self) : self.pdgs = frozenset([-11, -13])
class genIndiceslneg(genIndicesSherpa) :
    def __init__(self) : self.pdgs = frozenset([+11, +13])
class genIndicesb(genIndicesSherpa) :
    def __init__(self) : self.pdgs = frozenset([+5])
class genIndicesbbar(genIndicesSherpa) :
    def __init__(self) : self.pdgs = frozenset([-5])
class genIndicesv(genIndicesSherpa) :
    def __init__(self) : self.pdgs = frozenset([+12,+14])
class genIndicesvbar(genIndicesSherpa) :
    def __init__(self) : self.pdgs = frozenset([-12,-14])
