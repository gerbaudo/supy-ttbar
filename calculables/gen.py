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
        self.label = label
        self.fixes = collection
        self.stash(['pdgId','status','parent_index'])
        self.PDGs = frozenset(pdgs)
        self.statuses = frozenset(status)
        self.parentPdgs = frozenset(parentPdgs)
        self.parentIndexLabel = parentIndexLabel
        self.maxLen = maxLen
        self.moreName = "; ".join(filter(lambda x:x,
                                         ["pdgId in %s" %str(list(self.PDGs)) if self.PDGs else '',
                                          "status in %s"%str(list(self.statuses)) if self.statuses else '',
                                          "parentPdgs in %s"%str(list(self.parentPdgs)) if self.parentPdgs else '',
                                          "parentIndexLabel in %s"%self.parentIndexLabel if self.parentIndexLabel else '',
                                          "maxLen %s"%str(self.maxLen) if self.maxLen else '',
                                          ])
                                  )
    def update(self,_) :
        pdg = self.source[self.pdgId]
        status = self.source[self.status]
        parents = self.source[self.parent_index]
        parentsPdgs = [[pdg[p] for p in par] for par in parents]
        parentsIndices = self.source[self.parentIndexLabel] if self.parentIndexLabel else []
        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.statuses) or (status.at(i) in self.statuses) ) and \
                                 ( (not self.parentPdgs) or any(p in self.parentPdgs for p in parentsPdgs[i])) and \
                                 ( (not self.parentIndexLabel) or any(p in parents[i] for p in parentsIndices))
                                 ,
                             range(pdg.size()) )
        if self.maxLen : self.value = self.value[:self.maxLen] # trim if needed
#___________________________________________________________
class smTopIndex(wrappedChain.calculable) :
    def __init__(self, collection=defaultColl):
        self.fixes = collection
        self.stash(['pdgId', 'child_index'])
        self.PDGs = frozenset([-6,+6])
        self.childrenPdgs = frozenset([-5,5,-24,24])
        self.moreName = "; ".join(["pdgId in %s" %str(list(self.PDGs)),
                                   "childrenPdgs in %s"%str(list(self.childrenPdgs))])
    def update(self,_) :
        pdg = self.source[self.pdgId]
        childrenIndices = self.source[self.child_index]
        childrenPdgs = [[pdg[c] for c in chi] for chi in childrenIndices]
        self.value = filter( lambda i: ( (not self.PDGs) or (pdg.at(i) in self.PDGs) ) and \
                                 ( (not self.childrenPdgs) or frozenset(childrenPdgs[i]).issubset(self.childrenPdgs)),
                             range(pdg.size()) )
#___________________________________________________________
class genIndicesWqq(wrappedChain.calculable) :
    def __init__(self, collection=defaultColl):
        self.fixes = collection
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
    def __init__(self, collection=defaultColl, pdgs=[]) :
        self.fixes = collection
        self.stash(['pdgId', 'sherpaTtbarProductsIndices'])
        self.pdgs = frozenset(pdgs)
    def update(self,_) :
        pdg = self.source[self.pdgId]
        sherpaIndices = self.source[self.sherpaTtbarProductsIndices]
        self.value = filter( lambda i: pdg[i] in self.pdgs, sherpaIndices )

class genIndiceslpos(genIndicesSherpa) :
    def __init__(self, collection) : super(genIndiceslpos, self).__init__(collection, [-11, -13])
class genIndiceslneg(genIndicesSherpa) :
    def __init__(self, collection) : super(genIndiceslneg, self).__init__(collection, [+11, +13])
class genIndicesb(genIndicesSherpa) :
    def __init__(self, collection) : super(genIndicesb,    self).__init__(collection, [+5])
class genIndicesbbar(genIndicesSherpa) :
    def __init__(self, collection) : super(genIndicesbbar, self).__init__(collection, [-5])
class genIndicesv(genIndicesSherpa) :
    def __init__(self, collection) : super(genIndicesv,    self).__init__(collection, [+12,+14])
class genIndicesvbar(genIndicesSherpa) :
    def __init__(self, collection) : super(genIndicesvbar, self).__init__(collection, [-12,-14])
#___________________________________________________________
class higgsIndices(wrappedChain.calculable) :
    @property
    def name(self) : return 'higgsIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId','status','parent_index','child_index'])
        self.moreName = 'higgsIndices'
    def update(self, _) :
        pdg        = self.source[self.pdgId]
        status     = self.source[self.status]
        parents    = self.source[self.parent_index]
        children   = self.source[self.child_index]
        higgsPdg = 25
        iHiggs = [i for i,p in enumerate(pdg) if p==higgsPdg]
        higgsChildren = [[pdg[i] if i<len(pdg) else 0 for i in hhc] for hhc in [children[i] for i in iHiggs] ]
        higgsParents = [[pdg[i] if i<len(pdg) else 0 for i in hhp] for hhp in [parents[i] for i in iHiggs] ]
        # there can be several intermediate higgs
        boringIhiggs = [ih for i, ih in enumerate(iHiggs) if len(higgsChildren[i])<2 or higgsPdg in higgsChildren[i]]
        interestingIhiggs = [i for i in iHiggs if i not in boringIhiggs]
        higgsChildren = [hc for hc,ih in zip(higgsChildren, iHiggs) if ih in interestingIhiggs]
        higgsParents = [hp for hp,ih in zip(higgsParents, iHiggs) if ih in interestingIhiggs]
        self.value = interestingIhiggs
#___________________________________________________________
class higgsChildrenIndices(wrappedChain.calculable) :
    @property
    def name(self) : return 'higgsChildrenIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['higgsIndices','child_index'])
    def update(self, _) :
        higgsIndices = self.source[self.higgsIndices]
        children   = self.source[self.child_index]
        self.value = [i for chs in [children[h] for h in higgsIndices] for i in chs]
#___________________________________________________________
class higgsParentsIndices(wrappedChain.calculable) :
    @property
    def name(self) : return 'higgsParentsIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['higgsIndices','parent_index'])
    def update(self, _) :
        higgsIndices = self.source[self.higgsIndices]
        parents   = self.source[self.parent_index]
        self.value = [i for pas in [parents[h] for h in higgsIndices] for i in pas]
#___________________________________________________________
class higgsDecayType(wrappedChain.calculable) :
    @property
    def name(self) : return 'higgsDecayType'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId','higgsChildrenIndices'])
    def update(self, _) :
        children = self.source[self.higgsChildrenIndices]
        pdgs     = self.source[self.pdgId]
        ch = [pdgs[i] for i in children]
        decay = None
        if   any(p in ch for p in [-24, +24])  : decay = 0 # 'WW'
        elif any(p in ch for p in [23])        : decay = 1 # 'ZZ'
        elif any(p in ch for p in [-15, +15])  : decay = 2 # 'tautau'
        elif any(p in ch for p in [-5, +5])    : decay = 3 # 'bbbar'
        elif any(p in ch for p in [-13, +13])  : decay = 4 # 'mumu'
        else                                   : decay = -1# 'unknown'
        self.value = decay
#___________________________________________________________
class higgsIsDecayType(wrappedChain.calculable) :
    def __init__(self, decayType, collection=defaultColl) :
        self.fixes = collection
        self.key = 'higgsDecayType'.join(collection)
        self.stash([self.key])
        self.decayType = decayType # see higgsDecayType above
    def update(self, _) : self.value = self.source[self.key]==self.decayType
#___________________________________________________________
class higgsIsWw(higgsIsDecayType) :
    @property
    def name(self) : return 'higgsIsWw'.join(self.fixes)
    def __init__(self, collection) : super(higgsIsWw, self).__init__(0, collection)
#___________________________________________________________
class higgsIsTauTau(higgsIsDecayType) :
    @property
    def name(self) : return 'higgsIsTauTau'.join(self.fixes)
    def __init__(self, collection) : super(higgsIsTauTau, self).__init__(2, collection)
#___________________________________________________________
class wIndices(wrappedChain.calculable) :
    "Index of the W that comes from the chargino decay"
    @property
    def name(self) : return 'wIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId','parent_index'])
    def update(self,_) :
        ws        = frozenset([+24, -24])
        charginos = frozenset([+1000024, -1000024])
        pdg     = self.source[self.pdgId]
        parents = self.source[self.parent_index]
        self.value = filter( lambda i: pdg[i] in ws and len(parents[i])==1 and
                             all([pdg[p] in charginos for p in parents[i]]),
                             range(pdg.size()) )
#___________________________________________________________
def extractXchildrenIndices(xPdgs=[], xIndex=None, pdgs=[], childrenIndices=[[]]) :
    """look for the children's indices of particle X; sometimes there
    are intermediate Xs...need to follow the thread"""
    xPdgs     = frozenset(xPdgs)
    xI        = xIndex
    children  = childrenIndices
    childrenX = children[xI]
    while len(childrenX) and any([pdgs[c] in xPdgs for c in childrenX]) :
        xI = filter(lambda i: pdgs[i] in xPdgs, childrenX)[0]
        childrenX = [i for i in children[xI]]
    return childrenX
def extractWchildrenIndices(wIndex, pdgs, childrenIndices) :
    return extractXchildrenIndices(xPdgs=[+24, -24], xIndex=wIndex, pdgs=pdgs, childrenIndices=childrenIndices)
def extractTauchildrenIndices(tauIndex, pdgs, childrenIndices) :
    return extractXchildrenIndices(xPdgs=[+15, -15], xIndex=tauIndex, pdgs=pdgs, childrenIndices=childrenIndices)
#___________________________________________________________
class wChildrenIndices(wrappedChain.calculable) :
    @property
    def name(self) : return 'wChildrenIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId', 'child_index','wIndices'])
    def update(self,_) :
        ws        = frozenset([+24, -24])
        wI        = self.source[self.wIndices]
        assert len(wI)==1,"don't know how to deal with multiple W"
        wI        = wI[0]
        self.value = extractWchildrenIndices(wI, self.source[self.pdgId],
                                             self.source[self.child_index])
#___________________________________________________________
def extractWdecayType(pdgChildren=[]):
    pdgs = frozenset(pdgChildren)
    el,ve, mu, vm, ta, vt = 11, 12, 13, 14, 15, 16
    d, u, s, c, b = 1, 2, 3, 4, 5
    lv, qqbar, unknown = 0, 1, -1
    if (frozenset([-el, +ve]).issubset(pdgs) or  frozenset([+el, -ve]).issubset(pdgs) or
        frozenset([-mu, +vm]).issubset(pdgs) or  frozenset([+mu, -vm]).issubset(pdgs) or
        frozenset([-ta, +vt]).issubset(pdgs) or  frozenset([+ta, -vt]).issubset(pdgs)) :
        return lv
    elif (frozenset([-u, +d]).issubset(pdgs) or  frozenset([+u, -d]).issubset(pdgs) or
          frozenset([-u, +s]).issubset(pdgs) or  frozenset([+u, -s]).issubset(pdgs) or
          frozenset([-u, +b]).issubset(pdgs) or  frozenset([+u, -b]).issubset(pdgs) or
          frozenset([-c, +d]).issubset(pdgs) or  frozenset([+c, -d]).issubset(pdgs) or
          frozenset([-c, +s]).issubset(pdgs) or  frozenset([+c, -s]).issubset(pdgs) or
          frozenset([-c, +b]).issubset(pdgs) or  frozenset([+c, -b]).issubset(pdgs)) :
        return qqbar
    else :
        return unknown
#___________________________________________________________
def extractTaudecayType(pdgChildren=[]):
    pdgs = frozenset(pdgChildren)
    el,ve, mu, vm, ta, vt = 11, 12, 13, 14, 15, 16
    leptonic, hadronic, unknown = 0, 1, -1
    if (frozenset([-el, +ve]).issubset(pdgs) or  frozenset([+el, -ve]).issubset(pdgs) or
        frozenset([-mu, +vm]).issubset(pdgs) or  frozenset([+mu, -vm]).issubset(pdgs) or
        frozenset([-ta, +vt]).issubset(pdgs) or  frozenset([+ta, -vt]).issubset(pdgs)) :
        return leptonic
    else : return hadronic # we could try to distinguish b/w had and unknown...not now
#___________________________________________________________
class wDecayType(wrappedChain.calculable) :
    @property
    def name(self) : return 'wDecayType'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId', 'wChildrenIndices'])
    def update(self, _) :
        children = self.source[self.wChildrenIndices]
        pdgs     = self.source[self.pdgId]
        pdgs     = [pdgs[i] for i in children]
        self.value = extractWdecayType(pdgs)
#___________________________________________________________
class wIsLeptonic(wrappedChain.calculable) :
    @property
    def name(self) : return 'wIsLeptonic'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['wDecayType'])
    def update(self, _) :
        self.value = self.source[self.wDecayType]==0
#___________________________________________________________
class wIsHadronic(wrappedChain.calculable) :
    @property
    def name(self) : return 'wIsHadronic'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['wDecayType'])
    def update(self, _) :
        self.value = self.source[self.wDecayType]==1
#___________________________________________________________
class HtautauTauIndices(wrappedChain.calculable) :
    "Indices of the taus from H->tautau"
    @property
    def name(self) : return 'tauIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId', 'higgsChildrenIndices'])
    def update(self, _) :
        pdgs             = self.source[self.pdgId]
        hChildrenIndices = self.source[self.higgsChildrenIndices]
        tauPgds          = [-15, +15]
        self.value = [i for i in hChildrenIndices if pdgs[i] in tauPgds]
#___________________________________________________________
class tausDecayType(wrappedChain.calculable) :
    "determine the decay type for given taus"
    @property
    def name(self) : return 'tausDecayType'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['pdgId','child_index','tauIndices'])
    def update(self, _) :
        children     = self.source[self.child_index]
        pdgs         = self.source[self.pdgId]
        tauIndices   = self.source[self.tauIndices]
        tauPgds      = [-15, +15]
        assert all([pdgs[t] in tauPgds for t in tauIndices]),"not taus : %s"%str([pdgs[t] for t in tauIndices])
        tausChildren = [extractTauchildrenIndices(tI, pdgs, children) for tI in tauIndices]
        pdgsChildren = [[pdgs[i]  for i in tc] for tc in tausChildren]
        self.value = [extractTaudecayType(pdgs) for pdgs in pdgsChildren]
#___________________________________________________________
class HwwWchildrenIndices(wrappedChain.calculable) :
    "Indices of the children of the W from H->WW"
    @property
    def name(self) : return ("HwwW%dchildrenIndices"%self.wIndex()).join(self.fixes)
    def __init__(self, collection=defaultColl, firstW=False, secondW=False) :
        assert firstW!=secondW,"HwwWchildrenIndices: pick either first or second"
        for attr in ['firstW', 'secondW'] : setattr(self, attr, eval(attr))
        self.fixes = collection
        self.stash(['pdgId','child_index','higgsChildrenIndices','higgsIsWw'])
    def wIndex(self) : return 0 if self.firstW else 1
    def update(self, _) :
        pdgs      = self.source[self.pdgId]
        isWw      = self.source[self.higgsIsWw]
        children  = self.source[self.child_index]
        hChildren = self.source[self.higgsChildrenIndices]
        chPdgs = [pdgs[i] for i in hChildren]
        self.value = []
        Ws = [+24,-24]
        if not isWw : return
        wIndices = [i for i,p in zip(hChildren, chPdgs) if p in Ws]
        nonWindices = [i for i in hChildren if i not in wIndices] # one can either have the W or its decay products (lv)
        if not (len(nonWindices)/2 + len(wIndices))==2 : return # don't know how to interpret
        nWs = len(wIndices)
        if   nWs==2 :
            w0Indices = extractWchildrenIndices(wIndices[0], pdgs, children)
            w1Indices = extractWchildrenIndices(wIndices[1], pdgs, children)
        elif nWs==1 :
            w0Indices = extractWchildrenIndices(wIndices[0], pdgs, children)
            w1Indices = nonWindices
        else :
            w0Indices = nonWindices[:2]
            w1Indices = nonWindices[2:]
        self.value = w0Indices if self.firstW else w1Indices
class HwwW0childrenIndices(HwwWchildrenIndices) :
    def __init__(self, collection=defaultColl) :
        super(HwwW0childrenIndices, self).__init__(collection, firstW=True)
class HwwW1childrenIndices(HwwWchildrenIndices) :
    def __init__(self, collection=defaultColl) :
        super(HwwW1childrenIndices, self).__init__(collection, secondW=True)
#___________________________________________________________
class HwwWdecayType(wrappedChain.calculable) :
    "Decay type of the W from H->WW"
    @property
    def name(self) : return ("HwwW%ddecayType"%self.wIndex()).join(self.fixes)
    def __init__(self, collection=defaultColl, firstW=False, secondW=False) :
        assert firstW!=secondW,"HwwWdecayType: pick either first or second"
        for attr in ['firstW', 'secondW'] : setattr(self, attr, eval(attr))
        self.fixes = collection
        self.stash(['pdgId','HwwW0childrenIndices','HwwW1childrenIndices'])
    def wIndex(self) : return 0 if self.firstW else 1
    def update(self, _) :
        pdgs       = self.source[self.pdgId]
        childrenW0 = self.source[self.HwwW0childrenIndices]
        childrenW1 = self.source[self.HwwW1childrenIndices]
        childrenIndices = childrenW0 if self.firstW else childrenW1
        childrenPgds = [pdgs[i] for i in childrenIndices]
        self.value = extractWdecayType(childrenPgds)
class HwwW0decayType(HwwWdecayType) :
    def __init__(self, collection=defaultColl) :
        super(HwwW0decayType, self).__init__(collection, firstW=True)
class HwwW1decayType(HwwWdecayType) :
    def __init__(self, collection=defaultColl) :
        super(HwwW1decayType, self).__init__(collection, secondW=True)
#___________________________________________________________
class HwwHadronicIndices(wrappedChain.calculable) :
    "Indices of the quarks from the hadronic W's from H->WW"
    @property
    def name(self) : return 'HwwHadronicIndices'.join(self.fixes)
    def __init__(self, collection=defaultColl) :
        self.fixes = collection
        self.stash(['HwwW0childrenIndices', 'HwwW0decayType',
                    'HwwW1childrenIndices', 'HwwW1decayType'])
    def update(self, _) :
        w0chI = self.source[self.HwwW0childrenIndices]
        w1chI = self.source[self.HwwW1childrenIndices]
        w0dec = self.source[self.HwwW0decayType]
        w1dec = self.source[self.HwwW1decayType]
        def isWhad(decay) : return 1==decay
        w0isHad, w1isHad = isWhad(w0dec), isWhad(w1dec)
        self.value = (w0chI if isWhad(w0dec) else []) + (w1chI if isWhad(w1dec) else [])
#___________________________________________________________
class UnmatchedGenIndices(wrappedChain.calculable) :
    "Unmatched truth particles (generally partons without a matching reco jet"
    def __init__(self, genIndices='', otherP4Coll='', otherIndices='') :
        self.genIndices   = genIndices
        self.otherP4Coll  = otherP4Coll
        self.otherIndices = otherIndices
        self.maxDr = 0.4
        self.moreName = "%s w/out %.1f match %s"%(genIndices, self.maxDr, otherIndices)
    def update(self, _) :
        gP4s = self.source['genP4']
        gIds = self.source[self.genIndices]
        rP4s = self.source[self.otherP4Coll]
        rIds = self.source[self.otherIndices]
        unmatched = []
        for ig in gIds :
            jg = gP4s[ig]
            match = None
            for ir in rIds :
                jr =  rP4s[ir]
                if r.Math.VectorUtil.DeltaR(jg, jr) < self.maxDr :
                    match = True
                    break
            if not match : unmatched.append(ig)
        self.value = unmatched
#___________________________________________________________
