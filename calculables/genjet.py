from supy import wrappedChain,utils
#import ROOT as r



class genJetP4(wrappedChain.calculable) :
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

    # nTruthJets = tree.jet_AntiKt4TruthJets_n
    # truthJets = [tlv(0., 0., 0., 0.) for i in xrange(nTruthJets)]
    # pdg, parents, children = tree.mc_pdgId, tree.mc_parent_index, tree.mc_child_index
    # interestingIhiggs, higgsChildren, higgsParents = findHiggs(pdg, parents, children)
    # if len(interestingIhiggs)>1 :
    #     print "skip event with multiple (%d) higgs"%len(interestingIhiggs)
    #     continue
    # higgsChildren, higgsParents = higgsChildren[0], higgsParents[0]
    # hDecay = gen.guessHdecayLabel(higgsChildren)
    # for j, pt, eta, phi, en in zip(truthJets,
    #                                tree.jet_AntiKt4TruthJets_pt,
    #                                tree.jet_AntiKt4TruthJets_eta,
    #                                tree.jet_AntiKt4TruthJets_phi,
    #                                tree.jet_AntiKt4TruthJets_E) :
    #     j.SetPtEtaPhiE(pt*MeV2GeV, eta, phi, en*MeV2GeV)
    # truthJets = [j for j in truthJets if j.Pt()>minJetPt]
