# largely based on
# https://github.com/betchart/susycaf/blob/master/calculables/gen.py
import collections, ROOT as r
from supy import utils,analysisStep
#####################################
try:
    import pdgLookup
    pdgLookupExists = True
except ImportError:
    pdgLookupExists = False
#####################################


class PdgLookup:
    def __init__(self) :
        self.names = {
            -1:'/d', 1:'d',
             -2:'/u', 2:'u',
             -3:'/s', 3:'s',
             -4:'/c', 4:'c',
             -5:'/b', 5:'b',
             -6:'/t', 6:'t',
             -11:'e+', 11:'e-',
             -12:'ve', 12:'ve',
             -13:'mu+', 13:'mu-',
             -14:'vmu', 14:'vmu',
             -15:'tau+', 15:'tau-',
             -16:'vtau', 16:'vtau',
             21:'g',
             -22:'gamma', 22:'gamma',
             -23:'Z', 23:'Z',
             -24:'W-', 24:'W+',
             -25:'h', 25:'h',
             -511:'/B0',511:'B0',
             -513:'/B*0',513:'B*0',
             -521:'B-',521:'B+',
             -523:'B*-',523:'B*+',
             -531:'/B0s',531:'B0s',
             -533:'/B*0s',533:'B*0s',
        }
    def pdgid_to_name(self,id) : return self.names[id] if id in self.names else 'unknown'

#####################################
class ParticleCountFilter(analysisStep) :
    def __init__(self, reqDict) :
        self.reqDict = reqDict
    def select (self,eventVars) :
        for key,value in self.reqDict.iteritems() :
            if eventVars["GenParticleCategoryCounts"][key]!=value : return False
        return True
#####################################
#####################################
class particlePrinter(analysisStep) :

    def __init__(self,minPt=-1.0,minStatus=-1, indices=''):
        self.oneP4=utils.LorentzV()
        self.sumP4=utils.LorentzV()
        self.zeroP4=utils.LorentzV()
        self.minPt=minPt
        self.minStatus=minStatus
        self.indices=indices

    def uponAcceptance (self,eventVars) :
        pdgLookupExists = True
        pdgLookup = PdgLookup()
        self.sumP4.SetCoordinates(0.0,0.0,0.0,0.0)

        parents=set([p for pp in eventVars['mc_parent_index'] for p in pp])
        #print "parents: ",parents
        print "-----------------------------------------------------------------------------------"
        print " i  st   par         id            name        E        pt       eta    phi    mass"
        print "-----------------------------------------------------------------------------------"

        size=len(eventVars["genP4"])
        maxPrintSize=50
        indices = eventVars[self.indices] if self.indices else range(min([maxPrintSize,size]))
        for iGen in indices :

            p4=eventVars["genP4"][iGen]
            if p4.pt()<self.minPt : continue

            status=eventVars['mc_status'][iGen]
            if status<self.minStatus : continue

            pars = [i for i in eventVars['mc_parent_index'][iGen]]
            pdgId=eventVars['mc_pdgId'][iGen]
            outString=""
            outString+="%#2d"%iGen
            outString+=" %#3d"%status
            outString+= str(pars).rjust(6)
            outString+=" %#10d"%pdgId
            if pdgLookupExists : outString+=" "+pdgLookup.pdgid_to_name(pdgId).rjust(15)
            else :                 outString+="".rjust(16)
            outString+="  %#7.1f"%p4.E()
            outString+="  %#8.1f"%p4.pt()
            outString+="  %#8.1f"%p4.eta()
            outString+="  %#5.1f"%p4.phi()
            outString+="  %#6.1f"%p4.mass()
            if not (iGen in parents) : outString+="   non-mo"
            print outString
        print
#####################################
class ttbarPrinter(analysisStep) :

    def uponAcceptance (self,ev) :
        ids = list(ev['mc_pdgId'])
        if not (6 in ids and -6 in ids) : return
        iTop = ids.index(6)
        iTbar= ids.index(-6)
        parents = ev['mc_parent_index'] #['mc_parents']
        parentIds = [[ids[p] for p in par] for par in parents]
        p4s = ev['genP4']
        status = ev['mc_status']
        pdg = PdgLookup()

        iGs = filter(lambda i: ids[i]==21, range(max(iTop,iTbar)+1,len(ids)))
        iTopChildren = [i for i,pars in enumerate(parentIds) if 6 in pars or -6 in pars]
        iWChildren = [i for i,pars in enumerate(parentIds) if 24 in pars or -24 in pars]
        iHChildren = [i for i,pars in enumerate(parentIds) if 25 in pars]
        width=15
        fieldNames = ['item','parents','parentIds','pt','eta','phi','status']
        print '-'*width*(len(fieldNames))
        print ''.join([("%s"%n).rjust(width) for n in fieldNames])
        print
        for i in [iTop,iTbar]+iTopChildren+iWChildren+iHChildren :
            print ''.join((("%s"%d).rjust(width) for d in ["%d (%s)"%(i,pdg.pdgid_to_name(ids[i])),
                                                           list(parents[i]),
                                                           list(parentIds[i])]\
                               +["%.3f"%v for v in [p4s[i].pt(),p4s[i].eta(), p4s[i].phi()]]\
                               +[status[i]]))
        print
#___________________________________________________________
class particleDecay(analysisStep) :
    def __init__(self, var='', title='') :
        self.var = var
        self.title = title
    def makeLabels(self, eventVars) : print 'to be implemented in inheriting class'
    def uponAcceptance(self, eventVars) :
        if not hasattr(self, 'labels') : self.makeLabels(eventVars)
        iBin = self.labels.index(self.decays[eventVars[self.var]])
        self.book.fill(iBin, self.var, self.nBins, 0.0, self.nBins,
                       title = self.title, xAxisLabels = self.labels)
#___________________________________________________________
class higgsDecay(particleDecay) :
    def __init__(self, var='mc_higgsDecayType', title=';Higgs Decay;events / bin') :
        super(higgsDecay, self).__init__(var, title)
    def makeLabels(self, eventVars) :
        self.decays = { # see gen.higgsDecayType
            0  : 'WW',
            1  : 'ZZ',
            2  : 'tautau',
            3  : 'bbbar',
            4  : 'mumu',
            -1 : 'unknown',
            }
        self.labels = sorted(list(set(self.decays.values())))
        self.nBins  = len(self.labels)
#___________________________________________________________
class wDecay(particleDecay) :
    def __init__(self, var='mc_wDecayType', title=';W Decay;events / bin') :
        super(wDecay, self).__init__(var, title)
    def makeLabels(self, eventVars) :
        self.decays = { # see gen.wDecayType
            0  : 'l#nu',
            1  : "q#bar{q}",
            -1 : 'unknown',
            }
        self.labels = sorted(list(set(self.decays.values())))
        self.nBins  = len(self.labels)
#___________________________________________________________
class deltaR(analysisStep) :
    def __init__(self, var='genP4', indices='', title='', N=50, xmin=0.0, xmax=4.0) :
        for item in ['var','indices', 'title', 'N', 'xmin', 'xmax'] : setattr(self, item, eval(item))
        self.moreName = 'deltaR'+var+':'+indices
        self.title = ';#Delta R ['+indices+'];events / bin'
    def uponAcceptance(self,eventVars) :
        val = eventVars[self.var]
        indices = eventVars[self.indices]
        if val is None : return
        if len(indices)!=2 :
            print "skipping event (%d indices)"%len(indices)
            return
        p0, p1 = val[indices[0]], val[indices[1]]
        self.book.fill(r.Math.VectorUtil.DeltaR(p0, p1),
                       self.moreName, self.N, self.xmin, self.xmax, title=self.title)
#___________________________________________________________
class tauDecay(particleDecay) :
    def __init__(self, indices='', index=None, var='mc_tauDecayType', title=';W Decay;events / bin') :
        super(wDecay, self).__init__(var, title)
    def makeLabels(self, eventVars) :
        self.decays = { # see gen.tauDecayType
            0  : 'l#nu',
            1  : "q#bar{q}",
            -1 : 'unknown',
            }
        self.labels = sorted(list(set(self.decays.values())))
        self.nBins  = len(self.labels)
