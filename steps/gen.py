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

    def __init__(self,minPt=-1.0,minStatus=-1):
        self.oneP4=utils.LorentzV()
        self.sumP4=utils.LorentzV()
        self.zeroP4=utils.LorentzV()
        self.minPt=minPt
        self.minStatus=minStatus
        
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
        for iGen in range(min([maxPrintSize,size])) :

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
