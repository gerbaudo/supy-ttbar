#!/usr/bin/env python

import supy
import calculables,steps,samples, ROOT as r

MeV2GeV = 1.0e+3
GeV=1.0e+3
TeV=1.0e+3*GeV

mcColl = ('mc_','')

class acceptanceLook(supy.analysis) :
    def parameters(self) :
        objects = self.vary()
        leptons = self.vary()
        fields = [ "jet", "met", "muon", "electron",]
        objects['std'] =  dict(zip(fields, [("jet_AntiKt4LCTopo_",""),
                                            "MET_RefFinal_STVF_",
                                            ("mu_staco_",""),
                                            ("el_","")]))

        fieldsLepton =                           ['name', 'ptMin', 'etaMax',   'isoVar', 'isoType']
        leptons['muon'] = dict(zip(fieldsLepton, ['muon',      10,      2.4, 'ptcone30', 'relative']))
        jetPars = {'maxEta' : 2.5, 'minPt' : 20.0}
        return {
            'objects'  : objects,
            'lepton'   : leptons,
            'jetPars'  : jetPars
                }

    def listOfSteps(self,config) :
        obj = config["objects"]
        lepton = obj['muon']
        print lepton
        indices = ['genIndices'+p for p in ['lpos','lneg', 'b', 'bbar', 'v', 'vbar']]
        def dropInd(idxName) :
            return idxName.replace('Index','').replace('genIndices','').replace('Indices','')
        hi, wi = 'higgsIndices'.join(mcColl), 'wIndices'.join(mcColl)
        stepsList = [
            supy.steps.printer.progressPrinter(),
#              supy.steps.printer.printstuff([hi,
#                                             'higgsChildrenIndices'.join(mcColl),
#                                             wi,
#                                             'genIndicesW',
#                                             ]),
#             steps.gen.particlePrinter(), #indices='whIndices'.join(mcColl)),
#             steps.gen.particlePrinter(indices='whIndices'.join(mcColl)),
            supy.steps.histos.multiplicity("genP4", max=50),
            supy.steps.histos.multiplicity('genJetIndices', max=50),
            supy.steps.histos.multiplicity(hi, max=5),
            supy.steps.histos.multiplicity(wi, max=5),
            supy.steps.histos.multiplicity(wi, max=5),
            supy.steps.histos.multiplicity('wChildrenIndices'.join(mcColl), max=10),
            steps.gen.higgsDecay(),
            steps.gen.wDecay(),
            supy.steps.histos.pt    ('genP4',100, 0., 250*GeV, hi, xtitle='H_{truth}'),
            supy.steps.histos.absEta('genP4',100, 0.,       5, hi, xtitle='H_{truth}'),
            supy.steps.histos.pt    ('genP4',100, 0., 250*GeV, wi, xtitle='W_{truth}'),
            supy.steps.histos.absEta('genP4',100, 0.,       5, wi, xtitle='W_{truth}'),
            #supy.steps.filters.multiplicity('wChildrenIndices'.join(mcColl), min=2, max=2), #.invert(),

            supy.steps.filters.value('wDecayType'.join(mcColl), min=1, max=1),
            
            ]
        stepsList+=[supy.steps.histos.pt("genP4",
                                         20, 0.0, 300*GeV,
                                         indices = ii.join(mcColl), index = 0,
                                         xtitle = dropInd(ii))
                    for ii in indices]
        shv = supy.steps.histos.value
        isoVar = "%sRelativeIso%s"%lepton
        lepInd = 'mu_staco_Indices'
        # stepsList += [supy.steps.printer.printstuff(["%sRelativeIso%s"%lepton])]
        stepsList+= [shv(isoVar, 50, 0.0, 1.0, indices = lepInd, index=0),]

        return stepsList

    def listOfCalculables(self,config) :
        obj = config['objects']
        lepton = config['lepton']
        ptMin, etaMax = lepton['ptMin'], lepton['etaMax']
        jetPars = config['jetPars']
        listOfCalculables = supy.calculables.zeroArgs(supy.calculables)
        listOfCalculables += [calculables.gen.genIndices(mcColl, [+24,-24],'W')]
        listOfCalculables += [calculables.muon.Indices(obj['muon'], ptMin=ptMin),
                              calculables.genjet.genJetP4(),
                              calculables.genjet.genJetIndices(ptMin=jetPars['minPt'],
                                                               etaMax=jetPars['maxEta']),
                              ]
        print 'mcColl : ',mcColl
        listOfCalculables += supy.calculables.fromCollections(calculables.gen, [mcColl])
        listOfCalculables += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])

        return listOfCalculables

    def listOfSampleDictionaries(self) :
        exampleDict = supy.samples.SampleHolder()
        exampleDict.add('WH_2Lep_11',
                        '["/tmp/gerbaudo/wA_noslep_WH_2Lep_11/NTUP_SUSY.01176858._000001.root.1"'
                        ',"/tmp/gerbaudo/wA_noslep_WH_2Lep_11/NTUP_SUSY.01176858._000002.root.1"]',
                        xs = 1.140) #pb # 1.1402753294*0.30636*0.3348500000
        return [exampleDict]

    def listOfSamples(self,config) :
        test = True #False
        nEventsMax= 4 if test else None
        return (supy.samples.specify(names='WH_2Lep_11', color = r.kViolet, nEventsMax=nEventsMax)
                )

    def conclude(self,pars) :
        org = self.organizer(pars)
        org.scale(lumiToUseInAbsenceOfData=20.0)
        supy.plotter(org,
                     pdfFileName = self.pdfFileName(org.tag),
                     doLog=False,
                     ).plotAll()
