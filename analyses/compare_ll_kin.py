#!/usr/bin/env python

import supy
import calculables,steps,samples, ROOT as r

MeV2GeV = 1.0e+3
GeV=1.0e+3
TeV=1.0e+3*GeV

class compare_ll_kin(supy.analysis) :
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
        return {
            'objects'  : objects,
            'lepton'   : leptons,
            'minJetPt' : 10.0,
                }

    def listOfSteps(self,config) :
        indices = ['genIndices'+p for p in ['lpos','lneg', 'b', 'bbar', 'v', 'vbar']]
        def dropInd(idxName) :
            return idxName.replace('Index','').replace('genIndices','').replace('Indices','')
        stepsList = [
            supy.steps.printer.progressPrinter(),
            # sometimes multiple b or bbar?
            supy.steps.filters.multiplicity('genIndicesb',max=1),
            supy.steps.filters.multiplicity('genIndicesbbar',max=1),
            #steps.gen.particlePrinter(),
            supy.steps.histos.multiplicity("genP4", max=50),
            ]
        stepsList += [
            supy.steps.histos.multiplicity(ii, max=4) for ii in indices]
        stepsList+=[supy.steps.histos.pt("genP4",
                                         20, 0.0, 300*GeV,
                                         indices = ii, index = 0, xtitle = dropInd(ii))
                    for ii in indices]
#            supy.steps.histos.multiplicity(var = "jet_pt", max = 20),
#            supy.steps.histos.multiplicity(var = "jet_Indices", max = 20),
#            supy.steps.filters.multiplicity(min = 4, var = "jet_Indices"),
#            supy.steps.histos.multiplicity(var = "jet_Indices", max = 20),
#            supy.steps.histos.eta(var = "jet_P4", N = 20, low = -2., up = +2., indices = "jet_Indices"),
#            supy.steps.histos.value(var = "jet_M01" , N = 50, low = 0., up = 1.0e+3*GeV),
#            supy.steps.filters.value(var = "jet_M01", min = 1.0*TeV),
#            supy.steps.other.skimmer()
        return stepsList

    def listOfCalculables(self,config) :
        print config
        obj = config["objects"]
        listOfCalculables = supy.calculables.zeroArgs(supy.calculables)
        listOfCalculables += [calculables.gen.genP4(),
                              calculables.gen.sherpaTtbarProductsIndices(),
                              calculables.gen.genIndiceslpos(), calculables.gen.genIndiceslneg(),
                              calculables.gen.genIndicesb(),    calculables.gen.genIndicesbbar(),
                              calculables.gen.genIndicesv(),    calculables.gen.genIndicesvbar(),
                              ]
        listOfCalculables += [supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
                             ]

#        listOfCalculables += [calculables.Davide.Indices(collection = ("jet_",""), # prefix,suffix
#                                                         ptMin=30.*GeV,
#                                                         etaMax=1.5),
#                              calculables.Davide.P4(collection = ("jet_",""))]
#        listOfCalculables += [calculables.Davide.M01(collection = ("jet_",""))]
        return listOfCalculables

    def listOfSampleDictionaries(self) :
        protocol="root://xrootd-disk.pic.es/"
        basedir="/pnfs-disk/pic.es/at3/projects/TOPD3PD/2011/Skimming/DPD_prod01_02_October11"
        exampleDict = supy.samples.SampleHolder()
        exampleDict.add("ttbar_sherpa",
                        '["/tmp/gerbaudo/data/NTUP_SUSY.01116143._000136.root.1"]',
                        xs = 89.3615) #pb
        print "Fix cross sections"
        return [exampleDict]

    def listOfSamples(self,config) :
        test = True
        nEventsMax= 100 if test else None
        print 'nEventsMax :',nEventsMax
        return (supy.samples.specify(names = "ttbar_sherpa", color = r.kViolet,# effectiveLumi = 10.0,
                                     nEventsMax=nEventsMax)
                )

    def conclude(self,pars) :
        #make a pdf file with plots from the histograms created above
        org = self.organizer(pars)
        org.scale(lumiToUseInAbsenceOfData=20.0)
        supy.plotter( org,
                      pdfFileName = self.pdfFileName(org.tag),
                      #samplesForRatios = ("Example_Skimmed_900_GeV_Data","Example_Skimmed_900_GeV_MC"),
                      #sampleLabelsForRatios = ("data","sim"),
                      ).plotAll()
