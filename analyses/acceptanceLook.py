#!/usr/bin/env python

import supy
import calculables,steps,samples, ROOT as r

MeV2GeV = 1.0e+3
GeV=1.0e+3
TeV=1.0e+3*GeV


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
        return {
            'objects'  : objects,
            'lepton'   : leptons,
            'minJetPt' : 10.0,
                }

    def listOfSteps(self,config) :
        obj = config["objects"]
        lepton = obj['muon']
        print lepton
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
        shv = supy.steps.histos.value
        isoVar = "%sRelativeIso%s"%lepton
        lepInd = 'mu_staco_Indices'
        # stepsList += [supy.steps.printer.printstuff(["%sRelativeIso%s"%lepton])]
        stepsList+= [shv(isoVar, 50, 0.0, 1.0, indices = lepInd, index=0),]

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
                              calculables.muon.Indices(obj['muon'], ptMin=10.)
                              ]
        listOfCalculables += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])

        return listOfCalculables

    def listOfSampleDictionaries(self) :
        exampleDict = supy.samples.SampleHolder()
        # one file from mc12_8TeV.176584*_p1328/
        exampleDict.add("ttbar_sherpa",
                        '["/tmp/gerbaudo/wA_noslep_WH_2Lep_11/NTUP_SUSY.01176858._000001.root.1"'
                        ',"/tmp/gerbaudo/wA_noslep_WH_2Lep_11/NTUP_SUSY.01176858._000002.root.1"]',
                        xs = 1.140) #pb # 1.1402753294*0.30636*0.3348500000
        print "Fix cross sections"
        return [exampleDict]

    def listOfSamples(self,config) :
        test = False #True
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
