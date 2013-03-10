#!/usr/bin/env python

import supy
import calculables,steps,samples, ROOT as r

GeV=1.0e+3
TeV=1.0e+3*GeV

class compare_ll_kin(supy.analysis) :
    def parameters(self) :
        return {"minJetPt" : 10.0,
                }

    def listOfSteps(self,config) :
        
        outList=[
            supy.steps.printer.progressPrinter(),
            supy.steps.histos.multiplicity(var = "jet_pt", max = 20), 
            supy.steps.histos.multiplicity(var = "jet_Indices", max = 20), 
            supy.steps.filters.multiplicity(min = 4, var = "jet_Indices"),
            supy.steps.histos.multiplicity(var = "jet_Indices", max = 20),
            supy.steps.histos.eta(var = "jet_P4", N = 20, low = -2., up = +2., indices = "jet_Indices"),
            supy.steps.histos.value(var = "jet_M01" , N = 50, low = 0., up = 1.0e+3*GeV),
            supy.steps.filters.value(var = "jet_M01", min = 1.0*TeV),
            supy.steps.other.skimmer()
            ]
        return outList
    
    def listOfCalculables(self,config) :
        listOfCalculables = supy.calculables.zeroArgs(supy.calculables)
        listOfCalculables += [calculables.Davide.Indices(collection = ("jet_",""), # prefix,suffix
                                                         ptMin=30.*GeV,
                                                         etaMax=1.5),
                              calculables.Davide.P4(collection = ("jet_",""))]
        listOfCalculables += [calculables.Davide.M01(collection = ("jet_",""))]
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
        return (supy.samples.specify(names = "ttbar_sherpa", color = r.kViolet, effectiveLumi = 10.0e+3)
                )

    def conclude(self,pars) :
        #make a pdf file with plots from the histograms created above
        org = self.organizer(pars)
        org.scale()
        supy.plotter( org,
                      pdfFileName = self.pdfFileName(org.tag),
                      #samplesForRatios = ("Example_Skimmed_900_GeV_Data","Example_Skimmed_900_GeV_MC"),
                      #sampleLabelsForRatios = ("data","sim"),
                      ).plotAll()
