#!/usr/bin/env python

import supy
import calculables,steps,samples, ROOT as r

GeV=1.0e+3
TeV=1.0e+3*GeV

class example_ttbar(supy.analysis) :
    def mainTree(self) : return ("/","physics")
    def otherTreesToKeepWhenSkimming(self) : return []
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
        exampleDict.add("Data_skimMu",
                        '["%s%s/data11_7TeV.physics_Muons.merge.NTUP_TOP_PeriodJ2_part1_skimMu.root"]'%(protocol,basedir),
                        lumi = 1.0e+3 ) #/pb
        exampleDict.add("Zmumu_skimMu",
                        '["%s%s/mc10_7TeV.107665.AlpgenJimmyZmumuNp5_pt20.merge.NTUP_TOP_skimMu.root"]'%(protocol,basedir),
                        xs = 96.25*0.9003) #pb * skim eff
        exampleDict.add("ttbar_skimMu",
                        '["%s%s/mc10_7TeV.105892.AlpgenJimmyttbarlnlnNp2_baseline.merge.NTUP_TOP_skimMu.root"]'%(protocol,basedir),
                        xs = 89.3615 * 0.5247) #pb * skim eff
        print "Fix cross sections"
        return [exampleDict]

    def listOfSamples(self,config) :
        return (supy.samples.specify(names = "Data_skimMu", color = r.kBlack, markerStyle = 20) +
                supy.samples.specify(names = "Zmumu_skimMu", color = r.kRed, effectiveLumi = 10.0e+3) +
                supy.samples.specify(names = "ttbar_skimMu", color = r.kViolet, effectiveLumi = 10.0e+3)
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
