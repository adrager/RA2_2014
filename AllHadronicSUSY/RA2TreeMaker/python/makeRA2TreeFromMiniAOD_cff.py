# $Id: makeTreeFromPAT_cff.py,v 1.16 2013/01/24 15:42:53 mschrode Exp $
#

import FWCore.ParameterSet.Config as cms

def makeRA2TreeTreeFromMiniADO(process,
                    outFileName,
                    NJetsMin=2,
                    HTMin=350.,
                    MHTMin=0.,
                    reportEveryEvt=10,
                    testFileName="",
		    Global_Tag="",
		    MC=False,
		    debug = False,
		    QCD=False,
		    StoreAll=False,
                    numProcessedEvt=1000):

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag.globaltag = Global_Tag

    ## --- Log output ------------------------------------------------------
    process.load("FWCore.MessageService.MessageLogger_cfi")
    process.MessageLogger.cerr = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
        )
    process.MessageLogger.cout = cms.untracked.PSet(
        INFO = cms.untracked.PSet(reportEvery = cms.untracked.int32(reportEveryEvt))
        )
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        ) 


    ## --- Files to process ------------------------------------------------
    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(numProcessedEvt)
        )
    process.source = cms.Source("PoolSource",

        fileNames = cms.untracked.vstring(testFileName)
 #       fileNames = cms.untracked.vstring(
 #		'file:/nfs/dust/cms/user/csander/LHE/workdir/simulation_test/T1qqqqHV/output_66.root'
#		)
        )
        
    hltPath=['HLT_PFHT350_PFMET100_v*','HLT_PFNoPUHT350_PFMET100_v*']
    process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
    process.hltHighLevel.HLTPaths = cms.vstring(hltPath)
    process.hltHighLevel.andOr = cms.bool(True)
    process.hltHighLevel.throw = cms.bool(False)

    process.HLTSelection = cms.Sequence(
        process.hltHighLevel
        )
    if MC:
        print "Running on MC: removing HLT selection"
        process.HLTSelection.remove(process.hltHighLevel)
    elif not hltPath:
        print "Empty list of HLT paths: removing HLT selection"
        process.HLTSelection.remove(process.hltHighLevel)
    ## --- Output file -----------------------------------------------------
    process.TFileService = cms.Service(
        "TFileService",
        fileName = cms.string(outFileName+".root")
        )
	    
    ## --- Selection sequences ---------------------------------------------
    # leptons
    process.load("PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi")
    process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")
    process.selectedIDIsoMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>5. &&
    (pfIsolationR04().sumChargedHadronPt+
    max(0.,pfIsolationR04().sumNeutralHadronEt+
    pfIsolationR04().sumPhotonEt-
    0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedIDMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>5. && 
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedIDIsoElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>5. &&
    gsfTrack.isAvailable() &&
    gsfTrack.trackerExpectedHitsInner.numberOfLostHits<2 &&
    (pfIsolationVariables().sumChargedHadronPt+
    max(0.,pfIsolationVariables().sumNeutralHadronEt+
    pfIsolationVariables().sumPhotonEt-
    0.5*pfIsolationVariables().sumPUPt))/pt < 0.20'''))
    process.selectedIDElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>5. &&
    gsfTrack.isAvailable() &&
    gsfTrack.trackerExpectedHitsInner.numberOfLostHits<2'''))
    
    
       ## --- Setup of TreeMaker ----------------------------------------------
    FilterNames = cms.VInputTag()
    FilterNames.append(cms.InputTag("HBHENoiseFilterRA2","HBHENoiseFilterResult","PAT"))
    FilterNames.append(cms.InputTag("beamHaloFilter"))
    FilterNames.append(cms.InputTag("eeNoiseFilter"))
    FilterNames.append(cms.InputTag("trackingFailureFilter"))
    FilterNames.append(cms.InputTag("inconsistentMuons"))
    FilterNames.append(cms.InputTag("greedyMuons"))
    FilterNames.append(cms.InputTag("ra2EcalTPFilter"))
    FilterNames.append(cms.InputTag("ra2EcalBEFilter"))
    FilterNames.append(cms.InputTag("hcalLaserEventFilter"))
    FilterNames.append(cms.InputTag("ecalLaserCorrFilter"))
    FilterNames.append(cms.InputTag("eeBadScFilter"))
    FilterNames.append(cms.InputTag("PBNRFilter"))
    FilterNames.append(cms.InputTag("HCALLaserEvtFilterList2012"))
    FilterNames.append(cms.InputTag("manystripclus53X"))
    FilterNames.append(cms.InputTag("toomanystripclus53X"))
    FilterNames.append(cms.InputTag("logErrorTooManyClusters"))
    FilterNames.append(cms.InputTag("RA2HONoiseFilter"))
    
    
    ## --- Setup WeightProducer -------------------------------------------
    from AllHadronicSUSY.WeightProducer.getWeightProducer_cff import getWeightProducer
    process.WeightProducer = getWeightProducer(testFileName)
    process.WeightProducer.Lumi                       = cms.double(5000)
    process.WeightProducer.PU                         = cms.int32(0) # PU S10 3 for S10 2 for S7
    process.WeightProducer.FileNamePUDataDistribution = cms.string("NONE")
    print process.WeightProducer.PU
    from AllHadronicSUSY.RA2TreeMaker.ra2TreeMaker import RA2TreeMaker
    process.RA2TreeMaker2 = RA2TreeMaker.clone(
    	TreeName          = cms.string("RA2PreSelection"),
    	VertexCollection  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    	VarsDouble        = cms.VInputTag(cms.InputTag('WeightProducer:weight')),
    	VarsDoubleNamesInTree = cms.vstring('Weight'),
    	Filters           = cms.VInputTag(FilterNames), #FilterNames
    	MC = MC,
    	QCD = QCD,  # use this switch to store only information needed for qcd background estimation method
    	StoreAll = StoreAll, # override QCD only switch and stores in addition to non QCD selection also gen jet information use this for maximum information in tree
    	debug =debug,
    	LeptonTag = cms.VInputTag(cms.InputTag('selectedIDIsoMuons'),cms.InputTag('selectedIDMuons'),cms.InputTag('selectedIDIsoElectrons'),cms.InputTag('selectedIDElectrons')),
    	LeptonTagName = cms.vstring('RecoIsoMuon','RecoMuon','RecoIsoElec','RecoElec'),
    #RA2JetsTag = cms.InputTag("patJetsAK5PFCHS"),   
    	RA2DefaultJetsTag = cms.InputTag("slimmedJets"),  
    	ra2JetsCollectionInputTag = cms.VInputTag(cms.InputTag('slimmedJets')),
    	ra2JetsCollectionNameInTree = cms.vstring('ak4'),
    	ra2JetsBTagInputTag = cms.vstring('combinedSecondaryVertexBJetTags','combinedSecondaryVertexBJetTags'),
    	ra2JetsBTagValueInput_ = cms.vdouble(0.679,0.679),
    )

    ## --- Final paths ----------------------------------------------------

    process.dump = cms.EDAnalyzer("EventContentAnalyzer")
    process.WriteTree = cms.Path(
    	process.selectedIDIsoMuons *
    	process.selectedIDMuons *
    	process.selectedIDIsoElectrons *
    	process.selectedIDElectrons *
    	process.WeightProducer *
 #   	process.PrintDecay *
    	process.RA2TreeMaker2

        )
