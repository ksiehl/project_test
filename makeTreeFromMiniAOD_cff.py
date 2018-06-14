# $Id: makeTreeFromPAT_cff.py,v 1.16 2013/01/24 15:42:53 mschrode Exp $
#

import FWCore.ParameterSet.Config as cms
import os

def makeTreeTreeFromMiniAOD(process,
outFileName,
NJetsMin=2,
HTMin=350.,
MHTMin=0.,
reportEveryEvt=10,
testFileName="",
Global_Tag="",
METFiltersProcess="",
MC=False,
debug = False,
QCD=False,
LostLepton=False,
numProcessedEvt=1000,
doJECCorrection=False,
doPuppi=True,
leptonFilter=True,
doAK8Reclustering=False,
doAK10Reclustering=False,
doAK12Reclustering=False,
genJetsAK8Reclustering=True,
genJetsAK10Reclustering=False,
genJetsAK12Reclustering=False,
customizeHBHENoiseForEarlyData=False,
customizeHBHENoiseForRun2015D=True,
jsonFileName="",
isCrab=False,
reDoPruningAndSoftdrop=False,
reDoPruningAndSoftdropPuppi=True,
isRun2016BCD=False,
isRun2016E=False,
isRun2016F=False,
isRun2016GH=False
):

    if (MC):
#        process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
        process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
        customizeHBHENoiseForRun2015D=False
    else:
        process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

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

 ## ----------------------------------------------------------------------------------------------
## Triggers
## ----------------------------------------------------------------------------------------------
# The trigger results are saved to the tree as a vector
# Three vectors are saved:
# 1) names of the triggers
# 2) trigger results
# 3) trigger prescales
# the indexing of these vectors must match
# If the version number of the input trigger name is omitted,
# any matching trigger will be included (default behavior)
    from AllHadronicSUSY.Utils.triggerproducer_cfi import triggerProducer
    process.TriggerProducer = triggerProducer.clone( 
        trigTagArg1 = cms.string('TriggerResults::HLT'),
        trigTagArg2 = cms.string(''),
        trigTagArg3 = cms.string('HLT'),
        prescaleTagArg1 = cms.string('patTrigger'),
        prescaleTagArg2 = cms.string(''),
        prescaleTagArg3 = cms.string(''),
        triggerNameList = cms.vstring( # list of trigger names
            'HLT_Ele105_CaloIdVT_GsfTrkIdT_v',
            'HLT_Mu45_eta2p1_v',
            'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
            'HLT_Mu50_v',
            'HLT_Ele27_WPLoose_Gsf_v',
            'HLT_Mu24_eta2p1_v',
            'HLT_Mu24_v',
            'HLT_PFJet450_v',
            'HLT_PFJet500_v',
            'HLT_AK8PFJet360_TrimMass30_v',
            'HLT_Ele22_eta2p1_WPLoose_Gsf_v',
            'HLT_Ele22_eta2p1_WPTight_Gsf_v',
            'HLT_Ele23_WPLoose_Gsf_v',
            'HLT_Ele27_WPLoose_Gsf_v',
            'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
            'HLT_Ele32_eta2p1_WPTight_Gsf_v',
            'HLT_Ele27_eta2p1_WPTight_Gsf_v',
            'HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v',
            'HLT_IsoMu27_v',
            'HLT_IsoMu24_v',
            'HLT_IsoMu22_v',
            'HLT_IsoMu20_v',
            'HLT_IsoTkMu20_v',
            'HLT_Ele22_eta2p1_WP75_Gsf_v1',
            'HLT_Ele27_WP85_Gsf_v1',
            'HLT_Ele27_eta2p1_WP75_Gsf_v1',
            'HLT_Ele32_eta2p1_WP75_Gsf_v1'
            )
        )

## this is here if you want to apply the trigger selection in this step, instead of saving the trigger result        
    hltPath=['HLT_Mu45_eta2p1_v*','HLT_Ele105_CaloIdVT_GsfTrkIdT_v*']
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

############### JSON Filter            
    import FWCore.PythonUtilities.LumiList as LumiList
    import sys

    if not MC:
        if(len(jsonFileName)>0):
            import FWCore.PythonUtilities.LumiList as LumiList
#luca            process.source.lumisToProcess = LumiList.LumiList(filename = jsonFileName).getVLuminosityBlockRange()
        else:
            print "ERROR!! running on data with no json file applied!"
            sys.exit()

####### some gen infos

#    process.GenEventInfo = cms.Sequence()

    from AllHadronicSUSY.Utils.geneventinfo_cfi import geneventinfo

#    if MC:
    process.GenEventInfo = geneventinfo.clone()
    
	    
    ## --- Selection sequences ---------------------------------------------
    # leptons
    process.load("PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi")
    process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")
    process.selectedIDIsoMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    (pfIsolationR04().sumChargedHadronPt+
    max(0.,pfIsolationR04().sumNeutralHadronEt+
    pfIsolationR04().sumPhotonEt-
    0.50*pfIsolationR04().sumPUPt))/pt < 0.20 &&
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedIDMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    (isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedIDIsoElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    gsfTrack.isAvailable() &&
    gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS')<2 &&
    (pfIsolationVariables().sumChargedHadronPt+
    max(0.,pfIsolationVariables().sumNeutralHadronEt+
    pfIsolationVariables().sumPhotonEt-
    0.5*pfIsolationVariables().sumPUPt))/pt < 0.20'''))
    process.selectedIDElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
    gsfTrack.isAvailable() &&
    gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS')<2'''))


###Lepton Filter
    process.filterSeq = cms.Sequence ()

#    from AllHadronicSUSY.Utils.leptonfilter_cfi import leptonFilter
    process.load('AllHadronicSUSY.Utils.leptonfilter_cfi')

    process.leptonFilter.electronsInputTag = cms.InputTag("slimmedElectrons")
    process.leptonFilter.muonsInputTag = cms.InputTag("slimmedMuons")
    process.leptonFilter.eleFilterPtCut = cms.double(20.0)
    process.leptonFilter.muFilterPtCut = cms.double(20.0)

    if (leptonFilter):
        process.filterSeq = cms.Sequence (process.leptonFilter)
    
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

    from RecoBTag.Configuration.RecoBTag_cff import *
    from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
    process.slimmedJetsPFJetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
      j2tParametersVX,
      jets = cms.InputTag("iterativeCone5PFJets")
    )
    process.slimmedJetsPFJetTracksAssociatorAtVertex.jets = "slimmedJets"
    process.slimmedJetsPFJetTracksAssociatorAtVertex.tracks = "generalTracks"
    
    process.slimmedJetsPFImpactParameterTagInfos = impactParameterTagInfos.clone()
    process.slimmedJetsPFImpactParameterTagInfos.jetTracks = "slimmedJetsPFJetTracksAssociatorAtVertex"
    process.slimmedJetsPFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
    process.slimmedJetsPFSecondaryVertexTagInfos.trackIPTagInfos = "slimmedJetsPFImpactParameterTagInfos"
    #slimmedJetsPFSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
    #slimmedJetsPFSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("slimmedJetsPFSecondaryVertexTagInfos") )
    process.slimmedJetsPFCombinedSecondaryVertexBJetTags = combinedInclusiveSecondaryVertexV2BJetTags.clone()
 #   process.slimmedJetsPFStandardCombinedSecondaryVertex = combinedSecondaryVertex.clone()
 #   process.slimmedJetsPFCombinedSecondaryVertexBJetTags.jetTagComputer = cms.string('slimmedJetsPFStandardCombinedSecondaryVertex')
 #   process. slimmedJetsPFCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("slimmedJetsPFImpactParameterTagInfos"), cms.InputTag("slimmedJetsPFSecondaryVertexTagInfos") )
    
    process.slimmedJetsPFJetBtaggingSV = cms.Sequence(
    	process.slimmedJetsPFImpactParameterTagInfos *
    process.slimmedJetsPFSecondaryVertexTagInfos 
    # slimmedJetsPFStandardCombinedSecondaryVertex *
#    process.slimmedJetsPFCombinedSecondaryVertexBJetTags
    )
    process.slimmedJetsPFJetsBtag = cms.Sequence(
    process.slimmedJetsPFJetTracksAssociatorAtVertex *
    process.slimmedJetsPFJetBtaggingSV
    )
    
    ## isotrack producer
    from AllHadronicSUSY.Utils.trackIsolationMaker_cfi import trackIsolationFilter
    from AllHadronicSUSY.Utils.trackIsolationMaker_cfi import trackIsolationCounter
    ## default
    process.IsolatedTracks = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(15.0),
      isoCut              = cms.double(0.1),
      )
    #study
    process.IsolatedTracksPT10 = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(10.0),
      isoCut              = cms.double(0.1),
      )
    process.IsolatedTracksPT10IsoCut08 = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(10.0),
      isoCut              = cms.double(0.08),
      )
    process.IsolatedTracksPT10IsoCut12 = trackIsolationFilter.clone(
      doTrkIsoVeto= False,
      vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
      pfCandidatesTag = cms.InputTag("packedPFCandidates"),
      dR_ConeSize         = cms.double(0.3),
      dz_CutValue         = cms.double(0.05),
      minPt_PFCandidate   = cms.double(10.0),
      isoCut              = cms.double(0.12),
      )
    process.CountIsoTracks = trackIsolationCounter.clone(
      src = cms.InputTag("IsolatedTracks"),
      minNumber = 1,
      )

    '''
    jecLevelsAK8 = [
        "Spring16_25nsV3_DATA_L1FastJet_AK8PFchs.txt",
        "Spring16_25nsV3_DATA_L2Relative_AK8PFchs.txt",
        "Spring16_25nsV3_DATA_L3Absolute_AK8PFchs.txt"
        ]
    jecLevelsAK4 = [
#        "AAA/src/AllHadronicSUSY/PHYS14_25_V2_All_L1FastJet_AK4PFchs.txt",
        "Spring16_25nsV3_DATA_L1FastJet_AK4PFchs.txt",
        "Spring16_25nsV3_DATA_L2Relative_AK4PFchs.txt",
        "Spring16_25nsV3_DATA_L3Absolute_AK4PFchs.txt"
        ]
    if (MC):
        jecLevelsAK8 = [
            "Spring16_25nsV3_MC_L1FastJet_AK8PFchs.txt",
            "Spring16_25nsV3_MC_L2Relative_AK8PFchs.txt",
            "Spring16_25nsV3_MC_L3Absolute_AK8PFchs.txt"
            ]
        jecLevelsAK4 = [
            "Spring16_25nsV3_MC_L1FastJet_AK4PFchs.txt",
            "Spring16_25nsV3_MC_L2Relative_AK4PFchs.txt",
            "Spring16_25nsV3_MC_L3Absolute_AK4PFchs.txt"
            ]
    '''
#    print cms.string(os.environ['CMSSW_BASE']+"/src/AllHadronicSUSY/PHYS14_25_V2_All_L1FastJet_AK8PFchs.txt")

    process.substructureSequence = cms.Sequence()
    process.softdrop_onMiniAOD = cms.Sequence()
    process.pruning_onMiniAOD = cms.Sequence()
    process.redoPatJets = cms.Sequence()
    process.puppi_onMiniAOD = cms.Sequence()
    process.redoPuppiAK8Jets = cms.Sequence()

    if (doAK8Reclustering):
        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.chs = cms.EDFilter("CandPtrSelector",
                               src = cms.InputTag('packedPFCandidates'),
                               cut = cms.string('fromPV')
                               )

        process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak8PFJetsCHS"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.ak4PFJetsCHS = ak4PFJetsCHS.clone(src = 'chs', rParam = cms.double(0.4))
        process.ak8PFJetsCHS = ak8PFJetsCHS.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(0.8) )

        process.ak8PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'chs', jetPtMin = 100.0 )
        process.ak8PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone()
        process.ak8PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'chs', jetPtMin = 100.0 )
        process.ak8PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone()

        process.substructureSequence+=process.chs
        process.substructureSequence+=process.ak8PFJetsCHS
        process.substructureSequence+=process.NjettinessAK8

        process.softdrop_onMiniAOD += process.ak8PFJetsCHSSoftDrop + process.ak8PFJetsCHSSoftDropMass
        process.pruning_onMiniAOD += process.ak8PFJetsCHSPruned + process.ak8PFJetsCHSPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from ak8PFJetsCHS

        process.patJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHS' )
        process.patJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsCHS' )
        process.patJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8PFJetsCHSPrunedMass"), cms.InputTag("ak8PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK8:tau1"), cms.InputTag("NjettinessAK8:tau2"), cms.InputTag("NjettinessAK8:tau3")]
        process.patJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8") )
        process.selectedPatJetsAK8 = selectedPatJetsAK8.clone( cut = cms.string('pt > 20') )

        process.redoPatJets+=process.patJetCorrFactorsAK8
        process.redoPatJets+=process.patJetsAK8
        process.redoPatJets+=process.selectedPatJetsAK8

        if (reDoPruningAndSoftdrop):
            process.patJetCorrFactorsAK8Pruned = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSPruned' )
            process.patJetsAK8Pruned = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSPruned' )
            process.patJetsAK8Pruned.userData.userFloats.src = [ "" ]
            process.patJetsAK8Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Pruned") )
            process.selectedPatJetsAK8Pruned = selectedPatJetsAK8.clone( 
                src = cms.InputTag('patJetsAK8Pruned'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK8Pruned
            process.redoPatJets+=process.patJetsAK8Pruned
            process.redoPatJets+=process.selectedPatJetsAK8Pruned

            process.patJetCorrFactorsAK8Softdrop = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSSoftDrop' )
            process.patJetsAK8Softdrop = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSSoftDrop' )
            process.patJetsAK8Softdrop.userData.userFloats.src = [ "" ]
            process.patJetsAK8Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Softdrop") )
            process.selectedPatJetsAK8Softdrop = selectedPatJetsAK8.clone(
                src = cms.InputTag('patJetsAK8Softdrop'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK8Softdrop
            process.redoPatJets+=process.patJetsAK8Softdrop
            process.redoPatJets+=process.selectedPatJetsAK8Softdrop

    if (doAK10Reclustering):
        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.chs = cms.EDFilter("CandPtrSelector",
                               src = cms.InputTag('packedPFCandidates'),
                               cut = cms.string('fromPV')
                               )

        process.NjettinessAK10 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak10PFJetsCHS"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 1.0 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.ak10PFJetsCHS = ak8PFJetsCHS.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(1.0) )

        process.ak10PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(1.0) )
        process.ak10PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            distMax = cms.double(1.0),
            matched = cms.InputTag("ak10PFJetsCHSPruned"),
            src = cms.InputTag("ak10PFJetsCHS")
            )
        process.ak10PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(1.0) )
        process.ak10PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            distMax = cms.double(1.0),
            matched = cms.InputTag("ak10PFJetsCHSSoftDrop"),
            src = cms.InputTag("ak10PFJetsCHS")
            )

        process.substructureSequence+=process.chs
        process.substructureSequence+=process.ak10PFJetsCHS
        process.substructureSequence+=process.NjettinessAK10

        process.softdrop_onMiniAOD += process.ak10PFJetsCHSSoftDrop + process.ak10PFJetsCHSSoftDropMass
        process.pruning_onMiniAOD += process.ak10PFJetsCHSPruned + process.ak10PFJetsCHSPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from ak8PFJetsCHS

        process.patJetCorrFactorsAK10 = patJetCorrFactorsAK8.clone( src = 'ak10PFJetsCHS' )
        process.patJetsAK10 = patJetsAK8.clone( jetSource = 'ak10PFJetsCHS' )
        process.patJetsAK10.userData.userFloats.src = [ cms.InputTag("ak10PFJetsCHSPrunedMass"), cms.InputTag("ak10PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK10:tau1"), cms.InputTag("NjettinessAK10:tau2"), cms.InputTag("NjettinessAK10:tau3")]
        process.patJetsAK10.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK10") )
        process.selectedPatJetsAK10 = selectedPatJetsAK8.clone( cut = cms.string('pt > 20'), src = cms.InputTag("patJetsAK10") ) #LUCA

        process.redoPatJets+=process.patJetCorrFactorsAK10
        process.redoPatJets+=process.patJetsAK10
        process.redoPatJets+=process.selectedPatJetsAK10

        if (reDoPruningAndSoftdrop):
            process.patJetCorrFactorsAK10Pruned = patJetCorrFactorsAK8.clone( src = 'ak10PFJetsCHSPruned' )
            process.patJetsAK10Pruned = patJetsAK8.clone( jetSource = 'ak10PFJetsCHSPruned' )
            process.patJetsAK10Pruned.userData.userFloats.src = [ "" ]
            process.patJetsAK10Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK10Pruned") )
            process.selectedPatJetsAK10Pruned = selectedPatJetsAK8.clone( 
                src = cms.InputTag('patJetsAK10Pruned'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK10Pruned
            process.redoPatJets+=process.patJetsAK10Pruned
            process.redoPatJets+=process.selectedPatJetsAK10Pruned

            process.patJetCorrFactorsAK10Softdrop = patJetCorrFactorsAK8.clone( src = 'ak10PFJetsCHSSoftDrop' )
            process.patJetsAK10Softdrop = patJetsAK8.clone( jetSource = 'ak10PFJetsCHSSoftDrop' )
            process.patJetsAK10Softdrop.userData.userFloats.src = [ "" ]
            process.patJetsAK10Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK10Softdrop") )
            process.selectedPatJetsAK10Softdrop = selectedPatJetsAK8.clone(
                src = cms.InputTag('patJetsAK10Softdrop'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK10Softdrop
            process.redoPatJets+=process.patJetsAK10Softdrop
            process.redoPatJets+=process.selectedPatJetsAK10Softdrop



    if (doAK12Reclustering):
        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.chs = cms.EDFilter("CandPtrSelector",
                               src = cms.InputTag('packedPFCandidates'),
                               cut = cms.string('fromPV')
                               )

        process.NjettinessAK12 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak12PFJetsCHS"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 1.2 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.ak12PFJetsCHS = ak8PFJetsCHS.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(1.2) )

        process.ak12PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(1.2) )
        process.ak12PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            distMax = cms.double(1.2),
            matched = cms.InputTag("ak12PFJetsCHSPruned"),
            src = cms.InputTag("ak12PFJetsCHS")
            )
        process.ak12PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'chs', jetPtMin = 100.0,     rParam = cms.double(1.2) )
        process.ak12PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            distMax = cms.double(1.2),
            matched = cms.InputTag("ak12PFJetsCHSSoftDrop"),
            src = cms.InputTag("ak12PFJetsCHS")
            )

        process.substructureSequence+=process.chs
        process.substructureSequence+=process.ak12PFJetsCHS
        process.substructureSequence+=process.NjettinessAK12

        process.softdrop_onMiniAOD += process.ak12PFJetsCHSSoftDrop + process.ak12PFJetsCHSSoftDropMass
        process.pruning_onMiniAOD += process.ak12PFJetsCHSPruned + process.ak12PFJetsCHSPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from ak8PFJetsCHS

        process.patJetCorrFactorsAK12 = patJetCorrFactorsAK8.clone( src = 'ak12PFJetsCHS' )
        process.patJetsAK12 = patJetsAK8.clone( jetSource = 'ak12PFJetsCHS' )
        process.patJetsAK12.userData.userFloats.src = [ cms.InputTag("ak12PFJetsCHSPrunedMass"), cms.InputTag("ak12PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK12:tau1"), cms.InputTag("NjettinessAK12:tau2"), cms.InputTag("NjettinessAK12:tau3")]
        process.patJetsAK12.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK12") )
        process.selectedPatJetsAK12 = selectedPatJetsAK8.clone( cut = cms.string('pt > 20'), src = cms.InputTag("patJetsAK12") ) #LUCA

        process.redoPatJets+=process.patJetCorrFactorsAK12
        process.redoPatJets+=process.patJetsAK12
        process.redoPatJets+=process.selectedPatJetsAK12

        if (reDoPruningAndSoftdrop):
            process.patJetCorrFactorsAK12Pruned = patJetCorrFactorsAK8.clone( src = 'ak12PFJetsCHSPruned' )
            process.patJetsAK12Pruned = patJetsAK8.clone( jetSource = 'ak12PFJetsCHSPruned' )
            process.patJetsAK12Pruned.userData.userFloats.src = [ "" ]
            process.patJetsAK12Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK12Pruned") )
            process.selectedPatJetsAK12Pruned = selectedPatJetsAK8.clone( 
                src = cms.InputTag('patJetsAK12Pruned'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK12Pruned
            process.redoPatJets+=process.patJetsAK12Pruned
            process.redoPatJets+=process.selectedPatJetsAK12Pruned

            process.patJetCorrFactorsAK12Softdrop = patJetCorrFactorsAK8.clone( src = 'ak12PFJetsCHSSoftDrop' )
            process.patJetsAK12Softdrop = patJetsAK8.clone( jetSource = 'ak12PFJetsCHSSoftDrop' )
            process.patJetsAK12Softdrop.userData.userFloats.src = [ "" ]
            process.patJetsAK12Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK12Softdrop") )
            process.selectedPatJetsAK12Softdrop = selectedPatJetsAK8.clone(
                src = cms.InputTag('patJetsAK12Softdrop'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsAK12Softdrop
            process.redoPatJets+=process.patJetsAK12Softdrop
            process.redoPatJets+=process.selectedPatJetsAK12Softdrop

    if (doPuppi):

        from CommonTools.PileupAlgos.Puppi_cff import puppi
        from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi

        process.ak8PFJetsPuppi = ak4PFJetsPuppi.clone( rParam = 0.8 )
        process.puppi = puppi.clone( candName = cms.InputTag('packedPFCandidates'),
                             vertexName = cms.InputTag('offlineSlimmedPrimaryVertices'))

        process.puppi_onMiniAOD = cms.Sequence(process.puppi + process.ak8PFJetsPuppi)

        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.NjettinessAK8Puppi = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak8PFJetsPuppi"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.ak4PFJetsPuppi = ak4PFJetsPuppi.clone(src = 'puppi')
        process.ak8PFJetsPuppi = process.ak8PFJetsPuppi.clone( src = 'puppi', jetPtMin = 100.0 )

        process.ak8PFJetsPuppiPruned = ak8PFJetsCHSPruned.clone( src = 'puppi', jetPtMin = 100.0 )
        process.ak8PFJetsPuppiPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            matched = cms.InputTag("ak8PFJetsPuppiPruned"),
            src = cms.InputTag("ak8PFJetsPuppi")
            )
        process.ak8PFJetsPuppiSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'puppi', jetPtMin = 100.0 )
        process.ak8PFJetsPuppiSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            matched = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
            src = cms.InputTag("ak8PFJetsPuppi")
            )

        process.substructureSequence+=process.puppi
        process.substructureSequence+=process.ak8PFJetsPuppi
        process.substructureSequence+=process.NjettinessAK8Puppi

        process.softdrop_onMiniAOD += process.ak8PFJetsPuppiSoftDrop + process.ak8PFJetsPuppiSoftDropMass
        process.pruning_onMiniAOD += process.ak8PFJetsPuppiPruned + process.ak8PFJetsPuppiPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from puppi AK8

        process.puppiJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppi',
                                                                     levels = cms.vstring('L2Relative',
                                                                                          'L3Absolute')
                                                                     )
        process.puppiJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsPuppi' )
        process.puppiJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8PFJetsPuppiPrunedMass"), cms.InputTag("ak8PFJetsPuppiSoftDropMass"), cms.InputTag("NjettinessAK8Puppi:tau1"), cms.InputTag("NjettinessAK8Puppi:tau2"), cms.InputTag("NjettinessAK8Puppi:tau3")]
        process.puppiJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("puppiJetCorrFactorsAK8") )
        process.selectedPuppiJetsAK8 = selectedPatJetsAK8.clone( src = 'puppiJetsAK8', cut = cms.string('pt > 20') )

        process.redoPuppiAK8Jets+=process.puppiJetCorrFactorsAK8
        process.redoPuppiAK8Jets+=process.puppiJetsAK8
        process.redoPuppiAK8Jets+=process.selectedPuppiJetsAK8

        if (reDoPruningAndSoftdropPuppi):
            process.patJetCorrFactorsPuppiAK8Pruned = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppiPruned' )
            process.patJetsPuppiAK8Pruned = patJetsAK8.clone( jetSource = 'ak8PFJetsPuppiPruned' )
            process.patJetsPuppiAK8Pruned.userData.userFloats.src = [ "" ]
            process.patJetsPuppiAK8Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsPuppiAK8Pruned") )
            process.selectedPatJetsPuppiAK8Pruned = selectedPatJetsAK8.clone( 
                src = cms.InputTag('patJetsPuppiAK8Pruned'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsPuppiAK8Pruned
            process.redoPatJets+=process.patJetsPuppiAK8Pruned
            process.redoPatJets+=process.selectedPatJetsPuppiAK8Pruned

            process.patJetCorrFactorsPuppiAK8Softdrop = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppiSoftDrop' )
            process.patJetsPuppiAK8Softdrop = patJetsAK8.clone( jetSource = 'ak8PFJetsPuppiSoftDrop' )
            process.patJetsPuppiAK8Softdrop.userData.userFloats.src = [ "" ]
            process.patJetsPuppiAK8Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsPuppiAK8Softdrop") )
            process.selectedPatJetsPuppiAK8Softdrop = selectedPatJetsAK8.clone(
                src = cms.InputTag('patJetsPuppiAK8Softdrop'),
                cut = cms.string('pt > 20') 
                )

            process.redoPatJets+=process.patJetCorrFactorsPuppiAK8Softdrop
            process.redoPatJets+=process.patJetsPuppiAK8Softdrop
            process.redoPatJets+=process.selectedPatJetsPuppiAK8Softdrop

#######AK8 GEN JETS################

    process.substructureSequenceGen = cms.Sequence()
    process.softdropGen_onMiniAOD = cms.Sequence()
    process.pruningGen_onMiniAOD = cms.Sequence()
    process.redoGenJets = cms.Sequence()
#    process.puppi_onMiniAOD = cms.Sequence()

    if (genJetsAK8Reclustering and MC):    
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        
        process.ak8GenJets = ak4GenJets.clone(src = cms.InputTag('packedGenParticles'),
                                          rParam = cms.double(0.8)
                                          )

        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.NjettinessGenAK8 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak8GenJets"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
                                             src = cms.InputTag("packedGenParticles"),
                                             ignoreParticleIDs = cms.vuint32(
                1000022,
                1000012, 1000014, 1000016,
                2000012, 2000014, 2000016,
                1000039, 5100039,
                4000012, 4000014, 4000016,
                9900012, 9900014, 9900016,
                39),
                                             partonicFinalState = cms.bool(False),
                                             excludeResonances = cms.bool(False),
                                             excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                                             tausAsJets = cms.bool(False)
                                             )

        from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

        process.ak8GenJetsPruned = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(0.8),
            src = cms.InputTag("genParticlesForJets"),
            usePruning = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak8GenJetsPruned = ak8PFJetsCHSPruned.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak8GenJetsPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            matched = cms.InputTag("ak8GenJetsPruned"),
            src = cms.InputTag("ak8GenJets")
            )

        process.ak8GenJetsSoftDrop = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(0.8),
            src = cms.InputTag("genParticlesForJets"),
            useSoftDrop = cms.bool(True),
            R0 = cms.double(0.8),
            beta = cms.double(0.0),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak8GenJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak8GenJetsSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            matched = cms.InputTag("ak8GenJetsSoftDrop"),
            src = cms.InputTag("ak8GenJets")
            )

#        process.substructureSequence+=process.puppi
        process.substructureSequenceGen+=process.genParticlesForJets
        process.substructureSequenceGen+=process.ak8GenJets
        process.substructureSequenceGen+=process.NjettinessGenAK8

        process.softdropGen_onMiniAOD += process.ak8GenJetsSoftDrop + process.ak8GenJetsSoftDropMass
        process.pruningGen_onMiniAOD += process.ak8GenJetsPruned + process.ak8GenJetsPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from gen AK8

#        process.GenJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppi' )
        process.genJetsAK8 = patJetsAK8.clone( jetSource = 'ak8GenJets' )
        process.genJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8GenJetsPrunedMass"), cms.InputTag("ak8GenJetsSoftDropMass"), cms.InputTag("NjettinessGenAK8:tau1"), cms.InputTag("NjettinessGenAK8:tau2"), cms.InputTag("NjettinessGenAK8:tau3")]
        process.genJetsAK8.addJetCorrFactors = cms.bool(False)
        process.genJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("") )
        process.selectedGenJetsAK8 = selectedPatJetsAK8.clone( src = 'genJetsAK8', cut = cms.string('pt > 20') )

#        process.redoGenJets+=process.genJetCorrFactorsAK8
        process.redoGenJets+=process.genJetsAK8
        process.redoGenJets+=process.selectedGenJetsAK8


    if (genJetsAK10Reclustering and MC):    
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        
        process.ak10GenJets = ak4GenJets.clone(src = cms.InputTag('packedGenParticles'),
                                          rParam = cms.double(1.0)
                                          )

        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.NjettinessGenAK10 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak10GenJets"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 1.0 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
                                             src = cms.InputTag("packedGenParticles"),
                                             ignoreParticleIDs = cms.vuint32(
                1000022,
                1000012, 1000014, 1000016,
                2000012, 2000014, 2000016,
                1000039, 5100039,
                4000012, 4000014, 4000016,
                9900012, 9900014, 9900016,
                39),
                                             partonicFinalState = cms.bool(False),
                                             excludeResonances = cms.bool(False),
                                             excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                                             tausAsJets = cms.bool(False)
                                             )

        from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

        process.ak10GenJetsPruned = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(1.0),
            src = cms.InputTag("genParticlesForJets"),
            usePruning = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak10GenJetsPruned = ak8PFJetsCHSPruned.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak10GenJetsPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            matched = cms.InputTag("ak10GenJetsPruned"),
            src = cms.InputTag("ak10GenJets")
            )

        process.ak10GenJetsSoftDrop = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(1.0),
            src = cms.InputTag("genParticlesForJets"),
            useSoftDrop = cms.bool(True),
            R0 = cms.double(1.0),
            beta = cms.double(0.0),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak10GenJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak10GenJetsSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            matched = cms.InputTag("ak10GenJetsSoftDrop"),
            src = cms.InputTag("ak10GenJets")
            )

#        process.substructureSequence+=process.puppi
        process.substructureSequenceGen+=process.genParticlesForJets
        process.substructureSequenceGen+=process.ak10GenJets
        process.substructureSequenceGen+=process.NjettinessGenAK10

        process.softdropGen_onMiniAOD += process.ak10GenJetsSoftDrop + process.ak10GenJetsSoftDropMass
        process.pruningGen_onMiniAOD += process.ak10GenJetsPruned + process.ak10GenJetsPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from gen AK10

#        process.GenJetCorrFactorsAK10 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppi' )
        process.genJetsAK10 = patJetsAK8.clone( jetSource = 'ak10GenJets' )
        process.genJetsAK10.userData.userFloats.src = [ cms.InputTag("ak10GenJetsPrunedMass"), cms.InputTag("ak10GenJetsSoftDropMass"), cms.InputTag("NjettinessGenAK10:tau1"), cms.InputTag("NjettinessGenAK10:tau2"), cms.InputTag("NjettinessGenAK10:tau3")]
        process.genJetsAK10.addJetCorrFactors = cms.bool(False)
        process.genJetsAK10.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("") )
        process.selectedGenJetsAK10 = selectedPatJetsAK8.clone( src = 'genJetsAK10', cut = cms.string('pt > 20') )

#        process.redoGenJets+=process.genJetCorrFactorsAK10
        process.redoGenJets+=process.genJetsAK10
        process.redoGenJets+=process.selectedGenJetsAK10



    if (genJetsAK12Reclustering and MC):    
        from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
        
        process.ak12GenJets = ak4GenJets.clone(src = cms.InputTag('packedGenParticles'),
                                          rParam = cms.double(1.2)
                                          )

        from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

        process.NjettinessGenAK12 = cms.EDProducer("NjettinessAdder",
                            src=cms.InputTag("ak12GenJets"),
                            Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                            # variables for measure definition : 
                            measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                            beta = cms.double(1.0),              # CMS default is 1
                            R0 = cms.double( 1.2 ),              # CMS default is jet cone size
                            Rcutoff = cms.double( 999.0),       # not used by default
                            # variables for axes definition :
                            axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                            nPass = cms.int32(999),             # not used by default
                            akAxesR0 = cms.double(999.0)        # not used by default
                            )

        process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
                                             src = cms.InputTag("packedGenParticles"),
                                             ignoreParticleIDs = cms.vuint32(
                1000022,
                1000012, 1000014, 1000016,
                2000012, 2000014, 2000016,
                1000039, 5100039,
                4000012, 4000014, 4000016,
                9900012, 9900014, 9900016,
                39),
                                             partonicFinalState = cms.bool(False),
                                             excludeResonances = cms.bool(False),
                                             excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                                             tausAsJets = cms.bool(False)
                                             )

        from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

        process.ak12GenJetsPruned = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(1.2),
            src = cms.InputTag("genParticlesForJets"),
            usePruning = cms.bool(True),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak12GenJetsPruned = ak8PFJetsCHSPruned.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak12GenJetsPrunedMass = ak8PFJetsCHSPrunedMass.clone(    
            matched = cms.InputTag("ak12GenJetsPruned"),
            src = cms.InputTag("ak12GenJets")
            )

        process.ak12GenJetsSoftDrop = ak4GenJets.clone(
            SubJetParameters,
            rParam = cms.double(1.2),
            src = cms.InputTag("genParticlesForJets"),
            useSoftDrop = cms.bool(True),
            R0 = cms.double(1.2),
            beta = cms.double(0.0),
            writeCompound = cms.bool(True),
            jetCollInstanceName=cms.string("SubJets")
            )

#        process.ak12GenJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'packedGenParticles', jetPtMin = 100.0 )
        process.ak12GenJetsSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
            matched = cms.InputTag("ak12GenJetsSoftDrop"),
            src = cms.InputTag("ak12GenJets")
            )

#        process.substructureSequence+=process.puppi
        process.substructureSequenceGen+=process.genParticlesForJets
        process.substructureSequenceGen+=process.ak12GenJets
        process.substructureSequenceGen+=process.NjettinessGenAK12

        process.softdropGen_onMiniAOD += process.ak12GenJetsSoftDrop + process.ak12GenJetsSoftDropMass
        process.pruningGen_onMiniAOD += process.ak12GenJetsPruned + process.ak12GenJetsPrunedMass

        ####### Redo pat jets sequence ##########
        from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

        # Redo pat jets from gen AK12

#        process.GenJetCorrFactorsAK12 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsPuppi' )
        process.genJetsAK12 = patJetsAK8.clone( jetSource = 'ak12GenJets' )
        process.genJetsAK12.userData.userFloats.src = [ cms.InputTag("ak12GenJetsPrunedMass"), cms.InputTag("ak12GenJetsSoftDropMass"), cms.InputTag("NjettinessGenAK12:tau1"), cms.InputTag("NjettinessGenAK12:tau2"), cms.InputTag("NjettinessGenAK12:tau3")]
        process.genJetsAK12.addJetCorrFactors = cms.bool(False)
        process.genJetsAK12.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("") )
        process.selectedGenJetsAK12 = selectedPatJetsAK8.clone( src = 'genJetsAK12', cut = cms.string('pt > 20') )

#        process.redoGenJets+=process.genJetCorrFactorsAK12
        process.redoGenJets+=process.genJetsAK12
        process.redoGenJets+=process.selectedGenJetsAK12


######### A4PF-nonCHS jets ###########

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

    process.ak4PFJets = ak4PFJets.clone(src = "packedPFCandidates")

    from AllHadronicSUSY.Utils.ak4pfjets_cfi import patJetCorrFactorsAK4, patJetsAK4

    process.patJetCorrFactorsAK4 = patJetCorrFactorsAK4.clone( src = 'ak4PFJets' )
    process.patJetsAK4 = patJetsAK4.clone( jetSource = 'ak4PFJets' )

    #
    # Set up electron ID (VID framework)
    #
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    # turn on VID producer, indicate data format to be
    # DataFormat.AOD or DataFormat.MiniAOD, as appropriate
#    if useAOD == True :
#        dataFormat = DataFormat.AOD
#    else :
#        dataFormat = DataFormat.MiniAOD 
    dataFormat=DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process,dataFormat)

    process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

    # define which IDs we want to produce
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
    #add them to the VID producer
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

    # Producers
    from AllHadronicSUSY.Utils.electron_cfi import electron
    process.Electrons = electron.clone(
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        EleTag = cms.InputTag("slimmedElectrons"),
        MinPt = cms.double(-1),
        RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
        eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
        eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
#        eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51")
        eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60")
    )
    from AllHadronicSUSY.Utils.muon_cfi import muon
    process.Muons = muon.clone(
        VertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
        MuTag = cms.InputTag("slimmedMuons"),
        MinPt = cms.double(-1),
        RhoTag = cms.InputTag("fixedGridRhoFastjetAll")
    )
    from AllHadronicSUSY.Utils.subJetSelection_cfi import SubJetSelection
    process.HTJets = SubJetSelection.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt								  = cms.double(50),
    MaxEta								  = cms.double(2.5),
    )
    from AllHadronicSUSY.Utils.htdouble_cfi import htdouble
    process.HT = htdouble.clone(
    JetTag  = cms.InputTag('HTJets'),
    )
    from AllHadronicSUSY.Utils.njetint_cfi import njetint
    process.NJets = njetint.clone(
#    JetTag  = cms.InputTag('HTJets'),
    JetTag  = cms.InputTag('slimmedJets'),
    )
    from AllHadronicSUSY.Utils.btagint_cfi import btagint
    process.BTags = btagint.clone(
#    JetTag  = cms.InputTag('HTJets'),
    JetTag  = cms.InputTag('slimmedJets'),
    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    BTagCutValue					= cms.double(0.679)
    )


#    JECstring="Spring16_25nsV10"
    JECstring="Spring16_25nsV3"
    if (isRun2016BCD):
        JECstring="Spring16_25nsV10BCD"
    elif (isRun2016E):
        JECstring="Spring16_25nsV10E"
    elif (isRun2016F):
        JECstring="Spring16_25nsV10F"
    elif (isRun2016GH):
        JECstring="Spring16_25nsV10p2"
    

    L1string_AK4_data=JECstring+"_DATA_L1FastJet_AK4PFchs.txt"
    L2string_AK4_data=JECstring+"_DATA_L2Relative_AK4PFchs.txt"
    L3string_AK4_data=JECstring+"_DATA_L3Absolute_AK4PFchs.txt"
    L2L3string_AK4_data=JECstring+"_DATA_L2L3Residual_AK4PFchs.txt"
    Uncstring_AK4_data=JECstring+"_DATA_Uncertainty_AK4PFchs.txt"

    L1string_AK8_data=JECstring+"_DATA_L1FastJet_AK8PFchs.txt"
    L2string_AK8_data=JECstring+"_DATA_L2Relative_AK8PFchs.txt"
    L3string_AK8_data=JECstring+"_DATA_L3Absolute_AK8PFchs.txt"
    L2L3string_AK8_data=JECstring+"_DATA_L2L3Residual_AK8PFchs.txt"
    Uncstring_AK8_data=JECstring+"_DATA_Uncertainty_AK8PFchs.txt"

    L1string_AK4_mc=JECstring+"_MC_L1FastJet_AK4PFchs.txt"
    L2string_AK4_mc=JECstring+"_MC_L2Relative_AK4PFchs.txt"
    L3string_AK4_mc=JECstring+"_MC_L3Absolute_AK4PFchs.txt"
    L2L3string_AK4_mc=JECstring+"_MC_L2L3Residual_AK4PFchs.txt"
    Uncstring_AK4_mc=JECstring+"_MC_Uncertainty_AK4PFchs.txt"

    L1string_AK8_mc=JECstring+"_MC_L1FastJet_AK8PFchs.txt"
    L2string_AK8_mc=JECstring+"_MC_L2Relative_AK8PFchs.txt"
    L3string_AK8_mc=JECstring+"_MC_L3Absolute_AK8PFchs.txt"
    L2L3string_AK8_mc=JECstring+"_MC_L2L3Residual_AK8PFchs.txt"
    Uncstring_AK8_mc=JECstring+"_MC_Uncertainty_AK8PFchs.txt"

    L1string_puppi_AK4_data=JECstring+"_DATA_L1FastJet_AK4PFPuppi.txt"
    L2string_puppi_AK4_data=JECstring+"_DATA_L2Relative_AK4PFPuppi.txt"
    L3string_puppi_AK4_data=JECstring+"_DATA_L3Absolute_AK4PFPuppi.txt"
    L2L3string_puppi_AK4_data=JECstring+"_DATA_L2L3Residual_AK4PFPuppi.txt"
    Uncstring_puppi_AK4_data=JECstring+"_DATA_Uncertainty_AK4PFPuppi.txt"

    L1string_puppi_AK8_data=JECstring+"_DATA_L1FastJet_AK8PFPuppi.txt"
    L2string_puppi_AK8_data=JECstring+"_DATA_L2Relative_AK8PFPuppi.txt"
    L3string_puppi_AK8_data=JECstring+"_DATA_L3Absolute_AK8PFPuppi.txt"
    L2L3string_puppi_AK8_data=JECstring+"_DATA_L2L3Residual_AK8PFPuppi.txt"
    Uncstring_puppi_AK8_data=JECstring+"_DATA_Uncertainty_AK8PFPuppi.txt"

    L1string_puppi_AK4_mc=JECstring+"_MC_L1FastJet_AK4PFPuppi.txt"
    L2string_puppi_AK4_mc=JECstring+"_MC_L2Relative_AK4PFPuppi.txt"
    L3string_puppi_AK4_mc=JECstring+"_MC_L3Absolute_AK4PFPuppi.txt"
    L2L3string_puppi_AK4_mc=JECstring+"_MC_L2L3Residual_AK4PFPuppi.txt"
    Uncstring_puppi_AK4_mc=JECstring+"_MC_Uncertainty_AK4PFPuppi.txt"

    L1string_puppi_AK8_mc=JECstring+"_MC_L1FastJet_AK8PFPuppi.txt"
    L2string_puppi_AK8_mc=JECstring+"_MC_L2Relative_AK8PFPuppi.txt"
    L3string_puppi_AK8_mc=JECstring+"_MC_L3Absolute_AK8PFPuppi.txt"
    L2L3string_puppi_AK8_mc=JECstring+"_MC_L2L3Residual_AK8PFPuppi.txt"
    Uncstring_puppi_AK8_mc=JECstring+"_MC_Uncertainty_AK8PFPuppi.txt"


    from AllHadronicSUSY.Utils.jetproperties_cfi import jetproperties
    process.JetsProperties = jetproperties.clone(
    JetTag  = cms.InputTag('slimmedJets'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string(L1string_AK4_data),
    L2File = cms.string(L2string_AK4_data),
    L3File = cms.string(L3string_AK4_data),
    L2L3File = cms.string(L2L3string_AK4_data),
    uncFile = cms.string(Uncstring_AK4_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK4),
#    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    )
    from AllHadronicSUSY.Utils.jetpropertiesAK8_cfi import jetpropertiesAK8
    process.JetsPropertiesAK8 = jetpropertiesAK8.clone(
    JetTag  = cms.InputTag('slimmedJetsAK8'),
#    puppiJetTag = cms.InputTag('selectedPuppiJetsAK8'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    doReclusteringForPrunedAndSoftdrop = cms.bool(reDoPruningAndSoftdrop),
    L1File = cms.string(L1string_AK8_data),
    L2File = cms.string(L2string_AK8_data),
    L3File = cms.string(L3string_AK8_data),
    L2L3File = cms.string(L2L3string_AK8_data),
    uncFile = cms.string(Uncstring_AK8_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK8),
#    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    )
    from AllHadronicSUSY.Utils.jetpropertiesAK10_cfi import jetpropertiesAK10
    process.JetsPropertiesAK10 = jetpropertiesAK10.clone(
    JetTag  = cms.InputTag('selectedPatJetsAK10'),
#    puppiJetTag = cms.InputTag('selectedPuppiJetsAK8'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    doReclusteringForPrunedAndSoftdrop = cms.bool(reDoPruningAndSoftdrop),
    L1File = cms.string(L1string_AK8_data),
    L2File = cms.string(L2string_AK8_data),
    L3File = cms.string(L3string_AK8_data),
    L2L3File = cms.string(L2L3string_AK8_data),
    uncFile = cms.string(Uncstring_AK8_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK8),
#    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    )
    from AllHadronicSUSY.Utils.jetpropertiesAK12_cfi import jetpropertiesAK12
    process.JetsPropertiesAK12 = jetpropertiesAK12.clone(
    JetTag  = cms.InputTag('selectedPatJetsAK12'),
#    puppiJetTag = cms.InputTag('selectedPuppiJetsAK8'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    doReclusteringForPrunedAndSoftdrop = cms.bool(reDoPruningAndSoftdrop),
    L1File = cms.string(L1string_AK8_data),
    L2File = cms.string(L2string_AK8_data),
    L3File = cms.string(L3string_AK8_data),
    L2L3File = cms.string(L2L3string_AK8_data),
    uncFile = cms.string(Uncstring_AK8_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK8),
#    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    )


    from AllHadronicSUSY.Utils.jetpropertiesPuppiAK8_cfi import jetpropertiesPuppiAK8
#    process.JetsPropertiesPuppiAK8 = cms.Sequence()
#    if doPuppi:
    process.JetsPropertiesPuppiAK8 = jetpropertiesPuppiAK8.clone(
        JetTag  = cms.InputTag('selectedPuppiJetsAK8'),
        #    puppiJetTag = cms.InputTag('selectedPuppiJetsAK8'),
        MinPt = cms.double(-1),
        doJEC  = cms.bool(doJECCorrection),
        doReclusteringForPrunedAndSoftdrop = cms.bool(reDoPruningAndSoftdropPuppi),
        L1File = cms.string(L1string_puppi_AK8_data),
        L2File = cms.string(L2string_puppi_AK8_data),
        L3File = cms.string(L3string_puppi_AK8_data),
        L2L3File = cms.string(L2L3string_puppi_AK8_data),
        uncFile = cms.string(Uncstring_puppi_AK8_data),
        #    jecPayloadNames      = cms.vstring(jecLevelsAK8),
        #    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
        )
    #        process.JetsPropertiesAK8.JetTag = cms.InputTag('selectedPuppiJetsAK8')
    from AllHadronicSUSY.Utils.jetpropertiesPuppi_cfi import jetpropertiesPuppi
    process.JetsPropertiesPuppi = jetpropertiesPuppi.clone(
    JetTag  = cms.InputTag('slimmedJetsPuppi'),
    MinPt = cms.double(-1),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string(L1string_puppi_AK4_data),
    L2File = cms.string(L2string_puppi_AK4_data),
    L3File = cms.string(L3string_puppi_AK4_data),
    L2L3File = cms.string(L2L3string_puppi_AK4_data),
    uncFile = cms.string(Uncstring_puppi_AK4_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK4),
#    BTagInputTag	        = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
    )

    if (reDoPruningAndSoftdrop):
        process.JetsPropertiesAK8.prunedJetTag  = cms.InputTag('selectedPatJetsAK8Pruned')
        process.JetsPropertiesAK8.softdropJetTag  = cms.InputTag('selectedPatJetsAK8Softdrop')
        process.JetsPropertiesAK10.prunedJetTag  = cms.InputTag('selectedPatJetsAK10Pruned')
        process.JetsPropertiesAK10.softdropJetTag  = cms.InputTag('selectedPatJetsAK10Softdrop')
        process.JetsPropertiesAK12.prunedJetTag  = cms.InputTag('selectedPatJetsAK12Pruned')
        process.JetsPropertiesAK12.softdropJetTag  = cms.InputTag('selectedPatJetsAK12Softdrop')
    else:
        process.JetsPropertiesAK8.prunedJetTag  = cms.InputTag('slimmedJetsAK8')
        process.JetsPropertiesAK8.softdropJetTag  = cms.InputTag('slimmedJetsAK8')
        process.JetsPropertiesAK10.prunedJetTag  = cms.InputTag('slimmedJetsAK8')
        process.JetsPropertiesAK10.softdropJetTag  = cms.InputTag('slimmedJetsAK8')
        process.JetsPropertiesAK12.prunedJetTag  = cms.InputTag('slimmedJetsAK8')
        process.JetsPropertiesAK12.softdropJetTag  = cms.InputTag('slimmedJetsAK8')

    if (reDoPruningAndSoftdropPuppi):
        process.JetsPropertiesPuppiAK8.prunedJetTag  = cms.InputTag('selectedPatJetsPuppiAK8Pruned')
        process.JetsPropertiesPuppiAK8.softdropJetTag  = cms.InputTag('selectedPatJetsPuppiAK8Softdrop')
    else:
        process.JetsPropertiesPuppiAK8.prunedJetTag  = cms.InputTag('slimmedJetsPuppiAK8')
        process.JetsPropertiesPuppiAK8.softdropJetTag  = cms.InputTag('slimmedJetsPuppiAK8')

    if (MC):
        process.JetsProperties.L1File = cms.string(L1string_AK4_mc)
        process.JetsProperties.L2File = cms.string(L2string_AK4_mc)
        process.JetsProperties.L3File = cms.string(L3string_AK4_mc)
        process.JetsProperties.L2L3File = cms.string("NONE")
        process.JetsProperties.uncFile = cms.string(Uncstring_AK4_mc)
        process.JetsPropertiesAK8.L1File = cms.string(L1string_AK8_mc)
        process.JetsPropertiesAK8.L2File = cms.string(L2string_AK8_mc)
        process.JetsPropertiesAK8.L3File = cms.string(L3string_AK8_mc)
        process.JetsPropertiesAK8.L2L3File = cms.string("NONE")
        process.JetsPropertiesAK8.uncFile = cms.string(Uncstring_AK8_mc)
        process.JetsPropertiesAK10.L1File = cms.string(L1string_AK8_mc)
        process.JetsPropertiesAK10.L2File = cms.string(L2string_AK8_mc)
        process.JetsPropertiesAK10.L3File = cms.string(L3string_AK8_mc)
        process.JetsPropertiesAK10.L2L3File = cms.string("NONE")
        process.JetsPropertiesAK10.uncFile = cms.string(Uncstring_AK8_mc)
        process.JetsPropertiesAK12.L1File = cms.string(L1string_AK8_mc)
        process.JetsPropertiesAK12.L2File = cms.string(L2string_AK8_mc)
        process.JetsPropertiesAK12.L3File = cms.string(L3string_AK8_mc)
        process.JetsPropertiesAK12.L2L3File = cms.string("NONE")
        process.JetsPropertiesAK12.uncFile = cms.string(Uncstring_AK8_mc)
        process.JetsPropertiesPuppiAK8.L1File = cms.string(L1string_puppi_AK8_mc)
        process.JetsPropertiesPuppiAK8.L2File = cms.string(L2string_puppi_AK8_mc)
        process.JetsPropertiesPuppiAK8.L3File = cms.string(L3string_puppi_AK8_mc)
        process.JetsPropertiesPuppiAK8.L2L3File = cms.string("NONE")
        process.JetsPropertiesPuppiAK8.uncFile = cms.string(Uncstring_puppi_AK8_mc)
        process.JetsPropertiesPuppi.L1File = cms.string(L1string_puppi_AK4_mc)
        process.JetsPropertiesPuppi.L2File = cms.string(L2string_puppi_AK4_mc)
        process.JetsPropertiesPuppi.L3File = cms.string(L3string_puppi_AK4_mc)
        process.JetsPropertiesPuppi.L2L3File = cms.string("NONE")
        process.JetsPropertiesPuppi.uncFile = cms.string(Uncstring_puppi_AK4_mc)

    if doAK8Reclustering:
        process.JetsPropertiesAK8.JetTag = cms.InputTag('selectedPatJetsAK8')

    from AllHadronicSUSY.Utils.metdouble_cfi import metdouble
    process.MET = metdouble.clone(
    METTag  = cms.InputTag("slimmedMETs"),
    JetTag  = cms.InputTag('slimmedJets'),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string(L1string_AK4_data),
    L2File = cms.string(L2string_AK4_data),
    L3File = cms.string(L3string_AK4_data),
    L2L3File = cms.string(L2L3string_AK4_data),
    uncFile = cms.string(Uncstring_AK4_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK4),
    MuTag = cms.InputTag("slimmedMuons"),
    RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    corrMet = cms.bool(doJECCorrection),
#    corrMet = cms.bool(True),
    )

    from AllHadronicSUSY.Utils.metpuppidouble_cfi import metpuppidouble
    process.METpuppi = metpuppidouble.clone(
    METTag  = cms.InputTag("slimmedMETsPuppi"),
    JetTag  = cms.InputTag('slimmedJetsPuppi'),
    doJEC  = cms.bool(doJECCorrection),
    L1File = cms.string(L1string_puppi_AK4_data),
    L2File = cms.string(L2string_puppi_AK4_data),
    L3File = cms.string(L3string_puppi_AK4_data),
    L2L3File = cms.string(L2L3string_puppi_AK4_data),
    uncFile = cms.string(Uncstring_puppi_AK4_data),
#    jecPayloadNames      = cms.vstring(jecLevelsAK4),
    MuTag = cms.InputTag("slimmedMuons"),
    RhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
    corrMet = cms.bool(doJECCorrection),
#    corrMet = cms.bool(True),
    )

    if (MC):
        process.MET.L1File = cms.string(L1string_AK4_mc)
        process.MET.L2File = cms.string(L2string_AK4_mc)
        process.MET.L3File = cms.string(L3string_AK4_mc)
        process.MET.L2L3File = cms.string("NONE")
        process.MET.uncFile = cms.string(Uncstring_AK4_mc)
        process.METpuppi.L1File = cms.string(L1string_puppi_AK4_mc)
        process.METpuppi.L2File = cms.string(L2string_puppi_AK4_mc)
        process.METpuppi.L3File = cms.string(L3string_puppi_AK4_mc)
        process.METpuppi.L2L3File = cms.string("NONE")
        process.METpuppi.uncFile = cms.string(Uncstring_puppi_AK4_mc)


    from AllHadronicSUSY.Utils.leptonint_cfi import leptonint
    process.Leptons = leptonint.clone(
    LeptonTag = cms.VInputTag(cms.InputTag('selectedIDIsoMuons'),cms.InputTag('selectedIDIsoElectrons')),
    )
    from AllHadronicSUSY.Utils.primaryverticies_cfi import primaryverticies
    process.NVtx = primaryverticies.clone(
    VertexCollection  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    )
    from AllHadronicSUSY.Utils.genLeptonRecoCand_cfi import genLeptonRecoCand
    process.GenLeptons = genLeptonRecoCand.clone(
    PrunedGenParticleTag  = cms.InputTag("prunedGenParticles"),
    )
    from AllHadronicSUSY.Utils.genJet_cfi import genJet
    process.GenJets = genJet.clone(
                            GenJetCollTag  = cms.InputTag("slimmedGenJets"),
                            )
    from AllHadronicSUSY.Utils.genJetAK8_cfi import genJetAK8
    process.GenJetsAK8 = genJetAK8.clone(
                            GenJetCollTag  = cms.InputTag("selectedGenJetsAK8"),
                            )
    from AllHadronicSUSY.Utils.genJetAK10_cfi import genJetAK10
    process.GenJetsAK10 = cms.Sequence()
    process.GenJetsAK12 = cms.Sequence()
    if genJetsAK10Reclustering:
        process.GenJetsAK10 = genJetAK10.clone(
            GenJetCollTag  = cms.InputTag("selectedGenJetsAK10"),
            )
    from AllHadronicSUSY.Utils.genJetAK12_cfi import genJetAK12
    if genJetsAK12Reclustering:
        process.GenJetsAK12 = genJetAK12.clone(
            GenJetCollTag  = cms.InputTag("selectedGenJetsAK12"),
            )
    if not MC:
        process.GenLeptons = cms.Sequence()
        process.GenJets = cms.Sequence()
        process.GenJetsAK8 = cms.Sequence()
        process.GenJetsAK10 = cms.Sequence()
        process.GenJetsAK12 = cms.Sequence()



    ##### MET filters #####

    #### -----> MET Filter Flags from MiniAOD/TWiki <----- ####
    import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
    process.metBits_miniAOD = hlt.triggerResultsFilter.clone()
    # default is to use the latest process (but can set different process through Commandline Args)
    process.metBits_miniAOD.hltResults = cms.InputTag('TriggerResults::%s'%METFiltersProcess) 
    process.metBits_miniAOD.l1tResults = cms.InputTag('')
    #currently configured for CSCTightHaloFilter + GoodVertices
    met_bits = ['(Flag_CSCTightHaloFilter)','(Flag_goodVertices)','(Flag_eeBadScFilter)']
    bitsexpr = ' AND '.join(met_bits)
    process.metBits_miniAOD.triggerConditions = cms.vstring(bitsexpr)

    #### -----> HBHE noise filter <----- ####
##___________________________HCAL_Noise_Filter________________________________||
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#    if customizeHBHENoiseForRun2015D:
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
#                                                        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#                                                        reverseDecision = cms.bool(False)
#                                                        )
    
#    process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
#                                                           inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
#                                                           reverseDecision = cms.bool(False)
#                                                           )

#    noiseFilterSelection_HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
#    noiseFilterSelection_HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
#    noiseFilterSelection_HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),    

#### obsolete - for early run 2 data (runB)
#    if customizeHBHENoiseForEarlyData:
#        process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    
#    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter(
#        'BooleanFlagFilter',
#        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResultRun2Loose'), #for 25ns data-taking
#        inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'), #for 50ns data-taking
#        reverseDecision = cms.bool(False)
#        )

#    process.metFilters = cms.Sequence()
#    if customizeHBHENoiseForRun2015D and not MC:
#    process.metFilters = cms.Sequence(process.metBits_miniAOD
#                                          *process.HBHENoiseFilterResultProducer
#                                          *process.ApplyBaselineHBHENoiseFilter
#                                          *process.ApplyBaselineHBHEIsoNoiseFilter
#                                      )

########## save flags for filters
    from AllHadronicSUSY.Utils.filterproducer_cfi import filterProducer
    from RecoMET.METFilters.BadChargedCandidateFilter_cfi import BadChargedCandidateFilter
    from RecoMET.METFilters.BadPFMuonFilter_cfi import BadPFMuonFilter
    process.FilterProducer = filterProducer.clone(
                                 #triggerNameList = cms.vstring(),
                                 noiseFilterTag = cms.InputTag("TriggerResults"),
#                                 noiseFilterTag = cms.InputTag('TriggerResults','', 'PAT'),
                                 BadChCandNoiseFilterTag = cms.InputTag("badChargedCandidateFilter"),
                                 BadPFMuonNoiseFilterTag = cms.InputTag("badPFMuonFilter"),
                                 HBHENoiseFilter_Selector_ = cms.string("Flag_HBHENoiseFilter"),
                                 HBHENoiseIsoFilter_Selector_ = cms.string("Flag_HBHENoiseIsoFilter"),
                                 CSCHaloNoiseFilter_Selector_ = cms.string("Flag_CSCTightHaloFilter"),
                                 GoodVtxNoiseFilter_Selector_ = cms.string("Flag_goodVertices"),
                                 EEBadScNoiseFilter_Selector_ = cms.string("Flag_eeBadScFilter"),
                                 EcalDeadCellTriggerPrimitiveFilter_Selector_ = cms.string("Flag_EcalDeadCellTriggerPrimitiveFilter"),
                                 GlobalTightHalo2016Filter_Selector_ = cms.string("Flag_globalTightHalo2016Filter"),
                                 HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
                                 HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
                                 HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult")
                                 )


    process.badChargedCandidateFilter = BadChargedCandidateFilter.clone(
        PFCandidates  = cms.InputTag("packedPFCandidates"),   # Collection to test
        muons         = cms.InputTag("slimmedMuons"),   # Collection to test
        taggingMode   = cms.bool(False),
        debug         = cms.bool(False),
        maxDR         = cms.double(0.01),               # Maximum DR between reco::muon->innerTrack and pfCandidate 
        minPtErrorRel = cms.double(-0.5),               # lower threshold on difference between pt of reco::muon->innerTrack and pfCandidate
        minMuonPt     = cms.double(20),                 # minimum muon pt 
        )

    process.badPFMuonFilter = BadPFMuonFilter.clone(
        PFCandidates  = cms.InputTag("packedPFCandidates"),   # Collection to test 
        muons         = cms.InputTag("slimmedMuons"),   # Collection to test
        taggingMode   = cms.bool(False),
        debug         = cms.bool(False),
        minDZ         = cms.double(0.1),              # dz threshold on PF muons to consider
        minMuPt       = cms.double(100),               # pt threshold on PF muons 
        minTrkPtError  = cms.double(-1),               # threshold on inner track pt Error 
        )

    RecoCandVector = cms.vstring()
#    RecoCandVector.extend(['IsolatedTracks']) # basic muons electrons and isoalted tracks
#    RecoCandVector.extend(['selectedIDIsoMuons','selectedIDIsoElectrons','IsolatedTracks']) # basic muons electrons and isoalted tracks
#    RecoCandVector.extend(['selectedIDMuons','selectedIDElectrons']) # mu and e no isolation cuts
    RecoCandVector.extend(['GenLeptons:Boson(GenBoson)|GenLeptons:BosonPDGId(I_GenBosonPDGId)|GenLeptons:isBosonLeptonic(I_isBosonLeptonic)','GenLeptons:Muon(GenMu)|GenLeptons:MuonTauDecay(I_GenMuFromTau)' ,'GenLeptons:Electron(GenElec)|GenLeptons:ElectronTauDecay(I_GenElecFromTau)','GenLeptons:Tau(GenTau)|GenLeptons:TauHadronic(I_GenTauHad)','GenLeptons:Neutrino(GenNu)','GenLeptons:Top(GenTop)|GenLeptons:TopPDGId(I_GenTopPDGId)',] ) # gen information on leptons
    RecoCandVector.extend(['GenJets:GenJet(GenJets)'] ) # gen information on jets
    RecoCandVector.extend(['GenJetsAK8:GenJetAK8(GenJetsAK8)|GenJetsAK8:GenAK8prunedMass(F_prunedMass)|GenJetsAK8:GenAK8softdropMass(F_softdropMass)|GenJetsAK8:GenAK8softdropPt(F_softdropPt)|GenJetsAK8:GenAK8tau1(F_tau1)|GenJetsAK8:GenAK8tau2(F_tau2)|GenJetsAK8:GenAK8tau3(F_tau3)'] ) # gen information on AK8 jets
#    RecoCandVector.extend(['GenJetsAK10:GenJetAK10(GenJetsAK10)|GenJetsAK10:GenAK10prunedMass(F_prunedMass)|GenJetsAK10:GenAK10softdropMass(F_softdropMass)|GenJetsAK10:GenAK10softdropPt(F_softdropPt)|GenJetsAK10:GenAK10tau1(F_tau1)|GenJetsAK10:GenAK10tau2(F_tau2)|GenJetsAK10:GenAK10tau3(F_tau3)'] ) # gen information on AK10 jets
#    RecoCandVector.extend(['GenJetsAK12:GenJetAK12(GenJetsAK12)|GenJetsAK12:GenAK12prunedMass(F_prunedMass)|GenJetsAK12:GenAK12softdropMass(F_softdropMass)|GenJetsAK12:GenAK12softdropPt(F_softdropPt)|GenJetsAK12:GenAK12tau1(F_tau1)|GenJetsAK12:GenAK12tau2(F_tau2)|GenJetsAK12:GenAK12tau3(F_tau3)'] ) # gen information on AK12 jets
    RecoCandVector.extend(['JetsProperties(Jets)|JetsProperties:bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsProperties:bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsProperties:chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsProperties:chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsProperties:chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsProperties:electronMultiplicity(I_electronMultiplicity)|JetsProperties:jetArea(F_jetArea)|JetsProperties:muonEnergyFraction(F_muonEnergyFraction)|JetsProperties:muonMultiplicity(I_muonMultiplicity)|JetsProperties:neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsProperties:neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsProperties:photonEnergyFraction(F_photonEnergyFraction)|JetsProperties:photonMultiplicity(I_photonMultiplicity)|JetsProperties:isLooseJetId(b_isLooseJetId)|JetsProperties:isTightJetId(b_isTightJetId)|JetsProperties:PtCorr(F_PtCorr)|JetsProperties:EtaCorr(F_EtaCorr)|JetsProperties:PhiCorr(F_PhiCorr)|JetsProperties:ECorr(F_ECorr)|JetsProperties:AK4correction(F_AK4correction)|JetsProperties:AK4correctionUp(F_AK4correctionUp)|JetsProperties:AK4correctionDown(F_AK4correctionDown)'] ) # jet information on various variables
    RecoCandVector.extend(['JetsPropertiesAK8(AK8Jets)|JetsPropertiesAK8:AK8bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsPropertiesAK8:AK8bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsPropertiesAK8:AK8chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsPropertiesAK8:AK8chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsPropertiesAK8:AK8chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsPropertiesAK8:AK8electronMultiplicity(I_electronMultiplicity)|JetsPropertiesAK8:AK8jetArea(F_jetArea)|JetsPropertiesAK8:AK8muonEnergyFraction(F_muonEnergyFraction)|JetsPropertiesAK8:AK8muonMultiplicity(I_muonMultiplicity)|JetsPropertiesAK8:AK8neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsPropertiesAK8:AK8neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsPropertiesAK8:AK8photonEnergyFraction(F_photonEnergyFraction)|JetsPropertiesAK8:AK8photonMultiplicity(I_photonMultiplicity)|JetsPropertiesAK8:AK8prunedMass(F_prunedMass)|JetsPropertiesAK8:AK8softDropMass(F_softDropMass)|JetsPropertiesAK8:AK8softDropPt(F_softDropPt)|JetsPropertiesAK8:AK8trimmedMass(F_trimmedMass)|JetsPropertiesAK8:AK8filteredMass(F_filteredMass)|JetsPropertiesAK8:AK8tau1(F_tau1)|JetsPropertiesAK8:AK8tau2(F_tau2)|JetsPropertiesAK8:AK8tau3(F_tau3)|JetsPropertiesAK8:AK8isLooseJetId(b_AK8isLooseJetId)|JetsPropertiesAK8:AK8isTightJetId(b_AK8isTightJetId)|JetsPropertiesAK8:PtCorr(F_PtCorr)|JetsPropertiesAK8:EtaCorr(F_EtaCorr)|JetsPropertiesAK8:PhiCorr(F_PhiCorr)|JetsPropertiesAK8:ECorr(F_ECorr)|JetsPropertiesAK8:mass(F_mass)|JetsPropertiesAK8:AK8correction(F_AK8correction)|JetsPropertiesAK8:AK8correctionUp(F_AK8correctionUp)|JetsPropertiesAK8:AK8correctionDown(F_AK8correctionDown)|JetsPropertiesAK8:AK8massCorrection(F_AK8massCorrection)|JetsPropertiesAK8:AK8massCorrectionUp(F_AK8massCorrectionUp)|JetsPropertiesAK8:AK8massCorrectionDown(F_AK8massCorrectionDown)'] ) # AK8 jet information on various variables
# 
#    RecoCandVector.extend(['JetsPropertiesAK10(AK10Jets)|JetsPropertiesAK10:AK10bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsPropertiesAK10:AK10bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsPropertiesAK10:AK10chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsPropertiesAK10:AK10chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsPropertiesAK10:AK10chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsPropertiesAK10:AK10electronMultiplicity(I_electronMultiplicity)|JetsPropertiesAK10:AK10jetArea(F_jetArea)|JetsPropertiesAK10:AK10muonEnergyFraction(F_muonEnergyFraction)|JetsPropertiesAK10:AK10muonMultiplicity(I_muonMultiplicity)|JetsPropertiesAK10:AK10neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsPropertiesAK10:AK10neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsPropertiesAK10:AK10photonEnergyFraction(F_photonEnergyFraction)|JetsPropertiesAK10:AK10photonMultiplicity(I_photonMultiplicity)|JetsPropertiesAK10:AK10prunedMass(F_prunedMass)|JetsPropertiesAK10:AK10softDropMass(F_softDropMass)|JetsPropertiesAK10:AK10softDropPt(F_softDropPt)|JetsPropertiesAK10:AK10trimmedMass(F_trimmedMass)|JetsPropertiesAK10:AK10filteredMass(F_filteredMass)|JetsPropertiesAK10:AK10tau1(F_tau1)|JetsPropertiesAK10:AK10tau2(F_tau2)|JetsPropertiesAK10:AK10tau3(F_tau3)|JetsPropertiesAK10:AK10isLooseJetId(b_AK10isLooseJetId)|JetsPropertiesAK10:AK10isTightJetId(b_AK10isTightJetId)|JetsPropertiesAK10:PtCorr(F_PtCorr)|JetsPropertiesAK10:EtaCorr(F_EtaCorr)|JetsPropertiesAK10:PhiCorr(F_PhiCorr)|JetsPropertiesAK10:ECorr(F_ECorr)|JetsPropertiesAK10:mass(F_mass)|JetsPropertiesAK10:AK10correction(F_AK10correction)|JetsPropertiesAK10:AK10correctionUp(F_AK10correctionUp)|JetsPropertiesAK10:AK10correctionDown(F_AK10correctionDown)|JetsPropertiesAK10:AK10massCorrection(F_AK10massCorrection)|JetsPropertiesAK10:AK10massCorrectionUp(F_AK10massCorrectionUp)|JetsPropertiesAK10:AK10massCorrectionDown(F_AK10massCorrectionDown)'] ) # AK10 jet information on various variables
# 
#    RecoCandVector.extend(['JetsPropertiesAK12(AK12Jets)|JetsPropertiesAK12:AK12bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsPropertiesAK12:AK12bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsPropertiesAK12:AK12chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsPropertiesAK12:AK12chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsPropertiesAK12:AK12chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsPropertiesAK12:AK12electronMultiplicity(I_electronMultiplicity)|JetsPropertiesAK12:AK12jetArea(F_jetArea)|JetsPropertiesAK12:AK12muonEnergyFraction(F_muonEnergyFraction)|JetsPropertiesAK12:AK12muonMultiplicity(I_muonMultiplicity)|JetsPropertiesAK12:AK12neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsPropertiesAK12:AK12neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsPropertiesAK12:AK12photonEnergyFraction(F_photonEnergyFraction)|JetsPropertiesAK12:AK12photonMultiplicity(I_photonMultiplicity)|JetsPropertiesAK12:AK12prunedMass(F_prunedMass)|JetsPropertiesAK12:AK12softDropMass(F_softDropMass)|JetsPropertiesAK12:AK12softDropPt(F_softDropPt)|JetsPropertiesAK12:AK12trimmedMass(F_trimmedMass)|JetsPropertiesAK12:AK12filteredMass(F_filteredMass)|JetsPropertiesAK12:AK12tau1(F_tau1)|JetsPropertiesAK12:AK12tau2(F_tau2)|JetsPropertiesAK12:AK12tau3(F_tau3)|JetsPropertiesAK12:AK12isLooseJetId(b_AK12isLooseJetId)|JetsPropertiesAK12:AK12isTightJetId(b_AK12isTightJetId)|JetsPropertiesAK12:PtCorr(F_PtCorr)|JetsPropertiesAK12:EtaCorr(F_EtaCorr)|JetsPropertiesAK12:PhiCorr(F_PhiCorr)|JetsPropertiesAK12:ECorr(F_ECorr)|JetsPropertiesAK12:mass(F_mass)|JetsPropertiesAK12:AK12correction(F_AK12correction)|JetsPropertiesAK12:AK12correctionUp(F_AK12correctionUp)|JetsPropertiesAK12:AK12correctionDown(F_AK12correctionDown)|JetsPropertiesAK12:AK12massCorrection(F_AK12massCorrection)|JetsPropertiesAK12:AK12massCorrectionUp(F_AK12massCorrectionUp)|JetsPropertiesAK12:AK12massCorrectionDown(F_AK12massCorrectionDown)'] ) # AK12 jet information on various variables

    RecoCandVector.extend(['JetsPropertiesPuppiAK8(PuppiAK8Jets)|JetsPropertiesPuppiAK8:PuppiAK8bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsPropertiesPuppiAK8:PuppiAK8bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsPropertiesPuppiAK8:PuppiAK8chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsPropertiesPuppiAK8:PuppiAK8chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsPropertiesPuppiAK8:PuppiAK8chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsPropertiesPuppiAK8:PuppiAK8electronMultiplicity(I_electronMultiplicity)|JetsPropertiesPuppiAK8:PuppiAK8jetArea(F_jetArea)|JetsPropertiesPuppiAK8:PuppiAK8muonEnergyFraction(F_muonEnergyFraction)|JetsPropertiesPuppiAK8:PuppiAK8muonMultiplicity(I_muonMultiplicity)|JetsPropertiesPuppiAK8:PuppiAK8neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsPropertiesPuppiAK8:PuppiAK8neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsPropertiesPuppiAK8:PuppiAK8photonEnergyFraction(F_photonEnergyFraction)|JetsPropertiesPuppiAK8:PuppiAK8photonMultiplicity(I_photonMultiplicity)|JetsPropertiesPuppiAK8:PuppiAK8prunedMass(F_prunedMass)|JetsPropertiesPuppiAK8:PuppiAK8softDropMass(F_softDropMass)|JetsPropertiesPuppiAK8:PuppiAK8softDropPt(F_softDropPt)|JetsPropertiesPuppiAK8:PuppiAK8trimmedMass(F_trimmedMass)|JetsPropertiesPuppiAK8:PuppiAK8filteredMass(F_filteredMass)|JetsPropertiesPuppiAK8:PuppiAK8tau1(F_tau1)|JetsPropertiesPuppiAK8:PuppiAK8tau2(F_tau2)|JetsPropertiesPuppiAK8:PuppiAK8tau3(F_tau3)|JetsPropertiesPuppiAK8:PuppiAK8isLooseJetId(b_PuppiAK8isLooseJetId)|JetsPropertiesPuppiAK8:PuppiAK8isTightJetId(b_PuppiAK8isTightJetId)|JetsPropertiesPuppiAK8:PtCorr(F_PtCorr)|JetsPropertiesPuppiAK8:EtaCorr(F_EtaCorr)|JetsPropertiesPuppiAK8:PhiCorr(F_PhiCorr)|JetsPropertiesPuppiAK8:ECorr(F_ECorr)|JetsPropertiesPuppiAK8:mass(F_mass)|JetsPropertiesPuppiAK8:PuppiAK8correction(F_PuppiAK8correction)|JetsPropertiesPuppiAK8:PuppiAK8correctionUp(F_PuppiAK8correctionUp)|JetsPropertiesPuppiAK8:PuppiAK8correctionDown(F_PuppiAK8correctionDown)|JetsPropertiesPuppiAK8:PuppiAK8massCorrection(F_PuppiAK8massCorrection)|JetsPropertiesPuppiAK8:PuppiAK8massCorrectionUp(F_PuppiAK8massCorrectionUp)|JetsPropertiesPuppiAK8:PuppiAK8massCorrectionDown(F_PuppiAK8massCorrectionDown)'] ) # PuppiAK8 jet information on various variables

    RecoCandVector.extend(['JetsPropertiesPuppi(JetsPuppi)|JetsPropertiesPuppi:bDiscriminatorCSV(F_bDiscriminatorCSV)|JetsPropertiesPuppi:bDiscriminatorICSV(F_bDiscriminatorICSV)|JetsPropertiesPuppi:chargedEmEnergyFraction(F_chargedEmEnergyFraction)|JetsPropertiesPuppi:chargedHadronEnergyFraction(F_chargedHadronEnergyFraction)|JetsPropertiesPuppi:chargedHadronMultiplicity(I_chargedHadronMultiplicity)|JetsPropertiesPuppi:electronMultiplicity(I_electronMultiplicity)|JetsPropertiesPuppi:jetArea(F_jetArea)|JetsPropertiesPuppi:muonEnergyFraction(F_muonEnergyFraction)|JetsPropertiesPuppi:muonMultiplicity(I_muonMultiplicity)|JetsPropertiesPuppi:neutralEmEnergyFraction(F_neutralEmEnergyFraction)|JetsPropertiesPuppi:neutralHadronMultiplicity(I_neutralHadronMultiplicity)|JetsPropertiesPuppi:photonEnergyFraction(F_photonEnergyFraction)|JetsPropertiesPuppi:photonMultiplicity(I_photonMultiplicity)|JetsPropertiesPuppi:isLooseJetId(b_isLooseJetId)|JetsPropertiesPuppi:isTightJetId(b_isTightJetId)|JetsPropertiesPuppi:PtCorr(F_PtCorr)|JetsPropertiesPuppi:EtaCorr(F_EtaCorr)|JetsPropertiesPuppi:PhiCorr(F_PhiCorr)|JetsPropertiesPuppi:ECorr(F_ECorr)|JetsPropertiesPuppi:AK4correction(F_AK4correction)|JetsPropertiesPuppi:AK4correctionUp(F_AK4correctionUp)|JetsPropertiesPuppi:AK4correctionDown(F_AK4correctionDown)'] ) # jet information on various variables

    RecoCandVector.extend(['Electrons(Electrons)|Electrons:charge(I_charge)|Electrons:isHEEP(b_isHEEP)|Electrons:type(I_type)|Electrons:mass(F_mass)|Electrons:pfDeltaCorrRelIso(F_pfDeltaCorrRelIso)|Electrons:pfRhoCorrRelIso04(F_pfRhoCorrRelIso04)|Electrons:pfRhoCorrRelIso03(F_pfRhoCorrRelIso03)|Electrons:pfRelIso(F_pfRelIso)|Electrons:photonIso(F_photonIso)|Electrons:neutralHadIso(F_neutralHadIso)|Electrons:chargedHadIso(F_chargedHadIso)|Electrons:trackIso(F_trackIso)|Electrons:isLoose(b_isLoose)|Electrons:isMedium(b_isMedium)|Electrons:isTight(b_isTight)|Electrons:SCEnergy(F_SCEnergy)|Electrons:deltaEtaSCTracker(F_deltaEtaSCTracker)|Electrons:deltaPhiSCTracker(F_deltaPhiSCTracker)|Electrons:sigmaIetaIeta(F_sigmaIetaIeta)|Electrons:sigmaIphiIphi(F_sigmaIphiIphi)'] ) # electron information on various variables
    RecoCandVector.extend(['Muons(Muons)|Muons:charge(I_charge)|Muons:isHighPt(b_isHighPt)|Muons:type(I_type)|Muons:mass(F_mass)|Muons:pfDeltaCorrRelIso(F_pfDeltaCorrRelIso)|Muons:pfRelIso(F_pfRelIso)|Muons:photonIso(F_photonIso)|Muons:neutralHadIso(F_neutralHadIso)|Muons:chargedHadIso(F_chargedHadIso)|Muons:trackIso(F_trackIso)|Muons:isLoose(b_isLoose)|Muons:isMedium(b_isMedium)|Muons:isTight(b_isTight)|Muons:isPFMuon(b_isPFMuon)|'] ) # muon information on various variables
    
    from AllHadronicSUSY.TreeMaker.treeMaker import TreeMaker
    process.TreeMaker2 = TreeMaker.clone(
    	TreeName          = cms.string("PreSelection"),
    	VarsRecoCand = RecoCandVector,
    	#VarsRecoCand = cms.vstring('selectedIDIsoMuons','selectedIDIsoElectrons','IsolatedTracks','HTJets'),
    	VarsDouble  	  = cms.vstring('WeightProducer:weight(Weight)','MET:Pt(METPt)','MET:Phi(METPhi)','MET:PtUp(METPtUp)','MET:PhiUp(METPhiUp)','MET:PtDown(METPtDown)','MET:PhiDown(METPhiDown)','MET:PtRaw(METPtRaw)','MET:PhiRaw(METPhiRaw)','MET:CaloMetPt(CaloMetPt)','MET:CaloMetPhi(CaloMetPhi)','GenEventInfo:genEventWeight(genEventWeight)','GenEventInfo:PUWeight(PUWeight)','METpuppi:Pt(METpuppiPt)','METpuppi:Phi(METpuppiPhi)','METpuppi:PtUp(METpuppiPtUp)','METpuppi:PhiUp(METpuppiPhiUp)','METpuppi:PtDown(METpuppiPtDown)','METpuppi:PhiDown(METpuppiPhiDown)','METpuppi:PtRaw(METpuppiPtRaw)','METpuppi:PhiRaw(METpuppiPhiRaw)','METpuppi:CaloMetPt(METpuppiCaloMetPt)','METpuppi:CaloMetPhi(METpuppiCaloMetPhi)'), #'MHT','HT','DeltaPhi:DeltaPhi1(DeltaPhi1)','DeltaPhi:DeltaPhi2(DeltaPhi2)','DeltaPhi:DeltaPhi3(DeltaPhi3)',
    	VarsInt = cms.vstring('NJets','NVtx','GenEventInfo:npT(npT)','FilterProducer:passFilterHBHE(passFilterHBHE)','FilterProducer:passFilterHBHEIso(passFilterHBHEIso)','FilterProducer:passFilterCSCHalo(passFilterCSCHalo)','FilterProducer:passFilterGoodVtx(passFilterGoodVtx)','FilterProducer:passFilterEEBadSC(passFilterEEBadSC)','FilterProducer:passFilterEcalDeadCellTriggerPrimitive(passFilterEcalDeadCellTriggerPrimitive)','FilterProducer:passFilterGlobalTightHalo2016(passFilterGlobalTightHalo2016)','FilterProducer:passFilterBadChCand(passFilterBadChCand)','FilterProducer:passFilterBadPFMuon(passFilterBadPFMuon)','FilterProducer:passFilterHBHELooseRerun(passFilterHBHELooseRerun)','FilterProducer:passFilterHBHETightRerun(passFilterHBHETightRerun)','FilterProducer:passFilterHBHEIsoRerun(passFilterHBHEIsoRerun)'),#,'Leptons'),
    #	VarsDoubleNamesInTree = cms.vstring('WeightProducer'),
#        VarsBool = cms.vstring('FilterProducer:passFilterEEBadSC(passFilterEEBadSC)'),
        debug = debug,
    	)

    process.TreeMaker2.VectorBool.extend(['TriggerProducer:TriggerPass'])
    process.TreeMaker2.VectorInt.extend(['TriggerProducer:TriggerPrescales'])
    process.TreeMaker2.VectorString.extend(['TriggerProducer:TriggerNames'])



    ## --- Final paths ----------------------------------------------------
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string("output.root"),
                                   )

    process.additionalJets = cms.Sequence() 
    if doAK10Reclustering:
        process.additionalJets += process.JetsPropertiesAK10;
    if doAK12Reclustering:
        process.additionalJets += process.JetsPropertiesAK12;
    if doPuppi:
        process.additionalJets += process.JetsPropertiesPuppiAK8;
        process.additionalJets += process.METpuppi;
        
    process.dump = cms.EDAnalyzer("EventContentAnalyzer")
    process.WriteTree = cms.Path(
#        process.HLTSelection*
        process.TriggerProducer* #NOT PRESENT IN 80x MINIAOD
        ### MET Filter Bits
#        process.HBHENoiseFilterResultProducer*
        process.badChargedCandidateFilter*
        process.badPFMuonFilter*
        process.FilterProducer* #this now contains all the met filters
#        process.metFilters* #this now contains all the met filters
#        process.metBits_miniAOD*
        ### HBHE noise filter
#        process.HBHENoiseFilterResultProducer* #produces HBHE baseline bools
#        process.ApplyBaselineHBHENoiseFilter*  #reject events based 
        #process.ApplyBaselineHBHEIsoNoiseFilter*   #reject events based - disable the module, performance is being investigated further
        ### rest of ntupling starts after here
        process.filterSeq *
        process.GenEventInfo *
        process.Muons *
        process.egmGsfElectronIDSequence*
        process.Electrons *
#    	process.selectedIDIsoMuons *
#    	process.selectedIDMuons *
#    	process.selectedIDIsoElectrons *
#    	process.selectedIDElectrons *
    	process.WeightProducer *
#    	process.IsolatedTracks *
        process.substructureSequenceGen *
        process.softdropGen_onMiniAOD *
        process.pruningGen_onMiniAOD *
        process.redoGenJets*
        process.GenJetsAK8 *
        process.GenJetsAK10 *
        process.GenJetsAK12 *
 #   	process.IsolatedTracksPT10 *
 #   	process.IsolatedTracksPT10IsoCut08 *
 #   	process.IsolatedTracksPT10IsoCut12 *
  #  	process.slimmedJetsPFCombinedSecondaryVertexBJetTags *
        process.puppi_onMiniAOD *
        process.substructureSequence *
        process.softdrop_onMiniAOD *
        process.pruning_onMiniAOD *
        process.redoPatJets*
        process.redoPuppiAK8Jets*
#      process.HTJets *
#      process.HT *
      process.NJets *
#      process.BTags *
#      process.ak4PFJets *
#      process.patJetCorrFactorsAK4 *
#      process.patJetsAK4 *
      process.JetsProperties *
      process.JetsPropertiesAK8 *
      process.JetsPropertiesPuppi *
      process.additionalJets *
#      process.Leptons *
      process.MET *
      process.NVtx *
      process.GenLeptons *
      process.GenJets *
    	#process.dump *
 #   	process.CountIsoTracks *
 #   	process.PrintDecay *
    	process.TreeMaker2

        )

#    process.outPath = cms.EndPath(process.out)
