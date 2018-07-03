import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

### ADDITIONS

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

#############

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource",# fileNames = cms.untracked.vstring())
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/uscms/home/ksiehl/nobackup/CMSSW_8_0_8_patch3/src/output/miniAOD-prod_PAT.root'
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/ECBB1C69-57D8-E611-ADA8-0025904A9430.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/D0D0C3B9-C0C4-E611-9CEE-0023AEEEB559.root'#,
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/8CC312BC-C0C4-E611-BB61-BC305B390AB4.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/461F0E7B-D2C4-E611-B465-BC305B390AB4.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36F3E28A-C0C4-E611-834A-0023AEEEB799.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/1C1FD740-C8C4-E611-A2D9-009C02AAB554.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/1476F163-D2C4-E611-9307-008CFAF558EE.root',
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/062203BC-C0C4-E611-A14A-BC305B390AB4.root'
    )
)

# MET filters
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

# for miniaod running
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.debug = cms.bool(False)
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.debug = cms.bool(False)

### Puppi (https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI)
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.useExistingWeights = cms.bool(True)
#'''
#### Toolbox (https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox)
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

#ak8Cut = 'pt > 200 && abs(eta) < 2.5'
#isMC   = True

#listBTagInfos = [
#     'pfInclusiveSecondaryVertexFinderTagInfos',
#]
#listBtagDiscriminatorsAK8 = [ 
#    'pfJetProbabilityBJetTags',
#    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
#    'pfCombinedMVAV2BJetTags',
#    # 'pfCombinedCvsLJetTags',
#    # 'pfCombinedCvsBJetTags',
#    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
#    # 'pfBoostedDoubleSecondaryVertexCA15BJetTags',
#]
#############################

#jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
#  runOnMC = isMC, 
#  PUMethod='CHS', 
#  addSoftDropSubjets = True, 
#  addTrimming = True, rFiltTrim=0.2, ptFrac=0.05,
#  addPruning = True, 
#  addFiltering = True, 
#  addSoftDrop = True, 
#  addNsub = True, 
#  bTagInfos = listBTagInfos, 
#  bTagDiscriminators = listBtagDiscriminatorsAK8, 
#  addCMSTopTagger = False, 
#  Cut = ak8Cut , 
#  addNsubSubjets = True, 
#  subjetMaxTau = 4 )

##############################


#jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
#  runOnMC = isMC, 
#  PUMethod='Puppi', 
#  addSoftDropSubjets = True, 
#  addTrimming = True,  rFiltTrim=0.2, ptFrac=0.05,
#  addPruning = True, 
#  addFiltering = True, 
#  addSoftDrop = True, 
#  addNsub = True, 
#  bTagInfos = listBTagInfos, 
#  bTagDiscriminators = listBtagDiscriminatorsAK8, 
#  addCMSTopTagger = False, 
#  Cut = ak8Cut , 
#  addNsubSubjets = True, 
#  subjetMaxTau = 4 )
#'''

######### NTUPLIZER
process.demo = cms.EDAnalyzer('ntuplizer',
    useToolbox = cms.bool(False),
    verbose = cms.bool(False),
    verboseGen = cms.bool(False),
    runGenLoop = cms.bool(True),
    isZprime = cms.bool(False),
    isttbar = cms.bool(False),
    isRSG = cms.bool(True),
    ak8chsInput = cms.InputTag("slimmedJetsAK8"),
    #ak8chsInput = cms.InputTag("selectedPatJetsAK8PFCHS"),
    #ak8puppiInput = cms.InputTag("slimmedJetsAK8"),
    ak8puppiInput = cms.InputTag("selectedPatJetsAK8PFPuppi"),
    ak8chsSubjetsInput   = cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked","SubJets"),
    #ak8chsSubjetsInput   = cms.InputTag("selectedPatJetsAK8PFCHSSoftDropPacked","Subjets"),
    triggerBits = cms.InputTag("TriggerResults", "", "HLT"),
    #ak8puppiSubjetsInput = cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked","SubJets"),
    ak8puppiSubjetsInput = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropPacked","SubJets"),
    theSrc = cms.InputTag("externalLHEProducer", "", "LHE")#,
    #jecPayloadsAK8chs = cms.vstring([
    #                                'JECs/Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2L3Residual_AK8PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt'
    #                                ]),
    #jecPayloadsAK4chs = cms.vstring([
    #                                'JECs/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt',
    #                                'JECs/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt'
    #                                ]),
    #jecPayloadsAK8pup = cms.vstring([
    #                                'JECs/Spring16_25nsV6_MC_L1FastJet_AK8PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2L3Residual_AK8PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_Uncertainty_AK8PFPuppi.txt'
    #                                ]),
    #jecPayloadsAK4pup = cms.vstring([
    #                                'JECs/Spring16_25nsV6_MC_L1FastJet_AK4PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2Relative_AK4PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_L3Absolute_AK4PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_L2L3Residual_AK4PFPuppi.txt',
    #                                'JECs/Spring16_25nsV6_MC_Uncertainty_AK4PFPuppi.txt'
    #                                ]),
    #jerSFtext         = cms.string('JERs/Spring16_25nsV6_MC_SF_AK8PFchs.txt'
    #                                )
)

process.TFileService = cms.Service("TFileService",
	#fileName = cms.string("signal_t.root"),
        #fileName = cms.string("background_t.root"),
        #fileName = cms.string("signal_WW600to800.root"),
        #fileName = cms.string("signal_WW800.root"),
        #fileName = cms.string("signal_WZ600to800.root"),
        #fileName = cms.string("signal_WZ800.root"),
        #fileName = cms.string("background_Wjet100to200.root"),
        #fileName = cms.string("background_Wjet200to400.root"),
        #fileName = cms.string("background_Wjet400to600.root"),
        #fileName = cms.string("background_Wjet600to800.root"),
        fileName = cms.string("background_Wjet800to1200.root"),
        #fileName = cms.string("background_Wjet1200to2500.root"),
        #fileName = cms.string("background_Wjet2500toinf.root"),
        #fileName = cms.string("background_WW.root"),
        #fileName = cms.string("background_WZ.root"),
        #fileName = cms.string("background_ttbar.root"),
        #fileName = cms.string("background_t_ch_top.root"),
        #fileName = cms.string("background_t_ch_tbar.root"),
        #fileName = cms.string("background_s_ch.root"),
        #fileName = cms.string("background_tW_ch_top.root"),
        #fileName = cms.string("background_tW_ch_tbar.root"),
	closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.BadChargedCandidateFilter*
    process.BadPFMuonFilter*
    process.demo
)
