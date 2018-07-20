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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",# fileNames = cms.untracked.vstring())
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:/uscms/home/ksiehl/nobackup/CMSSW_8_0_8_patch3/src/output/miniAOD-prod_PAT.root'
	'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_1.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_2.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_3.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_4.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_5.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_6.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_7.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_8.root',
        #'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_9.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_10.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_11.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_12.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_13.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_14.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_15.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_16.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_17.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_18.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_19.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_20.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_21.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_22.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_23.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_24.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_25.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_26.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_27.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_28.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_29.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_30.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_31.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_32.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_33.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_34.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_35.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_36.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_37.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_38.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_39.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_40.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_41.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_42.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_43.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_44.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_45.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_46.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_47.root',
        'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry/multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_48.root'
       	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WWToLNuQQ_MWW-800_PtW-180_aTGC_ShowerReconfig_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/FCB4520F-B1EE-E611-9678-1CC1DE782F02.root'
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
	fileName = cms.string("ntuple_output.root"),
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
        #fileName = cms.string("background_Wjet800to1200.root"),
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
