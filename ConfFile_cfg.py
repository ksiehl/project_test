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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = 'file:/uscms/home/ksiehl/nobackup/CMSSW_8_0_8_patch3/src/output/miniAOD-prod_PAT.root'
    fileNames = cms.untracked.vstring(
	'root://cms-xrd-global.cern.ch//store/user/ksiehl/multiweight_rereretry//multiweights-nextstep/180515_234552/0000/miniAOD-prod_PAT_1.root'
  
#        'file:/uscms/home/ksiehl/nobackup/CMSSW_8_0_8_patch3/src/output/miniAOD-prod_PAT.root'
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
'''
### Toolbox (https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox)
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

ak8Cut = 'pt > 200 && abs(eta) < 2.5'
isMC   = True

listBTagInfos = [
     'pfInclusiveSecondaryVertexFinderTagInfos',
]
listBtagDiscriminatorsAK8 = [ 
    'pfJetProbabilityBJetTags',
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfCombinedMVAV2BJetTags',
    # 'pfCombinedCvsLJetTags',
    # 'pfCombinedCvsBJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    # 'pfBoostedDoubleSecondaryVertexCA15BJetTags',
]
############################

jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
  runOnMC = isMC, 
  PUMethod='CHS', 
  addSoftDropSubjets = True, 
  addTrimming = True, rFiltTrim=0.2, ptFrac=0.05,
  addPruning = True, 
  addFiltering = True, 
  addSoftDrop = True, 
  addNsub = True, 
  bTagInfos = listBTagInfos, 
  bTagDiscriminators = listBtagDiscriminatorsAK8, 
  addCMSTopTagger = False, 
  Cut = ak8Cut , 
  addNsubSubjets = True, 
  subjetMaxTau = 4 )

#############################


jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', 
  runOnMC = isMC, 
  PUMethod='Puppi', 
  addSoftDropSubjets = True, 
  addTrimming = True,  rFiltTrim=0.2, ptFrac=0.05,
  addPruning = True, 
  addFiltering = True, 
  addSoftDrop = True, 
  addNsub = True, 
  bTagInfos = listBTagInfos, 
  bTagDiscriminators = listBtagDiscriminatorsAK8, 
  addCMSTopTagger = False, 
  Cut = ak8Cut , 
  addNsubSubjets = True, 
  subjetMaxTau = 4 )
'''
'''
from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass

process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
	src = InputTag("ak8PFJetsCHS"),
	Njets = cms.vuint32(1,2,3,4),
	measureDefinition = cms.uint32(0),
	beta = cms.double(1.0),
	R0 = cms.double(0.8),
	Rcutoff = cms.double(999.0),
	axesDefinition = cms.uint32(6),
	nPass = cms.int32(999),
	akAxesR0 = cms.double(999)
)
'''
'''
from AllHadronicSUSY.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()

dataSetName = parameters.value("dataset","")
#dataSetName = parameters.value("dataset","file:/pnfs/desy.de/cms/tier2/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/F2742E0D-F603-E411-A246-0025905A60BE.root")
#dataSetName = parameters.value("dataset","file:/pnfs/desy.de/cms/tier2/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F452BBD7-BE76-E411-B1D7-002590DB928E.root")
#dataSetName = parameters.value("dataset","file:/pnfs/desy.de/cms/tier2/store/mc/Spring14miniaod/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/B6E695EA-DE18-E411-B4D9-002590596498.root")
global_tag = parameters.value("global_tag","")
MC= parameters.value("MC", True)
QCD= parameters.value("QCD", False)
LostLepton= parameters.value("LostLepton", True)
debug= parameters.value("debug", False)
nJetsMin    = parameters.value("njets_min",0)
htMin       = parameters.value("ht_min",0)
mhtMin      = parameters.value("mht_min",0)
NumProcessedEvt=parameters.value("NumProcessedEvt",-1)
METFiltersProcess=parameters.value("METFiltersProcess","")
DoAK8Reclustering=parameters.value("DoAK8Reclustering",False)
DoAK10Reclustering=parameters.value("DoAK10Reclustering",False)
DoAK12Reclustering=parameters.value("DoAK12Reclustering",False)
DoJECCorrection=parameters.value("DoJECCorrection",False)
DoPuppi=parameters.value("DoPuppi",True)
LeptonFilter=parameters.value("leptonFilter",True)
GenJetsAK8Reclustering=parameters.value("genJetsAK8Reclustering",True)
GenJetsAK10Reclustering=parameters.value("genJetsAK10Reclustering",False)
GenJetsAK12Reclustering=parameters.value("genJetsAK12Reclustering",False)
isHBHEEarlyData = parameters.value("isHBHEEarlyData",False)
isHBHERun2015D=parameters.value("isHBHERun2015D",True)
JsonFileName=parameters.value("jsonFileName","/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt")
IsCrab=parameters.value("isCrab",False)
ReDoPruningAndSoftdrop=parameters.value("ReDoPruningAndSoftdrop",False)
ReDoPruningAndSoftdropPuppi=parameters.value("ReDoPruningAndSoftdropPuppi",True)
#IsRun2016BCD=parameters.Value("IsRun2016BCD",False)
#IsRun2016E=parameters.Value("IsRun2016E",False)
#IsRun2016F=parameters.Value("IsRun2016F",False)
#IsRun2016GH=parameters.Value("IsRun2016GH",False)

processName      = parameters.value("name","RSGraviton1000")


print "***** SETUP ************************************"
print "  dataSetName : "+dataSetName
print " global_tag : "+global_tag
print " runningOnMC : "+str(MC)
print " runningOnQCD : "+str(QCD)
print " LostLepton(MC) : "+str(LostLepton)
print "     nJetsMin : "+str(nJetsMin)
print "        htMin : "+str(htMin)
print "       mhtMin : "+str(mhtMin)
print "       debug : "+str(debug)
print "       num of events : "+str(NumProcessedEvt)
print "       doAK8Reclustering : "+str(DoAK8Reclustering)
print "       doAK10Reclustering : "+str(DoAK10Reclustering)
print "       doAK12Reclustering : "+str(DoAK12Reclustering)
print "       doJECCorrection : "+str(DoJECCorrection)
print "       ReDoPruningAndSoftdrop :"+str(ReDoPruningAndSoftdrop)
print "       ReDoPruningAndSoftdropPuppi :"+str(ReDoPruningAndSoftdropPuppi)
print "       doPuppi : "+str(DoPuppi)
print "       leptonFilter : "+str(LeptonFilter)
print "       genJetsAK8Reclustering : "+str(GenJetsAK8Reclustering)
print "       genJetsAK10Reclustering : "+str(GenJetsAK10Reclustering)
print "       genJetsAK12Reclustering : "+str(GenJetsAK12Reclustering)
print "       isHBHEEarlyData : "+str(isHBHEEarlyData)
print "       isHBHERun2015D : "+str(isHBHERun2015D)
print "       jsonFileName : "+str(JsonFileName)
#print "       isRun2016BCD : "+str(IsRun2016BCD)
#print "       isRun2016E : "+str(IsRun2016E)
#print "       isRun2016F : "+str(IsRun2016F)
#print "       isRun2016GH : "+str(IsRun2016GH)
print "       isCrab : "+str(False)
print "************************************************"

# The process needs to be defined AFTER reading sys.argv,
# otherwise edmConfigHash fails
import FWCore.ParameterSet.Config as cms
#process = cms.Process("RA2EventSelection")
process = cms.Process(processName)

from AllHadronicSUSY.TreeMaker.makeTreeFromMiniAOD_cff import makeTreeTreeFromMiniAOD
makeTreeTreeFromMiniAOD(process,
                outFileName="ReducedSelection",
                NJetsMin=nJetsMin,
                HTMin=htMin,
                MHTMin=mhtMin,
                reportEveryEvt=1000,
                testFileName=dataSetName,
		Global_Tag=global_tag,
                METFiltersProcess=METFiltersProcess,
		MC=MC,
		QCD=QCD,
		LostLepton=LostLepton,
		debug = debug,
                numProcessedEvt=NumProcessedEvt,
                doAK8Reclustering=DoAK8Reclustering,
                doAK10Reclustering=DoAK10Reclustering,
                doAK12Reclustering=DoAK12Reclustering,
                doJECCorrection=DoJECCorrection,
                doPuppi=DoPuppi,
                leptonFilter=LeptonFilter,
                genJetsAK8Reclustering=GenJetsAK8Reclustering,
                genJetsAK10Reclustering=GenJetsAK10Reclustering,
                genJetsAK12Reclustering=GenJetsAK12Reclustering,
                customizeHBHENoiseForEarlyData=isHBHEEarlyData,
                customizeHBHENoiseForRun2015D=isHBHERun2015D,
                jsonFileName=JsonFileName,
                isCrab=IsCrab,
                reDoPruningAndSoftdrop=ReDoPruningAndSoftdrop,
                reDoPruningAndSoftdropPuppi=ReDoPruningAndSoftdropPuppi
#                isRun2016BCD=IsRun2016BCD,
#                isRun2016E=IsRun2016E,
#                isRun2016F=IsRun2016F,
#                isRun2016GH=IsRun2016GH
                )
'''
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
    theSrc = cms.InputTag("externalLHEProducer", "", "LHE"),
    jecPayloadsAK8chs = cms.vstring([
                                    'JECs/Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_L2L3Residual_AK8PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt'
                                    ]),
    jecPayloadsAK4chs = cms.vstring([
                                    'JECs/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt',
                                    'JECs/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt'
                                    ]),
    jecPayloadsAK8pup = cms.vstring([
                                    'JECs/Spring16_25nsV6_MC_L1FastJet_AK8PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_L2L3Residual_AK8PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_Uncertainty_AK8PFPuppi.txt'
                                    ]),
    jecPayloadsAK4pup = cms.vstring([
                                    'JECs/Spring16_25nsV6_MC_L1FastJet_AK4PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_L2Relative_AK4PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_L3Absolute_AK4PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_L2L3Residual_AK4PFPuppi.txt',
                                    'JECs/Spring16_25nsV6_MC_Uncertainty_AK4PFPuppi.txt'
                                    ]),
    jerSFtext         = cms.string('JERs/Spring16_25nsV6_MC_SF_AK8PFchs.txt'
                                    )
)

process.TFileService = cms.Service("TFileService",
	#fileName = cms.string("background_t.root")
	fileName = cms.string("signal_t.root"),
	closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.BadChargedCandidateFilter*
    process.BadPFMuonFilter*
    process.demo
)
