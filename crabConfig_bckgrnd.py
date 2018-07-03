from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'background_studies_tryagain'
#config.General.workArea
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'

config.Data.inputDataset = '/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic' #'FileBased'
#config.Data.unitsPerJob ?
#config.DatatotalUnits ?
#config.Data.outLFNDirBase
config.Data.publication = True
config.Data.outputDatasetTag = 'Analysis_bckgrnd'
config.Data.ignoreLocality = True

config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US*']
