# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('OutputDAS-Data.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('InclusiveJet',
                                     # filename        = cms.vstring('file:/nfs/dust/cms/user/gunnep/DASInclusiveJets/Ntuples-Data-DAS2016.root'),
                                     filename        = cms.vstring('file:/nfs/dust/cms/user/connorpa/SMPJ/Ntuples-Data-RunD-ReReco16Dec2015-76X-part1.root'),

                                     treename        = cms.string('ProcessedTree'),
                                     dirname         = cms.string('ak4'),
                                     
                                     minPt           = cms.double(50.0),
                                     ymax            = cms.double(4.7),
                                     JetID           = cms.int32(2),
                                     
                                     printOk         = cms.int32(-10),
                                     
                                     isMCarlo        = cms.untracked.bool(False),
                                     
                                     # uncomment to apply JEC on the fly                                     
                                     # pseudoglobaltag = cms.string('Summer15_50nsV5'), # Make sure that the name is the same as in data directory   
                                     # jettype         = cms.string('AK4PFchs'),        # AKXPF or AKXPFchs  -  this is only for JECs                   
                                     
                                     # jecUncSrc  = cms.string('/nfs/dust/cms/user/gunnep/DASInclusiveJets/SMPJ/AnalysisFW/data/Fall15_25nsV2/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt'),
                                     # jecUncSrcNames  = cms.vstring('')

                                     #'AbsoluteStat,AbsoluteFlavMap', 'AbsoluteMPFBias', 'Fragmentation', 'SinglePionECAL',
                                     #'SinglePionHCAL', 'FlavorQCD', 'TimeEta', 'TimePt', 'RelativeJEREC1', 'RelativeJEREC2', 'RelativeJERHF',
                                     #'RelativePtBB', 'RelativePtEC1', 'RelativePtEC2', 'RelativePtHF', 'RelativeFSR', 'RelativeStatFSR',
                                     #'RelativeStatEC2', 'RelativeStatHF', 'PileUpDataMC', 'PileUpPtRef', 'PileUpPtBB', 'PileUpPtEC1', 
                                     #'PileUpPtEC2', 'PileUpPtHF', 'PileUpMuZero', 'PileUpEnvelope', 
                                     #'SubTotalPileUp', 'SubTotalRelative', 'SubTotalPt', 'SubTotalScale', 'SubTotalAbsolute', 'SubTotalMC', 'Total','TotalNoFlavor',
                                     #'TotalNoTime','TotalNoFlavorNoTime',
                                     #'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark','FlavorPureCharm',
                                     #'FlavorPureBottom','TimeRunA','TimeRunB','TimeRunC','TimeRunD','CorrelationGroupMPFInSitu',
                                     #'CorrelationGroupIntercalibration','CorrelationGroupbJES','CorrelationGroupFlavor',
                                     #'CorrelationGroupUncorrelated'),
                                     
                                     )

process.p = cms.Path(process.efficiency)

