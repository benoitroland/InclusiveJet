# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string("OutputDAS-MC.root"))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")


##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('InclusiveJet',
                                     #filename        = cms.vstring('file:/nfs/dust/cms/user/gunnep/DASInclusiveJets/Ntuples-MC-DAS2016.root'),
                                     filename        = cms.vstring('file:/nfs/dust/cms/user/gunnep/Ntuples/DASBenoit.root'),
                                     treename        = cms.string('ProcessedTree'),
                                     dirname         = cms.string('ak4'),
                                  
                                     minPt           = cms.double(50.0),   
                                     ymax            = cms.double(4.7),    
                                     JetID           = cms.int32(2),       

                                     printOk         = cms.int32(-10),                                                
                                     
                                     isMCarlo        = cms.untracked.bool(True),                                     
                                     PUReweighting   = cms.untracked.bool(False),                                    
                                     LowPileUp  = cms.untracked.bool(False),                                         
                                     
                                     # uncomment to apply JEC on the fly
                                     # pseudoglobaltag = cms.string('Summer15_50nsV5'), # Make sure that the name is the same as in data direcotry
                                     # jettype         = cms.string('AK4PFchs'),        # AKXPF or AKXPFchs  -  this is only for JECs
                                     
                                     # jecUncSrc  = cms.string('/nfs/dust/cms/user/gunnep/DASInclusiveJets/SMPJ/AnalysisFW/data/Fall15_25nsV2/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt'),
                                     # jecUncSrcNames  = cms.vstring(''),                                                             

                                     )

process.p = cms.Path(process.efficiency)
