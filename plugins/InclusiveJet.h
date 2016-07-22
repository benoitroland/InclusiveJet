#ifndef InclusiveJet_h
#define InclusiveJet_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDEvent.h"
#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJet.h"
#include "SMPJ/AnalysisFW/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "PhysicsTools/Utilities/interface/LumiReweighting.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include <boost/shared_ptr.hpp>
using namespace edm;
using namespace std;

class InclusiveJet : public edm::EDAnalyzer
 {
  public:
    explicit InclusiveJet(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~InclusiveJet();

  private:
    //---- configurable parameters --------    
    std::string mjettype,mGlobalTag,mTreeName,mDirName;
    std::vector<std::string> mFileName;
    double mMinPt, mYMax;
    int    mJetID;        // looseID==1 tightID==2
    int    mprintOk;       // noPrint=0  Print=1
    int    mMCSlice;       // noPrint=0  Print=1
    bool mIsMCarlo;
    bool mPUReweighting;
    bool mLowPileUp;
    std::vector<std::string> mJECUncSrcNames;
    string mJECUncSrc;    

    edm::Service<TFileService> fs;
    //std::vector<TTree> *mTree;
    //std::vector<TFile> *mInf;
    TFile *mPuf,*mInf;
    TTree *mTree;
    //std::vector<TDirectoryFile> *mDir;
    TDirectoryFile *mDir;
   
    //-- trigger									 
    static const int Ntrigger = 9;	
    static const int Ntriggeredge = 10;
    double triggeredge[Ntriggeredge] = {0,114,133,220,300,430,507,638,737,3000};	 
    
    //-- pt bin																		  
    static const int Nptbin = 80;															  
    double Ptbinning[Nptbin+1] = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 	  
			    330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 	  
			    1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961,   
			    5220, 5492, 5777, 6076, 6389, 6717, 7000};                                                                                    

    //-- y bin	
    static const int Nybin = 7;
    static const int Nyedge = 9;								 
    double yedge[Nyedge] = {0,0.5,1,1.5,2,2.5,3,3.2,4.7};                              
    
    //- tree
    QCDEvent *Event;

    //-- pt hat
    TH1F *hmc_pthat;
    TH1F *hmc_pthat_weighted;   

    //-- PU and vertex		  
    TH1F *hnum_of_Vtx;	  	  
    TH1F *hnum_of_VtxGood;	  
				  
    TH2F *hPileUpVSVertex;	  
				  
    TH1F *hTruePileUpMCInteger;	  
    TH1F *hTruePileUpMC;          

    //-- generator level jets
    TH1F *hpt0_GENJet;
    TH1F *hy0_GENJet;
    TH1F *hphi0_GENJet;

    TH1F *hMultiplicity_GENJet;
    
    TH1F *hpt_GENJet[Nybin];

    TH1F *hpt_GENJetReso[Nptbin];
    TH1F *hy_GENJetReso[Nybin];     

    //-- detector level jets                                                                                                                                            
    TH1F *hpt0_DETJet;	   
    TH1F *hpt0_DETJetUncor;	   
    TH1F *hy0_DETJet;	   
    TH1F *hphi0_DETJet;        

    TH1F *hMultiplicity_DETJet;

    TH1F *hpt_DETJet[Nybin];     
    TH1F *hpt_DETJetUP[Nybin];   
    TH1F *hpt_DETJetDOWN[Nybin]; 
    TH1F *hpt_DETJetUncor[Nybin];

    //-- matching					  
    TH1F *hDeltaRpt;	      		  
    TH1F *hDeltaRy;  	      		  
			      		  
    TH1F *hDeltaRptmin;           	  
    TH1F *hDeltaRymin;        		  
			      		  
    //-- resolution
    TH1F *hResolutionPtvsY[Nybin];    
    TH1F *hResolutionPtvsPt[Nptbin];

    //-- purity, stability, acceptance, fake, miss
    TH1F *hpt_DETJetMatched[Nybin]; 
    TH1F *hpt_DETJetMatchedSame[Nybin]; 
    TH1F *hpt_DETJetUnMatched[Nybin];

    TH1F *hpt_GENJetMatched[Nybin];
    TH1F *hpt_GENJetMatchedSame[Nybin];
    TH1F *hpt_GENJetUnMatched[Nybin];

    TH2F *hpt_JetMatched[Nybin];

    TH1F *hpt_GENJetStability[Nybin];
    TH1F *hpt_GENJetAcceptance[Nybin];
    TH1F *hpt_GENJetMiss[Nybin];

    TH1F *hpt_DETJetPurity[Nybin];
    TH1F *hpt_DETJetFake[Nybin];

    //-- trigger
    TH1F *hleading_pt_all_Jet[Ntrigger];
    TH1F *hleading_pt_emulated_Jet[Ntrigger];
    TH1F *hleading_pt_HLTeffi_Jet[Ntrigger];
    
    TH1F *hleading_eta_all_Jet[Ntrigger];
    TH1F *hleading_eta_emulated_Jet[Ntrigger];
    TH1F *hleading_eta_HLTeffi_Jet[Ntrigger];
    
    //-- MET
    TH1F *hMET_DET;
    TH1F *hMETPhi_DET;
    TH1F *hFractionMET_DET;
 };

#endif

