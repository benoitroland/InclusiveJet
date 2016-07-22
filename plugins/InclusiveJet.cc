#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TMath.h"
#include "TRandom.h"

#include "SMPJ/AnalysisFW/plugins/InclusiveJet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
InclusiveJet::InclusiveJet(edm::ParameterSet const& cfg)
{
  TH1::SetDefaultSumw2(true);
  
  mFileName       = cfg.getParameter<std::vector<std::string>>            ("filename");
  mTreeName       = cfg.getParameter<std::string>               ("treename");
  mDirName        = cfg.getParameter<std::string>               ("dirname");
  
  mMinPt     = cfg.getParameter<double> ("minPt");
  mYMax      = cfg.getParameter<double> ("ymax");
  mJetID     = cfg.getParameter<int>  ("JetID");
  
  mprintOk   = cfg.getParameter<int>  ("printOk");

  //-- uncomment to apply JEC on the fly
  //-- mGlobalTag        = cfg.getParameter<std::string>             ("pseudoglobaltag");
  //-- mjettype        = cfg.getParameter<std::string>               ("jettype");
  //-- mJECUncSrc = cfg.getParameter<std::string> ("jecUncSrc");
  //-- mJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");

  
  mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");  
}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void InclusiveJet::beginJob()
 {

   triggeredge[0] = mMinPt;
 
   //-- pt hat
   hmc_pthat = fs->make<TH1F>("mc_pthat","mc_pthat",200,0.,2000.);			     
   hmc_pthat_weighted = fs->make<TH1F>("mc_pthat_weighted","mc_pthat_weighted",100,0.,2000.);

   //-- PU and vertex
   hnum_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
   hnum_of_VtxGood = fs->make<TH1F>("num_of_VtxGood","num_of_VtxGood",100,0.,100.);
   
   hPileUpVSVertex  = fs->make<TH2F>("PileUpVSVertex","PileUpVSVertex",31,-0.5,30.5,31,-0.5,30.5); 

   hTruePileUpMC = fs->make<TH1F>("TruePileUpMC","TruePileUpMC",100,0,100); 			
   hTruePileUpMCInteger = fs->make<TH1F>("TruePileUpMCInteger","TruePileUpMCInteger",100,0,100);
  
   //-- generator level jets
   hpt0_GENJet  = fs->make<TH1F>("pt0_GENJet","pt0_GENJet",int(Nptbin),Ptbinning); 		  
   hy0_GENJet = fs->make<TH1F>("y0_GENJet","y0_GENJet",100,-5,5); 				  
   hphi0_GENJet = fs->make<TH1F>("phi0_GENJet","phi0_GENJet",60, -TMath::Pi(),TMath::Pi());   

   hMultiplicity_GENJet = fs->make<TH1F>("Multiplicity_GENJet","Multiplicity_GENJet",21, -0.5,20.5); 

   TString histo_label;
   TString histo_title;

   for(int i = 0; i < Nybin; i++) {
     histo_label = "pt_GENJet_" + TString::Format("%d",i+1) + "bin";
     histo_title = "pt_GENJet_" + TString::Format("%d",i+1) + "bin";
     hpt_GENJet[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);

     histo_label = "y_GENJetReso_" + TString::Format("%d",i+1) + "bin";	 
     histo_title = "y_GENJetReso_" + TString::Format("%d",i+1) + "bin";	 
     hy_GENJetReso[i] = fs->make<TH1F>(histo_label,histo_title,500,0,5);
   }

   for(int i = 0; i < Nptbin; i++) {							
     histo_label = "pt_GENJetReso_" + TString::Format("%d",i+1) + "bin";		
     histo_title = "pt_GENJetReso_" + TString::Format("%d",i+1) + "bin";		
     hpt_GENJetReso[i]  = fs->make<TH1F>(histo_label,histo_title,70000,0,7000);
   }

   //-- detector level jets
   hpt0_DETJet  = fs->make<TH1F>("pt0_DETJet","pt0_DETJet",int(Nptbin),Ptbinning); 		  
   hpt0_DETJetUncor  = fs->make<TH1F>("pt0_DETJetUncor","pt0_DETJetUncor",int(Nptbin),Ptbinning);  
   hy0_DETJet = fs->make<TH1F>("y0_DETJet","y0_DETJet",100,-5,5); 				  
   hphi0_DETJet = fs->make<TH1F>("phi0_DETJet","phi0_DETJet",60, -TMath::Pi(),TMath::Pi());   

   hMultiplicity_DETJet = fs->make<TH1F>("Multiplicity_DETJet","Multiplicity_DETJet",21, -0.5,20.5); 

   for(int i = 0; i < Nybin; i++) {
    
     histo_label = "pt_DETJet_" + TString::Format("%d",i+1) + "bin";
     histo_title = "pt_DETJet_" + TString::Format("%d",i+1) + "bin";
     hpt_DETJet[i] = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);

     histo_label = "pt_DETJetUP_" + TString::Format("%d",i+1) + "bin";
     histo_title = "pt_DETJetUP_" + TString::Format("%d",i+1) + "bin";
     hpt_DETJetUP[i] = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);

     histo_label = "pt_DETJetDOWN_" + TString::Format("%d",i+1) + "bin";
     histo_title = "pt_DETJetDOWN_" + TString::Format("%d",i+1) + "bin";
     hpt_DETJetDOWN[i] = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);

     histo_label = "pt_DETJetUncor_" + TString::Format("%d",i+1) + "bin";	  
     histo_title = "pt_DETJetUncor_" + TString::Format("%d",i+1) + "bin";	  
     hpt_DETJetUncor[i] = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);                                                                
   }

   //-- matching
   hDeltaRpt = fs->make<TH1F>("DeltaRpt","DeltaRpt",20000,0,10); 			       
   hDeltaRy = fs->make<TH1F>("DeltaRy","DeltaRy",20000,0,10);    			       
											       
   hDeltaRptmin = fs->make<TH1F>("DeltaRptmin","DeltaRptmin",20000,0,10); 		       
   hDeltaRymin = fs->make<TH1F>("DeltaRymin","DeltaRymin",20000,0,10);    		       
											       
   //-- resolution
   for(int i = 0; i < Nybin; i++) {    
     histo_label = "ResolutionPtvsY_" + TString::Format("%d",i+1) + "bin";
     histo_title = "ResolutionPtvsY_" + TString::Format("%d",i+1) + "bin";
     hResolutionPtvsY[i] = fs->make<TH1F>(histo_label,histo_title,1000,-5,5);
   }

   for(int i = 0; i < Nptbin; i++) {
     histo_label = "ResolutionPtvsPt_" + TString::Format("%d",i+1) + "bin";
     histo_title = "ResolutionPtvsPt_" + TString::Format("%d",i+1) + "bin";
     hResolutionPtvsPt[i] = fs->make<TH1F>(histo_label,histo_title,1000,-5,5);
   }

   //-- purity, stability, acceptance, fake, miss
   for(int i = 0; i < Nybin; i++) {

     //-- detector level 
     histo_label = "pt_DETJetMatched_" + TString::Format("%d",i+1) + "bin";
     histo_title = "pt_DETJetMatched_" + TString::Format("%d",i+1) + "bin";
     hpt_DETJetMatched[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);

     histo_label = "pt_DETJetMatchedSame_" + TString::Format("%d",i+1) + "bin";		   
     histo_title = "pt_DETJetMatchedSame_" + TString::Format("%d",i+1) + "bin";		   
     hpt_DETJetMatchedSame[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);

     histo_label = "pt_DETJetUnMatched_" + TString::Format("%d",i+1) + "bin";		   
     histo_title = "pt_DETJetUnMatched_" + TString::Format("%d",i+1) + "bin";		   
     hpt_DETJetUnMatched[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      

     //-- generator level
     histo_label = "pt_GENJetMatched_" + TString::Format("%d",i+1) + "bin";			   
     histo_title = "pt_GENJetMatched_" + TString::Format("%d",i+1) + "bin";			   
     hpt_GENJetMatched[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);	   
												   
     histo_label = "pt_GENJetMatchedSame_" + TString::Format("%d",i+1) + "bin";		   	   
     histo_title = "pt_GENJetMatchedSame_" + TString::Format("%d",i+1) + "bin";		   	   
     hpt_GENJetMatchedSame[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);	   
												   
     histo_label = "pt_GENJetUnMatched_" + TString::Format("%d",i+1) + "bin";		   	   
     histo_title = "pt_GENJetUnMatched_" + TString::Format("%d",i+1) + "bin";		   	   
     hpt_GENJetUnMatched[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      
     
     //-- correlation
     histo_label = "pt_JetMatched_" + TString::Format("%d",i+1) + "bin";		   
     histo_title = "pt_JetMatched_" + TString::Format("%d",i+1) + "bin";		   
     hpt_JetMatched[i]  = fs->make<TH2F>(histo_label,histo_title,int(Nptbin),Ptbinning,int(Nptbin),Ptbinning); 
     
     //-- stability
     histo_label = "pt_GENJetStability_" + TString::Format("%d",i+1) + "bin";		   	   
     histo_title = "pt_GENJetStability_" + TString::Format("%d",i+1) + "bin";		   	   
     hpt_GENJetStability[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      

     //-- acceptance										   
     histo_label = "pt_GENJetAcceptance_" + TString::Format("%d",i+1) + "bin";		   	   
     histo_title = "pt_GENJetAcceptance_" + TString::Format("%d",i+1) + "bin";		   	   
     hpt_GENJetAcceptance[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      

     //-- miss										    
     histo_label = "pt_GENJetMiss_" + TString::Format("%d",i+1) + "bin";		   	    
     histo_title = "pt_GENJetMiss_" + TString::Format("%d",i+1) + "bin";		   	    
     hpt_GENJetMiss[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      

     //-- purity										      
     histo_label = "pt_DETJetPurity_" + TString::Format("%d",i+1) + "bin";		      
     histo_title = "pt_DETJetPurity_" + TString::Format("%d",i+1) + "bin";		      
     hpt_DETJetPurity[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      

     //-- fake										
     histo_label = "pt_DETJetFake_" + TString::Format("%d",i+1) + "bin";		      	
     histo_title = "pt_DETJetFake_" + TString::Format("%d",i+1) + "bin";		      	
     hpt_DETJetFake[i]  = fs->make<TH1F>(histo_label,histo_title,int(Nptbin),Ptbinning);      
   }

   //-- trigger efficiency
   TString HLT[Ntrigger] = {"60","80","140","200","260","320","400","450","500"};
   
   for(int itrig = 0; itrig < Ntrigger; itrig++) {
     
     histo_label = "leading_pt_all_Jet" + HLT[itrig];
     histo_title = "leading_pt_all_Jet" + HLT[itrig];
     hleading_pt_all_Jet[itrig] = fs->make<TH1F>(histo_label,histo_title,300,0,3000); // int(Nptbin),Ptbinning); 

     histo_label = "leading_pt_emulated_Jet" + HLT[itrig];
     histo_title = "leading_pt_emulated_Jet" + HLT[itrig];
     hleading_pt_emulated_Jet[itrig] = fs->make<TH1F>(histo_label,histo_title,300,0,3000); // int(Nptbin),Ptbinning); 

     histo_label = "leading_pt_HLTeffi_Jet" + HLT[itrig];					 
     histo_title = "leading_pt_HLTeffi_Jet" + HLT[itrig];
     hleading_pt_HLTeffi_Jet[itrig] = fs->make<TH1F>(histo_label,histo_title,300,0,3000); // int(Nptbin),Ptbinning); 

     histo_label = "leading_eta_all_Jet" + HLT[itrig];					    
     histo_title = "leading_eta_all_Jet" + HLT[itrig];					    
     hleading_eta_all_Jet[itrig] = fs->make<TH1F>(histo_label,histo_title,24,-5.2,5.2);

     histo_label = "leading_eta_emulated_Jet" + HLT[itrig];				       
     histo_title = "leading_eta_emulated_Jet" + HLT[itrig];				       
     hleading_eta_emulated_Jet[itrig] = fs->make<TH1F>(histo_label,histo_title,24,-5.2,5.2);

     histo_label = "leading_eta_HLTeffi_Jet" + HLT[itrig];				       
     histo_title = "leading_eta_HLTeffi_Jet" + HLT[itrig];				       
     hleading_eta_HLTeffi_Jet[itrig] = fs->make<TH1F>(histo_label,histo_title,24,-5.2,5.2);   
   }
   
   //-- MET
   hMET_DET  = fs->make<TH1F>("MET_DET","MET_DET",500,0,2000); 
   hMETPhi_DET  = fs->make<TH1F>("METPhi_DET","METPhi_DET",35,-7,7); 
   hFractionMET_DET  = fs->make<TH1F>("FractionMET_DET","FractionMET_DET",200,0,1); 

 } // end of function beginJob()


 //------------------------ endjob() function declaration ---------------------- //
 void InclusiveJet::endJob()
 {

   mInf->Close();

 } // closing endJob()


 //--------------------------- analyze() fuction declaration ------------------ //
void InclusiveJet::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

  for(int itrig = 0; itrig < Ntriggeredge-1; itrig++) {
    cout<<"range trigger "<<itrig+1<<" = "<<triggeredge[itrig]<<" - "<<triggeredge[itrig+1]<<endl;
  }

  cout<<endl;

  for(int iy = 0; iy < Nyedge-1; iy++) {
    if(iy == 6) continue;
    
    int ibin = 0;
    if(iy <= 5) ibin = iy;
    else ibin = iy - 1;

    cout<<"y bin "<<ibin+1<<" = "<<yedge[iy]<<" - "<<yedge[iy+1]<<endl;
  }
  
  cout<<endl;

  cout<<"Number of input files = "<<mFileName.size()<<endl<<endl;
  
  //-- loop over files
  for(unsigned ifile=0;ifile<mFileName.size();ifile++) {
     
    mInf = TFile::Open(mFileName[ifile].c_str());
    mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());  
    mTree=(TTree*)mDir->Get(mTreeName.c_str());
     
    Event = new QCDEvent();
    TBranch *branch = mTree->GetBranch("events");
    branch->SetAddress(&Event);  

    unsigned NEntries = mTree->GetEntries();
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl<<endl;  
     
    int decade = 0 ;
     
    float hweight = 1.;  

    //-- L1 and HLT thresholds
    int HLTth[Ntrigger] = {60,80,140,200,260,320,400,450,500};
    int L1Th[Ntrigger] = {10,10,10,10,10,10,10,10,10};
    
    //-- trigger information
    bool hltPassj[Ntrigger]; 
    bool l1cut[Ntrigger];
    bool hltcut[Ntrigger];
    float prescalej[Ntrigger];    
   
    //-- Define trigger
    TString HLTJetName[Ntrigger];
    int HLTJetThName[Ntrigger] = {40,60,80,140,200,260,320,400,450};

    for (int itrig = 0; itrig < Ntrigger; itrig++)
      HLTJetName[itrig] = "HLT_PFJet" + TString::Format("%d",HLTJetThName[itrig]) + "_v3";

    //-- Retrieve requested trigger index 
    int ihltj[Ntrigger];
    for (int itrig = 0; itrig < Ntrigger; itrig++) ihltj[itrig] = -1;
    
    TH1F *hTrigNames = (TH1F*)mInf->Get("ak4/TriggerNames");
    cout<<"Finding trigger mapping: "<<endl;
    
     for(int ibin = 0; ibin < hTrigNames->GetNbinsX(); ibin++) {
      
       TString TrigName(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));

       for (int itrig = 0; itrig < Ntrigger; itrig++){
	 if (TrigName == HLTJetName[itrig]) {
	   ihltj[itrig] = ibin;
	   continue;
	 }
       } 
     }
     
     for (int itrig = 0; itrig < Ntrigger; itrig++){
       if (ihltj[itrig] == -1) 
	 cout<<"The requested trigger ("<<HLTJetName[itrig]<<") is not found "<<endl;
       else 
	 cout<<HLTJetName[itrig]<<" --> "<<ihltj[itrig]<<endl;
     }

     //-- loop over tree entries
     for (unsigned ientry = 0; ientry < NEntries; ientry++) {
     
       hweight = 1;

       double WeightMC = 1;
       unsigned NEntriesNorm = NEntries;

       mTree->GetEntry(ientry);

       //-- pileup - weight
       double TruePileUpInteger = 0;
       if(mIsMCarlo) {
	 TruePileUpInteger = Event->evtHdr().pu()/Event->evtHdr().nbx();
	 WeightMC = 2022100000*Event->evtHdr().weight()/NEntriesNorm; //-- flat MC
	 hweight = WeightMC;	 
       }

       if(hweight < 0) continue; 

       //-- progress report
       double progress = 10.0*ientry/(1.0*NEntriesNorm);
       int k = TMath::FloorNint(progress);
       if (k > decade) cout<<10*k<<" %"<<endl;
       decade = k;

       //-- trigger information
       if(!mIsMCarlo){
	 
	 for (int itrig = 0; itrig < Ntrigger; itrig++){

	   hltPassj[itrig] = false ;
	   prescalej[itrig] = 1;
	 
	   if (Event->fired(ihltj[itrig]) > 0) {
	     hltPassj[itrig] = true;

	     if(Event->preL1(ihltj[itrig]) >= 0) 
	       prescalej[itrig] = Event->preL1(ihltj[itrig]) * Event->preHLT(ihltj[itrig]);

	     else 
	       prescalej[itrig] = Event->preHLT(ihltj[itrig]);

	     if(mprintOk != 4) continue;
	     cout<<itrig<<" "<<"HLT trigger = "<<HLTJetName[itrig]<<" with prescale = "<<prescalej[itrig]<<"has fired"<<endl;
	     cout<<"L1 prescale = "<<Event->preL1(ihltj[itrig])<<", HLT prescale = "<< Event->preHLT(ihltj[itrig])<<endl;	     
	   }
	 }
       }
    
       int DetJets = 0;				         
       int DETjet_ok[100]; 			  
       for(int i = 0;i < 100; ++i) DETjet_ok[i] = 0;

       if(mIsMCarlo){

	 //-- pt hat - MC wieght
	 double pthat = Event->evtHdr().pthat();
	 double mc_weight = Event->evtHdr().weight();

	 if(mprintOk==10) cout<<"pt hat = "<<pthat<<", MC weight = "<<mc_weight<<endl;	 

	 hmc_pthat->Fill(pthat);
	 hmc_pthat_weighted->Fill(pthat,hweight);

	 //-- vertex information
	 hnum_of_Vtx->Fill(Event->evtHdr().nVtx(),hweight);

	 if(mprintOk == 1) 
	   cout<<"N Vtx = "<<Event->evtHdr().nVtx()<<", N good Vtx ="<<Event->evtHdr().nVtxGood()<<", PV good = "<<Event->evtHdr().isPVgood()<<endl;
	
	 //-- good vertex filter
	 if (Event->evtHdr().isPVgood() != 1) continue;
	 
	 //-- generated jets

	 int SelGenJets = 0 ;	 
	 unsigned n_genJets = Event->nGenJets(); 
	
	 if(mprintOk == 2) cout<<"Number of GENJets = "<<n_genJets<<endl;
	 
	 //-- loop over generated jets
	 for(unsigned ijet = 0; ijet < n_genJets; ++ijet) {
	   
	   double ptjet = Event->genjet(ijet).pt();
	   double yjet = Event->genjet(ijet).Rapidity(); 
	   double phijet = Event->genjet(ijet).phi();

	   if(mprintOk == 2) cout<<ijet<<": pT = "<<ptjet<<", y = "<<yjet<<", phi = "<<phijet<<endl;	   

	   if(ptjet < mMinPt) continue;

	   if(fabs(yjet ) > mYMax) continue;
	  	   
	   SelGenJets++;

	   //-- full y							     
	   hpt0_GENJet->Fill(ptjet,hweight);			     
	   hy0_GENJet->Fill(yjet,hweight);				     
	   hphi0_GENJet->Fill(phijet,hweight);				     
  
	   //-- loop y bin
	   for(int iy = 0; iy < Nyedge-1; iy++) {

	     if(iy == 6) continue;

	     int ihisto = 0;
	     if(iy <= 5) ihisto = iy;
	     else ihisto = iy - 1;

	     //-- y bin condition
	     if(fabs(yjet) > yedge[iy] && fabs(yjet) <= yedge[iy+1]) {
	       hpt_GENJet[ihisto]->Fill(ptjet,hweight);                    
	       hy_GENJetReso[ihisto]->Fill(fabs(yjet),hweight);
	     } //-- end y bin condition	     	   

	   } //-- end loop y bin

	   //-- loop pt bin 
	   for(int ipt = 0; ipt < Nptbin; ipt++) {							   
	     if(ptjet > Ptbinning[ipt] && ptjet < Ptbinning[ipt+1]) {	   
	       hpt_GENJetReso[ipt]->Fill(ptjet,hweight);					   
	     }											   
	   }                                                                                             

	 } //-- end loop over generated jets
	 
	 hMultiplicity_GENJet->Fill(SelGenJets,hweight);
       
	 //-- Vertex selection
	 hTruePileUpMC->Fill(Event->evtHdr().trpu(),hweight);
	 hPileUpVSVertex->Fill(TruePileUpInteger,Event->evtHdr().nVtx(),hweight);

	 //-- Vertex selection
	 if(mprintOk == 1) {
	   cout<<"N Vtx = "<<Event->evtHdr().nVtx()<<", N good Vtx = "<<Event->evtHdr().nVtxGood()<<", PV good = "<<Event->evtHdr().isPVgood()<<endl;
	   cout<<"PF rho = "<<Event->evtHdr().pfRho()<<endl;
	 }

	 //-- Dump all jets no cuts only good PV
	 if(mprintOk==1){
	   cout<<"Number of PFJets = "<<Event->nPFJetsCHS()<<endl;
	   for(unsigned ijet = 0; ijet < Event->nPFJetsCHS(); ++ijet) {
	     cout<<ijet<<": pT raw = "<<Event->pfjetchs(ijet).pt()<<", pT corr = "<<Event->pfjetchs(ijet).ptCor()<<endl;
	     cout<<"y = "<<Event->pfjetchs(ijet).y()<<", phi = "<<Event->pfjetchs(ijet).phi()<<endl;
	     cout<<"correction = "<<Event->pfjetchs(ijet).cor()<<", tight id = "<<Event->pfjetchs(ijet).tightID()<<endl;
	   }
	 }
       
	 //-- vector of selected jets at detector level
	 vector<double> RecoJet_y;
	 vector<double> RecoJet_phi;
	 vector<double> RecoJet_pt;

	 for(unsigned ijet = 0; ijet < Event->nPFJetsCHS(); ++ijet){

	   double yjet = Event->pfjetchs(ijet).y();
	   double phijet = Event->pfjetchs(ijet).phi();
	   double ptjet = Event->pfjetchs(ijet).ptCor();
	   double tightid = Event->pfjetchs(ijet).tightID();

	   if(ptjet < mMinPt) continue;
	   
	   if(fabs(yjet) > mYMax) continue;

	   bool jetid_cen = tightid;
	   bool jetid_fwd = (Event->pfjetchs(ijet).nemf() < 0.90) && (Event->pfjetchs(ijet).ncand() > 10);
		
	   if((fabs(yjet) <= 3.0) && !jetid_cen) continue;
	   if((fabs(yjet) > 3.2) && (fabs(yjet) <= 4.7) && !jetid_fwd) continue;
 
	   RecoJet_y.push_back(yjet);
	   RecoJet_phi.push_back(phijet);
	   RecoJet_pt.push_back(ptjet);  
	 }

	 //-- vector of selected jets at generator level
	 vector<double> GenJet_y;	   
	 vector<double> GenJet_phi;   
	 vector<double> GenJet_pt;

	 for(unsigned ijet = 0; ijet < n_genJets; ++ijet){

	   double yjet = Event->genjet(ijet).Rapidity();
	   double phijet = Event->genjet(ijet).phi();
	   double ptjet = Event->genjet(ijet).pt();
	   
	   if(ptjet < mMinPt) continue;
	   
	   if(fabs(yjet) > mYMax) continue;

	   GenJet_y.push_back(yjet);     
	   GenJet_phi.push_back(phijet); 
	   GenJet_pt.push_back(ptjet);   
	 }

	 //-- distribution without matching for pt and y resolution
	 //-- pt resolution: Rpt (phi,y) - Rpt (phi,y) min
	 //-- y resolution:  Ry (phi,pt) - Ry (phi,pt) min

	 //-- loop over selected reco jets
	 for(unsigned ireco = 0; ireco < RecoJet_y.size(); ++ireco) {
	   
	   double dRptmin = 1000;
	   double dRymin = 1000;

	   //-- loop over selected gen jets
	   for(unsigned igen = 0; igen< GenJet_y.size(); ++igen) {
	     
	     double dphi = GenJet_phi.at(igen) - RecoJet_phi.at(ireco);
	     if(dphi < -TMath::Pi()) dphi+=2*TMath::Pi();
	     if(dphi > TMath::Pi()) dphi-=2*TMath::Pi();
	     dphi=fabs(dphi);
	     
	     double dy = GenJet_y.at(igen) - RecoJet_y.at(ireco);
	     double dpt = (GenJet_pt.at(igen) - RecoJet_pt.at(ireco))/GenJet_pt.at(igen);
	     
	     double dRpt = sqrt(pow(dy,2) + pow(dphi,2));
	     double dRy = sqrt(pow(dpt,2) + pow(dphi,2));

	     if(dRpt < dRptmin) dRptmin = dRpt;
	     if(dRy < dRymin) dRymin = dRy;

	     hDeltaRpt->Fill(dRpt,hweight); 
	     hDeltaRy->Fill(dRy,hweight);   

	   } //-- end loop over gen jets

	   hDeltaRptmin->Fill(dRptmin,hweight); 
	   hDeltaRymin->Fill(dRymin,hweight);   
	   
	 } //-- end loop over reco jets

	 //-- pt resolution 
	 //-- jets ordered in pt: FIRST jet with Rpt (phi,y) < Rpt CUT 

	 double dRptCUT = 0.2;

	 //-- loop over selected reco jets
	 for(unsigned ireco = 0; ireco < RecoJet_y.size(); ++ireco) {

	   bool firstjet = true;

	   double dRptmin = 1000;  //-- can be different from the true min
	   int igenmin = -1;     //-- can be different from the true min
 
	   //-- loop over selected gen jets
	   for(unsigned igen = 0; igen< GenJet_y.size(); ++igen) {
	     
	     double dphi = GenJet_phi.at(igen) - RecoJet_phi.at(ireco);
	     if(dphi < -TMath::Pi()) dphi+=2*TMath::Pi();
	     if(dphi > TMath::Pi()) dphi-=2*TMath::Pi();
	     dphi=fabs(dphi);
	     
	     double dy = GenJet_y.at(igen) - RecoJet_y.at(ireco);

	     double dRpt = sqrt(pow(dy,2) + pow(dphi,2));
	     
	     //-- matching condition
	     if(!(dRpt < dRptCUT && firstjet)) continue;

	     dRptmin  = dRpt;
	     igenmin = igen;
	     firstjet = false;

	   } //-- end loop over selected gen jets  
	   
	   //-- loop over y bin 												   
	   for(int iy = 0; iy < Nyedge-1; iy++) {						     				   
	     
	     if(iy == 6) continue;								     				   
	     
	     int ihisto = 0;									     				   
	     if(iy <= 5) ihisto = iy;								     				   
	     else ihisto = iy - 1;                 						     				   
	     
	     //-- y bin condition								     				   
	     if(fabs(RecoJet_y.at(ireco)) > yedge[iy] && fabs(RecoJet_y.at(ireco)) <= yedge[iy+1]) {                             
	 
	       //-- no gen jet associated to the reco jet  			     
	       if (dRptmin > dRptCUT) 	     				     
		 hpt_DETJetUnMatched[ihisto]->Fill(RecoJet_pt.at(ireco),hweight);	     	       

	       //-- matching condition (gen jet associated to the reco jet)
	       if (dRptmin < dRptCUT) {						       
		 hpt_DETJetMatched[ihisto]->Fill(RecoJet_pt.at(ireco),hweight);							 		 
		 hpt_JetMatched[ihisto]->Fill(RecoJet_pt.at(ireco),GenJet_pt.at(igenmin),hweight);                               
		 
		 double resolution = (RecoJet_pt.at(ireco)-GenJet_pt.at(igenmin))/GenJet_pt.at(igenmin);		 
		 hResolutionPtvsY[ihisto]->Fill(resolution,hweight);

		 //-- loop over pt bin
		 for(int ipt = 0; ipt < Nptbin; ipt++) {
		   
		   //-- gen pt condition
		   if(GenJet_pt.at(igenmin) > Ptbinning[ipt] && GenJet_pt.at(igenmin) < Ptbinning[ipt+1]) {
		     hResolutionPtvsPt[ipt]->Fill(resolution,hweight);					 	     
		     
		     //-- reco pt condition
		     if(RecoJet_pt.at(ireco) > Ptbinning[ipt] && RecoJet_pt.at(ireco) < Ptbinning[ipt+1]) {
		       hpt_DETJetMatchedSame[ihisto]->Fill(RecoJet_pt.at(ireco),hweight);                         

		     } //-- end reco pt condition 		   
		   } //-- end gen pt condition
		 } //-- end loop over pt bin
	       
	       } //-- end matching condition	       
	     } //-- end y bin condition	     
	   } //-- end loop over y bin
  
	   //-- erase gen jet if matched to a reco jet
	   if (dRptmin < dRptCUT) {						       
	     GenJet_y.erase(GenJet_y.begin()+igenmin);	   
	     GenJet_phi.erase(GenJet_phi.begin()+igenmin);	   
	     GenJet_pt.erase(GenJet_pt.begin()+igenmin);	   
	   } 
	   	   
	 } //-- end loop over selected reco jets

	 //-- vector of selected gen jets RE-INITIALIZED
	 GenJet_y.clear();	   
	 GenJet_phi.clear();   
	 GenJet_pt.clear();
       
	 for(unsigned ijet = 0; ijet < n_genJets; ++ijet){
	   
	   double yjet = Event->genjet(ijet).Rapidity();
	   double phijet = Event->genjet(ijet).phi();
	   double ptjet = Event->genjet(ijet).pt();
	   
	   if(ptjet < mMinPt) continue;
	   
	   if(fabs(yjet) > mYMax) continue;
	   
	   GenJet_y.push_back(yjet);     
	   GenJet_phi.push_back(phijet); 
	   GenJet_pt.push_back(ptjet);   
	 }
	 
	 //-- gen miss (gen without reco)
	 //-- jet with Rpt (phi,y) min
	 
	 //-- loop over selected gen jets
	 for(unsigned igen = 0; igen < GenJet_y.size(); ++igen) {
	   
	   double dRptmin = 1000;
	   unsigned irecomin = -1;
	   
	   //-- loop over selected reco jets
	   for(unsigned ireco = 0; ireco < RecoJet_y.size(); ++ireco) {
	     
	     double dphi = GenJet_phi.at(igen) - RecoJet_phi.at(ireco);
	     if(dphi < -TMath::Pi()) dphi+=2*TMath::Pi();
	     if(dphi > TMath::Pi()) dphi-=2*TMath::Pi();
	     dphi=fabs(dphi);
	     
	     double dy = GenJet_y.at(igen) - RecoJet_y.at(ireco);
	     double dRpt = sqrt(pow(dy,2) + pow(dphi,2));
	     
	     if(!(dRpt < dRptmin)) continue;
	     dRptmin = dRpt;
	     irecomin = ireco;
	     
	   } //-- end loop over selected reco jets
	   
	   //-- loop over y bin 												 
	   for(int iy = 0; iy < Nyedge-1; iy++) {						     				 
	     
	     if(iy == 6) continue;								     				 
	     
	     int ihisto = 0;									     				 
	     if(iy <= 5) ihisto = iy;								     				 
	     else ihisto = iy - 1;                 						     				 
	     
	     //-- y bin condition								     				 
	     if(fabs(GenJet_y.at(igen)) > yedge[iy] && fabs(GenJet_y.at(igen)) <= yedge[iy+1]) {				 
	       
	       //-- matching condition: reco jet associated to the gen jet 	       
	       if(dRptmin < dRptCUT) {
		 hpt_GENJetMatched[ihisto]->Fill(GenJet_pt.at(igen),hweight);							 

		 //-- loop over pt bin										  
		 for(int ipt = 0; ipt < Nptbin; ipt++) {							  
		   
		   //-- reco pt condition									  
		   if(RecoJet_pt.at(irecomin) > Ptbinning[ipt] && RecoJet_pt.at(irecomin) < Ptbinning[ipt+1]) {	  	  

		     //-- gen pt condition									  
		     if(GenJet_pt.at(igen) > Ptbinning[ipt] && GenJet_pt.at(igen) < Ptbinning[ipt+1]) {	  
		       hpt_GENJetMatchedSame[ihisto]->Fill(GenJet_pt.at(igen),hweight);							 
		     }
		   }
		 } //-- end loop over pt bin

	       } //-- end matching condition

	       //-- no reco jet associated to the gen jet - gen miss		 
	       if(dRptmin > dRptCUT)													
		 hpt_GENJetUnMatched[ihisto]->Fill(GenJet_pt.at(igen),hweight);					

	     } //-- end y bin condition	     
	   } //-- end loop over y bin
	   
	   //-- erase reco jet matched to a gen jet
	   if(dRptmin < dRptCUT) {
	     RecoJet_y.erase(RecoJet_y.begin()+irecomin);	   
	     RecoJet_phi.erase(RecoJet_phi.begin()+irecomin);	   
	     RecoJet_pt.erase(RecoJet_pt.begin()+irecomin);	   
	   }
	   
	 } //-- end loop over selected gen jets
	 
       } //-- end isMC
	 
       //--  Trigger efficiency

       for(int itrig = 0; itrig < Ntrigger; itrig++) {
	 hltcut[itrig] = false; 
	 l1cut[itrig] = false; 
       }

       if(!mIsMCarlo){
	 
	 double leading_pt = -10.0;
	 double leading_eta = -10.0;

	 hnum_of_Vtx->Fill(Event->evtHdr().nVtx(),hweight);

	 ///-- fillter on good primary vertex
	 if(Event->evtHdr().isPVgood() != 1) continue;
	 
	 for(int itrig = 0; itrig < Ntrigger; itrig++) {
	   
	   if (hltPassj[itrig]) {

	     //-- L1 theshold check
	     for(unsigned l1iobj = 0; l1iobj < Event->nL1Obj(ihltj[itrig]); l1iobj++) {	       
	       //-- cout<<"L1 object id: "<<l1iobj<<", pT = "<<Event->l1obj(ihltj[itrig],l1iobj).pt()<<", pT th = "<<L1Th[itrig]<<endl;
	       if(Event->l1obj(ihltj[itrig],l1iobj).pt() > L1Th[itrig]) l1cut[itrig] = true;
	     } 
	       
	     //-- HLT theshold check
	     for(unsigned hltiobj = 0; hltiobj < Event->nHLTObj(ihltj[itrig]); hltiobj++) {
	       //-- cout<<"HLT object id: "<<hltiobj<<", pT = "<< Event->hltobj(ihltj[itrig],hltiobj).pt()<<", pT th = "<< HLTth[itrig]<<endl;
	       if(Event->hltobj(ihltj[itrig],hltiobj).pt() > HLTth[itrig]) hltcut[itrig] = true;
	     } 	       
	
	   }
	 }
       	 
	 //-- find leading pT jet
	 for(unsigned int ijet = 0; ijet < Event->nPFJetsCHS(); ijet++) {

	   double pt = Event->pfjetchs(ijet).ptCor();
	   double eta = Event->pfjetchs(ijet).eta();
	   bool tightid = Event->pfjetchs(ijet).tightID();
	   
	   bool jetid_cen = tightid;
	   bool jetid_fwd = (Event->pfjetchs(ijet).nemf() < 0.90) && (Event->pfjetchs(ijet).ncand() > 10);

	   if(pt < 20) continue;

	   if((fabs(eta) <= 3.0) && !jetid_cen) continue;
	   if((fabs(eta) > 3.2) && (fabs(eta) <= 4.7) && !jetid_fwd) continue;

	   if (pt < leading_pt) continue; 

	   leading_pt = pt; 
	   leading_eta = eta; 
	 }
	   
  	 hweight = 1.;
	 l1cut[0] = true;	 

	 if(fabs(leading_eta)<1.5){	   

	   for(int itrig = 0; itrig < Ntrigger; itrig++) {
	     
	     if(hltPassj[itrig]){
	       hleading_pt_all_Jet[itrig]->Fill(leading_pt,hweight);
	       hleading_eta_all_Jet[itrig]->Fill(leading_eta,hweight);	  
	       
	       if(!(l1cut[itrig] && hltcut[itrig])) continue;
	       
	       hleading_pt_emulated_Jet[itrig]->Fill(leading_pt,hweight);
	       hleading_eta_emulated_Jet[itrig]->Fill(leading_eta,hweight);
	     }
	   }
	 }
       } //-- end trigger information


       //-- PF jets - data and MC 
       
       unsigned n_PFJets = Event->nPFJetsCHS();		              
       if(mprintOk == 1) cout<<"Number of PFJets = "<<n_PFJets<<endl;
       
       //-- loop over jets
       for(unsigned int ijet=0; ijet < n_PFJets; ijet++) {

	 double yjet = Event->pfjetchs(ijet).y();
	 double phijet = Event->pfjetchs(ijet).phi();
	 
	 bool jetid_cen = Event->pfjetchs(ijet).tightID();
	 bool jetid_fwd = (Event->pfjetchs(ijet).nemf() < 0.90) && (Event->pfjetchs(ijet).ncand() > 10);
	 
	 double ptjetcorr = Event->pfjetchs(ijet).ptCor();
	 double ptjetraw = Event->pfjetchs(ijet).pt();
	 double jetptunc = Event->pfjetchs(ijet).unc();
	 
	 if(mprintOk == 1)  cout<<ijet<<": pT = "<<ptjetcorr<<", y = "<<yjet<<", phi = "<<phijet<<endl;
	 
	 if(ptjetcorr<mMinPt) continue;
	 
	 if(fabs(yjet)>mYMax) continue;
	 
	 //-- loop over trigger
	 for(int itrig = 0; itrig < Ntriggeredge-1; itrig++) {
	   
	   //-- cout<<"range trigger "<<itrig+1<<" = "<<triggeredge[itrig]<<" - "<<triggeredge[itrig+1]<<": trigger = "<<HLTJetName[itrig]<<endl;
	   //-- getchar();
	   
	   //-- MC - no trigger condition
	   if(mIsMCarlo) hltPassj[itrig] = true; 
	   if(!mIsMCarlo) hweight = 1;
	   
	   //-- trigger condition
	   if(Event->pfjetchs(0).ptCor()>= triggeredge[itrig] && Event->pfjetchs(0).ptCor() < triggeredge[itrig+1] && hltPassj[itrig]) {
	     
	     //-- prescale only applied to data
	     if(!mIsMCarlo) hweight*=TMath::Abs(prescalej[itrig]);

	     //-- full y
	     if((fabs(yjet) <= 3.0 && jetid_cen) ||(fabs(yjet) > 3.2 && fabs(yjet) <= 4.7 && jetid_fwd)) {
	       
	       DetJets++;
	       DETjet_ok[ijet] = 1;       
	       
	       hpt0_DETJet->Fill(ptjetcorr,hweight);
	       hpt0_DETJetUncor->Fill(ptjetraw,hweight);
	       hy0_DETJet->Fill(yjet,hweight);
	       hphi0_DETJet->Fill(phijet,hweight);
	       
	     } //-- end full y                                                     
	     
	       //-- loop y bin
	     for(int iy = 0; iy < Nyedge-1; iy++) {
	       
	       if(iy == 6) continue;
	       
	       int ihisto = 0;
	       if(iy <= 5) ihisto = iy;
	       else ihisto = iy - 1;
	       
	       bool jetid = false;
	       if(iy <= 5) jetid = jetid_cen;
	       else jetid = jetid_fwd;
	       
	       //-- y bin condition
	       if(fabs(yjet) > yedge[iy] && fabs(yjet) <= yedge[iy+1] && jetid) {
		 
		 //-- cout<<"y bin = "<<yedge[iy]<<" - "<<yedge[iy+1]<<endl;
		 
		 hpt_DETJet[ihisto]->Fill(ptjetcorr,hweight);
		 hpt_DETJetUncor[ihisto]->Fill(ptjetraw,hweight);
		 hpt_DETJetUP[ihisto]->Fill(ptjetcorr*(1+jetptunc),hweight);
		 hpt_DETJetDOWN[ihisto]->Fill(ptjetcorr*(1-jetptunc),hweight);
	   	 
	       } //-- end y bin condition
	     } //-- end loop y bin
	     
	   } //-- end trigger condition
	 } //-- end loop over trigger
	 
       } //-- end loop over jets
       
       if(DETjet_ok[0] == 1) {
	 
	 if(mIsMCarlo)  hTruePileUpMCInteger->Fill(TruePileUpInteger,hweight);
	 hMultiplicity_DETJet->Fill(DetJets,hweight);
	 hnum_of_VtxGood->Fill(Event->evtHdr().nVtxGood(),hweight);
	 
	 if(fabs(Event->pfjetchs(0).y()) <= 3.0 && Event->pfjetchs(0).tightID()) {
	   hMET_DET->Fill(Event->pfmet().met(),hweight); 
	   hMETPhi_DET->Fill(fabs(Event->pfmet().phi()),hweight); 
	   hFractionMET_DET->Fill(Event->pfmet().met_o_sumet(),hweight); 
	 }
	 
	 if(fabs(Event->pfjetchs(0).y()) > 3.2 && fabs(Event->pfjetchs(0).y()) <= 4.7){
	   if(Event->pfjetchs(0).nemf()<0.90 && Event->pfjetchs(0).ncand()>10) {
	     hMET_DET->Fill(Event->pfmet().met(),hweight); 
	     hMETPhi_DET->Fill(fabs(Event->pfmet().phi()),hweight); 
	     hFractionMET_DET->Fill(Event->pfmet().met_o_sumet(),hweight); 
	   }
	 }
	 
       }

     } //-- end loop over tree entries
     
     mInf->Close();
     
  } //-- end file loop
  
  for(int itrig = 0; itrig < Ntrigger; itrig++) {
    hleading_pt_HLTeffi_Jet[itrig]->Divide(hleading_pt_emulated_Jet[itrig],hleading_pt_all_Jet[itrig],1.,1.,"B");
    hleading_eta_HLTeffi_Jet[itrig]->Divide(hleading_eta_emulated_Jet[itrig],hleading_eta_all_Jet[itrig],1.,1.,"B");
  }
  
  //-- loop y bin						   
  for(int iy = 0; iy < Nybin; iy++) {			   
    hpt_GENJetStability[iy]->Divide(hpt_GENJetMatchedSame[iy],hpt_GENJetMatched[iy],1.,1.,"B");
    hpt_GENJetAcceptance[iy]->Divide(hpt_GENJetMatched[iy],hpt_GENJet[iy],1.,1.,"B");
    hpt_GENJetMiss[iy]->Divide(hpt_GENJetUnMatched[iy],hpt_GENJet[iy],1.,1.,"B");

    hpt_DETJetPurity[iy]->Divide(hpt_DETJetMatchedSame[iy],hpt_DETJetMatched[iy],1.,1.,"B");
    hpt_DETJetFake[iy]->Divide(hpt_DETJetUnMatched[iy],hpt_DETJet[iy],1.,1.,"B");
  }

} //-- end analyze()

InclusiveJet::~InclusiveJet() { }

DEFINE_FWK_MODULE(InclusiveJet);



