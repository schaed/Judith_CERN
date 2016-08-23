#include "optimizeCuts.h"

int optimizeCuts(const char *fileNameIn,
		 const string outputFolder,
		 const string flag,
		 const bool holdPlots){

  gErrorIgnoreLevel = kWarning;
 
  // opening input file
  TFile *fileIn = TFile::Open(fileNameIn);
  if(fileIn == NULL){
    cout << " - ERROR!!! - cannot open file " << fileNameIn << endl;
    return 1;
  }

  // getting tree
  const char *treeName = "Plane0/Hits";
  TTree *hits = (TTree *) fileIn -> Get(treeName);
  if(hits == NULL){
    cout << " - ERROR!!! - cannot get tree " << treeName << endl;
    return 1;
  }

  // defining variables
  Int_t    numHits;
  Int_t    hitPixX[MAX_HITS];
  Int_t    hitPixY[MAX_HITS];
  Double_t hitPosX[MAX_HITS];
  Double_t hitPosY[MAX_HITS];
  Double_t hitPosZ[MAX_HITS];
  Double_t hitValue[MAX_HITS];
  Double_t hitT0[MAX_HITS];  
  Double_t hitTiming[MAX_HITS];
  Double_t hitChi2[MAX_HITS];
  Double_t hitIsHit[MAX_HITS];
  Double_t hitValidFit[MAX_HITS];

  // defining branches
  TBranch* bNumHits;
  TBranch* bHitPixX;
  TBranch* bHitPixY;
  TBranch* bHitPosX;
  TBranch* bHitPosY;
  TBranch* bHitPosZ;
  TBranch* bHitValue;
  TBranch* bHitT0;
  TBranch* bHitTiming;
  TBranch* bHitChi2;
  TBranch* bHitIsHit;
  TBranch* bHitValidFit;

  // setting branch addresses
  hits -> SetBranchAddress("NHits", &numHits, &bNumHits);
  hits -> SetBranchAddress("PixX", hitPixX, &bHitPixX);
  hits -> SetBranchAddress("PixY", hitPixY, &bHitPixY);
  hits -> SetBranchAddress("Value", hitValue, &bHitValue);
  hits -> SetBranchAddress("T0", hitT0, &bHitT0);
  hits -> SetBranchAddress("Timing", hitTiming, &bHitTiming);
  hits -> SetBranchAddress("PosX", hitPosX, &bHitPosX);
  hits -> SetBranchAddress("PosY", hitPosY, &bHitPosY);
  hits -> SetBranchAddress("PosZ", hitPosZ, &bHitPosZ);
  hits -> SetBranchAddress("IsHit", hitIsHit, &bHitIsHit);
  hits -> SetBranchAddress("ValidFit", hitValidFit, &bHitValidFit);
  hits -> SetBranchAddress("Chi2", hitChi2, &bHitChi2);

  // allocating 1D histograms
  TH1F **h1_T0 = allocateH1Array(MAX_HITS, "h1_T0", "channel %d",
				 T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
				 5);
  TH1F **h1_T0_TimingCut = allocateH1Array(MAX_HITS, "h1_T0_TimingCut", "channel %d",
					   T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
					   3);
  TH1F **h1_exclusionLeft_T0 = allocateH1Array(MAX_HITS, "h1_exclusionLeft_T0", "channel %d",
					       (T0NBINS+1)*10., T0MIN, T0MAX, "T_{0} [ns]",
					       1);
  TH1F **h1_exclusionRight_T0 = allocateH1Array(MAX_HITS, "h1_exclusionRight_T0", "channel %d",
						(T0NBINS+1)*10., T0MIN, T0MAX, "T_{0} [ns]",
						1);
  TH1F **h1_Charge = allocateH1Array(MAX_HITS, "h1_Charge", "channel %d",
				     CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]",
				     5);
  TH1F **h1_Charge_T0Cut_bg = allocateH1Array(MAX_HITS, "h1_Charge_T0Cut_bg", "channel %d",
					   CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]",
					   3);
  TH1F **h1_Charge_bgSubtracted = allocateH1Array(MAX_HITS, "h1_Charge_bgSubtracted", "channel %d",
						  CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]",
						  0);
  TH1F **h1_exclusion_Charge = allocateH1Array(MAX_HITS, "h1_exclusionLeft_CHARGE", "channel %d",
					       (CHARGENBINS+1)*10., CHARGEMIN, CHARGEMAX, "Charge [mV]",
					       1);
  TH1F **h1_Timing = allocateH1Array(MAX_HITS, "h1_Timing", "channel %d",
				     TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]",
				     5);
  TH1F **h1_Timing_T0Cut_ChargeCut = allocateH1Array(MAX_HITS, "h1_Timing_T0Cut_ChargeCut", "channel %d",
						     TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]",
						     3);

  // allocating 2D histograms
  TH2F **h2_Timing_vs_T0 = allocateH2Array(MAX_HITS, "h2_Timing_vs_T0_%d", "channel %d",
					   T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
					   TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]");
  
  TH2F **h2_Charge_vs_T0 = allocateH2Array(MAX_HITS, "h2_Charge_vs_T0_%d", "channel %d",
					   T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
					   CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [V]");

  TH2F **h2_Charge_vs_T0_bgSubtracted = allocateH2Array(MAX_HITS, "h2_Charge_vs_T0_bgSubtracted_%d", "channel %d",
							T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
							CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [V]");
  
  TH2F **h2_Charge_vs_Timing = allocateH2Array(MAX_HITS, "h2_Charge_vs_Timing_%d", "channel %d",
					       TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]",
					       CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [V]");
  

  // setting preliminary timing cuts
  double cut_Timing[MAX_HITS];
  setTimingCuts(cut_Timing);
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "timing cut = " << cut_Timing[iHit] << " ns"
	 << endl;
  }

  // first loop
  const unsigned int nEntries = hits -> GetEntries();
  cout << " - nEntries = " << nEntries << endl;
  for(unsigned int iEntry=0; iEntry<nEntries; iEntry++){
    hits -> GetEntry(iEntry);
    for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
      h2_Timing_vs_T0[iHit] -> Fill(hitT0[iHit], hitTiming[iHit]);
      h2_Charge_vs_T0[iHit] -> Fill(hitT0[iHit], hitValue[iHit]);
      h2_Charge_vs_Timing[iHit] -> Fill(hitTiming[iHit], hitValue[iHit]);
      h1_T0[iHit] -> Fill(hitT0[iHit]);
      if(hitTiming[iHit] > cut_Timing[iHit] && hitTiming[iHit] > TIMINGMIN && hitTiming[iHit] < TIMINGMAX) h1_T0_TimingCut[iHit] -> Fill(hitT0[iHit]);
      h1_Charge[iHit] -> Fill(hitValue[iHit]);
    }
  }

  // optimizing T0
  double cutLow_T0[MAX_HITS];
  double cutHigh_T0[MAX_HITS];
  double cutNEventsTotal_T0[MAX_HITS];
  double cutNEventsTotalErr_T0[MAX_HITS];
  double cutNEventsBg_T0[MAX_HITS];
  double cutNEventsBgErr_T0[MAX_HITS];
  double cutNEventsSignal_T0[MAX_HITS];
  double cutNEventsSignalErr_T0[MAX_HITS];
  double cutBgRejection_T0[MAX_HITS];
  double cutBgRejectionErr_T0[MAX_HITS];
  TF1 **f_T0 = optimizeCutT0(MAX_HITS,
			     h1_T0_TimingCut,
			     cutLow_T0,
			     cutHigh_T0,
			     cutNEventsTotal_T0,
			     cutNEventsTotalErr_T0,
			     cutNEventsBg_T0,
			     cutNEventsBgErr_T0,
			     cutNEventsSignal_T0,
			     cutNEventsSignalErr_T0,
			     cutBgRejection_T0,
			     cutBgRejectionErr_T0);
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "T0 cut low = " << cutLow_T0[iHit] << ", "
	 << "high = " << cutHigh_T0[iHit] << ", "
	 << "nEvents total = (" << cutNEventsTotal_T0[iHit] << "+-" << cutNEventsTotalErr_T0[iHit] << "), "
	 << "nEvents background = (" << cutNEventsBg_T0[iHit] << "+-" << cutNEventsBgErr_T0[iHit] << "), "
	 << "nEvents signal = (" << cutNEventsSignal_T0[iHit] << "+-" << cutNEventsSignalErr_T0[iHit] << "), "
	 << "bgRejection = (" << cutBgRejection_T0[iHit] << "+-" << cutBgRejectionErr_T0[iHit] << ")"
	 << endl;
  }

  // second loop
  for(unsigned int iEntry=0; iEntry<nEntries; iEntry++){
    hits -> GetEntry(iEntry);
    for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
      h1_Charge[iHit] -> Fill(hitValue[iHit]);
      if(!(hitT0[iHit] > cutLow_T0[iHit] && hitT0[iHit] < cutHigh_T0[iHit])) h1_Charge_T0Cut_bg[iHit] -> Fill(hitValue[iHit]);
    }
  }

  // optiminzing charge
  double cutChargeSharing_Charge[MAX_HITS];
  double cutNEventsBg_Charge[MAX_HITS];
  double cutNEventsBgErr_Charge[MAX_HITS];
  double cutNEventsSignal_Charge[MAX_HITS];
  double cutNEventsSignalErr_Charge[MAX_HITS];
  double cutBgRejection_Charge[MAX_HITS];
  double cutBgRejectionErr_Charge[MAX_HITS];
  TF1 **f_Charge = optimizeChargeCut(MAX_HITS,
				     h1_Charge_T0Cut_bg,
				     h2_Charge_vs_T0,
				     h2_Charge_vs_T0_bgSubtracted,
				     h1_Charge_bgSubtracted,
				     cutLow_T0,
				     cutHigh_T0,
				     cutChargeSharing_Charge,
				     cutNEventsBg_Charge,
				     cutNEventsBgErr_Charge,
				     cutNEventsSignal_Charge,
				     cutNEventsSignalErr_Charge,
				     cutBgRejection_Charge,
				     cutBgRejectionErr_Charge);
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "Charge cut charge sharing = " << cutChargeSharing_Charge[iHit] << ", "
	 << "nEvents background = (" << cutNEventsBg_Charge[iHit] << "+-" << cutNEventsBgErr_Charge[iHit] << "), "
	 << "nEvents signal = (" << cutNEventsSignal_Charge[iHit] << "+-" << cutNEventsSignalErr_Charge[iHit] << "), "
	 << "bgRejection = (" << cutBgRejection_Charge[iHit] << "+-" << cutBgRejectionErr_Charge[iHit] << ")"
	 << endl;
  }
  double cut_Charge[MAX_HITS];
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cut_Charge[iHit] = cutChargeSharing_Charge[iHit];
  }

  // third loop
  for(unsigned int iEntry=0; iEntry<nEntries; iEntry++){
    hits -> GetEntry(iEntry);
    for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
      if(hitTiming[iHit] <= 0.) continue;
      h1_Timing[iHit] -> Fill(hitTiming[iHit]);
      if((hitT0[iHit] > cutLow_T0[iHit] && hitT0[iHit] < cutHigh_T0[iHit]) && hitValue[iHit] > cut_Charge[iHit]) h1_Timing_T0Cut_ChargeCut[iHit] -> Fill(hitTiming[iHit]);
    }
  }

  // closing input file
  delete hits;
  fileIn -> Close();

  // drawing

  draw(MAX_HITS,
       h2_Timing_vs_T0,
       h2_Charge_vs_T0,
       h2_Charge_vs_Timing,
       h1_T0,
       h1_T0_TimingCut,
       h1_exclusionLeft_T0,
       h1_exclusionRight_T0,
       cut_Timing,
       T0NPAR,
       f_T0,
       cutLow_T0,
       cutHigh_T0,
       cutBgRejection_T0,
       cutBgRejectionErr_T0,
       h1_Charge,
       h1_Charge_T0Cut_bg,
       h2_Charge_vs_T0_bgSubtracted,
       h1_Charge_bgSubtracted,
       f_Charge,
       cutChargeSharing_Charge,
       h1_exclusion_Charge,
       h1_Timing,
       h1_Timing_T0Cut_ChargeCut,
       outputFolder,
       flag,
       holdPlots);

  if(holdPlots) return 0;

  // cleaning memory
  deleteFArray(MAX_HITS, f_T0);
  //  deleteFArray(MAX_HITS, f_Charge); // memory leak here, must investigate
  deleteH1Array(MAX_HITS, h1_T0);
  deleteH1Array(MAX_HITS, h1_T0_TimingCut);
  deleteH1Array(MAX_HITS, h1_exclusionLeft_T0);
  deleteH1Array(MAX_HITS, h1_exclusionRight_T0);
  deleteH1Array(MAX_HITS, h1_Charge);
  deleteH1Array(MAX_HITS, h1_Charge_T0Cut_bg);
  deleteH1Array(MAX_HITS, h1_Charge_bgSubtracted);
  deleteH1Array(MAX_HITS, h1_exclusion_Charge);
  deleteH1Array(MAX_HITS, h1_Timing);
  deleteH1Array(MAX_HITS, h1_Timing_T0Cut_ChargeCut);
  deleteH2Array(MAX_HITS, h2_Timing_vs_T0);
  deleteH2Array(MAX_HITS, h2_Charge_vs_T0);
  deleteH2Array(MAX_HITS, h2_Charge_vs_T0_bgSubtracted);
  deleteH2Array(MAX_HITS, h2_Charge_vs_Timing);
  delete fileIn;

  return 0;
}
