#include "optimizeCuts.h"

int optimizeCuts(const char *fileNameIn, // name of the input root file
		 const string outputFolder, // name of the output folder (must exist already)
		 const string flag, // flag for the output plot names
		 const bool drawPlots = false, // if true, plots are also individually saved as .png, .eps and .pdf
		 const bool holdPlots = false){ // if true leaves the TCanvas open at the end of the script

  gErrorIgnoreLevel = kWarning;

  // opening logfile
  string logfileName = outputFolder + flag + "-cuts.log";
  ofstream logfile(logfileName.c_str());
  if(logfile == NULL){
    cout << "- ERROR!!! cannot open logfile " << logfileName << endl;
    return 1;
  }
  logfile << " - input file = " << fileNameIn << endl;
  logfile << " - output folder = " << outputFolder << endl;

  // opening input file
  logfile << "+++++++  opening input file " << fileNameIn << endl;
  TFile *fileIn = TFile::Open(fileNameIn);
  if(fileIn == NULL){
    cout << " - ERROR!!! - cannot open file " << fileNameIn << endl;
    logfile << " -  - ERROR!!! - cannot open file " << fileNameIn << endl;
    return 1;
  }
  logfile << " - input file correctly open" << endl;

  // getting tree
  const char *treeName = "Plane0/Hits";
  logfile << "+++++++  getting tree " << treeName << endl;
  TTree *hits = (TTree *) fileIn -> Get(treeName);
  if(hits == NULL){
    cout << " - ERROR!!! - cannot get tree " << treeName << endl;
    logfile << " -  - ERROR!!! - cannot get tree " << treeName << endl;
    return 1;
  }
  logfile << " - tree correctly loaded" << endl;
  const unsigned int nEntries = hits -> GetEntries();
  cout << " - nEntries = " << nEntries << endl;
  logfile << " - nEntries = " << nEntries << endl;
  if(nEntries < HIGHSTATISTICSLIMIT){
    logfile << "WARNING: low statistics run" << endl;
  }

  // defining variables
  logfile << "+++++++  defining tree variables" << endl;
  Int_t    numHits = 0;
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
  logfile << "+++++++  defining tree branches" << endl;
  TBranch* bNumHits = NULL;
  TBranch* bHitPixX = NULL;
  TBranch* bHitPixY = NULL;
  TBranch* bHitPosX = NULL;
  TBranch* bHitPosY = NULL;
  TBranch* bHitPosZ = NULL;
  TBranch* bHitValue = NULL;
  TBranch* bHitT0 = NULL;
  TBranch* bHitTiming = NULL;
  TBranch* bHitChi2 = NULL;
  TBranch* bHitIsHit = NULL;
  TBranch* bHitValidFit = NULL;

  // setting branch addresses
  logfile << "+++++++  setting branch addresses" << endl;
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
  logfile << "+++++++  allocating 1D histograms" << endl;
  TH1F **h1_T0 = allocateH1Array(MAX_HITS, "h1_T0", "channel %d",
				 T0NBINS, T0MIN, T0MAX, "T_{0} [ns]", "#left[ ns^{-1} #right]",
				 5,
				 logfile);
  TH1F **h1_T0_TimingCut = allocateH1Array(MAX_HITS, "h1_T0_TimingCut", "channel %d",
					   T0NBINS, T0MIN, T0MAX, "T_{0} [ns]", "#left[ ns^{-1} #right]",
					   3,
					   logfile);
  TH1F **h1_exclusionLeft_T0 = allocateH1Array(MAX_HITS, "h1_exclusionLeft_T0", "channel %d",
					       (T0NBINS+1)*10., T0MIN, T0MAX, "T_{0} [ns]", "#left[ ns^{-1} #right]",
					       1,
					       logfile);
  TH1F **h1_exclusionRight_T0 = allocateH1Array(MAX_HITS, "h1_exclusionRight_T0", "channel %d",
						(T0NBINS+1)*10., T0MIN, T0MAX, "T_{0} [ns]", "#left[ ns^{-1} #right]",
						1,
						logfile);
  TH1F **h1_Charge = allocateH1Array(MAX_HITS, "h1_Charge", "channel %d",
				     CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]", "#left[ mV^{-1} #right]",
				     5,
				     logfile);
  TH1F **h1_Charge_T0Cut_bg = allocateH1Array(MAX_HITS, "h1_Charge_T0Cut_bg", "channel %d",
					      CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]", "#left[ mV^{-1} #right]",
					      3,
					      logfile);
  TH1F **h1_Charge_bgSubtracted = allocateH1Array(MAX_HITS, "h1_Charge_bgSubtracted", "channel %d",
						  CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]", "#left[ mV^{-1} #right]",
						  5,
						  logfile);
  TH1F **h1_Charge_bgSubtracted_filtered = allocateH1Array(MAX_HITS, "h1_Charge_bgSubtracted_filtered", "channel %d",
							      CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [mV]", "#left[ mV^{-1} #right]",
							      5,
							      logfile);
  TH1F **h1_exclusion_Charge = allocateH1Array(MAX_HITS, "h1_exclusionLeft_CHARGE", "channel %d", 
					       (CHARGENBINS+1)*10., CHARGEMIN, CHARGEMAX, "Charge [mV]",  "#left[ mV^{-1} #right]",
					       1,
					       logfile);
  TH1F **h1_Timing = allocateH1Array(MAX_HITS, "h1_Timing", "channel %d",
				     TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]", "#left[ ns^{-1} #right]",
				     5,
				     logfile);
  TH1F **h1_Timing_T0Cut_ChargeCut = allocateH1Array(MAX_HITS, "h1_Timing_T0Cut_ChargeCut", "channel %d",
						     TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]", "#left[ ns^{-1} #right]",
						     3,
						     logfile);
  TH1F **h1_exclusionLeft_Timing = allocateH1Array(MAX_HITS, "h1_exclusionLeft_Timing", "channel %d",
						   (TIMINGNBINS+1)*10., TIMINGMIN, TIMINGMAX, "Timing [ns]", "#left[ ns^{-1} #right]",
						   1,
						   logfile);
  TH1F **h1_exclusionRight_Timing = allocateH1Array(MAX_HITS, "h1_exclusionRight_Timing", "channel %d",
						    (TIMINGNBINS+1)*10., TIMINGMIN, TIMINGMAX, "Timing [ns]", "#left[ ns^{-1} #right]",
						    1,
						    logfile);

  // allocating 2D histograms
  logfile << "+++++++  allocating 2D histograms" << endl;
  TH2F **h2_Timing_vs_T0 = allocateH2Array(MAX_HITS, "h2_Timing_vs_T0_%d", "channel %d",
					   T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
					   TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]",
					   logfile);
  
  TH2F **h2_Charge_vs_T0 = allocateH2Array(MAX_HITS, "h2_Charge_vs_T0_%d", "channel %d",
					   T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
					   CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [V]",
					   logfile);

  TH2F **h2_Charge_vs_T0_bgSubtracted = allocateH2Array(MAX_HITS, "h2_Charge_vs_T0_bgSubtracted_%d", "channel %d",
							T0NBINS, T0MIN, T0MAX, "T_{0} [ns]",
							CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [V]",
							logfile);
  
  TH2F **h2_Charge_vs_Timing = allocateH2Array(MAX_HITS, "h2_Charge_vs_Timing_%d", "channel %d",
					       TIMINGNBINS, TIMINGMIN, TIMINGMAX, "Timing [ns]",
					       CHARGENBINS, CHARGEMIN, CHARGEMAX, "Charge [V]",
					       logfile);
  

  // setting preliminary timing cuts
  logfile << "+++++++  setting preliminary timing cuts" << endl;
  double cut_Timing[MAX_HITS];
  setTimingCuts(cut_Timing, 
		logfile);
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "timing cut = " << cut_Timing[iHit] << " ns"
	 << endl;
    logfile << " - channel " << iHit << ": "
	    << "timing cut = " << cut_Timing[iHit] << " ns"
	    << endl;
  }

  // first loop
  logfile << "+++++++  first loop" << endl;
  logfile << " - will fill following plots:" << endl;
  logfile << "\t- correlation plot timing vs T0" << endl;
  logfile << "\t- correlation plot charge vs T0" << endl;
  logfile << "\t- correlation plot charge vs timing" << endl;
  logfile << "\t- T0 distribution" << endl;
  logfile << "\t- T0 distribution after preliminary timing cut" << endl;
  logfile << "\t- charge distribution" << endl;
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

  // optimizing T0 cut
  logfile << "+++++++  optimizing T0 cut" << endl;
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
			     cutBgRejectionErr_T0,
			     logfile);
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "T0 cut low = " << cutLow_T0[iHit] << ", "
	 << "high = " << cutHigh_T0[iHit] << ", "
	 << "nEvents total = (" << cutNEventsTotal_T0[iHit] << "+-" << cutNEventsTotalErr_T0[iHit] << "), "
	 << "nEvents background = (" << cutNEventsBg_T0[iHit] << "+-" << cutNEventsBgErr_T0[iHit] << "), "
	 << "nEvents signal = (" << cutNEventsSignal_T0[iHit] << "+-" << cutNEventsSignalErr_T0[iHit] << "), "
	 << "bgRejection = (" << cutBgRejection_T0[iHit] << "+-" << cutBgRejectionErr_T0[iHit] << ")"
	 << endl;
    logfile << " - channel " << iHit << ": "
	    << "T0 cut low = " << cutLow_T0[iHit] << ", "
	    << "high = " << cutHigh_T0[iHit] << ", "
	    << "nEvents total = (" << cutNEventsTotal_T0[iHit] << "+-" << cutNEventsTotalErr_T0[iHit] << "), "
	    << "nEvents background = (" << cutNEventsBg_T0[iHit] << "+-" << cutNEventsBgErr_T0[iHit] << "), "
	    << "nEvents signal = (" << cutNEventsSignal_T0[iHit] << "+-" << cutNEventsSignalErr_T0[iHit] << "), "
	    << "bgRejection = (" << cutBgRejection_T0[iHit] << "+-" << cutBgRejectionErr_T0[iHit] << ")"
	    << endl;
  }

  // second loop
  logfile << "+++++++  second loop" << endl;
  logfile << " - will fill following plots:" << endl;
  logfile << "\t- charge distribution after T0 cut" << endl;
  for(unsigned int iEntry=0; iEntry<nEntries; iEntry++){
    hits -> GetEntry(iEntry);
    for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
      if(!(hitT0[iHit] > cutLow_T0[iHit] && hitT0[iHit] < cutHigh_T0[iHit])) h1_Charge_T0Cut_bg[iHit] -> Fill(hitValue[iHit]);
    }
  }

  // optimizing charge cut
  logfile << "+++++++  optimizing charge cut" << endl;
  double cutChargeSharing_Charge[MAX_HITS];
  double cutNEventsBg_Charge[MAX_HITS];
  double cutNEventsBgErr_Charge[MAX_HITS];
  double cutNEventsSignal_Charge[MAX_HITS];
  double cutNEventsSignalErr_Charge[MAX_HITS];
  double cutBgRejection_Charge[MAX_HITS];
  double cutBgRejectionErr_Charge[MAX_HITS];
  TF1 **f_Charge = optimizeCutCharge(MAX_HITS,
				     h1_Charge_T0Cut_bg,
				     h2_Charge_vs_T0,
				     h2_Charge_vs_T0_bgSubtracted,
				     h1_Charge_bgSubtracted,
				     h1_Charge_bgSubtracted_filtered,
				     cutLow_T0,
				     cutHigh_T0,
				     cutChargeSharing_Charge,
				     cutNEventsBg_Charge,
				     cutNEventsBgErr_Charge,
				     cutNEventsSignal_Charge,
				     cutNEventsSignalErr_Charge,
				     cutBgRejection_Charge,
				     cutBgRejectionErr_Charge,
				     logfile);
  if(f_Charge == NULL){
    cout << " - ERROR!!! cannot optimize charge cut" << endl;
    return 1;
  }
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "Charge cut charge sharing = " << cutChargeSharing_Charge[iHit] << ", "
	 << "nEvents background = (" << cutNEventsBg_Charge[iHit] << "+-" << cutNEventsBgErr_Charge[iHit] << "), "
	 << "nEvents signal = (" << cutNEventsSignal_Charge[iHit] << "+-" << cutNEventsSignalErr_Charge[iHit] << "), "
	 << "bgRejection = (" << cutBgRejection_Charge[iHit] << "+-" << cutBgRejectionErr_Charge[iHit] << ")"
	 << endl;
    logfile << " - channel " << iHit << ": "
	    << "Charge cut charge sharing = " << cutChargeSharing_Charge[iHit] << ", "
	    << "nEvents background = (" << cutNEventsBg_Charge[iHit] << "+-" << cutNEventsBgErr_Charge[iHit] << "), "
	    << "nEvents signal = (" << cutNEventsSignal_Charge[iHit] << "+-" << cutNEventsSignalErr_Charge[iHit] << "), "
	    << "bgRejection = (" << cutBgRejection_Charge[iHit] << "+-" << cutBgRejectionErr_Charge[iHit] << ")"
	    << endl;
  }

  // assigning charge cut
  logfile << "+++++++  assigning charge cut" << endl;
  double cut_Charge[MAX_HITS];
  if(CHARGECUTTYPE == 0){
    logfile << " - charge cut type: crossing point between charge sharing and signal distributions" << endl;
    for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
      cut_Charge[iHit] = cutChargeSharing_Charge[iHit];
    }
  }
  else{
    cout << " - ERROR!!! charge cut type undefined: " << CHARGECUTTYPE << endl;
    logfile << " - ERROR!!! charge cut type undefined: " << CHARGECUTTYPE << endl;
    return 1;
  }

  // third loop
  logfile << "+++++++  third loop" << endl;
  logfile << " - will fill following plots:" << endl;
  logfile << "\t- timing distribution" << endl;
  logfile << "\t- timing distribution after T0 and charge cuts" << endl;
  for(unsigned int iEntry=0; iEntry<nEntries; iEntry++){
    hits -> GetEntry(iEntry);
    for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
      h1_Timing[iHit] -> Fill(hitTiming[iHit]);
      if((hitT0[iHit] > cutLow_T0[iHit] && hitT0[iHit] < cutHigh_T0[iHit]) && hitValue[iHit] > cut_Charge[iHit]) h1_Timing_T0Cut_ChargeCut[iHit] -> Fill(hitTiming[iHit]);
    }
  }

  // closing input file
  logfile << "+++++++  closing input file" << endl;
  delete hits;
  fileIn -> Close();
  delete fileIn;

  // optimizing timing cut
  logfile << "+++++++  optimizing timing cut" << endl;
  double cutLow_Timing[MAX_HITS];
  double cutHigh_Timing[MAX_HITS];
  logfile << " - allocating quantile graphs" << endl;
  TGraph **gr_Timing_T0Cut_ChargeCut_quantiles = new TGraph*[MAX_HITS];
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    gr_Timing_T0Cut_ChargeCut_quantiles[iHit] = new TGraph();
  }
  optimizeCutTiming(MAX_HITS,
		    h1_Timing_T0Cut_ChargeCut, gr_Timing_T0Cut_ChargeCut_quantiles,
		    cutLow_Timing, cutHigh_Timing,
		    logfile);
  for(unsigned int iHit=0; iHit<MAX_HITS; iHit++){
    cout << " - channel " << iHit << ": "
	 << "Timing cut low = " << cutLow_Timing[iHit] << ", "
	 << "high = " << cutHigh_Timing[iHit]
	 << endl;
    logfile << " - channel " << iHit << ": "
	    << "Timing cut low = " << cutLow_Timing[iHit] << ", "
	    << "high = " << cutHigh_Timing[iHit]
	    << endl;
  }

  // saving plots
  logfile << "+++++++  saving plots" << endl;
  savePlots(MAX_HITS,
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
	    h1_Charge_bgSubtracted_filtered,
	    f_Charge,
	    cutChargeSharing_Charge,
	    h1_exclusion_Charge,
	    h1_Timing,
	    h1_Timing_T0Cut_ChargeCut,
	    gr_Timing_T0Cut_ChargeCut_quantiles,
	    cutLow_Timing,
	    cutHigh_Timing,
	    h1_exclusionLeft_Timing,
	    h1_exclusionRight_Timing,
	    outputFolder,
	    flag,
	    drawPlots,
	    holdPlots,
	    logfile);

  // cleaning memory
  logfile << "+++++++  cleaning memory" << endl;
  if(holdPlots) return 0;
  for(unsigned int i=0; i<MAX_HITS; i++){
    delete gr_Timing_T0Cut_ChargeCut_quantiles[i];
  }
  delete[] gr_Timing_T0Cut_ChargeCut_quantiles;
  deleteFArray(MAX_HITS, f_T0, logfile);
  //  deleteFArray(MAX_HITS, f_Charge, logfile); // memory leak here, must investigate
  deleteH1Array(MAX_HITS, h1_T0, logfile);
  deleteH1Array(MAX_HITS, h1_T0_TimingCut, logfile);
  deleteH1Array(MAX_HITS, h1_exclusionLeft_T0, logfile);
  deleteH1Array(MAX_HITS, h1_exclusionRight_T0, logfile);
  deleteH1Array(MAX_HITS, h1_Charge, logfile);
  deleteH1Array(MAX_HITS, h1_Charge_T0Cut_bg, logfile);
  deleteH1Array(MAX_HITS, h1_Charge_bgSubtracted, logfile);
  deleteH1Array(MAX_HITS, h1_Charge_bgSubtracted_filtered, logfile);
  deleteH1Array(MAX_HITS, h1_exclusion_Charge, logfile);
  deleteH1Array(MAX_HITS, h1_Timing, logfile);
  deleteH1Array(MAX_HITS, h1_Timing_T0Cut_ChargeCut, logfile);
  deleteH1Array(MAX_HITS, h1_exclusionLeft_Timing, logfile);
  deleteH1Array(MAX_HITS, h1_exclusionRight_Timing, logfile);
  deleteH2Array(MAX_HITS, h2_Timing_vs_T0, logfile);
  deleteH2Array(MAX_HITS, h2_Charge_vs_T0, logfile);
  deleteH2Array(MAX_HITS, h2_Charge_vs_T0_bgSubtracted, logfile);
  deleteH2Array(MAX_HITS, h2_Charge_vs_Timing, logfile);

  // closing logfile
  logfile.close();

  return 0;
}

int run904(){
  optimizeCuts("../../tmp/run_904_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_904_tj_W3R13_50um_6V_1DRS-cuts");
  return 0;
}

int run944(){
  optimizeCuts("../../tmp/run_944_tj_W3R13_20um_6V_1DRS.root", "plots/", "run_944_tj_W3R13_20um_6V_1DRS-cuts");
  return 0;
}

int loop(){

  optimizeCuts("../../tmp/run_802_tj_W3R15_50um_6V.root", "plots/", "run_802_tj_W3R15_50um_6V.root-cuts");
  //  optimizeCuts("../../tmp/run_803_tj_W3R15_50um_6V.root", "plots/", "run_803_tj_W3R15_50um_6V-cuts"); // low stat: 740
  //  optimizeCuts("../../tmp/run_801_tj_W3R15_50um_6V.root", "plots/", "run_801_tj_W3R15_50um_6V-cuts"); // low stat: 5071
  //  optimizeCuts("../../tmp/run_785_tj_W3R15_50um_6V_trigHits.root", "plots/", "run_785_tj_W3R15_50um_6V_trigHits-cuts"); // low stat: 68
  optimizeCuts("../../tmp/run_868_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_868_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_866_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_866_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_861_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_861_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_860_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_860_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_859_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_859_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_843_tj_W3R13_50um_6V.root", "plots/", "run_843_tj_W3R13_50um_6V-cuts");
  optimizeCuts("../../tmp/run_831_tj_W3R13_50um.root", "plots/", "run_831_tj_W3R13_50um-cuts");
  optimizeCuts("../../tmp/run_830_tj_W3R13_50um.root", "plots/", "run_830_tj_W3R13_50um-cuts");
  optimizeCuts("../../tmp/run_821_tj_W3R12_50um_6V_ROI_DUT.root", "plots/", "run_821_tj_W3R12_50um_6V_ROI_DUT-cuts");
  optimizeCuts("../../tmp/run_820_tj_W3R12_50um_6V_alignment.root", "plots/", "run_820_tj_W3R12_50um_6V_alignment-cuts");
  optimizeCuts("../../tmp/run_805_tj_W3R15_50um_3V.root", "plots/", "run_805_tj_W3R15_50um_3V-cuts");
  optimizeCuts("../../tmp/run_804_tj_W3R15_50um_6V.root", "plots/", "run_804_tj_W3R15_50um_6V-cuts");
  optimizeCuts("../../tmp/run_903_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_903_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_904_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_904_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_892_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_892_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_891_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_891_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_883_tj_W3R13_50um_6V_1DRS_alignment_ROI.root", "plots/", "run_883_tj_W3R13_50um_6V_1DRS_alignment_ROI-cuts");
  optimizeCuts("../../tmp/run_876_tj_W3R13_50um_6V_1DRS_alignment_ROI.root", "plots/", "run_876_tj_W3R13_50um_6V_1DRS_alignment_ROI-cuts");
  optimizeCuts("../../tmp/run_875_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_875_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_925_tj_W3R13_50um_3V_1DRS_noResetTriggerLogik.root", "plots/", "run_925_tj_W3R13_50um_3V_1DRS_noResetTriggerLogik-cuts");
  optimizeCuts("../../tmp/run_924_tj_W3R13_50um_3V_1DRS_noResetTriggerLogik.root", "plots/", "run_924_tj_W3R13_50um_3V_1DRS_noResetTriggerLogik-cuts");
  optimizeCuts("../../tmp/run_907_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_907_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_906_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_906_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_905_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_905_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_904_tj_W3R13_50um_6V_1DRS.root", "plots/", "run_904_tj_W3R13_50um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_942_tj_W3R13_20um_6V_1DRS.root", "plots/", "run_942_tj_W3R13_20um_6V_1DRS-cuts");
  optimizeCuts("../../tmp/run_930_tj_W3R13_20um_3V_1DRS_ROIFinding.root", "plots/", "run_930_tj_W3R13_20um_3V_1DRS_ROIFinding-cuts");
  optimizeCuts("../../tmp/run_926_tj_W3R13_50um_3V_1DRS.root", "plots/", "run_926_tj_W3R13_50um_3V_1DRS-cuts");
  //  optimizeCuts("../../tmp/run_951_tj_W3R13_20um_3V_1DRS.root", "plots/", "run_951_tj_W3R13_20um_3V_1DRS-cuts"); // low stat: 81
  optimizeCuts("../../tmp/run_949_tj_W3R13_20um_3V_1DRS.root", "plots/", "run_949_tj_W3R13_20um_3V_1DRS-cuts");
  optimizeCuts("../../tmp/run_938_tj_W3R13_20um_1V_1DRS.root", "plots/", "run_938_tj_W3R13_20um_1V_1DRS-cuts");
  optimizeCuts("../../tmp/run_944_tj_W3R13_20um_6V_1DRS.root", "plots/", "run_944_tj_W3R13_20um_6V_1DRS-cuts");

  return 0;

}
