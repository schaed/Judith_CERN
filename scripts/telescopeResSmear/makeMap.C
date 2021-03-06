// all sizes in microns

#include "makeMap.h"

int makeMap(const double pixX = 50., // pixel size X
	    const double pixY = 50., // pixel size Y
	    const double resX = 8., // telescope resolution X
	    const double resY = 8., // telescope resolution Y
	    const double resErrX = 1., // uncertainty on telescope resolution X
	    const double resErrY = 1., // uncertainty on telescope resolution Y
	    const double gridX = 1., // size of the divisions of the grid in X
	    const double gridY = 1., // size of the divisions of the grid in Y
	    const unsigned int clusterType = CLUSTERTYPE_SINGLE, // cluster type: 1 = single pixel, 2 = 2 adjacent pixels, 3 = 3 pixels (L shape), 4 = 2x2 pixels block
	    const unsigned int integratorType = INTEGRATORTYPE_ADAPTIVE, // integrator type: 0 = adaptive, 1 = GSLMC (slow)
	    const  bool holdPlot = false){ // leaves the TCanvas open at the end of the script

  // opening output file
  char fileName[1000];
  sprintf(fileName,
	  "maps/map_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.root",
	  clusterType, pixX, pixY, resX, resY, gridX, gridY, integratorType);
  TFile *file = TFile::Open(fileName, "RECREATE");
  if(file == NULL){
    cout << " - ERROR!!! - cannot open file " << fileName << endl;
    return 1;
  }

  // defining sensor's surface
  double xMin = DEFAULT_DOUBLE;
  double xMax = DEFAULT_DOUBLE;
  double yMin = DEFAULT_DOUBLE;
  double yMax = DEFAULT_DOUBLE;
  if(defineSurface(xMin, xMax, yMin, yMax,
		   pixX, pixY,
		   clusterType)){
    cout << " - ERROR!!! - cannot define sensor surface" << endl;
    return 1;
  }
  cout << " - xMin = " << xMin << endl;
  cout << " - xMax = " << xMax << endl;
  cout << " - yMin = " << yMin << endl;
  cout << " - yMax = " << yMax << endl;

  // defining map's binning
  unsigned int nX = DEFAULT_UINT;
  unsigned int nY = DEFAULT_UINT;
  defineBinning(nX, nY,
		xMax, xMin, gridX,
		yMax, yMin, gridY);
  cout << " - nX = " << nX << endl;
  cout << " - nY = " << nY << endl;

  // defining pdf integrand
  TF2 *pdfArg = new TF2("pdfArg", 
			"exp( - ( (x-[0]) * (x-[0]) / (2. * [1] * [1]) + (y-[2]) * (y-[2]) / (2. * [3] * [3]) ) ) / (2. * TMath::Pi() * [1] * [3])",
			xMin, 
			xMax, 
			yMin, 
			yMax);

  // defining error propagation formulas
  TF2 *pdfArgDerivativeX = new TF2("pdfArgDerivativeX", 
				   "exp( - ( (x-[0]) * (x-[0]) / (2. * [1] * [1]) + (y-[2]) * (y-[2]) / (2. * [3] * [3]) ) ) * ( (x-[0]) * (x-[0]) / ([1] * [1]) - 1.) / (2. * TMath::Pi() * [1] * [3] * [3])",
				   xMin, 
				   xMax, 
				   yMin, 
				   yMax);
  TF2 *pdfArgDerivativeY = new TF2("pdfArgDerivativeY", 
				   "exp( - ( (x-[0]) * (x-[0]) / (2. * [1] * [1]) + (y-[2]) * (y-[2]) / (2. * [3] * [3]) ) ) * ( (y-[2]) * (y-[2]) / ([3] * [3]) - 1.) / (2. * TMath::Pi() * [1] * [1] * [3])",
				   xMin, 
				   xMax, 
				   yMin, 
				   yMax);

  // allocating map
  TH2F *map = new TH2F("map", 
		       "map",
		       nX, xMin-gridX/2., xMax-gridX/2.,
		       nY, yMin-gridY/2., yMax-gridY/2.);
  map -> SetDirectory(0);

  // allocating map of uncertainties
  TH2F *mapErr = new TH2F("mapErr", 
			  "mapErr",
			  nX, xMin, xMax,
			  nY, yMin, yMax);
  mapErr -> SetDirectory(0);

  // building maps
  if(buildMap(pixX, pixY, resX, resY, resErrX, resErrY, nX, nY, xMin, xMax, yMin, yMax,
	      clusterType,
   	      integratorType, 
   	      pdfArg, 
	      pdfArgDerivativeX, pdfArgDerivativeY,
	      map, mapErr)){
    cout << " - ERROR!!! - cannot build maps" << endl;
    return 1;
  }

  //////////
  // drawing
  //////////

  // style
  style();

  // pixel edges
  vector<TLine *> edges;
  if(getPixelEdges(edges,
		   clusterType,
		   pixX, pixY)){
    cout << " - ERROR!!! - cannot get pixel edges" << endl;
    return 1;
  }
  for(unsigned int i=0; i<edges.size(); i++){
    edges[i] -> SetLineWidth(3);
  }

  // axes
  vector<TLine *> axes;
  if(getAxes(axes,
	     clusterType,
	     pixX, pixY,
	     gridX, gridY)){
    cout << " - ERROR!!! - cannot get axes" << endl;
    return 1;
  }

  // map
  TCanvas *cMap = getCanvas(clusterType, 
			    "cMap");
  map -> GetXaxis() -> SetTitle("x [#mu]");
  map -> GetYaxis() -> SetTitle("y [#mu]");
  map -> GetXaxis() -> SetTitleOffset(1.2);
  map -> GetYaxis() -> SetTitleOffset(1.2);
  map -> SetMaximum(1.);
  map -> Draw("colz");
  for(unsigned int i=0; i<edges.size(); i++){
    edges[i] -> SetLineWidth(3);
    edges[i] -> Draw();
  }
  for(unsigned int i=0; i<axes.size(); i++){
    axes[i] -> Draw();
  }
  char plotName[1000];
  sprintf(plotName,
	  "maps/map_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.png",
	  clusterType, pixX, pixY, resX, resY, gridX, gridY, integratorType);
  cMap -> SaveAs(plotName);
  sprintf(plotName,
	  "maps/map_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.eps",
	  clusterType, pixX, pixY, resX, resY, gridX, gridY, integratorType);
  cMap -> SaveAs(plotName);
  sprintf(plotName,
	  "maps/map_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.pdf",
	  clusterType, pixX, pixY, resX, resY, gridX, gridY, integratorType);
  cMap -> SaveAs(plotName);

  // map of uncertainties
  TCanvas *cMapErr = getCanvas(clusterType, 
			       "cMapErr");
  cMapErr -> SetLogz();
  mapErr -> GetXaxis() -> SetTitle("x [#mu]");
  mapErr -> GetYaxis() -> SetTitle("y [#mu]");
  mapErr -> GetXaxis() -> SetTitleOffset(1.2);
  mapErr -> GetYaxis() -> SetTitleOffset(1.2);
  mapErr -> SetMaximum(0.1);
  mapErr -> SetMinimum(0.001);
  mapErr -> Draw("colz");
  for(unsigned int i=0; i<edges.size(); i++){
    edges[i] -> SetLineWidth(3);
    edges[i] -> Draw();
  }
  for(unsigned int i=0; i<axes.size(); i++){
    axes[i] -> Draw();
  }
  sprintf(plotName,
	  "maps/map_err_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_resErrX_%lf_resErrY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.png",
	  clusterType, pixX, pixY, resX, resY, resErrX, resErrY, gridX, gridY, integratorType);
  cMapErr -> SaveAs(plotName);
  sprintf(plotName,
	  "maps/map_err_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_resErrX_%lf_resErrY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.pdf",
	  clusterType, pixX, pixY, resX, resY, resErrX, resErrY, gridX, gridY, integratorType);
  cMapErr -> SaveAs(plotName);
  sprintf(plotName,
	  "maps/map_err_clusterType_%d_val_pixX_%lf_pixY_%lf_resX_%lf_resY_%lf_resErrX_%lf_resErrY_%lf_gridX_%lf_gridY_%lf_integratorType_%d.eps",
	  clusterType, pixX, pixY, resX, resY, resErrX, resErrY, gridX, gridY, integratorType);
  cMapErr -> SaveAs(plotName);

  // saving to file
  map -> Write();
  mapErr -> Write();
  cMap -> Write();
  cMapErr -> Write();
  file -> Close();

  if(holdPlot) return 0;

  // cleaning memory
  delete file;
  delete pdfArg;
  delete pdfArgDerivativeX;
  delete pdfArgDerivativeY;
  delete map;
  delete mapErr;
  for(unsigned int i=0; i<edges.size(); i++){
    delete edges[i];
  }
  edges.clear();
  for(unsigned int i=0; i<axes.size(); i++){
    delete axes[i];
  }
  axes.clear();
  delete cMap;
  delete cMapErr;
  return 0;
}

int loop(){

  const double pixX = 50.;
  const double pixY = 50.;
  const double resX = 8.;
  const double resY = 8.;
  const double resErrX = 1.;
  const double resErrY = 1.;
  const double gridX = 1.;
  const double gridY = 1.;
  const unsigned int integratorType = 0;
  for(unsigned int clusterType=1; clusterType<=4; clusterType++){
    cout << " +++++++ computing map for clusterType = " << clusterType << endl;
    if(makeMap(pixX, pixY, resX, resY, resErrX, resErrY, gridX, gridY, clusterType, integratorType)) return 1;
  }

  return 0;
}
