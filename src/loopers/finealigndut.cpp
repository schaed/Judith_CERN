#include "finealigndut.h"

#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>

#include <Rtypes.h>

#include "../storage/storageio.h"
#include "../storage/event.h"
#include "../mechanics/device.h"
#include "../mechanics/sensor.h"
#include "../mechanics/alignment.h"
#include "../processors/processors.h"
#include "../processors/clustermaker.h"
#include "../processors/trackmaker.h"
#include "../analyzers/singleanalyzer.h"
#include "../analyzers/dualanalyzer.h"
#include "../analyzers/cuts.h"
#include "../analyzers/dutresiduals.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

using std::cout;
using std::endl;

namespace Loopers {

void FineAlignDut::loop()
{
  for (unsigned int niter = 0; niter < _numIterations; niter++)
  {
    cout << "Iteration " << niter << " of " << _numIterations - 1 << endl;

    // Use broad residual resolution for the first iteration
    const unsigned int numPixX = (niter == 0) ? _numPixXBroad : _numPixX;
    const double binxPerPix = (niter == 0) ? _binsPerPixBroad : _binsPerPix;

    std::stringstream name; // Build name strings for each histo
      double my_rotation=0.0;
      std::vector<std::vector<double> > rot_residuals;
      for(unsigned rot=0; rot<1; ++rot){
	std::cout << "rot: " << rot << std::endl;

	// setting up the rotations
	if(rot>0){
	  my_rotation = (0.1 * (6.0 - double(rot)) / (1.0+double(niter)*3.0));
	  for (unsigned int nsens = 0; nsens < _dutDevice->getNumSensors(); nsens++){
	    Mechanics::Sensor* tmp_sensor = _dutDevice->getSensor(nsens);
	    tmp_sensor->setRotZ(tmp_sensor->getRotZ() + my_rotation);
	  }
	}
	// Makes or gets a directory called from inside _dir with this name
	name.str("");
	name << "_perm" << niter; 	
	
	float ndiv = 10000.0;
	int n_residuals = int(float(_endEvent-_startEvent)/ndiv)+1;
	std::vector<Analyzers::DUTResiduals> v_residual;
	std::vector<double> v_timeStamp;
	// Residuals of DUT clusters to ref. tracks
	Analyzers::DUTResiduals residuals(_refDevice, _dutDevice, (rot>0) ? 0 : _dir, name.str().c_str(), numPixX, binxPerPix, _numBinsY);

	for(unsigned i=0; i<n_residuals; ++i){
	  Analyzers::DUTResiduals one_residual(_refDevice, _dutDevice, (rot>0) ? 0 : _dir, name.str().c_str(), numPixX, 10, _numBinsY);
	  v_residual.push_back(one_residual);
	  //long unsigned s=0_dutDevice->getTimeStart();
	  v_timeStamp.push_back(0.0);
	}
	
	// determine binning versus time stamp.
	// fill in here
	
	// Use events with only 1 track
	Analyzers::Cuts::EventTracks* cut1 =
	  new Analyzers::Cuts::EventTracks(1, Analyzers::EventCut::EQ);

	// Use tracks with one hit in each plane
	const unsigned int numClusters = _refDevice->getNumSensors();
	Analyzers::Cuts::TrackClusters* cut2 =
	  new Analyzers::Cuts::TrackClusters(numClusters, Analyzers::TrackCut::EQ);
	
	// Use tracks with low chi2
	//    Analyzers::Cuts::TrackClusters* cut3 =
	//        new Analyzers::Cuts::TrackChi2(1, Analyzers::TrackCut::LT);

	// Note: the analyzer will delete the cuts
	residuals.addCut(cut1);
	residuals.addCut(cut2);

	for(unsigned i=0; i<n_residuals; ++i){

	  // Use events with only 1 track
	  Analyzers::Cuts::EventTracks* cut1b =
	    new Analyzers::Cuts::EventTracks(1, Analyzers::EventCut::EQ);

	  // Use tracks with one hit in each plane
	  //const unsigned int numClusters = _refDevice->getNumSensors();
	  Analyzers::Cuts::TrackClusters* cut2b =
	    new Analyzers::Cuts::TrackClusters(numClusters, Analyzers::TrackCut::EQ);
	  
	  v_residual.at(i).addCut(cut1b);
	  v_residual.at(i).addCut(cut2b);
	}
	std::vector<std::string> tmp_funcX,tmp_funcY,tmp_funcTimeX,tmp_funcTimeY;
	unsigned last_res=0,my_evt_res=0;
	for (ULong64_t nevent = _startEvent; nevent <= _endEvent; nevent++)
	  {
	    Storage::Event* refEvent = _refStorage->readEvent(nevent);
	    Storage::Event* dutEvent = _dutStorage->readEvent(nevent);

	    if (refEvent->getNumClusters() || dutEvent->getNumClusters())
	      throw "FineAlignDut: can't recluster an event, mask the tree in the input";
	    for (unsigned int nplane = 0; nplane < refEvent->getNumPlanes(); nplane++)
	      _clusterMaker->generateClusters(refEvent, nplane);
	    for (unsigned int nplane = 0; nplane < dutEvent->getNumPlanes(); nplane++)
	      _clusterMaker->generateClusters(dutEvent, nplane);

	    // apply the existing position corrections
	    
	    for (unsigned int nsens = 0; nsens < _dutDevice->getNumSensors(); ++nsens){
	      Mechanics::Sensor* sensor = _dutDevice->getSensor(nsens);
	      sensor->setFrameNumber(refEvent->getFrameNumber());
	      sensor->setTimeStamp(refEvent->getTimeStamp());

	      if(!_DOITERATIVEPOSCORR){		
		sensor->setFrameNumberFuncX("0.0");
		sensor->setFrameNumberFuncY("0.0");
		sensor->setTimeStampFuncX("0.0");
		sensor->setTimeStampFuncY("0.0");
	      }
	      // save the alignment dependency for each sensor. Only want to fill for the first event
	      if(nevent==_startEvent){
		tmp_funcX.push_back(sensor->getFrameNumberFuncX());
		tmp_funcY.push_back(sensor->getFrameNumberFuncY());
		tmp_funcTimeX.push_back(sensor->timeStampFuncX());
		tmp_funcTimeY.push_back(sensor->timeStampFuncY());
	      }
	    }
	    
	    Processors::applyAlignment(refEvent, _refDevice);
	    Processors::applyAlignment(dutEvent, _dutDevice);

	    if (refEvent->getNumTracks())
	      throw "FineAlign: can't re-track an event, mask the tree in the input";
	    _trackMaker->generateTracks(refEvent,
					_refDevice->getBeamSlopeX(),
					_refDevice->getBeamSlopeY());
	    
	    residuals.processEvent(refEvent, dutEvent);
	    int this_residual = int(float(nevent-_startEvent)/ndiv);
	    v_residual.at(this_residual).processEvent(refEvent, dutEvent);
	    v_timeStamp.at(this_residual) += double(refEvent->getTimeStamp()/1.0e8/double(ndiv));
	    ++my_evt_res;
	    if(this_residual!=last_res){
	      my_evt_res=0;
	      last_res=this_residual;
	    }
	    
	    progressBar(nevent);

	    delete refEvent;
	    delete dutEvent;
	  }
	int end_residual = int(float(_endEvent-1-_startEvent)/ndiv);
	if(my_evt_res<ndiv){
	  v_timeStamp.at(end_residual)*=double(ndiv)/double(my_evt_res);
	}

	std::vector<TH1D*> x_hists,y_hists;
	for (unsigned int nsens = 0; nsens < _dutDevice->getNumSensors(); ++nsens)
	  {
	    std::cout << "sensors: " << nsens << " out of " << _dutDevice->getNumSensors() << std::endl;
	    Mechanics::Sensor* sensor = _dutDevice->getSensor(nsens);
	    
	    // list the X and Y residual per frameNumber
	    for(unsigned i=0; i<n_residuals; ++i){
	      std::cout << "X int: " << i << " "  << v_residual.at(i).getResidualX(nsens)->GetMean() << " " << (v_timeStamp.at(i)/1.0e6) << std::endl;
	      std::cout << "Y int: " << i << " "  << v_residual.at(i).getResidualY(nsens)->GetMean()  << " " << (v_timeStamp.at(i)/1.0e6) << std::endl;
	      x_hists.push_back(v_residual.at(i).getResidualX(nsens));
	      y_hists.push_back(v_residual.at(i).getResidualY(nsens));	      
	    }
	    bool fitGaus=true;
	    std::string funx = Processors::fitPosition(x_hists,unsigned(ndiv),_displayFits, fitGaus);

	    if(_DOFRAMENUMBER) sensor->setFrameNumberFuncX(tmp_funcX.at(nsens)=="0.0"?funx :funx+"+"+tmp_funcX.at(nsens));
	    else sensor->setFrameNumberFuncX(tmp_funcX.at(nsens)=="0.0"?"0.0" :tmp_funcX.at(nsens));
	    std::string funy = Processors::fitPosition(y_hists,unsigned(ndiv),_displayFits,fitGaus);
	    if(_DOFRAMENUMBER) sensor->setFrameNumberFuncY(tmp_funcY.at(nsens)=="0.0"?funy :funy+"+"+tmp_funcY.at(nsens));	    
	    else sensor->setFrameNumberFuncY(tmp_funcY.at(nsens)=="0.0"?"0.0" :tmp_funcY.at(nsens));	    

	    // versus time stamp
	    std::string tfunx = Processors::fitPosition(x_hists,v_timeStamp,_displayFits);
	    if(!_DOFRAMENUMBER) sensor->setTimeStampFuncX(tmp_funcTimeX.at(nsens)=="0.0"?tfunx :tfunx+"+"+tmp_funcTimeX.at(nsens));
	    else sensor->setTimeStampFuncX(tmp_funcTimeX.at(nsens)=="0.0"?"0.0": tmp_funcTimeX.at(nsens));	    
	    std::string tfuny = Processors::fitPosition(y_hists,v_timeStamp,_displayFits);
	    if(!_DOFRAMENUMBER) sensor->setTimeStampFuncY(tmp_funcTimeY.at(nsens)=="0.0"?tfuny :tfuny+"+"+tmp_funcTimeY.at(nsens));
	    else sensor->setTimeStampFuncY(tmp_funcTimeY.at(nsens)=="0.0"?"0.0": tmp_funcTimeY.at(nsens));
	    for(unsigned it=0;it<x_hists.size(); ++it){
	      delete x_hists.at(it);
	      delete y_hists.at(it);
	      
	    }
	    x_hists.clear();
	    y_hists.clear();

	    // Framenumber
	    //residuals
	    /*
	    for (unsigned int nsens = 0; nsens < _dutDevice->getNumSensors(); ++nsens)
	      {
		std::cout << "sensors: " << nsens << " out of " << _dutDevice->getNumSensors() << std::endl;
		Mechanics::Sensor* sensor = _dutDevice->getSensor(nsens);
		TH2D* hrestmp = residuals.getResidualFrameX(nsens);
		Processors::fitPosition(hrestmp,_displayFits);
		TH2D* hrestmpy = residuals.getResidualFrameY(nsens);
		Processors::fitPosition(hrestmpy,_displayFits);		
	      }
	    */
	    //std::cout << "end fit pos " << std::endl;
	    //double offsetX = 0, offsetY = 0, rotation = 0;
	    /*
	    if(rot==0){
	      Processors::residualAlignment(residuals.getResidualXY(nsens),
					    residuals.getResidualYX(nsens),
					    offsetX, offsetY, rotation, _displayFits);
	      
	      sensor->setOffX(sensor->getOffX() + offsetX);
	      sensor->setOffY(sensor->getOffY() + offsetY);
	      sensor->setRotZ(sensor->getRotZ() + rotation);
	    } else{
	      double offsetX=0.0, sigmaX=0.0;
	      double offsetY=0.0, sigmaY=0.0;	      
	      
	      Processors::fitGaussian( residuals.getResidualY(nsens), offsetY, sigmaY,_displayFits);
	      Processors::fitGaussian( residuals.getResidualX(nsens), offsetX, sigmaX,_displayFits);	      
	      std::cout << "Checking the rotation of " << my_rotation << " and the total residual is " << residuals.GetTotalResidual(nsens)
			<< " xRMS: " << residuals.getResidualX(nsens)->GetRMS()
			<< " yRMS: " << residuals.getResidualY(nsens)->GetRMS()
			<< " xFit: " << sigmaX
			<< " yFit: " << sigmaY
			<< std::endl;
	      if(rot_residuals.size()<=nsens) {
		std::vector<double> dummy_vec;
		rot_residuals.push_back(dummy_vec);
	      }
	      //rot_residuals.at(nsens).push_back(residuals.GetTotalResidual(nsens));
	      //rot_residuals.at(nsens).push_back(residuals.getResidualX(nsens)->GetRMS());
	      rot_residuals.at(nsens).push_back(sigmaX);
	      //un-do the rotations
	      sensor->setRotZ(sensor->getRotZ() - my_rotation);
	    }
	    */
	    
	  } // end loop over sensors
	//std::cout << "end SENSOR loop: " << rot << std::endl;
      } // end loop over rotations
      //std::cout << "end ROT loop" << std::endl;
      //std::cout << std::flush;
      /*
      for (unsigned int nsens = 0; nsens < _dutDevice->getNumSensors(); nsens++){
	Mechanics::Sensor* sensor = _dutDevice->getSensor(nsens);
	// find the minimum of residuals
	double min_rot_residuals=-1.0;
	int iter_rot_residuals=-1;      
	for(unsigned hh=0; hh<rot_residuals.at(nsens).size(); ++hh){
	  if(hh==0) min_rot_residuals = rot_residuals.at(nsens).at(hh);
	  if(rot_residuals.at(nsens).at(hh)<min_rot_residuals){
	    min_rot_residuals=rot_residuals.at(nsens).at(hh);
	    iter_rot_residuals=hh;
	  }
	}
	if(iter_rot_residuals>0){
	  std::cout << "best: " << (0.1 * (6.0 - double(iter_rot_residuals+1)) / (1.0+double(niter)*3.0))
		    << " found: " << iter_rot_residuals
		    << std::endl;
	  // apply the best rotation
	  sensor->setRotZ(sensor->getRotZ() + (0.1 * (6.0 - double(iter_rot_residuals+1)) / (1.0+double(niter)*3.0)));
	}
      } // end loop over sensors for rotations
      */
  }// end loop over iterations
  //std::cout << "_dutDevice->getAlignment(): " << _dutDevice->getAlignment() << std::endl;
  _dutDevice->getAlignment()->writeFile();
}

void FineAlignDut::setNumIteratioins(unsigned int value) { _numIterations = value; }
void FineAlignDut::setNumBinsY(unsigned int value) { _numBinsY = value; }
void FineAlignDut::setNumPixX(unsigned int value) { _numPixX = value; }
void FineAlignDut::setBinsPerPix(double value) { _binsPerPix = value; }
void FineAlignDut::setNumPixXBroad(unsigned int value) { _numPixXBroad = value; }
void FineAlignDut::setBinsPerPixBroad(double value) { _binsPerPixBroad = value; }
void FineAlignDut::setDisplayFits(bool value) { _displayFits = value; }
void FineAlignDut::setFrameNumberCorr(bool value) { _DOFRAMENUMBER = value; }  
void FineAlignDut::setIterativePosCorr(bool value) { _DOITERATIVEPOSCORR = value; }  
void FineAlignDut::setRelaxation(double value) { _relaxation = value; }

FineAlignDut::FineAlignDut(Mechanics::Device* refDevice,
                           Mechanics::Device* dutDevice,
                           Processors::ClusterMaker* clusterMaker,
                           Processors::TrackMaker* trackMaker,
                           Storage::StorageIO* refInput,
                           Storage::StorageIO* dutInput,
                           ULong64_t startEvent,
                           ULong64_t numEvents,
                           Long64_t eventSkip,
			   TDirectory* dir) :
			   
  Looper(refInput, dutInput, startEvent, numEvents, eventSkip),
  _refDevice(refDevice),
  _dutDevice(dutDevice),
  _clusterMaker(clusterMaker),
  _trackMaker(trackMaker),
  _numIterations(1),
  _numBinsY(15),
  _numPixX(5),
  _binsPerPix(10),
  _numPixXBroad(20),
  _binsPerPixBroad(1),
  _displayFits(true),
  _DOFRAMENUMBER(true),
  _DOITERATIVEPOSCORR(false),
  _relaxation(0.8),
  _dir(dir)
{
  assert(refInput && dutInput && refDevice && dutDevice && clusterMaker && trackMaker &&
         "Looper: initialized with null object(s)");
  assert(refInput->getNumPlanes() == refDevice->getNumSensors() &&
         dutInput->getNumPlanes() == dutDevice->getNumSensors() &&
         "Loopers: number of planes / sensors mis-match");
}

}
