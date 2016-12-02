#include "efficiency.h"

#include <cassert>
#include <sstream>
#include <math.h>
#include <vector>
#include <float.h>

#include <TDirectory.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TError.h>
#include <TMath.h>

// Access to the device being analyzed and its sensors
#include "../mechanics/device.h"
#include "../mechanics/sensor.h"
// Access to the data stored in the event
#include "../storage/hit.h"
#include "../storage/cluster.h"
#include "../storage/plane.h"
#include "../storage/track.h"
#include "../storage/event.h"
// Some generic processors to calcualte typical event related things
#include "../processors/processors.h"
// This header defines all the cuts
#include "cuts.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

//using namespace std;
using std::cout;
using std::endl;

namespace Analyzers {

void Efficiency::processEvent(const Storage::Event* refEvent,
                              const Storage::Event* dutEvent)
{
  assert(refEvent && dutEvent && "Analyzer: can't process null events");

  // Throw an error for sensor / plane mismatch
  eventDeivceAgree(refEvent, dutEvent);

  // Check if the event passes the cuts
  for (unsigned int ncut = 0; ncut < _numEventCuts; ncut++)
    if (!_eventCuts.at(ncut)->check(refEvent)) return;
  
  //std::cout << "EFFFrameNumber DUT: " << dutEvent->getFrameNumber() << " Telescope: " << refEvent->getFrameNumber() << std::endl;
  //fill in the amplitude distribution histogram
  for (unsigned int nsensor = 0; nsensor < _dutDevice->getNumSensors(); nsensor++)
  {
    //std::cout << "F: " << dutEvent->getFrameNumber() << std::endl;
    _dutDevice->getSensor(nsensor)->setFrameNumber(dutEvent->getFrameNumber());
    _dutDevice->getSensor(nsensor)->setTimeStamp(refEvent->getTimeStamp());     // used the RCE time stamp for fitting
    for (unsigned int nhit = 0; nhit < dutEvent->getNumHits(); nhit++)
    {
      const Storage::Hit* hit = dutEvent->getHit(nhit);

      _amplDist.at(nsensor)->Fill( hit->getValue() );
    }
  }

  //
  // Fill occupancy based
  //
  bool pass_track_selection = false;
  for (unsigned int ntrack = 0; ntrack < refEvent->getNumTracks(); ntrack++)
  {
    Storage::Track* track = refEvent->getTrack(ntrack);

    // Check if the track passes the cuts
    bool pass = true;
    for (unsigned int ncut = 0; ncut < _numTrackCuts; ncut++)
      if (!_trackCuts.at(ncut)->check(track)) { pass = false; break; }
    if (!pass) continue;

    // Hack cut on the residuals for strange events.
    if(false){
      Storage::Plane* hack_plane = refEvent->getPlane(3);
      Mechanics::Sensor* hack_sensor = _refDevice->getSensor(3);
      double hack_tx = 0, hack_ty = 0, hack_tz = 0;
      Processors::trackSensorIntercept(track, hack_sensor, hack_tx, hack_ty, hack_tz);
      // 
      for (unsigned int ncluster = 0; ncluster < hack_plane->getNumClusters(); ncluster++)
      {
        Storage::Cluster* cluster = hack_plane->getCluster(ncluster);

        // Check if the cluster passes the cuts
        for (unsigned int ncut = 0; ncut < _numClusterCuts; ncut++)
          if (!_clusterCuts.at(ncut)->check(cluster)) continue;

        const double rx = hack_tx - cluster->getPosX();
	//if(rx<10.0) return; // smaller peak
	if(rx>100.0) return; // This is the larger peak
      }
    }// end hack cut

    pass_track_selection=true;
    

    for (unsigned int nsensor = 0; nsensor < _dutDevice->getNumSensors(); nsensor++)
    {
      Mechanics::Sensor* sensor = _dutDevice->getSensor(nsensor);

      double tx = -9999.0, ty = -9999.0, tz = -9999.0;
      double ex = -9999.0, ey = -9999.0;
      double dx = -9999.0, dy = -9999.0;
      Processors::trackSensorIntercept(track, sensor, tx, ty, tz);
      // Marco Edit
      Processors::trackError(track, tz, ex, ey);
      _posErrorX.at(nsensor)->Fill(ex);
      _posErrorY.at(nsensor)->Fill(ey);
      //std::cout << "    EFFTrack_to_sensor x: " << (tx - sensor->getOffX()) << " y: " << (ty - sensor->getOffY()) << std::endl;      
      // Fill Track occupancy at DUT
      _trackOcc.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY());
      //_trackResFine.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY());      
      _trackRes.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY());
      _trackResFine.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY());
      fout << dutEvent->getFrameNumber() << "," << (tx - sensor->getOffX()) << "," << (ty - sensor->getOffY()) << "\n";
       //-50 to -30 x and 110 to 130 y
      //-60 to -40 y and -40 to -60 x
      //if((tx - sensor->getOffX())>-50.0 && (tx - sensor->getOffX())<-30.0
      //	 && (ty - sensor->getOffY())>110.0 && (ty - sensor->getOffY())<130.0)
      //if((tx - sensor->getOffX())>-105.0 && (tx - sensor->getOffX())<-95.0
      //	 && (ty - sensor->getOffY())>-30.0 && (ty - sensor->getOffY())<-50.0)
      //	std::cout << "inside sensor: " << dutEvent->getFrameNumber()
      //		  << " " << dutEvent->getTimeStamp()
      //		  << std::endl;
      bool isInsideHit=false;
      // Fill track residual
      Storage::Plane* plane = dutEvent->getPlane(nsensor);

      for (unsigned int ncluster = 0; ncluster < plane->getNumClusters(); ncluster++){
	Storage::Cluster* cluster = plane->getCluster(ncluster);	  
	//Marco 2edit

	// cout << "Clsuter pixx: " << cluster->getPixX() << endl;
	// cout << "Clsuter pixy: " << cluster->getPixY() << endl;
	// cout << "Clsuter posx: " << cluster->getPosX() << endl;
	// cout << "Clsuter posy: " << cluster->getPosY() << endl;
        Processors::trackClusterDistance( track, cluster, sensor, dx, dy, ex, ey);
	 _dutTrackResX.at(nsensor)->Fill(dx);
	 _dutTrackResY.at(nsensor)->Fill(dy);
	
	 _projXvsFN.at(nsensor)->Fill(dutEvent->getFrameNumber(),dx);
	 _projYvsFN.at(nsensor)->Fill(dutEvent->getFrameNumber(),dy);
	// hit residual 1 D
	//_hitResidualNoCut.at(nsensor)->Fill(sqrt((tx - cluster->getPosX())*(tx - cluster->getPosX()) + (ty - cluster->getPosY())*(ty - cluster->getPosY())));  
	_hitResidualxNoCut.at(nsensor)->Fill(tx - cluster->getPosX());
	_hitResidualyNoCut.at(nsensor)->Fill(ty - cluster->getPosY());
     
	// Check if the cluster passes the cuts
	bool cluster_pass = true;
	for (unsigned int ncut = 0; ncut < _numClusterCuts; ncut++)
	  if (!_clusterCuts.at(ncut)->check(cluster)) { cluster_pass = false; break; }
	if (!cluster_pass) continue;
	//std::cout << "    EFFTrack_to_cluster x: " << (tx - cluster->getPosX()) << " y: " << (ty - cluster->getPosY()) << " charge: " << cluster->getValue()
	// << " T0: " << cluster->getT0() << " NHits: " << cluster->getNumHits()
	//	  << " collection_time: " << cluster->getTiming() << std::endl;      	
	//std::cout << "    EFFCluster_to_sensor x: " << (sensor->getOffX()) << " cluster: " <<  (cluster->getPosX()) << " y: " << sensor->getOffY() << " cluster: " << cluster->getPosY() << std::endl; 
	//std::cout << "getclustX: " << cluster->getPosX() << " tx: " << tx << std::endl;
	//std::cout << "getclustY: " << cluster->getPosY() << " ty: " << ty << std::endl;	
	//_trackResHit.at(nsensor)->Fill(tx - cluster->getPosX(), ty - cluster->getPosY());
	//_trackResHitFine.at(nsensor)->Fill(tx - cluster->getPosX(), ty - cluster->getPosY());	
	_trackResHit.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY());
	_trackResHitFine.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY());	
      if((tx - sensor->getOffX())>-100.0 && (tx - sensor->getOffX())<-70.0
	  && (ty - sensor->getOffY())>-20.0 && (ty - sensor->getOffY())<10.0){
	// && (ty - sensor->getOffY())>130.0 && (ty - sensor->getOffY())<200.0){
	//std::cout << "     HIT sensor: " << dutEvent->getFrameNumber() << " "
	//	    << dutEvent->getTimeStamp()
	//	    << " " << cluster->getValue() << std::endl;
	  isInsideHit=true;
      }

      //if((tx - sensor->getOffX())>220.0 && (tx - sensor->getOffX())<2000.0
      // && (ty - sensor->getOffY())>-4000.0 && (ty - sensor->getOffY())<4000.0){
	//std::cout << "     OUTSIDE HIT sensor: " << dutEvent->getFrameNumber() << " "
	//<< dutEvent->getTimeStamp()
	//<< " " << cluster->getValue() << std::endl;
      // }
      /*
      std::cout << "cluster pix X: " << cluster->getPixX() << " Y: " <<  cluster->getPixY()
		<< " PoxX: " <<  cluster->getPosX() << " PoxY: " <<  cluster->getPosY() << " Nhit: " << cluster->getNumHits() << std::endl;
      for(unsigned s=0; s<cluster->getNumHits(); ++s)
	std::cout << "    hit: " << s << " pixX: " << cluster->getPixX() << " Y: " <<  cluster->getPixY()
		  << " PosX " <<  cluster->getHit(s)->getPosX() << " PoxY: " <<  cluster->getHit(s)->getPosY() << std::endl;
      */
      //std::cout << "cluster->getNumHits(): " << cluster->getNumHits() <<std::endl;
      _trackResT0.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY(), cluster->getT0());
      _trackResClusterSize.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY(), float(cluster->getNumHits()));
      _trackResCharge.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY(), cluster->getValue());
      _trackResCharge3D.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY(), cluster->getValue());
      for(unsigned mm=0; mm<cluster->getNumHits(); ++mm){
	
	_trackResTime.at(nsensor)->Fill(tx - sensor->getOffX(), ty - sensor->getOffY(), cluster->getHit(mm)->getTiming()); // timing is the first hit. the biggest
	if((tx - sensor->getOffX())>-100.0 && (tx - sensor->getOffX())<-70.0
	   && (ty - sensor->getOffY())>-70.0 && (ty - sensor->getOffY())<70.0){
	  _hitTimeVsCharge.at(nsensor)->Fill(cluster->getHit(mm)->getValue(), cluster->getTiming());
	}
	_hitT0.at(nsensor)->Fill(cluster->getHit(mm)->getT0());
      }
	//_hitCharge.at(nsensor)->Fill(cluster->getCharge());	
      // hit residual 1 D
      //_hitResidualCut.at(nsensor)->Fill(sqrt((tx - cluster->getPosX())*(tx - cluster->getPosX()) + (ty - cluster->getPosY())*(ty - cluster->getPosY())));
	_hitResidualxCut.at(nsensor)->Fill(tx - cluster->getPosX());
	_hitResidualyCut.at(nsensor)->Fill(ty - cluster->getPosY());	
      } // end cluster loop

      // Print missing hits
      /*
      if(!isInsideHit){
	if((tx - sensor->getOffX())>-100.0 && (tx - sensor->getOffX())<-70.0
	   && (ty - sensor->getOffY())>-20.0 && (ty - sensor->getOffY())<10.0){
	  std::cout << "inside sensor: " << dutEvent->getFrameNumber()
	  << " " << dutEvent->getTimeStamp()
	  << std::endl;
	}
      }
      */
      //if((tx - sensor->getOffX())>105.0 && (tx - sensor->getOffX())<115.0
      // && (ty - sensor->getOffY())>-15.0 && (ty - sensor->getOffY())<-5.0)
      //hhit->processEvent(dutEvent);
      
    }
  }// end track loop
  
  // if there is a track, then plot the occupancy
  if(pass_track_selection){
    
    // draw the clusters
    hcluster->processEvent(dutEvent);
    hhit->processEvent(dutEvent);
          
    for (unsigned int nplane = 0; nplane < dutEvent->getNumPlanes(); nplane++){
      Storage::Plane* plane = dutEvent->getPlane(nplane);	  
      for (unsigned int ncluster = 0; ncluster < plane->getNumClusters(); ncluster++){
	Storage::Cluster* cluster = plane->getCluster(ncluster);	  

        // Check if the cluster passes the cuts
        bool pass = true;
        for (unsigned int ncut = 0; ncut < _numClusterCuts; ncut++)
          if (!_clusterCuts.at(ncut)->check(cluster)) { pass = false; break; }
        if (!pass) continue;
	
	// DUT hit occupancy

	_dutHitOcc.at(nplane)->Fill(cluster->getPosX() - _dutDevice->getSensor(nplane)->getOffX(), cluster->getPosY() - _dutDevice->getSensor(nplane)->getOffY());
      }
    }
  } // end has a track
  //
  // end occupancy filling
  //
  
  // Check if the event passes the cuts
  for (unsigned int ncut = 0; ncut < _numEventCuts; ncut++)
    if (!_eventCuts.at(ncut)->check(refEvent)) return;

  for (unsigned int ntrack = 0; ntrack < refEvent->getNumTracks(); ntrack++)
  {
    Storage::Track* track = refEvent->getTrack(ntrack);

    // Check if the track passes the cuts
    bool pass = true;
    for (unsigned int ncut = 0; ncut < _numTrackCuts; ncut++)
      if (!_trackCuts.at(ncut)->check(track)) { pass = false; break; }
    if (!pass) continue;

    // Make a list of the planes with their matched cluster
    std::vector<Storage::Cluster*> matches;
    for (unsigned int nplane = 0; nplane < dutEvent->getNumPlanes(); nplane++)
      matches.push_back(0); // No matches

    // Get the matches from the track
    for (unsigned int nmatch = 0; nmatch < track->getNumMatchedClusters(); nmatch++)
    {
      Storage::Cluster* cluster = track->getMatchedCluster(nmatch);

      // Check if this cluster passes the cuts
      bool pass = true;
      for (unsigned int ncut = 0; ncut < _numClusterCuts; ncut++)
        if (!_clusterCuts.at(ncut)->check(cluster)) { pass = false; break; }
      if (!pass) continue;

      matches.at(cluster->getPlane()->getPlaneNum()) = cluster;
    }
    //std::cout << "matches: " << matches.size() << " relativeToSensor: " << _relativeToSensor << std::endl;
    assert(matches.size() == _dutDevice->getNumSensors() &&
           _relativeToSensor < (int)_dutDevice->getNumSensors() &&
           "Efficiency: matches has the wrong size");

    // If asking for relative efficiency to a plane, look for a match in that one first
    if (_relativeToSensor >= 0 && !matches.at(_relativeToSensor)) continue;

    for (unsigned int nsensor = 0; nsensor < _dutDevice->getNumSensors(); nsensor++)
    {
      Mechanics::Sensor* sensor = _dutDevice->getSensor(nsensor);

      double tx = 0, ty = 0, tz = 0;
      Processors::trackSensorIntercept(track, sensor, tx, ty, tz);

      double px = 0, py = 0;
      sensor->spaceToPixel(tx, ty, tz, px, py);

      const Storage::Cluster* match = matches.at(nsensor);

      //MAtevz 20141202 added 9-pixel efficiency plot.
      // Check if the intercepted pixel is hit
      bool trackMatchesPixel = false;
      unsigned int pixelX = px;
      unsigned int pixelY = py;
      double ninepixelX = 5.0;
      double ninepixelY = 5.0;

      if (match)
      {
        for (unsigned int nhit = 0; nhit < match->getNumHits(); nhit++)
        {
          const Storage::Hit* hit = match->getHit(nhit);

          if (abs( abs(hit->getPixX()) - abs(pixelX) ) <= 1 && abs( abs(hit->getPixY()) - abs(pixelY) ) <= 1)
          {
            if (abs(hit->getPixX()) - abs(pixelX) == 0)
              ninepixelX = px - (int)px;
            else if ( abs(hit->getPixX()) - abs(pixelX) == 1)
              ninepixelX = (px) - (int)px - 1.0; //)*(-1.0);
            else if (abs(hit->getPixX()) - abs(pixelX) == -1)
              ninepixelX = px - (int)px + 1.0;
            else cout<<"ERR: unexpected hit on X in efficiency.cpp line 128\n";

            if (abs(hit->getPixY()) - abs(pixelY) == 0)
              ninepixelY = py - (int)py;
            else if (abs(hit->getPixY()) - abs(pixelY) == 1)
              ninepixelY = py - (int)py - 1.0; //)*(-1.0);
            else if (abs(hit->getPixY()) - abs(pixelY) == -1)
              ninepixelY = py - (int)py + 1.0;
            else cout<<"ERR: unexpected hit on Y in efficiency.cpp line 133\n";

            _inPixelEfficiencyExtended.at(nsensor)->Fill(
                              (bool)match,
                              ninepixelX * sensor->getPitchX(),
                              ninepixelY * sensor->getPitchY());
            _inPixelCCE.at(nsensor)->Fill(
                              ninepixelX * sensor->getPitchX(),
                              ninepixelY * sensor->getPitchY(),
                              hit->getValue()  ); //ampl. value as weight
          }


          // make a spatial cut to get rid of the noisy amplitudes.
          // for now hardcode the cut to be smaller than 4x4 mm.
          double cutX = 0.8*(sensor->getPitchX()/2); //80% of the Xaxis
          double cutY = 0.8*(sensor->getPitchY()/2); //80% of the Yaxis

          if (  abs(ninepixelX) < cutX && abs(ninepixelY) < cutY )
            _amplDistCuts.at(nsensor)->Fill( hit->getValue() );


          if (hit->getPixX() == pixelX && hit->getPixY() == pixelY)
          {
           trackMatchesPixel = true;
           break;
          }
        }
      } // end of 9-pixel efficiency plotting.


      // Get the location within the pixel where the track interception
      const double trackPosX = px - (int)px;
      const double trackPosY = py - (int)py;

      _inPixelEfficiency.at(nsensor)->Fill((bool)match,
                                           trackPosX * sensor->getPitchX(),
                                           trackPosY * sensor->getPitchY());

      if (match)
      {
        _efficiencyMap.at(nsensor)->Fill(true, tx, ty);
        _matchedTracks.at(nsensor)->Fill(1);
        if (_efficiencyTime.size())
          _efficiencyTime.at(nsensor)->Fill(true, _refDevice->tsToTime(refEvent->getTimeStamp()));
      }
      else
      {
        _efficiencyMap.at(nsensor)->Fill(false, tx, ty);
        _matchedTracks.at(nsensor)->Fill(0);
        if (_efficiencyTime.size())
          _efficiencyTime.at(nsensor)->Fill(false, _refDevice->tsToTime(refEvent->getTimeStamp()));
      }
    }
  }
}

void Efficiency::postProcessing()
{
  if (_postProcessed) return;
  fout.close();
  for (unsigned int nsensor = 0; nsensor < _dutDevice->getNumSensors(); nsensor++)
  {
    //divide accumulated amplitude histogram with
    //number of hits per bin.
    _inPixelCCE.at(nsensor)->Divide(
            _inPixelEfficiencyExtended.at(nsensor)->GetTotalHistogram() );

    //superimpose the two amplitude histograms onto a new one.
    _amplDistCommon.at(nsensor) = (TH1D*)_amplDistCuts.at(nsensor)->DrawNormalized();
    _amplDistCommon.at(nsensor) = (TH1D*)_amplDist.at(nsensor)->DrawNormalized("same");

    // Get efficiency per pixel
    TEfficiency* efficiency = _efficiencyMap.at(nsensor);
    const TH1* values = efficiency->GetTotalHistogram();
    TH1D* distribution = _efficiencyDistribution.at(nsensor);

    // Loop over all pixel groups
    for (Int_t binx = 1; binx <= values->GetNbinsX(); binx++)
    {
      for (Int_t biny = 1; biny <= values->GetNbinsY(); biny++)
      {
        const Int_t bin = values->GetBin(binx, biny);
        if (values->GetBinContent(binx, biny) < 1) continue;
        const double value = efficiency->GetEfficiency(bin);
        const double sigmaLow = efficiency->GetEfficiencyErrorLow(bin);
        const double sigmaHigh = efficiency->GetEfficiencyErrorUp(bin);

        // Find the probability of this pixel group being found in all bins of the distribution
        double normalization = 0;
        for (Int_t distBin = 1; distBin <= distribution->GetNbinsX(); distBin++)
        {
          const double evaluate = distribution->GetBinCenter(distBin);
          const double sigma = (evaluate < value) ? sigmaLow : sigmaHigh;
          const double weight = TMath::Gaus(evaluate, value, sigma);
          normalization += weight;
        }
        for (Int_t distBin = 1; distBin <= distribution->GetNbinsX(); distBin++)
        {
          const double evaluate = distribution->GetBinCenter(distBin);
          const double sigma = (evaluate < value) ? sigmaLow : sigmaHigh;
          const double weight = TMath::Gaus(evaluate, value, sigma);
          distribution->Fill(evaluate, weight / normalization);
        }
      }
    }
  }

  // Compute the Efficiency
  assert(_dutHitOcc.size()==_trackOcc.size() && "Analyzer::Efficiency: track occupancy different in size from hit occ.");
  std::stringstream name; // Build name strings for each histo
  std::stringstream title; // Build title strings for each histo
  
  for(unsigned nplane = 0; nplane<_trackOcc.size(); ++nplane){
    TH2D *htmp = static_cast<TH2D *>(_dutHitOcc.at(nplane)->Clone());
    htmp->Divide(_trackOcc.at(nplane));

    name.str(""); title.str("");
    name << "sensor" << nplane << "_"
         << "hitEff" << _nameSuffix;
    title << "sensor" << nplane 
          << " Hit Efficiency";    
    htmp->SetName(name.str().c_str());
    htmp->SetTitle(title.str().c_str());
    htmp->SetDirectory(_dutHitOcc.at(nplane)->GetDirectory());
    _dutHitEff.push_back(htmp);

    // Track residual
    TH2D *htmp_res = static_cast<TH2D *>(_trackResHit.at(nplane)->Clone());
    htmp_res->Divide(_trackRes.at(nplane));

    name.str(""); title.str("");
    name << "sensor" << nplane << "_"
         << "TrackResEff" << _nameSuffix;
    title << "sensor" << nplane 
          << " Track Residual Efficiency";    
    htmp_res->SetName(name.str().c_str());
    htmp_res->SetTitle(title.str().c_str());
    htmp_res->SetDirectory(_trackRes.at(nplane)->GetDirectory());
    _trackResEff.push_back(htmp_res);
    //printEfficiency(_trackResHit.at(nplane), _trackRes.at(nplane), _dutDevice->getSensor(nplane),0.9);
    //printEfficiency(_trackResHit.at(nplane), _trackRes.at(nplane), _dutDevice->getSensor(nplane),0.6);

    // Track residual
    TH2D *htmp_res_fine = static_cast<TH2D *>(_trackResHitFine.at(nplane)->Clone());
    htmp_res_fine->Divide(_trackResFine.at(nplane));

    name.str(""); title.str("");
    name << "sensor" << nplane << "_"
         << "TrackResEffFine" << _nameSuffix;
    title << "sensor_fine" << nplane 
          << " Track Residual Efficiency";    
    htmp_res_fine->SetName(name.str().c_str());
    htmp_res_fine->SetTitle(title.str().c_str());
    htmp_res_fine->SetDirectory(_trackResFine.at(nplane)->GetDirectory());
    _trackResEffFine.push_back(htmp_res_fine);
    printEfficiency(_trackResHitFine.at(nplane), _trackResFine.at(nplane), _dutDevice->getSensor(nplane),0.90);
    printEfficiency(_trackResHitFine.at(nplane), _trackResFine.at(nplane), _dutDevice->getSensor(nplane),0.60); 
    
  }

  for(unsigned nplane = 0; nplane<_trackResCharge.size(); ++nplane){
    for(unsigned i=0; i<_trackResCharge.at(nplane)->GetNbinsX(); ++i ){
      for(unsigned j=0; j<_trackResCharge.at(nplane)->GetNbinsY(); ++j ){      
	float entries = float(_trackResHit.at(nplane)->GetBinContent(i,j));
	if (entries>0.0){
	  _trackResCharge.at(nplane)->SetBinContent(i,j, _trackResCharge.at(nplane)->GetBinContent(i,j)/entries);
	  _trackResTime.at(nplane)->SetBinContent(i,j, _trackResTime.at(nplane)->GetBinContent(i,j)/entries);
	  _trackResT0.at(nplane)->SetBinContent(i,j, _trackResT0.at(nplane)->GetBinContent(i,j)/entries);
	  _trackResClusterSize.at(nplane)->SetBinContent(i,j, _trackResClusterSize.at(nplane)->GetBinContent(i,j)/entries);	  	  	  
	}
      }
    }
  }
  
  _postProcessed = true;
}

  // print the efficiency for the size of the DUT
  void Efficiency::printEfficiency(const TH2D* hnum,
				   const TH2D *hden,
				   const Mechanics::Sensor *sensor,
				   const double size_of_dut){

    //const double size_of_dut = 0.8; // 80% of the chip size
    const double lowX =  -1.0* size_of_dut*sensor->getSensitiveX() / 2.0;
    const double uppX =   size_of_dut*sensor->getSensitiveX() / 2.0;
    const double lowY =  -1.0* size_of_dut*sensor->getSensitiveY() / 2.0;
    const double uppY =   size_of_dut*sensor->getSensitiveY() / 2.0;

    unsigned lowXbin=0, uppXbin=0, lowYbin=0, uppYbin=0;

    for(unsigned i=0; i<hnum->GetNbinsX()+1; ++i ){
      if(hnum->GetXaxis()->GetBinLowEdge(i)>=lowX) { lowXbin = i; break; }
    }
    for(unsigned i=lowXbin; i<hnum->GetNbinsX()+1; ++i ){
      if(hnum->GetXaxis()->GetBinUpEdge(i)>=uppX) { uppXbin = i; break; }
    }
    for(unsigned i=0; i<hnum->GetNbinsY()+1; ++i ){
      if(hnum->GetYaxis()->GetBinLowEdge(i)>=lowY) { lowYbin = i; break; }
    }
    for(unsigned i=lowYbin; i<hnum->GetNbinsY()+1; ++i ){
      if(hnum->GetYaxis()->GetBinUpEdge(i)>=uppY) { uppYbin = i; break; }
    }
    double stat_error =0.0;
    double den = hden->IntegralAndError(lowXbin,uppXbin,lowYbin,uppYbin, stat_error);
    double num = hnum->Integral(lowXbin,uppXbin,lowYbin,uppYbin);    
    if(den>0.0) num/=den;
    double error = sqrt(num*(1-num));
    double sqrtN = den;
    if(stat_error>0.0) sqrtN /=stat_error;
    if(sqrtN>0.0) error = error/sqrtN;

    // Print the efficiency results
    std::cout << "=========================================" << std::endl;
    std::cout << "Integral X: " << lowXbin << "-" << uppXbin << "  "
	      << hnum->GetXaxis()->GetBinLowEdge(lowXbin) << "-" << hnum->GetXaxis()->GetBinUpEdge(uppXbin) << "  "
	      << " Y: " << lowYbin << "-" << uppYbin << "  "
	      << hnum->GetYaxis()->GetBinLowEdge(lowYbin) << "-" << hnum->GetYaxis()->GetBinUpEdge(uppYbin)      
	      << " eff: " << num
	      << "+/-" << error << std::endl;
    std::cout << "=========================================" << std::endl;

    hcluster->postProcessing();
    hhit->postProcessing();
  }
  
Efficiency::Efficiency(const Mechanics::Device* refDevice,
                       const Mechanics::Device* dutDevice,
                       TDirectory* dir,
                       const char* suffix,
                       int relativeToSensor,
                       unsigned int rebinX,
                       unsigned int rebinY,
                       unsigned int pixBinsX,
                       unsigned int pixBinsY) :
  // Base class is initialized here and manages directory / device
  DualAnalyzer(refDevice, dutDevice, dir, suffix),
  // Initialize processing parameters here
  _relativeToSensor(relativeToSensor)
{
  assert(refDevice && dutDevice && "Analyzer: can't initialize with null device");

  // Makes or gets a directory called from inside _dir with this name
  TDirectory* plotDir = makeGetDirectory("Efficiency");
  fout.open("forEnrico.txt"); 
  // initialize the Tdirectory
  hcluster = new ClusterInfo(dutDevice,plotDir,"DUT",16,5);
  hhit = new HitInfo(dutDevice,plotDir,"DUT",16,5);
  
  std::stringstream name; // Build name strings for each histo
  std::stringstream title; // Build title strings for each histo

  if (relativeToSensor >= (int)_dutDevice->getNumSensors())
    throw "Efficiency: relative sensor exceeds range";

  // Generate a histogram for each sensor in the device
  for (unsigned int nsens = 0; nsens < _dutDevice->getNumSensors(); nsens++)
  {
    Mechanics::Sensor* sensor = _dutDevice->getSensor(nsens);

    unsigned int nx = sensor->getNumX() / rebinX;
    unsigned int ny = sensor->getNumY() / rebinY;
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    const double lowX = sensor->getOffX() - sensor->getSensitiveX() / 2.0;
    const double uppX = sensor->getOffX() + sensor->getSensitiveX() / 2.0;
    const double lowY = sensor->getOffY() - sensor->getSensitiveY() / 2.0;
    const double uppY = sensor->getOffY() + sensor->getSensitiveY() / 2.0;

    // TODO: change to pixel space

    cout << "   DEVICE NAME: "<< sensor->getDevice()->getName() << endl
	 << "   sensor NAME: " << sensor->getName() << endl
	 << "   map: "<< _nameSuffix << endl;
    // Efficiency map initialization
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "Map" << _nameSuffix;
    // Special titel includes axis
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Efficiency Map"
          << ";X position" << " [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position" << " [" << _dutDevice->getSpaceUnit() << "]"
          << ";Efficiency";
    TEfficiency* map = new TEfficiency(name.str().c_str(), title.str().c_str(),
                                       nx, lowX, uppX,
                                       ny, lowY, uppY);
    map->SetDirectory(plotDir);
    map->SetStatisticOption(TEfficiency::kFWilson);
    _efficiencyMap.push_back(map);


    // Pixel grouped efficieny initialization
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "GroupedDistribution" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Grouped Efficiency";
    TH1D* pixels = new TH1D(name.str().c_str(), title.str().c_str(),
                            20, 0, 1.001);
    pixels->GetXaxis()->SetTitle("Pixel group efficiency");
    pixels->GetYaxis()->SetTitle("Pixel groups / 0.05");
    pixels->SetDirectory(plotDir);
    _efficiencyDistribution.push_back(pixels);



    //Matevz 20141202 9-pixel in-pixel efficiency histogram initialization
    //Not inPixelEfficiency - inPixelEfficiencyExtended
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "In9PixelEfficiency" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " In 9 Pixel Eff"
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TEfficiency* inPixelEfficiencyExtended = new TEfficiency(name.str().c_str(), title.str().c_str(),
                             3*pixBinsX, -sensor->getPitchX(), 2*sensor->getPitchX(),
                             3*pixBinsY, -sensor->getPitchY(), 2*sensor->getPitchY());

    inPixelEfficiencyExtended->SetDirectory(plotDir);
    _inPixelEfficiencyExtended.push_back(inPixelEfficiencyExtended);



    //Matevz 20141203 charge collection map histogram initialization
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "inPixelCCE" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " in-pixel charge collection map"
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* inPixelCCE = new TH2D(name.str().c_str(), title.str().c_str(),
                             3*pixBinsX, -sensor->getPitchX(), 2*sensor->getPitchX(),
                             3*pixBinsY, -sensor->getPitchY(), 2*sensor->getPitchY());

    inPixelCCE->SetDirectory(plotDir);
    _inPixelCCE.push_back(inPixelCCE);



    //Matevz 20141203 amplitude histograms initialisation
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "PulseAmpl" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Pulse amplitude distribution";

    TH1D* amplDist = new TH1D(name.str().c_str(), title.str().c_str(), 100, 0, 0.1);

    amplDist->SetDirectory(plotDir);
    _amplDist.push_back(amplDist);

    //Matevz 20141203 amplitude histograms with CUTS initialisation
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "PulseAmplCuts" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Pulse amplitude distribution with cuts";

    TH1D* amplDistCuts = new TH1D(name.str().c_str(), title.str().c_str(), 100, 0, 0.1);

    amplDistCuts->SetDirectory(plotDir);
    _amplDistCuts.push_back(amplDistCuts);

    //Matevz 20141203 amplitude histograms for SAMEing of both - initialisation
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "PulseAmplCommon" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Pulse amplitude distribution with cuts";

    TH1D* amplDistCommon = new TH1D(name.str().c_str(), title.str().c_str(), 100, 0, 0.1);

    amplDistCommon->SetDirectory(plotDir);
    _amplDistCommon.push_back(amplDistCommon);



    // Track matching initialization
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         << "MatchedTracks" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Matched Tracks";
    TH1D* matched = new TH1D(name.str().c_str(), title.str().c_str(),
                             2, 0 - 0.5, 2 - 0.5);
    matched->GetXaxis()->SetTitle("0 : unmatched tracks | 1 : matched tracks");
    matched->GetZaxis()->SetTitle("Tracks");
    matched->SetDirectory(plotDir);
    _matchedTracks.push_back(matched);

    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "InPixelEfficiency" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " In Pixel Efficiency"
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TEfficiency* inPixelEfficiency = new TEfficiency(name.str().c_str(), title.str().c_str(),
                                                     pixBinsX, 0, sensor->getPitchX(),
                                                     pixBinsY, 0, sensor->getPitchY());
    inPixelEfficiency->SetDirectory(plotDir);
    _inPixelEfficiency.push_back(inPixelEfficiency);

    // Track occupancy extrapolated to DUT position
    float num_pixels = _dutDevice->getNumPixels()==1 ? 5.0 : float(_dutDevice->getNumPixels());
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackOccupancy" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Track Occupancy "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackOcc = new TH2D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      4*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackOcc->SetDirectory(plotDir);
    _trackOcc.push_back(trackOcc);


    // Average Charge 
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "DUTCharge" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Charge "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackResCharge = new TH2D(name.str().c_str(), title.str().c_str(),
				    20*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
				    20*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackResCharge->SetDirectory(plotDir);
    _trackResCharge.push_back(trackResCharge);


    // Average Charge 
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "DUTCharge3D" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Charge "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << "; Signal Size [V] "      
          << ";Tracks";
    TH3D* trackResCharge3D = new TH3D(name.str().c_str(), title.str().c_str(),
				      5*2*pixBinsX, -1.0*num_pixels*sensor->getPitchX(), 1.0*num_pixels*sensor->getPitchX(),
				      5*2*pixBinsY, -1.0*num_pixels*sensor->getPitchY(), 1.0*num_pixels*sensor->getPitchY(),
				      25, 0.0, 0.02);

    trackResCharge3D->SetDirectory(plotDir);
    _trackResCharge3D.push_back(trackResCharge3D);    


    // Average T0 
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "DUTT0" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Hit T0 "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackResT0 = new TH2D(name.str().c_str(), title.str().c_str(),
			      20*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      20*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackResT0->SetDirectory(plotDir);
    _trackResT0.push_back(trackResT0);

    // Average clustersize 
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "DUTClusterSize" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Cluster Size "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackResClusterSize = new TH2D(name.str().c_str(), title.str().c_str(),
			      20*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      20*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackResClusterSize->SetDirectory(plotDir);
    _trackResClusterSize.push_back(trackResClusterSize);     

    // Average time to collect charge
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "DUTTime" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Time "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackResTime = new TH2D(name.str().c_str(), title.str().c_str(),
			      20*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      20*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackResTime->SetDirectory(plotDir);
    _trackResTime.push_back(trackResTime);

    // Average time to collect charge
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TimeVsCharge" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " TimeVsCharge "
          << "; Charge " 
          << "; Time " 
          << ";Tracks";
    TH2D* hitTimeVsCharge = new TH2D(name.str().c_str(), title.str().c_str(),
				     600, -0.1, 0.2,
				     1000, 0.0, 1000.0);

    hitTimeVsCharge->SetDirectory(plotDir);
    _hitTimeVsCharge.push_back(hitTimeVsCharge);      

    // Track occupancy extrapolated to DUT position
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackResidual" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Track TrackResidual "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackRes = new TH2D(name.str().c_str(), title.str().c_str(),
			      20*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      20*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackRes->SetDirectory(plotDir);
    _trackRes.push_back(trackRes);    

    // Track occupancy extrapolated to DUT position
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackResidualHit" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Track TrackResidual Hit "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* trackResHit = new TH2D(name.str().c_str(), title.str().c_str(),
			      20*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      20*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());

    trackResHit->SetDirectory(plotDir);
    _trackResHit.push_back(trackResHit);


    ///----
    // Track occupancy extrapolated to DUT position
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackResidualFine" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Track TrackResidual Fine "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    /*    TH2D* trackResFine = new TH2D(name.str().c_str(), title.str().c_str(),
			      5*4*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      5*4*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());*/
   
    TH2D* trackResFine = new TH2D(name.str().c_str(), title.str().c_str(),
				  800, -400.0, 400.0,
				  800, -400.0, 400.0);
   
    trackResFine->SetDirectory(plotDir);
    _trackResFine.push_back(trackResFine);    

    // Track occupancy extrapolated to DUT position
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackResidualHitFine" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Track TrackResidual Hit Fine "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    /*TH2D* trackResHitFine = new TH2D(name.str().c_str(), title.str().c_str(),
			      5*4*pixBinsX, -20.0*num_pixels*sensor->getPitchX(), 20.0*num_pixels*sensor->getPitchX(),
			      5*4*pixBinsY, -20.0*num_pixels*sensor->getPitchY(), 20.0*num_pixels*sensor->getPitchY());*/
    TH2D* trackResHitFine = new TH2D(name.str().c_str(), title.str().c_str(),
				     800, -400.0, 400.0,
				     800, -400.0, 400.0);
   
    trackResHitFine->SetDirectory(plotDir);
    _trackResHitFine.push_back(trackResHitFine);  


    // DUT track cluster residual
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackClusterResidualX" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " x Track to hit residual"
          << ";x DUT Hit [" << _dutDevice->getSpaceUnit() << "]" 
          << ";Events";
    TH1D* dutTrackResX = new TH1D(name.str().c_str(), title.str().c_str(),
    			      4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );
    
    dutTrackResX->SetDirectory(plotDir);
    _dutTrackResX.push_back(dutTrackResX);


    // DUT track cluster residual
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackClusterResidualY" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " x Track to hit residual"
          << ";x DUT Hit [" << _dutDevice->getSpaceUnit() << "]" 
          << ";Events";
    TH1D* dutTrackResY = new TH1D(name.str().c_str(), title.str().c_str(),
    			      4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );
    
    dutTrackResY->SetDirectory(plotDir);
    _dutTrackResY.push_back(dutTrackResY);


    // Residual X vs Framenumber
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "ResidualX_vs_Framenumber" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Residual X vs Framenumber "
          << ";Framenumber" 
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* projXvsFN = new TH2D(name.str().c_str(), title.str().c_str(),
			       1000,0,400000,
			       2.0*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );
    
    projXvsFN->SetDirectory(plotDir);
    _projXvsFN.push_back(projXvsFN);


    // Residual Y vs Framenumber
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
      <<  "ResidualY_vs_Framenumber" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " Residual Y vs Framenumber "
          << ";Framenumber" 
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* projYvsFN = new TH2D(name.str().c_str(), title.str().c_str(),
			       1000,0,400000,
			       2*pixBinsY, -2.0*num_pixels*sensor->getPitchY(), 2.0*num_pixels*sensor->getPitchY());
    
    projYvsFN->SetDirectory(plotDir);
    _projYvsFN.push_back(projYvsFN);

 
    // Error Track Position
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackPositionErrorX" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " x Hit Track Position Error "
          << ";x DUT Hit [" << _dutDevice->getSpaceUnit() << "]" 
          << ";Events";
    TH1D* posErrorX = new TH1D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );

    posErrorX->SetDirectory(plotDir);
    _posErrorX.push_back(posErrorX);
 
    // Error Track Position
        name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "TrackPositionErrorY" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " x Hit Track Position Error "
          << ";x DUT Hit [" << _dutDevice->getSpaceUnit() << "]" 
          << ";Events";
    TH1D* posErrorY = new TH1D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );

    posErrorY->SetDirectory(plotDir);
    _posErrorY.push_back(posErrorY);


    // Hit Residual
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "HitResidualxNoCut" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " x Hit Track Residual No Cut "
          << ";x DUT Hit Residual [" << _dutDevice->getSpaceUnit() << "]" 
          << ";Events";
    TH1D* hitResxNoCut = new TH1D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );

    hitResxNoCut->SetDirectory(plotDir);
    _hitResidualxNoCut.push_back(hitResxNoCut);

    // Hit Residual
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "HitResidualxCut" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " x Hit TrackResidual Cut "
          << ";x DUT Hit Track Residual"
          << ";Events";
    TH1D* hitResxCut = new TH1D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX() );

    hitResxCut->SetDirectory(plotDir);
    _hitResidualxCut.push_back(hitResxCut);

    // Hit Residual
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "HitResidualyNoCut" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " y Hit Track Residual No Cut "
          << ";y DUT Hit Residual [" << _dutDevice->getSpaceUnit() << "]" 
          << ";Events";
    TH1D* hitResyNoCut = new TH1D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsY, -2.0*num_pixels*sensor->getPitchY(), 2.0*num_pixels*sensor->getPitchY() );

    hitResyNoCut->SetDirectory(plotDir);
    _hitResidualyNoCut.push_back(hitResyNoCut);

    // Hit Residual
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "HitResidualyCut" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " y Hit TrackResidual Cut "
          << ";y DUT Hit Track Residual"
          << ";Events";
    TH1D* hitResyCut = new TH1D(name.str().c_str(), title.str().c_str(),
			      4*pixBinsY, -2.0*num_pixels*sensor->getPitchY(), 2.0*num_pixels*sensor->getPitchY() );

    hitResyCut->SetDirectory(plotDir);
    _hitResidualyCut.push_back(hitResyCut);      

    // Hit T0
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "HitT0" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " T0 "
          << ";T_{0} [ns]"
          << ";Hits";
    TH1D* hitT0 = new TH1D(name.str().c_str(), title.str().c_str(),
			   1000,0,1000.0 );

    hitT0->SetDirectory(plotDir);
    _hitT0.push_back(hitT0);   

    // DUT hit occupancy
    name.str(""); title.str("");
    name << sensor->getDevice()->getName() << sensor->getName()
         <<  "DUTHitOccupancy" << _nameSuffix;
    title << sensor->getDevice()->getName() << " " << sensor->getName()
          << " DUT Hit Occupancy "
          << ";X position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Y position [" << _dutDevice->getSpaceUnit() << "]"
          << ";Tracks";
    TH2D* dutHitOcc = new TH2D(name.str().c_str(), title.str().c_str(),
			       4*pixBinsX, -2.0*num_pixels*sensor->getPitchX(), 2.0*num_pixels*sensor->getPitchX(),
			       4*pixBinsY, -2.0*num_pixels*sensor->getPitchY(), 2.0*num_pixels*sensor->getPitchY());

    dutHitOcc->SetDirectory(plotDir);
    _dutHitOcc.push_back(dutHitOcc);    
    

    if (_refDevice->getTimeEnd() > _refDevice->getTimeStart()) // If not used, they are both == 0
    {
      // Prevent aliasing
      const unsigned int nTimeBins = 100;
      const ULong64_t timeSpan = _refDevice->getTimeEnd() - _refDevice->getTimeStart() + 1;
      const ULong64_t startTime = _refDevice->getTimeStart();
      const ULong64_t endTime = timeSpan - (timeSpan % nTimeBins) + startTime;

      name.str(""); title.str("");
      name << sensor->getDevice()->getName() << sensor->getName()
           << "EfficiencyVsTime" << _nameSuffix;
      title << sensor->getDevice()->getName() << " " << sensor->getName()
            << " Efficiency Vs. Time"
            << ";Time [" << _refDevice->getTimeUnit() << "]"
            << ";Average efficiency";
      TEfficiency* efficiencyTime = new TEfficiency(name.str().c_str(), title.str().c_str(),
                                                    nTimeBins,
                                                    _refDevice->tsToTime(startTime),
                                                    _refDevice->tsToTime(endTime + 1));
      efficiencyTime->SetDirectory(plotDir);
      efficiencyTime->SetStatisticOption(TEfficiency::kFWilson);
      _efficiencyTime.push_back(efficiencyTime);
    }
  }
}

}
