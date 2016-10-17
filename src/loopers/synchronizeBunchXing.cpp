#include "synchronizeBunchXing.h"

#include <cassert>
#include <vector>
#include <iostream>
#include <math.h>
#include <string>

#include <Rtypes.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "../storage/storageio.h"
#include "../storage/event.h"
#include "../mechanics/device.h"
#include "../mechanics/alignment.h"
#include "../analyzers/singleanalyzer.h"
#include "../analyzers/dualanalyzer.h"
#include "../processors/processors.h"

#ifndef VERBOSE
#define VERBOSE 1
#endif

using std::cout;
using std::endl;
using std::flush;
using std::cin;

namespace Loopers {

unsigned int SynchronizeBunchXing::syncRatioLoop(ULong64_t start, ULong64_t num,
                                        unsigned int refOffset, unsigned int dutOffset)
{
  //Storage::Event* refStart = _refStorage->readEvent(start + refOffset);
  //Storage::Event* refEnd = _refStorage->readEvent(start + num + refOffset);
  //Storage::Event* dutStart = _dutStorage->readEvent(start + dutOffset);
  //Storage::Event* dutEnd = _dutStorage->readEvent(start + num + dutOffset);
  return 0.;
}

void SynchronizeBunchXing::loop()
{
  if (_bufferSize - _preDiscards <= 2)
    throw "SynchronizeBunchXing: buffer size and pre-discards make it impossible to sync";

  unsigned int refShift = 0;
  unsigned int dutShift = 0;

  // Some statistics for reporting
  ULong64_t singleRefOffsets = 0;
  ULong64_t multipleRefOffsets = 0;
  ULong64_t singleDutOffsets = 0;
  ULong64_t multipleDutOffsets = 0;
  ULong64_t questionableEvents = 0;
  ULong64_t failedSyncs = 0;
  ULong64_t goodSyncs = 0;
  ULong64_t numConsecutiveFails = 0;
  ULong64_t invalidEvents = 0;
  ULong64_t invalidEventsDUT = 0;
  ULong64_t invalidEventsRef = 0;  
  ULong64_t totalReadEvents = 0;
  ULong64_t totalWrittenEvents = 0;
  ULong64_t numConsecutiveSyncs = 0;
  ULong64_t consecutiveFails = 0;
  ULong64_t consecutiveSyncs = 0;
  // next we have a starting point that is aligned
  for (unsigned int nevent = _startEvent; nevent <= (_endEvent-1); ++nevent)
    {
      ++totalReadEvents;
      
      // read in events
      Storage::Event* refEvent = _refStorage->readEvent(nevent+_eventSkip);
      Storage::Event* dutEvent = _dutStorage->readEvent(dutShift+nevent+_eventSkip);
      Storage::Event* refEventPrev = _refStorage->readEvent(nevent+_eventSkip-1);
      Storage::Event* dutEventPrev = _dutStorage->readEvent(dutShift+nevent+_eventSkip-1);	
      
      if(dutEvent->getInvalid()){
	//std::cout << "Invalid Event DUT - skipping" << std::endl;
	++invalidEventsDUT;
	continue;
      }
      
      if( refEvent->getInvalid()){
	//std::cout << "Invalid Event Reference - skipping " << std::endl;
	++invalidEventsRef;
	continue;
      }
      
      // 40MHz trigger
      double difft = (refEvent->getTimeStamp()-refEventPrev->getTimeStamp())/40.0e6;
      double diffd = (dutEvent->getTimeStamp()-dutEventPrev->getTimeStamp())/1.0e5;
      if(difft>1.0){// telescope has a new bunch
	if(diffd>1.0){
	  // all good 
	  //DUTLastBunchXing=nevent
	}
	else if( nevent!=0){
	  //print 'NOT a DUT first Bunch Xing ',diffd
	  // search backward for start
	  for(unsigned j=0; j<_maxOffset; ++j){
	    if (j>nevent-1){
	      ///print 'NOT DONE because search passes beginning of file'
	      break;
	    }
	    // do backward search
	    ++dutShift;
	    Storage::Event* dutSearch = _dutStorage->readEvent(dutShift+nevent+_eventSkip);
	    Storage::Event* dutSearchPrev = _dutStorage->readEvent(dutShift+nevent+_eventSkip-1);
	    
	    double diffd_TMP = (dutSearch->getTimeStamp()-dutSearchPrev->getTimeStamp())/1.0e5;
	    // clean
	    delete dutSearchPrev;
	    if(diffd_TMP>1.0){
	      std::cout << "diffd_TMP: " << diffd_TMP << std::endl;
	      //diffd = (edut_e[i+nDUT][0]-edut_e[nDUT+i-1][0])/1.0e5;
	      delete dutEvent;
	      dutEvent=dutSearch;	      
	      ++goodSyncs;
	      break;
	    }// finish sync
	    else{ delete dutSearch; }
	    if(j==_maxOffset-1) ++failedSyncs;
	  }
	}//end fix sync
      }
      _refOutput->writeEvent(refEvent);
      _dutOutput->writeEvent(dutEvent);
      totalWrittenEvents++;
      
      //progressBar(nevent);
      
      // clean up
      delete refEvent;
      delete dutEvent;
      delete refEventPrev;
      delete dutEventPrev;
    } // end block
  
  
  if (VERBOSE)
    {
    cout << "\nSYNCHRONIZE STATISTICS:\n";
    cout << "  dutShift             : " << dutShift << "\n";
    cout << "  REF offsets of 1     : " << singleRefOffsets << "\n";    
    cout << "  REF offsets > 1      : " << multipleRefOffsets << "\n";
    cout << "  DUT offsets of 1     : " << singleDutOffsets << "\n";
    cout << "  DUT offsets > 1      : " << multipleDutOffsets << "\n";
    cout << "  Questionable events  : " << questionableEvents << "\n";
    cout << "  Good syncs           : " << goodSyncs << "\n";
    cout << "  Consecutive syncs    : " << numConsecutiveSyncs << "\n";
    cout << "  Failed syncs         : " << failedSyncs << "\n";
    cout << "  Consecutive fails    : " << numConsecutiveFails << "\n";
    cout << "  Invalid events       : " << invalidEvents << "\n";
    cout << "  Invalid events DUT   : " << invalidEventsDUT << "\n";
    cout << "  Invalid events Ref   : " << invalidEventsRef << "\n";        
    cout << "  Total read events    : " << totalReadEvents << "\n";
    cout << "  Total written events : " << totalWrittenEvents << "\n";
    cout << flush;
  }
}

void SynchronizeBunchXing::setThreshold(double threshold) { _threshold = threshold; }
void SynchronizeBunchXing::setSyncSample(unsigned int value) { _syncRatioSample = value; }
void SynchronizeBunchXing::setMaxOffset(unsigned int value) { _maxOffset = value; }
void SynchronizeBunchXing::setMaxLargeSyncAttempts(unsigned int value) { _maxLargeSyncAttempts = value; }
void SynchronizeBunchXing::setBufferSize(unsigned int value) { _bufferSize = value; }
void SynchronizeBunchXing::setPreDiscards(unsigned int value) { _preDiscards = value; }
void SynchronizeBunchXing::setMaxConsecutiveFails(unsigned int value) { _maxConsecutiveFails = value; }
void SynchronizeBunchXing::setDisplayDistributions(bool value) { _displayDistributions = value; }
  
SynchronizeBunchXing::SynchronizeBunchXing(Mechanics::Device* refDevice,
                         Mechanics::Device* dutDevice,
                         Storage::StorageIO* refOutput,
                         Storage::StorageIO* dutOutput,
                         Storage::StorageIO* refInput,
                         Storage::StorageIO* dutInput,
                         ULong64_t startEvent,
                         ULong64_t numEvents,
                         Long64_t eventSkip) :
  Looper(refInput, dutInput, startEvent, numEvents, eventSkip),
  _refDevice(refDevice),
  _dutDevice(dutDevice),
  _refOutput(refOutput),
  _dutOutput(dutOutput),
  _maxOffset(100),
  _maxLargeSyncAttempts(20),
  _bufferSize(10),
  _preDiscards(2),
  _displayDistributions(false),
  _maxConsecutiveFails(3)
{
  assert(refInput && dutInput && refOutput && dutOutput &&
         "Looper: initialized with null object(s)");
  assert(refInput->getNumPlanes() == refDevice->getNumSensors() &&
         dutInput->getNumPlanes() == dutDevice->getNumSensors() &&
         "Loopers: number of planes / sensors mis-match");  
}

}
