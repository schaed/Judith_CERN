#include <iostream>
#include <stdexcept>

#include <TSystem.h>

#include "storage/storageo.h"
#include "storage/storagei.h"
#include "storage/storageio.h"
#include "storage/event.h"
#include "storage/track.h"
#include "storage/plane.h"
#include "storage/cluster.h"
#include "storage/hit.h"

int test_storageio() {
  Storage::StorageIO store("tmp.root", Storage::StorageIO::OUTPUT);

  if (store.getNumEvents() != 0 ||
      store.getNumPlanes() != 1 ||
      store.getFileMode() != Storage::StorageIO::OUTPUT ||
      store.getMaskMode() != Storage::StorageIO::REMOVE) {
    std::cerr << "Storage::StorageIO: default values not as expected" << std::endl;
    return -1;
  }

  Storage::Event& event = store.newEvent();

  if (&store.newEvent() != &event) {
    std::cerr << "Storage::StorageIO: no event caching" << std::endl;
    return -1;
  }

  // Make some objects which will be cahed when a new event is requested

  Storage::Hit& hit = event.newHit(0);
  hit.setPix(1, 2);
  hit.setPos(.1, .2, .3);
  hit.setValue(.1);
  hit.setTiming(.2);
  hit.setMasked(true);

  Storage::Cluster& cluster = event.newCluster(0);
  cluster.setPix(.1, .2);
  cluster.setPixErr(.1, .2);
  cluster.setPos(.1, .2, .3);
  cluster.setPosErr(.1, .2, .3);
  cluster.setTiming(.1);
  cluster.setValue(.2);

  Storage::Track& track = event.newTrack();
  track.setOrigin(.1, .2);
  track.setOriginErr(.1, .2);
  track.setSlope(.1, .2);
  track.setSlopeErr(.1, .2);
  track.setCovariance(.1, .2);
  track.setChi2(.1);

  cluster.addHit(hit);
  track.addCluster(cluster);

  // Make a new event and a new hit from the new event. Note that the actual
  // object is the same so we can keep using Event
  store.newEvent();

  // New hit should recycle the old one
  if (&event.newHit(0) != &hit) {
    std::cerr << "Storage::StorageIO: hits not caching" << std::endl;
    return -1;
  }

  // It shoudl be zeroed
  if (hit.getPixX() != 0 ||
      hit.getPixY() != 0 ||
      hit.getPosX() != 0 ||
      hit.getPosY() != 0 ||
      hit.getPosZ() != 0 ||
      hit.getValue() != 0 ||
      hit.getTiming() != 0 ||
      hit.getMasked() != false ||
      hit.fetchCluster() != 0) {
    std::cerr << "Storage::StorageIO: cached hit not zeroed" << std::endl;
    return -1;
  }

  if (&event.newCluster(0) != &cluster) {
    std::cerr << "Storage::StorageIO: clusters not caching" << std::endl;
    return -1;
  }

  if (cluster.getNumHits() != 0 ||
      cluster.getMatchDistance() != 0 ||
      cluster.getPixX() != 0 ||
      cluster.getPixY() != 0 ||
      cluster.getPixErrX() != 0 ||
      cluster.getPixErrY() != 0 ||
      cluster.getPosX() != 0 ||
      cluster.getPosY() != 0 ||
      cluster.getPosZ() != 0 ||
      cluster.getPosErrX() != 0 ||
      cluster.getPosErrY() != 0 ||
      cluster.getPosErrZ() != 0 ||
      cluster.getTiming() != 0 ||
      cluster.getValue() != 0 ||
      cluster.getIndex() != 0 ||
      cluster.fetchTrack() != 0 ||
      cluster.fetchMatchedTrack() != 0) {
    std::cerr << "Storage::StorageIO: cached cluster not zeroed" << std::endl;
    return -1;
  }

  if (&event.newTrack() != &track) {
    std::cerr << "Storage::StorageIO: tracks not caching" << std::endl;
    return -1;
  }

  if (track.getNumClusters() != 0 ||
      track.getNumMatchedClusters() != 0 ||
      track.getOriginX() != 0 ||
      track.getOriginY() != 0 ||
      track.getOriginErrX() != 0 ||
      track.getOriginErrY() != 0 ||
      track.getSlopeX() != 0 ||
      track.getSlopeY() != 0 ||
      track.getSlopeErrX() != 0 ||
      track.getSlopeErrY() != 0 ||
      track.getCovarianceX() != 0 ||
      track.getCovarianceY() != 0 ||
      track.getChi2() != 0 ||
      track.getIndex() != 0) {
    std::cerr << "Storage::StorageIO: cahced track not zeroed" << std::endl;
    return -1;
  }

  // Add another hit
  Storage::Hit& hit2 = event.newHit(0);
  // Ensure it is in fact a new hit
  if (&hit2 == &hit) {
    std::cerr << "Storage::StorageIO: cache re-using hits" << std::endl;
    return -1;
  }

  Storage::Cluster& cluster2 = event.newCluster(0);
  if (&cluster2 == &cluster) {
    std::cerr << "Storage::StorageIO: cache re-using clusters" << std::endl;
    return -1;
  }

  Storage::Track& track2 = event.newTrack();
  if (&track2 == &track) {
    std::cerr << "Storage::StorageIO: cache re-using tracks" << std::endl;
    return -1;
  }

  // Load up a new event and ensure both objects were cached. Once more, the
  // event object is the same so our reference is still valid
  store.newEvent();

  // Each call to newHit() should retrieve from the cahce until a the cache
  // is empty
  if (&event.newHit(0) != &hit2 || &event.newHit(0) != &hit) {
    std::cerr << "Storage::StorageIO: hits cache not working back" << std::endl;
    return -1;
  }

  if (&event.newCluster(0) != &cluster2 || &event.newCluster(0) != &cluster) {
    std::cerr << "Storage::StorageIO: clusters cache not working back" << std::endl;
    return -1;
  }

  if (&event.newTrack() != &track2 || &event.newTrack() != &track) {
    std::cerr << "Storage::StorageIO: tracks cache not working back" << std::endl;
    return -1;
  }

  gSystem->Exec("rm -f tmp.root");
  return 0;
}

int test_storageioReadWrite() {
  const size_t nplanes = 1;
  const Int_t nevents = 2;
  {  // Scope so that the output is closed
    Storage::StorageO store("tmp.root", nplanes);

    // Make some number of events
    for (size_t n = 0; n < nevents; n++) {
      Storage::Event& event = store.newEvent();

      // Alternate the plane populated at each event
      Storage::Hit& hit = event.newHit(n%nplanes);
      // Scale the values by the event
      hit.setPix(1*n, 2*n);
      hit.setPos(.1*n, .2*n, .3*n);
      hit.setValue(.1*n);
      hit.setTiming(.2*n);
      hit.setMasked(true);

      Storage::Cluster& cluster = event.newCluster(n%nplanes);
      cluster.setPix(.1*n, .2*n);
      cluster.setPixErr(.1*n, .2*n);
      cluster.setPos(.1*n, .2*n, .3*n);
      cluster.setPosErr(.1*n, .2*n, .3*n);
      cluster.setTiming(.1*n);
      cluster.setValue(.2*n);

      Storage::Track& track = event.newTrack();
      track.setOrigin(.1*n, .2*n);
      track.setOriginErr(.1*n, .2*n);
      track.setSlope(.1*n, .2*n);
      track.setSlopeErr(.1*n, .2*n);
      track.setCovariance(.1*n, .2*n);
      track.setChi2(.1*n);

      cluster.addHit(hit);
      track.addCluster(cluster);

      store.writeEvent(event);
    }
  }  // output scope

  {  // Scope the read-back separately
    Storage::StorageI store("tmp.root");

    if (nevents != store.getNumEvents()) {
      std::cerr << "Storage::StorageI: incorrect number of events" << std::endl;
      return -1;
    }

    if (nplanes != store.getNumPlanes()) {
      std::cerr << "Storage::StorageI: incorrect number of planes" << std::endl;
      return -1;
    }

    // Make some number of events
    for (Int_t n = 0; n < store.getNumEvents(); n++) {
      Storage::Event& event = store.readEvent(n);

      Storage::Track& track = event.getTrack(0);
      if (track.getOriginX() != .1*n ||
          track.getOriginY() != .2*n ||
          track.getOriginErrX() != .1*n ||
          track.getOriginErrY() != .2*n ||
          track.getSlopeX() != .1*n ||
          track.getSlopeY() != .2*n ||
          track.getSlopeErrX() != .1*n ||
          track.getSlopeErrY() != .2*n ||
          track.getCovarianceX() != .1*n ||
          track.getCovarianceY() != .2*n ||
          track.getChi2() != .2*n) {
        std::cerr << "Storage::StorageI: track read back incorrect" << std::endl;
        return -1;
      }

      //// Alternate the plane populated at each event
      //Storage::Hit& hit = event.newHit(n%nplanes);
      //// Scale the values by the event
      //hit.setPix(1*n, 2*n);
      //hit.setPos(.1*n, .2*n, .3*n);
      //hit.setValue(.1*n);
      //hit.setTiming(.2*n);
      //hit.setMasked(true);

      //Storage::Cluster& cluster = event.newCluster(n%nplanes);
      //cluster.setPix(.1*n, .2*n);
      //cluster.setPixErr(.1*n, .2*n);
      //cluster.setPos(.1*n, .2*n, .3*n);
      //cluster.setPosErr(.1*n, .2*n, .3*n);
      //cluster.setTiming(.1*n);
      //cluster.setValue(.2*n);
    }

  }  // input scope

  gSystem->Exec("rm -f tmp.root");
  return 0;
}

int main() {
  int retval = 0;

  try {
    if ((retval = test_storageio()) != 0) return retval;
    if ((retval = test_storageioReadWrite()) != 0) return retval;
  }
  
  catch (std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return -1;
  }

  return 0;
}
