#ifndef HIT_H
#define HIT_H

namespace Storage {

class Cluster;
class Plane;

class Hit
{
private:
  unsigned int _pixX; // X pixel of this hit
  unsigned int _pixY;
  double _posX; // X position of the hit on the sensor
  double _posY;
  double _posZ;
  double _value; // Time over threshold, typically
  int _valueInt; // Time over threshold, typically
  double _timing; // Level 1 accept, typically
  double _t0; // start of hit detection
  int _isHit; // is a hit
  int _isValidFit; // is a valid fit
  double _lowFreqFFT; // FFT magnitude
  double _lowFreqFFTPhase; // FFT phase
  double _Chi2; // is a valid fit
  
  Cluster* _cluster; // The cluster containing this hit

protected:
  Plane* _plane; // Plane in which the hit is found
  Hit(); // Hits memory is managed by the event class
  ~Hit() { ; } // The user can get Hit pointers but can't delete them

public:
  void print();

  void setCluster(Cluster* cluster);

  // These are in the cpp file so that the classes can be included
  Cluster* getCluster() const;
  Plane* getPlane() const;

  // Inline setters and getters since they will be used frequently
  inline void setPix(unsigned int x, unsigned int y) { _pixX = x; _pixY = y; }
  inline void setPos(double x, double y, double z) { _posX = x; _posY = y; _posZ = z; }
  inline void setValue(double value) { _value = value; }
  inline void setValueInt(int value) { _valueInt = value; }
  inline void setTiming(double timing) { _timing = timing; }
  inline void setT0(double t0) { _t0 = t0; }

  inline void setIsHit(int isHit) { _isHit = isHit; }
  inline void setValidFit(int isValidFit) { _isValidFit = isValidFit; }
  inline void setLowFreqFFT(double LowFreqFFT) { _lowFreqFFT = LowFreqFFT; }
  inline void setLowFreqFFTPhase(double LowFreqFFTPhase) { _lowFreqFFTPhase = LowFreqFFTPhase; }    
  inline void setChi2(double chi2) { _Chi2 = chi2; }      

  inline unsigned int getPixX() const { return _pixX; }
  inline unsigned int getPixY() const { return _pixY; }
  inline double getPosX() const { return _posX; }
  inline double getPosY() const { return _posY; }
  inline double getPosZ() const { return _posZ; }
  inline double getValue() const { return _value; }
  inline int getValueInt() const { return _valueInt; }
  inline double getTiming() const { return _timing; }
  inline double getT0() const { return _t0; }

  inline int getIsHit() const { return _isHit; }
  inline int getValidFit() const { return _isValidFit; }
  inline double getLowFreqFFT() const { return _lowFreqFFT; }
  inline double getLowFreqFFTPhase() const { return _lowFreqFFTPhase; }    
  inline double getChi2() const { return _Chi2; }

  friend class Plane;     // Access set plane method
  friend class Event;     // Access cluster index
  friend class StorageIO; // Constructor and destructor
};

}

#endif // HIT_H
