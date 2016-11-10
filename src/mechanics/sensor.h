#ifndef SENSOR_H
#define SENSOR_H

#include <vector>
#include <string>
#include "TF1.h"

namespace Mechanics {

class Device;

class Sensor
{
private:
  const unsigned int _numX;
  const unsigned int _numY;
  const double _pitchX;
  const double _pitchY;
  const double _depth;
  const Device* _device;
  long unsigned _frameNumber;
  double _timeStamp;    
  std::string _name;
  const double _xox0;
  double _offX;
  double _offY;
  double _offZ;
  double _rotX;
  double _rotY;
  double _rotZ;
  const double _sensitiveX;
  const double _sensitiveY;
  unsigned int _numNoisyPixels;
  bool** _noisyPixels;

  double _rotation[3][3]; // The rotation matrix for the plane
  double _unRotate[3][3]; // Invert the rotation
  double _normalX;
  double _normalY;
  double _normalZ;

  void calculateRotation(); // Calculates the rotation matricies and normal vector
  void applyRotation(double& x, double& y, double& z, bool invert = false) const;

public:
  Sensor(unsigned int numX, unsigned int numY, double pitchX, double pitchY,
         double depth, Device* device, std::string name, double xox0 = 0,
         double offX = 0, double offY = 0, double offZ = 0,
         double rotX = 0, double rotY = 0, double rotZ = 0);
  ~Sensor();

  void print();

  void addNoisyPixel(unsigned int x, unsigned int y);
  void clearNoisyPixels();
  inline bool isPixelNoisy(unsigned int x, unsigned int y) const
      { return _noisyPixels[x][y]; }

  void setOffX(double offset);
  void setOffY(double offset);
  void setOffZ(double offset);
  void setRotX(double rotation);
  void setRotY(double rotation);
  void setRotZ(double rotation);
  void setFrameNumberFuncX(std::string a) {
    m_frameFuncX = a;
    f_frameFuncX = new TF1("xpos",m_frameFuncX.c_str());
  }
  void setFrameNumberFuncY(std::string a) {
    m_frameFuncY = a;
    f_frameFuncY = new TF1("ypos",m_frameFuncY.c_str());
  }
  void setTimeStampFuncX(std::string a) {
    m_timeStampFuncX = a;
    f_timeStampFuncX = new TF1("xtpos",m_timeStampFuncX.c_str());
  }
  void setTimeStampFuncY(std::string a) {
    m_timeStampFuncY = a;
    f_timeStampFuncY = new TF1("ytpos",m_timeStampFuncY.c_str());
  }

  void setFrameNumber(long unsigned frameNumber){_frameNumber = frameNumber;} 
  void setTimeStamp(double timeStamp){_timeStamp = timeStamp/1.0e8;} 

  void rotateToGlobal(double& x, double& y, double& z) const;
  void rotateToSensor(double& x, double& y, double& z) const;
  void pixelToSpace(double pixX, double pixY,
                    double& x, double& y, double& z) const;
  void spaceToPixel(double x, double y, double z,
                    double& pixX, double& pixY) const;

  void getGlobalOrigin(double& x, double& y, double& z) const;
  void getNormalVector(double& x, double& y, double& z) const;

  bool** getNoiseMask() const;
  unsigned int getNumX() const;
  unsigned int getNumY() const;
  unsigned int getPosNumX() const;
  unsigned int getPosNumY() const;
  double getPitchX() const;
  double getPitchY() const;
  double getPosPitchX() const;
  double getPosPitchY() const;
  double getDepth() const;
  double getXox0() const;
  double getOffX() const;
  double getOffY() const;
  double getOffZ() const;
  double getRotX() const;
  double getRotY() const;
  double getRotZ() const;
  double getSensitiveX() const;
  double getSensitiveY() const;
  double getPosSensitiveX() const;
  double getPosSensitiveY() const;

  std::string getFrameNumberFuncX() const { return m_frameFuncX; }
  std::string getFrameNumberFuncY() const { return m_frameFuncY; }
  std::string timeStampFuncX() const { return m_timeStampFuncX; }
  std::string timeStampFuncY() const { return m_timeStampFuncY; }
  
  const Device* getDevice() const;
  //const char* getName() const;
  const std::string getName() const;  
  std::string m_frameFuncX;
  std::string m_frameFuncY;
  std::string m_timeStampFuncX;
  std::string m_timeStampFuncY;
  TF1* f_frameFuncX;
  TF1* f_frameFuncY;
  TF1* f_timeStampFuncX;
  TF1* f_timeStampFuncY;
  static bool sort(const Sensor* s1, const Sensor* s2);
};

}

#endif // SENSOR_H
