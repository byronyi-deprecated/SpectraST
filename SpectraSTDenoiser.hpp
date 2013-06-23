#ifndef SPECTRASTDENOISER_HPP_
#define SPECTRASTDENOISER_HPP_

#include "SpectraSTPeakList.hpp"
#include <vector>

using namespace std;

class SpectraSTDenoiser {

public:

  SpectraSTDenoiser();

  virtual ~SpectraSTDenoiser();

  void addTrainingSpectrum(SpectraSTPeakList* pl);
  void generateBayesianModel();
  double predictPrior(SpectraSTPeakList* peakList, bool useOverallPrior = false);
  
  void filter(SpectraSTPeakList* peakList, unsigned int maxNumPeaks, double minSignalProb = -1.0);
  
  void useDefault();
  
  void printModel();
  
  bool isFilterReady() { return (m_ready); }
  
private:

  bool m_ready;
  
  unsigned int m_numSignal;
  unsigned int m_numNoise;
  
  vector<unsigned int> m_signalWithCNI;
  vector<unsigned int> m_noiseWithCNI;
   
  vector<unsigned int> m_signalWithMzPos;
  vector<unsigned int> m_noiseWithMzPos;
  
  unsigned int m_signalWithComplement;
  unsigned int m_noiseWithComplement;
  
  vector<double> m_sisters;
  vector<unsigned int> m_signalWithSister;
  vector<unsigned int> m_noiseWithSister;

  vector<double> m_logOddsWithCNI;
  vector<double> m_logOddsWithMzPos;
  double m_logOddsWithComplement;
  double m_logOddsWithoutComplement;
  vector<double> m_logOddsWithSister;
  vector<double> m_logOddsWithoutSister;
  
  
  static bool sortPeaksByOddsDesc(pair<double, Peak> a, pair<double, Peak> b);
  
};

#endif
