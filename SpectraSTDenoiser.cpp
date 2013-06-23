#include "SpectraSTDenoiser.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>

using namespace std;

SpectraSTDenoiser::SpectraSTDenoiser() :
  m_signalWithCNI(101, 1),
  m_noiseWithCNI(101, 1),
  m_signalWithMzPos(101, 1),
  m_noiseWithMzPos(101, 1),
  m_logOddsWithCNI(101, 0.0),
  m_logOddsWithMzPos(101, 0.0),
  m_ready(false) {

  m_sisters.push_back(-2.016);
  m_sisters.push_back(-1.008);
  m_sisters.push_back(1.008);
  m_sisters.push_back(2.016);
  m_sisters.push_back(9.005285);
  m_sisters.push_back(17.0027);
  m_sisters.push_back(18.01057);
  m_sisters.push_back(27.9949);
  m_sisters.push_back(57.02146);
  m_sisters.push_back(71.03711);
  m_sisters.push_back(87.03202);
  m_sisters.push_back(97.05276);
  m_sisters.push_back(99.06841);
  m_sisters.push_back(101.0476);
  m_sisters.push_back(103.0091);
  m_sisters.push_back(113.084);
  m_sisters.push_back(114.0429);
  m_sisters.push_back(115.0269);
  m_sisters.push_back(128.0585);
  m_sisters.push_back(128.0949);
  m_sisters.push_back(129.0425);
  m_sisters.push_back(131.0404);
  m_sisters.push_back(137.0589);
  m_sisters.push_back(147.0684);
  m_sisters.push_back(156.1011);
  m_sisters.push_back(163.0633);
  m_sisters.push_back(186.0793);
    
  m_signalWithSister.assign(m_sisters.size(), 1);
  m_noiseWithSister.assign(m_sisters.size(), 1);
  m_logOddsWithSister.assign(m_sisters.size(), 0);
  m_logOddsWithoutSister.assign(m_sisters.size(), 0);

}

SpectraSTDenoiser::~SpectraSTDenoiser() {

  
}

void SpectraSTDenoiser::useDefault() {
  
  m_logOddsWithCNI[0]=-4.60273;
  m_logOddsWithCNI[1]=-4.11211;
  m_logOddsWithCNI[2]=-3.36983;
  m_logOddsWithCNI[3]=-3.05712;
  m_logOddsWithCNI[4]=-2.79465;
  m_logOddsWithCNI[5]=-2.57589;
  m_logOddsWithCNI[6]=-2.39153;
  m_logOddsWithCNI[7]=-2.22484;
  m_logOddsWithCNI[8]=-2.07675;
  m_logOddsWithCNI[9]=-1.9475;
  m_logOddsWithCNI[10]=-1.83024;
  m_logOddsWithCNI[11]=-1.71716;
  m_logOddsWithCNI[12]=-1.61148;
  m_logOddsWithCNI[13]=-1.51311;
  m_logOddsWithCNI[14]=-1.41839;
  m_logOddsWithCNI[15]=-1.33093;
  m_logOddsWithCNI[16]=-1.2478;
  m_logOddsWithCNI[17]=-1.16848;
  m_logOddsWithCNI[18]=-1.08818;
  m_logOddsWithCNI[19]=-1.01368;
  m_logOddsWithCNI[20]=-0.941014;
  m_logOddsWithCNI[21]=-0.869332;
  m_logOddsWithCNI[22]=-0.801826;
  m_logOddsWithCNI[23]=-0.739038;
  m_logOddsWithCNI[24]=-0.674558;
  m_logOddsWithCNI[25]=-0.607569;
  m_logOddsWithCNI[26]=-0.543946;
  m_logOddsWithCNI[27]=-0.48515;
  m_logOddsWithCNI[28]=-0.424446;
  m_logOddsWithCNI[29]=-0.365842;
  m_logOddsWithCNI[30]=-0.308674;
  m_logOddsWithCNI[31]=-0.252784;
  m_logOddsWithCNI[32]=-0.196686;
  m_logOddsWithCNI[33]=-0.139565;
  m_logOddsWithCNI[34]=-0.0848028;
  m_logOddsWithCNI[35]=-0.031083;
  m_logOddsWithCNI[36]=0.0184339;
  m_logOddsWithCNI[37]=0.0720594;
  m_logOddsWithCNI[38]=0.125749;
  m_logOddsWithCNI[39]=0.17719;
  m_logOddsWithCNI[40]=0.222605;
  m_logOddsWithCNI[41]=0.268806;
  m_logOddsWithCNI[42]=0.319009;
  m_logOddsWithCNI[43]=0.369988;
  m_logOddsWithCNI[44]=0.421481;
  m_logOddsWithCNI[45]=0.467459;
  m_logOddsWithCNI[46]=0.519351;
  m_logOddsWithCNI[47]=0.571818;
  m_logOddsWithCNI[48]=0.620392;
  m_logOddsWithCNI[49]=0.666308;
  m_logOddsWithCNI[50]=0.710019;
  m_logOddsWithCNI[51]=0.756857;
  m_logOddsWithCNI[52]=0.793204;
  m_logOddsWithCNI[53]=0.842943;
  m_logOddsWithCNI[54]=0.890919;
  m_logOddsWithCNI[55]=0.943121;
  m_logOddsWithCNI[56]=0.989365;
  m_logOddsWithCNI[57]=1.03854;
  m_logOddsWithCNI[58]=1.08265;
  m_logOddsWithCNI[59]=1.13028;
  m_logOddsWithCNI[60]=1.18009;
  m_logOddsWithCNI[61]=1.22109;
  m_logOddsWithCNI[62]=1.26771;
  m_logOddsWithCNI[63]=1.32197;
  m_logOddsWithCNI[64]=1.37783;
  m_logOddsWithCNI[65]=1.41825;
  m_logOddsWithCNI[66]=1.4635;
  m_logOddsWithCNI[67]=1.50981;
  m_logOddsWithCNI[68]=1.56577;
  m_logOddsWithCNI[69]=1.6135;
  m_logOddsWithCNI[70]=1.65517;
  m_logOddsWithCNI[71]=1.70645;
  m_logOddsWithCNI[72]=1.758;
  m_logOddsWithCNI[73]=1.80127;
  m_logOddsWithCNI[74]=1.84201;
  m_logOddsWithCNI[75]=1.89882;
  m_logOddsWithCNI[76]=1.94849;
  m_logOddsWithCNI[77]=2.00161;
  m_logOddsWithCNI[78]=2.06052;
  m_logOddsWithCNI[79]=2.13645;
  m_logOddsWithCNI[80]=2.19232;
  m_logOddsWithCNI[81]=2.23578;
  m_logOddsWithCNI[82]=2.2549;
  m_logOddsWithCNI[83]=2.30154;
  m_logOddsWithCNI[84]=2.34564;
  m_logOddsWithCNI[85]=2.41686;
  m_logOddsWithCNI[86]=2.43023;
  m_logOddsWithCNI[87]=2.49433;
  m_logOddsWithCNI[88]=2.49112;
  m_logOddsWithCNI[89]=2.55153;
  m_logOddsWithCNI[90]=2.60136;
  m_logOddsWithCNI[91]=2.68186;
  m_logOddsWithCNI[92]=2.67592;
  m_logOddsWithCNI[93]=2.65865;
  m_logOddsWithCNI[94]=2.69613;
  m_logOddsWithCNI[95]=2.66071;
  m_logOddsWithCNI[96]=2.33451;
  m_logOddsWithCNI[97]=2.70504;
  m_logOddsWithCNI[98]=3.46921;
  m_logOddsWithCNI[99]=4.04662;
  m_logOddsWithCNI[100]=3.89577;

  m_logOddsWithMzPos[0]=-1.3336;
  m_logOddsWithMzPos[1]=-1.3336;
  m_logOddsWithMzPos[2]=-1.3336;
  m_logOddsWithMzPos[3]=-1.3336;
  m_logOddsWithMzPos[4]=-1.3336;
  m_logOddsWithMzPos[5]=-1.3336;
  m_logOddsWithMzPos[6]=-1.3336;
  m_logOddsWithMzPos[7]=-1.3336;
  m_logOddsWithMzPos[8]=-1.30137;
  m_logOddsWithMzPos[9]=-1.4293;
  m_logOddsWithMzPos[10]=-1.61753;
  m_logOddsWithMzPos[11]=-1.73912;
  m_logOddsWithMzPos[12]=-1.72183;
  m_logOddsWithMzPos[13]=-1.75837;
  m_logOddsWithMzPos[14]=-1.83157;
  m_logOddsWithMzPos[15]=-1.8871;
  m_logOddsWithMzPos[16]=-1.91406;
  m_logOddsWithMzPos[17]=-1.88434;
  m_logOddsWithMzPos[18]=-1.87967;
  m_logOddsWithMzPos[19]=-1.89946;
  m_logOddsWithMzPos[20]=-1.93436;
  m_logOddsWithMzPos[21]=-1.90978;
  m_logOddsWithMzPos[22]=-1.8728;
  m_logOddsWithMzPos[23]=-1.86679;
  m_logOddsWithMzPos[24]=-1.9159;
  m_logOddsWithMzPos[25]=-1.91467;
  m_logOddsWithMzPos[26]=-1.90991;
  m_logOddsWithMzPos[27]=-1.88939;
  m_logOddsWithMzPos[28]=-1.88594;
  m_logOddsWithMzPos[29]=-1.81491;
  m_logOddsWithMzPos[30]=-1.76876;
  m_logOddsWithMzPos[31]=-1.73546;
  m_logOddsWithMzPos[32]=-1.7187;
  m_logOddsWithMzPos[33]=-1.67569;
  m_logOddsWithMzPos[34]=-1.61635;
  m_logOddsWithMzPos[35]=-1.54452;
  m_logOddsWithMzPos[36]=-1.48825;
  m_logOddsWithMzPos[37]=-1.45447;
  m_logOddsWithMzPos[38]=-1.43703;
  m_logOddsWithMzPos[39]=-1.40376;
  m_logOddsWithMzPos[40]=-1.37975;
  m_logOddsWithMzPos[41]=-1.34361;
  m_logOddsWithMzPos[42]=-1.30068;
  m_logOddsWithMzPos[43]=-1.26981;
  m_logOddsWithMzPos[44]=-1.25402;
  m_logOddsWithMzPos[45]=-1.23638;
  m_logOddsWithMzPos[46]=-1.22953;
  m_logOddsWithMzPos[47]=-1.2325;
  m_logOddsWithMzPos[48]=-1.16099;
  m_logOddsWithMzPos[49]=-1.11821;
  m_logOddsWithMzPos[50]=-1.15766;
  m_logOddsWithMzPos[51]=-1.23143;
  m_logOddsWithMzPos[52]=-1.17794;
  m_logOddsWithMzPos[53]=-1.08828;
  m_logOddsWithMzPos[54]=-1.07973;
  m_logOddsWithMzPos[55]=-1.13349;
  m_logOddsWithMzPos[56]=-1.1993;
  m_logOddsWithMzPos[57]=-1.21481;
  m_logOddsWithMzPos[58]=-1.23543;
  m_logOddsWithMzPos[59]=-1.24045;
  m_logOddsWithMzPos[60]=-1.25301;
  m_logOddsWithMzPos[61]=-1.23781;
  m_logOddsWithMzPos[62]=-1.23868;
  m_logOddsWithMzPos[63]=-1.22308;
  m_logOddsWithMzPos[64]=-1.20037;
  m_logOddsWithMzPos[65]=-1.18554;
  m_logOddsWithMzPos[66]=-1.20716;
  m_logOddsWithMzPos[67]=-1.2271;
  m_logOddsWithMzPos[68]=-1.24172;
  m_logOddsWithMzPos[69]=-1.23101;
  m_logOddsWithMzPos[70]=-1.21749;
  m_logOddsWithMzPos[71]=-1.21544;
  m_logOddsWithMzPos[72]=-1.2234;
  m_logOddsWithMzPos[73]=-1.20528;
  m_logOddsWithMzPos[74]=-1.20659;
  m_logOddsWithMzPos[75]=-1.22384;
  m_logOddsWithMzPos[76]=-1.23646;
  m_logOddsWithMzPos[77]=-1.24574;
  m_logOddsWithMzPos[78]=-1.25593;
  m_logOddsWithMzPos[79]=-1.24542;
  m_logOddsWithMzPos[80]=-1.23818;
  m_logOddsWithMzPos[81]=-1.23835;
  m_logOddsWithMzPos[82]=-1.2218;
  m_logOddsWithMzPos[83]=-1.18733;
  m_logOddsWithMzPos[84]=-1.20016;
  m_logOddsWithMzPos[85]=-1.24036;
  m_logOddsWithMzPos[86]=-1.27879;
  m_logOddsWithMzPos[87]=-1.28698;
  m_logOddsWithMzPos[88]=-1.28077;
  m_logOddsWithMzPos[89]=-1.25385;
  m_logOddsWithMzPos[90]=-1.25648;
  m_logOddsWithMzPos[91]=-1.2699;
  m_logOddsWithMzPos[92]=-1.3137;
  m_logOddsWithMzPos[93]=-1.39817;
  m_logOddsWithMzPos[94]=-1.4997;
  m_logOddsWithMzPos[95]=-1.52831;
  m_logOddsWithMzPos[96]=-1.58549;
  m_logOddsWithMzPos[97]=-1.55305;
  m_logOddsWithMzPos[98]=-1.1903;
  m_logOddsWithMzPos[99]=-1.65296;
  m_logOddsWithMzPos[100]=-2.18582;

  m_logOddsWithSister[0]=0.412005;
  m_logOddsWithSister[1]=0.792883;
  m_logOddsWithSister[2]=0.438032;
  m_logOddsWithSister[3]=0.105349;
  m_logOddsWithSister[4]=0.258077;
  m_logOddsWithSister[5]=0.51076;
  m_logOddsWithSister[6]=0.619873;
  m_logOddsWithSister[7]=0.157325;
  m_logOddsWithSister[8]=0.0985658;
  m_logOddsWithSister[9]=0.105978;
  m_logOddsWithSister[10]=0.0983341;
  m_logOddsWithSister[11]=0.122971;
  m_logOddsWithSister[12]=0.189773;
  m_logOddsWithSister[13]=0.0954715;
  m_logOddsWithSister[14]=0.0194936;
  m_logOddsWithSister[15]=0.286911;
  m_logOddsWithSister[16]=0.253486;
  m_logOddsWithSister[17]=0.212954;
  m_logOddsWithSister[18]=0.259886;
  m_logOddsWithSister[19]=0.261061;
  m_logOddsWithSister[20]=0.280494;
  m_logOddsWithSister[21]=0.236982;
  m_logOddsWithSister[22]=0.0315007;
  m_logOddsWithSister[23]=0.265411;
  m_logOddsWithSister[24]=0.0808984;
  m_logOddsWithSister[25]=0.165485;
  m_logOddsWithSister[26]=0.126737;
  m_logOddsWithoutSister[0]=-0.214757;
  m_logOddsWithoutSister[1]=-0.751031;
  m_logOddsWithoutSister[2]=-0.335278;
  m_logOddsWithoutSister[3]=-0.0481388;
  m_logOddsWithoutSister[4]=-0.11088;
  m_logOddsWithoutSister[5]=-0.28493;
  m_logOddsWithoutSister[6]=-0.405398;
  m_logOddsWithoutSister[7]=-0.061188;
  m_logOddsWithoutSister[8]=-0.0358326;
  m_logOddsWithoutSister[9]=-0.0385674;
  m_logOddsWithoutSister[10]=-0.0344257;
  m_logOddsWithoutSister[11]=-0.0448237;
  m_logOddsWithoutSister[12]=-0.0720518;
  m_logOddsWithoutSister[13]=-0.0334376;
  m_logOddsWithoutSister[14]=-0.00615507;
  m_logOddsWithoutSister[15]=-0.119239;
  m_logOddsWithoutSister[16]=-0.101724;
  m_logOddsWithoutSister[17]=-0.081231;
  m_logOddsWithoutSister[18]=-0.0998066;
  m_logOddsWithoutSister[19]=-0.100301;
  m_logOddsWithoutSister[20]=-0.112067;
  m_logOddsWithoutSister[21]=-0.0885752;
  m_logOddsWithoutSister[22]=-0.00950305;
  m_logOddsWithoutSister[23]=-0.098089;
  m_logOddsWithoutSister[24]=-0.0248924;
  m_logOddsWithoutSister[25]=-0.053336;
  m_logOddsWithoutSister[26]=-0.0385232;

  m_logOddsWithComplement=0.533765;
  m_logOddsWithoutComplement = -0.129599;
  
  
  m_ready = true;
  
}



void SpectraSTDenoiser::addTrainingSpectrum(SpectraSTPeakList* pl) {

  if (m_ready) return;
  
  float cumInten = 0.0;

  for (unsigned int rank = pl->getNumPeaks(); rank >= 1; rank--) {
    
    Peak p;
    pl->getNthLargestPeak(rank, p);
    
    // ignore all precursor peaks?
    if (pl->isNearPrecursor(p.mz)) {
      continue;
    }

    // calculate CNI
    cumInten += p.intensity;
    float cni = cumInten / pl->getTotalIonCurrent();
    int cniBin = (int)(cni * 100.0);
    
    // calculate m/z position
    int mzPosBin = (int)(p.mz / pl->getMaxMz() * 100.0);

    // find complement
    double precursorMass = pl->getParentMz() * (double)(pl->getParentCharge());
    double complement = pl->findPeak(precursorMass - p.mz, 0.25);
    bool hasComplement = false;
    if (complement >= 10.0) hasComplement = true;
    
    if (p.info == "S") {
      m_numSignal++;
      m_signalWithCNI[cniBin]++;
      m_signalWithMzPos[mzPosBin]++; 
      if (hasComplement) {
        m_signalWithComplement++;
      } 
	
    } else {
      m_numNoise++;
      m_noiseWithCNI[cniBin]++;
      m_noiseWithMzPos[mzPosBin]++;
      if (hasComplement) {
        m_noiseWithComplement++;
      } 
    }
    
    // find sisters
    for (unsigned int sisIndex = 0; sisIndex < (unsigned int)(m_sisters.size()); sisIndex++) {
      double mzDiff = m_sisters[sisIndex];
      if (pl->findPeak(p.mz - mzDiff, 0.25) >= 10.0) {
        if (p.info == "S") {
	  m_signalWithSister[sisIndex]++;
	} else {
	  m_noiseWithSister[sisIndex]++;
	}
      }
    }

    
    
  } // END for all peaks 
  
  
}

void SpectraSTDenoiser::generateBayesianModel() {
 
  if (m_ready) return;
  
  for (int index = 100; index >= 0; index--) {

    m_logOddsWithCNI[index] = log((double)m_signalWithCNI[index]) - log((double)m_noiseWithCNI[index]);
    if (m_signalWithCNI[index] <= 100 || m_noiseWithCNI[index] <= 100) {
      m_logOddsWithCNI[index] = m_logOddsWithCNI[index+1];

    }
  
    m_logOddsWithMzPos[index] = log((double)m_signalWithMzPos[index]) - log((double)m_noiseWithMzPos[index]);
    if (m_signalWithMzPos[index] <= 100 || m_noiseWithMzPos[index] <= 100) {
      m_logOddsWithMzPos[index] = m_logOddsWithMzPos[index+1];
    }

  }
  
  vector<double> smooth(101,0);
  
  smooth = m_logOddsWithCNI;
  
  m_logOddsWithCNI[0] = (17 * smooth[0] + 12 * smooth[1] - 3 * smooth[2]) / 26.0;
  m_logOddsWithCNI[1] = (17 * smooth[1] + 12 * smooth[2] + 12 * smooth[0] - 3 * smooth[3]) / 38.0;  
  for (unsigned int index = 2; index <= 98; index++) {
    m_logOddsWithCNI[index] = (17 * smooth[index] + 12 * smooth[index+1] + 12 * smooth[index-1] - 3 * smooth[index+2] - 3 * smooth[index-2]) / 35.0;
  }
  m_logOddsWithCNI[99] = (17 * smooth[99] + 12 * smooth[100] + 12 * smooth[98] - 3 * smooth[97]) / 38.0;
  m_logOddsWithCNI[100] = (17 * smooth[100] + 12 * smooth[99] - 3 * smooth[98]) / 26.0;
  
  smooth = m_logOddsWithMzPos;
  
  m_logOddsWithMzPos[0] = (17 * smooth[0] + 12 * smooth[1] - 3 * smooth[2]) / 26.0;
  m_logOddsWithMzPos[1] = (17 * smooth[1] + 12 * smooth[2] + 12 * smooth[0] - 3 * smooth[3]) / 38.0; 
  for (unsigned int index = 2; index <= 98; index++) {
    m_logOddsWithMzPos[index] = (17 * smooth[index] + 12 * smooth[index+1] + 12 * smooth[index-1] - 3 * smooth[index+2] - 3 * smooth[index-2]) / 35.0;
  }
  m_logOddsWithMzPos[99] = (17 * smooth[99] + 12 * smooth[100] + 12 * smooth[98] - 3 * smooth[97]) / 38.0;
  m_logOddsWithMzPos[100] = (17 * smooth[100] + 12 * smooth[99] - 3 * smooth[98]) / 26.0;
  
  for (unsigned int sisIndex = 0; sisIndex < (unsigned int)(m_sisters.size()); sisIndex++) {
    m_logOddsWithSister[sisIndex] = log((double)m_signalWithSister[sisIndex]) + log((double)m_numNoise) - log((double)m_noiseWithSister[sisIndex]) - log((double)m_numSignal);
    m_logOddsWithoutSister[sisIndex] = log((double)(m_numSignal - m_signalWithSister[sisIndex])) + log((double)m_numNoise) - log((double)(m_numNoise - m_noiseWithSister[sisIndex])) - log((double)m_numSignal);
  }

  m_logOddsWithComplement = log((double)m_signalWithComplement) + log((double)m_numNoise) - log((double)m_noiseWithComplement) - log((double)m_numSignal);

  m_logOddsWithoutComplement = log((double)(m_numSignal - m_signalWithComplement)) + log((double)m_numNoise) - log((double)(m_numNoise - m_noiseWithComplement)) - log((double)m_numSignal);

  m_ready = true;
  
  printModel();

}

void SpectraSTDenoiser::printModel() {
  
   cout << "MODEL:" << endl << endl;;
  
  cout << "CNI" << endl;
  for (unsigned int index = 0; index <= 100; index++) {
    cout << "m_logOddsWithCNI[" << index << "]=" << m_logOddsWithCNI[index] << ";" << endl;
  }
  
  cout << "===" << endl;
  
  cout << "MZ" << endl;
  for (unsigned int index = 0; index <= 100; index++) {
    cout << "m_logOddsWithMzPos[" << index << "]=" << m_logOddsWithMzPos[index] << ";" << endl;
  }
  
  cout << "===" << endl;
  
  cout << "Sister" << '\t' << "With" << '\t' << "Without" << endl;
  for (unsigned int sisIndex = 0; sisIndex < (unsigned int)(m_sisters.size()); sisIndex++) {
    cout << "m_logOddsWithSister[" << sisIndex << "]=" << m_logOddsWithSister[sisIndex] << ";" <<endl;
  }
  for (unsigned int sisIndex = 0; sisIndex < (unsigned int)(m_sisters.size()); sisIndex++) {
    cout << "m_logOddsWithoutSister[" << sisIndex << "]=" << m_logOddsWithoutSister[sisIndex] << ";" <<endl;
  }
  
  cout << "===" << endl;
  
  cout << "m_logOddsWithComplement=" << m_logOddsWithComplement << endl;
  cout << "m_logOddsWithoutComplement = " << m_logOddsWithoutComplement << endl;
  
  cout << "===" << endl;
  
  cout << "number of signal = " << m_numSignal << endl;
  cout << "number of noise = " << m_numNoise << endl;
  
}

void SpectraSTDenoiser::filter(SpectraSTPeakList* pl, unsigned int maxNumPeaks, double minSignalProb) {
  
  if (!m_ready) return;
  // calculate signal prob for all peaks 
  
  vector<pair<double, Peak> > peakOdds; // vector of (logOdds, Peak) pairs
  
  float cumInten = 0.0;
  for (unsigned int rank = pl->getNumPeaks(); rank >= 1; rank--) {
    
    Peak p;
    pl->getNthLargestPeak(rank, p);
    
    // ignore all precursor peaks?
    if (pl->isNearPrecursor(p.mz)) {
      continue;
    }

    // calculate CNI
    cumInten += p.intensity;
    float cni = cumInten / pl->getTotalIonCurrent();
    int cniBin = (int)(cni * 100.0);

    
    // calculate m/z position
    int mzPosBin = (int)(p.mz / pl->getMaxMz() * 100.0);

    // find complement
    double precursorMass = pl->getParentMz() * (double)(pl->getParentCharge());
    double complement = pl->findPeak(precursorMass - p.mz, 0.25);
    bool hasComplement = false;
    if (complement >= 10.0) hasComplement = true;
    
    double logOdds = 0.0;
    
    logOdds += m_logOddsWithCNI[cniBin];
    logOdds += m_logOddsWithMzPos[mzPosBin];

    if (hasComplement) {
      logOdds += m_logOddsWithComplement;
    } else {
      logOdds += m_logOddsWithoutComplement;
    }
      
    // find sisters
    for (unsigned int sisIndex = 0; sisIndex < (unsigned int)(m_sisters.size()); sisIndex++) {
      double mzDiff = m_sisters[sisIndex];
      if (pl->findPeak(p.mz - mzDiff, 0.25) >= 10.0) {
 	logOdds += m_logOddsWithSister[sisIndex];
      } else {
        logOdds += m_logOddsWithoutSister[sisIndex];
      }
    }

    pair<double, Peak> op;
    op.first = logOdds;
    op.second = p;
    
    peakOdds.push_back(op);
  } // END for all peaks 
  
  // sort all peaks by odds
  sort(peakOdds.begin(), peakOdds.end(), SpectraSTDenoiser::sortPeaksByOddsDesc);
  
  // insert peaks back into new peak list
  pl->m_peaks.clear();
  
  double prior = 0.5;  
  if (minSignalProb > 0.0) {
    prior = predictPrior(pl);
  }
  
  double logOddsPrior = log(prior / (1.0 - prior));
  double logOddsCutoff = log(minSignalProb / (1.0 - minSignalProb));
  
  unsigned int numPeaksKept = 0;
  
  for (vector<pair<double, Peak> >::iterator i = peakOdds.begin(); i != peakOdds.end(); i++) {
    if (numPeaksKept >= maxNumPeaks) break;
    if (i->first + logOddsPrior < logOddsCutoff) break;
    
    pl->insert(i->second.mz, i->second.intensity, i->second.annotation, i->second.info);
    numPeaksKept++;
  }
  
  sort(pl->m_peaks.begin(), pl->m_peaks.end(), SpectraSTPeakList::sortPeaksByMzAsc);
  pl->m_isSortedByMz = true;
  
  if (pl->m_intensityRanked) {
    delete (pl->m_intensityRanked);
    pl->m_intensityRanked = NULL;
  }

  if (pl->m_peakMap) {
    delete (pl->m_peakMap);
    pl->m_peakMap = NULL;
  }
  
}

double SpectraSTDenoiser::predictPrior(SpectraSTPeakList* pl, bool useOverallPrior) {
  
  if (useOverallPrior) {
    return ((double)m_numSignal / (double)(m_numSignal + m_numNoise));
  }
  
  double prior = 0.2; // default value, based on observation of a large dataset
  
  double precursorMass = pl->getParentMz() * pl->getParentCharge();
  double TIC = pl->getTotalIonCurrent();
  double xrea = pl->calcXrea();
  unsigned int numPeaks = pl->getNumPeaks();
  unsigned int charge = pl->getParentCharge();

  prior = 3.177e-01 + -9.404e-06 * precursorMass + 8.920e-08 * TIC + 1.058e-01 * xrea + - 4.189e-04 * numPeaks + -1.280e-02 * charge;

  if(prior >= 1.0 || prior <= 0) prior = 0.2;
  
  return (prior);
  
}

bool SpectraSTDenoiser::sortPeaksByOddsDesc(pair<double, Peak> a, pair<double, Peak> b) {
 
  return (a.first > b.first);
  
}
