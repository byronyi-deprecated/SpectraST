#ifndef SPECTRASTSEARCHPARAMS_HPP_
#define SPECTRASTSEARCHPARAMS_HPP_

#include <vector>
#include <string>
#include <map>



/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA

Henry Lam
Insitute for Systems Biology
1441 North 34th St. 
Seattle, WA  98103  USA
hlam@systemsbiology.org

*/

/* Class: SpectraSTSearchParams
 * 
 * Manages the search parameters.
 *  
 */



using namespace std;

class SpectraSTSearchParams {

public:
	SpectraSTSearchParams();
	
	SpectraSTSearchParams(SpectraSTSearchParams& s);
	SpectraSTSearchParams& operator=(SpectraSTSearchParams& s);


	virtual ~SpectraSTSearchParams();
	
	// the parameters. all are public for easy retrieval.
	// unfortunately, no measure to prevent them from improperly modified
	
        // GENERAL
        string paramsFileName;
	string libraryFile;
	string databaseFile;
	string databaseType;
        bool indexCacheAll;
        string filterSelectedListFileName; 
        
        // CANDIDATE SELECTION AND SCORING
	string expectedCysteineMod;
	double indexRetrievalMzTolerance;
        bool indexRetrievalUseAverage;
        unsigned int detectHomologs;
        bool ignoreChargeOneLibSpectra;
        bool ignoreAbnormalSpectra;
        bool ignoreSpectraWithUnmodCysteine;
	bool searchAllCharges;
	        
        // OUTPUT AND DISPLAY
        string outputExtension;   
        string outputDirectory;
	double hitListTopHitFvalThreshold;
	double hitListLowerHitsFvalThreshold;
        bool hitListOnlyTopHit;
        bool hitListExcludeNoMatch;
        bool hitListShowHomologs;

  //        bool saveSpectra;
  //      bool tgzSavedSpectra;

        // SPECTRUM FILTERING AND PROCESSING 
        double filterAllPeaksBelowMz;
        unsigned int filterMinPeakCount;
        double filterCountPeakIntensityThreshold;
        double filterRemovePeakIntensityThreshold;
        double filterRemoveHuge515Threshold;
         
        unsigned int filterMaxPeaksUsed;
        double filterMaxDynamicRange;
        double filterMaxIntensityBelow;
	unsigned int filterLibMaxPeaksUsed;
        
        double peakScalingMzPower;
	double peakScalingIntensityPower;
	double peakScalingUnassignedPeaks;
	unsigned int peakBinningNumBinsPerMzUnit;
	double peakBinningFractionToNeighbor;
	
	double fvalFractionDelta;
         
	
	
        
	// methods
	bool addOption(string option);
	void finalizeOptions();
        void finalizeAdvancedOption(string option);
        
	void readFromFile();
        
        void printPepXMLSearchParams(ofstream& fout);
        
	static void printUsage(ostream& out);
        static void printAdvancedOptions(ostream& out);

private:
	// m_options - stores the options as strings, wait until finalizeOptions() is called to actually set them
	vector<string> m_options;

	void setDefault();
        bool isExpectingArg(string option);
	
	

};

#endif /*SPECTRASTSEARCHPARAMS_HPP_*/
