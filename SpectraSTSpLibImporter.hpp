#ifndef SPECTRASTSPLIBIMPORTER_HPP_
#define SPECTRASTSPLIBIMPORTER_HPP_

#include "SpectraSTLibImporter.hpp"
#include "SpectraSTDenoiser.hpp"
#include <map>
#include <set>

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

/* Class: SpectraSTSpLibImporter
 * 
 * Implements a library importer for the .splib file format. This takes a SpectraST processed
 * .splib file (which is already searchable), and performs certain actions to the library,
 * and writes the resulting new library to a new .splib file (with the accompanying indices).
 * 
 * These actions include:
 * - Uniquify spectra by (a) taking the best replicate; (b) taking consensus of all replicates
 * - Union (based on peptide) multiple libraries
 * - Intersect (based on peptide) multiple libraries
 * - Subtract one library from another (based on peptide)
 * 
 * Some others -- still under MAJOR CONSTRUCTION!
 * 
 */

// a bunch of numbers for keeping track of how many spectra failed each quality filters
// e.g. Q1Q3Q4 is the number of spectra that failed 3 filters (1 = inquorate, 3 = unconfirmed, 4 = conflicting ID)
typedef struct _qfstats {
  unsigned int immuneProb, immuneEngine;
  unsigned int Q1, Q2, Q3, Q4, Q5;
  unsigned int Q1Q2, Q1Q3, Q1Q4, Q1Q5, Q2Q3, Q2Q4, Q2Q5, Q3Q4, Q3Q5, Q4Q5;
  unsigned int Q1Q2Q3, Q1Q2Q4, Q1Q2Q5, Q1Q3Q4, Q1Q3Q5, Q1Q4Q5, Q2Q3Q4, Q2Q3Q5, Q2Q4Q5, Q3Q4Q5;
  unsigned int Q1Q2Q3Q4, Q1Q2Q3Q5, Q1Q2Q4Q5, Q1Q3Q4Q5, Q2Q3Q4Q5;
  unsigned int Q1Q2Q3Q4Q5;
} QFStats;


class SpectraSTSpLibImporter : public SpectraSTLibImporter { 

public:

  SpectraSTSpLibImporter(vector<string>& impFileNames, SpectraSTLib* lib, SpectraSTCreateParams& params);
  virtual ~SpectraSTSpLibImporter();
  
  virtual void import();
  
  // override SpectraSTLibImporter::constructOutputFileName for splib-specific behavior 
  virtual string constructOutputFileName();

  // override SpectraSTLibImporter::passAllFilters for splib-specific behavior (namely, refresh delete)
  virtual bool passAllFilters(SpectraSTLibEntry* entry);
  
private:
  
  // ifstream's for each of the .splib files to be imported
  vector<ifstream*> m_splibFins;
  
  // SpectraSTPeptideLibIndex objects for each of the .splib files to be imported
  // these are properties of this class
  vector<SpectraSTPeptideLibIndex*> m_pepIndices;
  
  // SpectraSTMzLibIndex objects for each of the .splib files to be imported
  // these are properties of this class  
  vector<SpectraSTMzLibIndex*> m_mzIndices;

  // a map to store all peptide-protein mappings during database refresh
  map<string, vector<pair<string, string> >* >* m_ppMappings;

  // a library object of the imported .splib instantiated for retrieval. 
  // and a default search params object. 
  // these are used for quality filter - conflicting ID
  // both are properties of this class
  SpectraSTLib* m_QFSearchLib;
  SpectraSTSearchParams* m_QFSearchParams;
  
  string m_plotPath;
  unsigned int m_count;
  
  SpectraSTDenoiser* m_denoiser;
  vector<pair<string, string> > m_singletonPeptideIons;

  // method to parse preambles of the imported .splib files
  void parsePreamble(ifstream& splibFin, bool binary);

  // join methods  
  void doSubtractHomologs();
  
  // uniquify methods
  void doBuildAction(vector<SpectraSTLibEntry*>& entries);
  SpectraSTLibEntry* findBestReplicate(vector<SpectraSTLibEntry*>& entries);
  SpectraSTLibEntry* makeConsensusSpectrum(vector<SpectraSTLibEntry*>& entries);

  // helper methods
  void processEntry(SpectraSTLibEntry* entry);
  string constructFileListStr();
  void openSplibs(bool openMzIndex, double mzIndexCacheRange, 
		  bool openPepIndex, bool checkUniqueness, bool refresh);

  // plot method
  void plot(SpectraSTLibEntry* entry);
  
  // quality filter methods
  void doQualityFilter();  
  bool applyQualityFilter(SpectraSTLibEntry* entry, QFStats& qfstats);
  bool isImpure(SpectraSTLibEntry* entry, unsigned int numUsedReps, bool penalizeSingletons, string& msg);
  bool isBadConflictingID(SpectraSTLibEntry* entry, unsigned int numUsedReps, bool penalizeSingletons, string& msg);
  bool hasSharedSequence(SpectraSTLibEntry* entry);
  
  // decoy generation methods
  void doGenerateDecoy();
    
  // sort by number of replicates
  void doSortByNreps();

  // spectrum perturbation by user-specified modifications
  void doUserSpecifiedModifications();
  void parseAllowableTokensStr(string& allowableTokensStr, vector<map<char, set<string> > >& allowableTokens);
  
  // similarity clustering
  void doSimilarityClustering();
  void findSpectralNeighbors(SpectraSTLibEntry* entry, double rootPrecursorMz, unsigned int round, vector<SpectraSTLibEntry*>& isobaricEntries, set<fstream::off_type>* cluster);
  
  
  // refresh peptide-protein mappings
  void addSequencesForRefresh(vector<string>& seqs);
  void refresh();
  
  void reloadAndProcessSingletons();
  
   bool hackDeamidation(SpectraSTLibEntry* entry);

};

#endif /*SPECTRASTSPLIBIMPORTER_HPP_*/
