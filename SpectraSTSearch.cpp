#include "SpectraSTSearch.hpp"
#include "SpectraSTLibEntry.hpp"
#include "FileUtils.hpp"

#include <algorithm>
#include <math.h>



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

/* Class: SpectraSTSearch
 * 
 * Implements the search for one query spectrum.
 *  
 */


using namespace std;


extern bool g_verbose;

// constructor
SpectraSTSearch::SpectraSTSearch(SpectraSTQuery* query,	SpectraSTSearchParams& params, SpectraSTSearchOutput* output) :
  m_query(query),
  m_params(params), 
  m_output(output),
  m_candidates() {
  
  m_candidates.clear();	
  
}

// destructor
SpectraSTSearch::~SpectraSTSearch() {
  if (m_query) delete m_query;
	
  for (vector<SpectraSTCandidate*>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
    delete (*i);	
  }
}

// search - main function to perform one search
void SpectraSTSearch::search(SpectraSTLib* lib) {

  double precursorMz = m_query->getPrecursorMz();
  
  if (g_verbose) {
    cout << endl;
    cout << "Now searching query: " << m_query->getName() << " (PrecursorMZ = " << precursorMz;
    
    if (m_query->getDefaultCharge() != 0) {
      cout << "; PrecursorCharge = " << m_query->getDefaultCharge();
      for (int ch = 1; ch <= 7; ch++) {
        if (m_query->isPossibleCharge(ch) && ch != m_query->getDefaultCharge()) cout << "," << ch;
      }
    }
    
    cout << ")" << endl;
  }
 
  // retrieves all entries from the library within the tolerable m/z range
  vector<SpectraSTLibEntry*> entries;
  
  double lowMz = precursorMz - m_params.indexRetrievalMzTolerance;
  if (m_params.indexRetrievalUseAverage) lowMz -= 1.0;
  double highMz = precursorMz + m_params.indexRetrievalMzTolerance;
  
  lib->retrieve(entries, lowMz, highMz);
  
  if (g_verbose) {
    cout << "\tFound " << entries.size() << " candidate(s)... " << " Comparing... ";
    cout.flush();
  }
  
  // for all retrieved entries, do the necessary filtering, add the good ones to m_candidates
  for (vector<SpectraSTLibEntry*>::iterator i = entries.begin(); i != entries.end(); i++) {

    if (!(m_params.ignoreChargeOneLibSpectra && (*i)->getCharge() == 1) &&
        !(!m_params.expectedCysteineMod.empty() && m_params.expectedCysteineMod != "ICAT_cl" && (*i)->isCleavableICAT()) &&
        !(!m_params.expectedCysteineMod.empty() && m_params.expectedCysteineMod != "ICAT_uc" && (*i)->isUncleavableICAT()) &&
        !(!m_params.expectedCysteineMod.empty() && m_params.expectedCysteineMod != "CAM" && (*i)->isCAMCysteine()) &&
        !(m_params.ignoreAbnormalSpectra && (*i)->getStatus() != "Normal") &&
        !(m_params.ignoreSpectraWithUnmodCysteine && (*i)->hasUnmodifiedCysteine()) &&
        (m_params.searchAllCharges || m_query->isPossibleCharge((*i)->getCharge()))) {	
 
      // apply library spectrum simplification
      double retained = (*i)->getPeakList()->simplify(m_params.filterLibMaxPeaksUsed, 999999999);
      
      SpectraSTCandidate* newCandidate = new SpectraSTCandidate(*i, m_params);	
      m_candidates.push_back(newCandidate);
    
  
    } 
  }
  
  
  
  // compare query to each candidate by calculating the dot product
  for (vector<SpectraSTCandidate*>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {

    SpectraSTLibEntry* entry = (*i)->getEntry();
    int charge = entry->getCharge();

    double dot = m_query->getPeakList(charge)->compare(entry->getPeakList(), (*i)->getSimScoresRef(), m_params);

    if (m_params.indexRetrievalUseAverage) {
      ((*i)->getSimScoresRef()).precursorMzDiff = precursorMz - entry->getAveragePrecursorMz();          	 
    } else {
      ((*i)->getSimScoresRef()).precursorMzDiff = precursorMz - entry->getPrecursorMz();    
    }
    (*i)->setSortKey(dot);	
    
  }

  // sort the hits by the sort key 
  // (in this case, the value of "dot" returned by the SpectraSTPeakList::compare function)
  sort(m_candidates.begin(), m_candidates.end(), SpectraSTCandidate::sortPtrsDesc);
  
  if (m_params.detectHomologs > 1) {
    detectHomologs();
  } else {
    if (m_candidates.size() > 1) {
      (m_candidates[0]->getSimScoresRef()).firstNonHomolog = 2;
    }
  }
 
  calcHitsStats();
  
  // to save time, we can throw away all except the top N hits
  // this is assuming nothing beyond N will be worth rescuing to the top due to a favorable dot_bias
  // unsigned int numHitsToKeep = m_params.detectHomologs + 6;
  //while (m_candidates.size() > numHitsToKeep) {
  //  delete (m_candidates.back());
  //  m_candidates.pop_back();
  //}
  
  // now that they are sorted, can calculate delta dots.
  calcDeltaSimpleDots();
	
  // sort again by F value (calculated in calcDeltaSimpleDots)
  sort(m_candidates.begin(), m_candidates.end(), SpectraSTCandidate::sortPtrsDesc);


}

// calcDeltaSimpleDots - calculates the delta dots
void SpectraSTSearch::calcDeltaSimpleDots() {

  if (m_candidates.empty()) {
    return;
  }
  
  for (unsigned int rank = 0; rank < (unsigned int)(m_candidates.size()); rank++) {
    SpectraSTSimScores& s = m_candidates[rank]->getSimScoresRef();
    
    if (s.firstNonHomolog > 0) {
      // delta is dot - dot(first non-homologous hit lower than this one)
      s.delta = s.dot - (m_candidates[s.firstNonHomolog - 1]->getSimScoresRef()).dot;
    } else {
      s.delta = 0.0;
    }
    
    double fval = s.calcFval(m_params.fvalFractionDelta);
    
    m_candidates[rank]->setSortKey(fval);
  }
}

// calcHitsStats - calculates such things as the mean and stdev of the dots of all candidates
void SpectraSTSearch::calcHitsStats() {
  
  unsigned int numHits;
  double mean;
  double stdev;
  
  if (m_candidates.empty()) {
    numHits = 0;
    mean = 0.0;
    stdev = 0.0;
  } else {
    
    double totalDot = 0;
    double totalSqDot = 0;
    for (vector<SpectraSTCandidate*>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
      SpectraSTSimScores& s = (*i)->getSimScoresRef();
      totalDot += s.dot;
      totalSqDot += s.dot * s.dot;
    }
    numHits = (unsigned int)(m_candidates.size());
    mean = totalDot / (double)numHits;
    stdev = totalSqDot / (double)numHits - mean * mean;
    if (stdev > 0.000001) stdev = sqrt(stdev);
  }
  for (vector<SpectraSTCandidate*>::iterator i = m_candidates.begin(); i != m_candidates.end(); i++) {
    SpectraSTSimScores& s = (*i)->getSimScoresRef();
    s.hitsNum = numHits;
    s.hitsMean = mean;
    s.hitsStDev = stdev;
    
  }
  	
}

// detectHomologs - try to see if the lower hits are homologous (or identical) to the first one
void SpectraSTSearch::detectHomologs() {
  
  if (m_candidates.empty() || m_candidates.size() == 1) return;
	
  Peptide* topHit = m_candidates[0]->getEntry()->getPeptidePtr();
  
  if (!topHit) {
    // not a peptide. no notion of homology
    (m_candidates[0]->getSimScoresRef()).firstNonHomolog = 2;  
    return;
  }
    
  bool homologFound = false;
  unsigned int curRank = 0;
  
  // go down the hit lists until hitting a nonhomologous hit
  do {
    homologFound = false;
    curRank++;
    
    Peptide* thisHit = m_candidates[curRank]->getEntry()->getPeptidePtr();
    
    if (!thisHit) {
      // not a peptide. definitely nonhomologous
      break;
    }
    
    int identity = 0;
    if (*topHit == *thisHit) {
      // identical!
      homologFound = true;
    } else if ((topHit->stripped.length() > thisHit->stripped.length() && topHit->isSubsequence(*thisHit, true)) ||
	       (thisHit->stripped.length() < topHit->stripped.length() && thisHit->isSubsequence(*topHit, true))) {
      // one is subsequence of the other!
      homologFound = true;
    } else if (topHit->isHomolog(*thisHit, 0.7, identity)) {
      homologFound = true;
      
    } 
    
  } while (homologFound && curRank < m_params.detectHomologs - 1 && curRank < (unsigned int)(m_candidates.size()) - 1);
  
  // setting the field firstNonHomolog for all the homologs found
  for (unsigned int rank = 0; rank < curRank; rank++) {
    (m_candidates[rank]->getSimScoresRef()).firstNonHomolog = curRank + 1;
  }
			
}



// print - prints out the search result
void SpectraSTSearch::print() {
  
	
  if (m_candidates.empty() || (!(m_candidates[0]->passTopHitFvalThreshold()))) {
    // no hit (or no hit above threshold
    if (g_verbose) {
      cout << " DONE! Top hit: NO_MATCH";
      cout.flush();
    }
    if (m_params.hitListExcludeNoMatch) {
      // told to exclude all NO_MATCH's from the final output, so just return
      return;
    }
  }
  
  // set the assumedCharge to that of the top hit, if any
  int assumedCharge = 0;
  if (!m_candidates.empty() && m_candidates[0]->passTopHitFvalThreshold()) {
    assumedCharge = m_candidates[0]->getEntry()->getCharge();
  }
  
  
  string name = m_query->getName();
  double precursorMz = m_query->getPrecursorMz();
  double rt = m_query->getRetentionTime();
  // prints all the query information
  m_output->printStartQuery(name, precursorMz, assumedCharge, rt);
  
  
  if (m_candidates.empty() || !(m_candidates[0]->passTopHitFvalThreshold())) {
    // no hit or no hit above threshold
    SpectraSTSimScores emptyScore;
		
    m_output->printHit(name, 1, NULL, emptyScore);
		
  } else {
    // print the top hit!
    if (g_verbose) {
      cout << " DONE! Top hit: " << m_candidates[0]->getEntry()->getFullName() << " (Dot = " << m_candidates[0]->getSortKey() << ")";
      cout.flush();
    }
    
    m_output->printHit(name, 1, m_candidates[0]->getEntry(), m_candidates[0]->getSimScoresRef());

    /*    
    if (m_params.saveSpectra) {
      // also save the query spectrum and the spectrum for the top hit
      string querySpectrumFileName(m_output->getOutputPath() + name + ".0.msp");
      m_query->printQuerySpectrum(querySpectrumFileName);
      
      string topMatchSpectrumFileName(m_output->getOutputPath() + name + ".1.msp");		
      m_candidates[0]->printSpectrum(topMatchSpectrumFileName);
    }
    */    
    
    if (!m_params.hitListOnlyTopHit) { 
      // told to print the lower hits too
      
      unsigned int firstNonHomolog = (m_candidates[0]->getSimScoresRef()).firstNonHomolog;
      
      for (unsigned int rank = 2; rank <= (unsigned int)(m_candidates.size()); rank++) {		
        if (m_candidates[rank - 1]->passLowerHitsFvalThreshold() || (m_params.hitListShowHomologs && rank < firstNonHomolog)) {
          m_output->printHit(name, rank, m_candidates[rank - 1]->getEntry(), m_candidates[rank - 1]->getSimScoresRef());
        } else {
          // below threshold now, stop printing
          break;
        }				
      }
    }
    
  }
  
  // print the closing tags to finish up
  m_output->printEndQuery(name);
  
  
}

// isLikelyGood. Simply returns if the F value is above 0.5, indicating a likely good hit (not always!)
bool SpectraSTSearch::isLikelyGood() {
  if (m_candidates.size() > 0 && m_candidates[0]->getSortKey() >= 0.5) {
    return (true);
  } else {
    return (false);
  }
  return (false);
}

// isLikelyBad. Simply returns if the F value is below 0.2.
bool SpectraSTSearch::isLikelyBad() {
  if (m_candidates.size() > 0 && m_candidates[0]->getSortKey() < 0.2) {
    return (true);
  } else {
    return (false);
  }
  return (false);
}

// isDecoy. returns if the hit is a decoy entry
bool SpectraSTSearch::isDecoy(unsigned int rank) {
  if ((unsigned int)(m_candidates.size()) >= rank) {
    string spec("");
    string rem("");
    return ((m_candidates[rank - 1]->getEntry()->getOneComment("Spec", spec) && spec == "Decoy") ||
	    (m_candidates[rank - 1]->getEntry()->getOneComment("Remark", rem) && rem.substr(0,5) == "DECOY"));
  } else {
    return (false);
  }
  return (false);
}

bool SpectraSTSearch::isSingleton(unsigned int rank) {
  if ((unsigned int)(m_candidates.size()) >= rank) {
    string nreps("");
    return (m_candidates[rank - 1]->getEntry()->getNrepsUsed() == 1);
  } 
  return (false);
}
