#include "SpectraSTQuery.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include <sstream>

#include <iostream>


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

/* Class: SpectraSTQuery
 * 
 * Manages a query.
 *  
 *  
 */


using namespace std;

extern SpectraSTLog* g_log;

// constructor
SpectraSTQuery::SpectraSTQuery(string name, double precursorMz, int charge, string comments, SpectraSTPeakList* peakList) : 
  m_name(name),
  m_precursorMz(precursorMz),
  m_defaultCharge(charge),
  m_comments(comments),
  m_defaultPeakList(peakList),
  m_retentionTime(-1.0) {	

  if (!m_defaultPeakList) {
    m_defaultPeakList = new SpectraSTPeakList(precursorMz, charge);
  } 
  if (charge != 0) {
    m_peakLists[charge] = m_defaultPeakList;
  }
}

// destructor
SpectraSTQuery::~SpectraSTQuery() {
  
  if (m_defaultPeakList) delete (m_defaultPeakList);
  /*
  cout << "Deleting " << m_peakLists.size() << endl;
  
  for (map<int, SpectraSTPeakList*>::iterator i = m_peakLists.begin(); i != m_peakLists.end(); i++) {
    cout << "  " << i->first << endl;
    if (i->second) {
      cout << "  " << "DEL" << endl;
      delete (i->second);
      i->second = NULL;
    }
  }
  */	
}

// printQuerySpectrum - outputs the query spectrum (as an .msp file) for later spectrum plotting
void SpectraSTQuery::printQuerySpectrum(string querySpectrumFileName) {
	
  ofstream specFout;
  if (!myFileOpen(specFout, querySpectrumFileName)) {
    g_log->error("SEARCH", "Cannot open \"" + querySpectrumFileName + "\" for printing query spectrum.");
    return;
  }
  
  specFout << "Name: " << m_name << endl;
    
  if (m_defaultCharge != 0) {
    specFout << "MW: " << m_precursorMz * m_defaultCharge << endl;
  } 
  specFout << "PrecursorMZ: " << m_precursorMz << endl;
  specFout << "Comment: " << m_comments << endl;
  m_defaultPeakList->writeToFile(specFout);
  specFout << endl;
  	
}

// constructQueryName - construct a query name of the format <prefix>.<scanNum>.<scanNum>.<charge>
// the <scanNum> fields will be padded with preceding zeros up to 5 characters long
string SpectraSTQuery::constructQueryName(string prefix, int scanNum, int charge) {
		
  stringstream ss;
  ss << prefix << '.';
  ss.width(5);
  ss.fill('0');
  ss << right << scanNum;
  ss << '.';
  ss.width(5);
  ss.fill('0');
  ss << right << scanNum;
  ss << '.';
  ss << right << charge;
  return (ss.str());	
  
}

SpectraSTPeakList* SpectraSTQuery::getPeakList(int charge) {
  if (charge == 0) {
    return (m_defaultPeakList);
  } else {
    map<int, SpectraSTPeakList*>::iterator found = m_peakLists.find(charge);
    if (found != m_peakLists.end()) {
      return (found->second);
    } else {
      return (m_defaultPeakList);
    }
  }
}

void SpectraSTQuery::addPossibleCharge(int charge) {
  
  if (m_defaultCharge == 0) m_defaultCharge = charge;
  
  map<int, SpectraSTPeakList*>::iterator found = m_peakLists.find(charge);
  if (found != m_peakLists.end()) {
    return;
  } else {
    m_peakLists[charge] = m_defaultPeakList; // don't copy -- assume same binning/scaling for all charges

    //    m_peakLists[charge] = new SpectraSTPeakList(*m_defaultPeakList); 
//    if (m_defaultPeakList->getFragType() == "ETD") {
//      m_peakLists[charge]->setParentCharge(charge, true);
//    } else {
    m_peakLists[charge]->setParentCharge(charge, false);
//    }
  }
}

bool SpectraSTQuery::isPossibleCharge(int charge) {
  return (charge == 0 || m_defaultCharge == 0 || m_peakLists.find(charge) != m_peakLists.end());
}
