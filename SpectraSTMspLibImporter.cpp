#include "SpectraSTMspLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"
#include <iostream>
#include <sstream>
#include <stdlib.h>


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

/* Class: SpectraSTMspLibImporter
 * 
 * Implements a library importer for the .msp file format (used by NIST).
 * 
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog* g_log;

// constructor - will open the first file
SpectraSTMspLibImporter::SpectraSTMspLibImporter(vector<string>& impFileNames, SpectraSTLib* lib, SpectraSTCreateParams& params) :
  SpectraSTLibImporter(impFileNames, lib, params),
  m_count(0) {
  
  
}

// destructor 
SpectraSTMspLibImporter::~SpectraSTMspLibImporter() {
  
}

// import - prints the preamble, then loops over all files and import them one by one
void SpectraSTMspLibImporter::import() {
  
  for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
    string fullName(*i);
    makeFullPath(fullName);
    string quoted("\"" + fullName + "\"");
    string desc = m_params.constructDescrStr(quoted, ".msp");
    m_preamble.push_back(desc);
  }
  
  m_lib->writePreamble(m_preamble);
  
  for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
    readFromFile(*i);
  }
  
}

// readFromCurFile - reads one .msp file
void SpectraSTMspLibImporter::readFromFile(string& impFileName) {
  
  ifstream fin;

  if (!myFileOpen(fin, impFileName)) {
    g_log->error("MSP IMPORT", "Cannot open " + impFileName + " for reading msp-format (NIST) library");
    return;
  }
  
  g_log->log("MSP IMPORT", "Importing .msp file \"" + impFileName + "\"."); 
  
  string line("");
  string::size_type pos = 0;
  
  if (g_verbose) {
    cout << "\nImporting spectra from .msp library file..." << endl;
  } 
  
  // start the progress count
  ProgressCount pc(!g_quiet && !g_verbose, 1000);
  pc.start("\nImporting spectra from .msp library file");
  
  while (true) {
    
    if (line.compare(0, 6, "Name: ") != 0) {
      // first record, skips over all lines until the line with Name: 
      while (nextLine(fin, line, "Name: ", ""));
      if (line == "_EOF_") {
        // empty file
	break;
      }
    }
    
    // line now starts with "Name: "
    string name = nextToken(line, 5, pos, "\r\n");
    double mw = 0.0;
    double precursorMz = 0.0;
    string comments;
    while (nextLine(fin, line, "Num peaks: ", "")) {
      // here are all the headers
      if (line.compare(0, 3, "MW:") == 0) {
        // only used for non-peptides. For peptide this will be calculated based on sequence.
	mw = atof(nextToken(line, 3, pos, "\r\n").c_str());
      } else if (line.compare(0, 12, "PrecursorMZ:") == 0) {
	precursorMz = atof(nextToken(line, 12, pos, "\r\n").c_str());
      } else if (line.compare(0, 8, "Comment:") == 0) {
        comments = nextToken(line, 8, pos, "\r\n");
      } else if (line.compare(0, 9, "NumPeaks:") == 0) {
	// reach here because this .msp probably is converted from a .sptxt,
	// hack the line (pad one char) so that the number of peaks will be correctly parsed later
	line = " " + line;
	break;	
      } else if (line == "_EOF_" || line.compare(0, 5, "Name:") == 0) {
        // reach the end unexpectedly, or see another name field before the Num peaks field. 
        // ignore this incomplete record, and return
        g_log->error("MSP IMPORT", "Badly formatted .msp file! Library creation truncated at point of error.");
        pc.done();
        return;
      } else {
        //        g_log->log("MSP IMPORT", "Unrecognized header field \"" + line + "\" ignored.");
      }
    }
    if (line == "_EOF_") {
      // no "Num peaks:" field. ignore this incomplete record, and return
      g_log->error("MSP IMPORT", "Badly formatted .msp file. Library creation truncated at point of error.");
      pc.done();
      return;
    }
    
    
    
    if (g_verbose) {
      cout << "Importing record #" << m_count << ": " << name << endl;
    }
    
    string mods("");
    string::size_type modsPos = 0;
    if ((modsPos = comments.find("Mods=", 0)) != string::npos) {
      if (modsPos + 5 > comments.length() - 1) {
        g_log->error("MSP IMPORT", "Badly formatted .msp file. Library creation truncated at point of error.");
        pc.done();
        return;
      }
      mods = nextToken(comments, modsPos + 5, modsPos, " \t\r\n");		
    }
    
    // parse out full name
    string fullName("");
    string::size_type fullNamePos = 0;
    if ((fullNamePos = comments.find("Fullname=", 0)) != string::npos) {
      if (fullNamePos + 9 > comments.length() - 1) {
        g_log->error("MSP IMPORT", "Badly formatted .msp file. Library creation truncated at point of error.");
        pc.done();
        return;
      }
      fullName = nextToken(comments, fullNamePos + 9, fullNamePos, " \t\r\n");		
    } else {
      // can't find it. just use name
      fullName = name;
    }
    
    // to avoid trouble, removes the '(O)' piece for methionine oxidation (mods should take care of it anyway)
    string::size_type openParenPos = fullName.find("(O)");
    
    if (openParenPos != string::npos) {
      string fixedFullName("");
      do {
        fixedFullName += fullName.substr(0, openParenPos);
        fullName = fullName.substr(openParenPos + 3);
        openParenPos = fullName.find("(O)");
      } while (openParenPos != string::npos);
      
      fullName = fixedFullName + fullName;
    }
    
    SpectraSTLibEntry* entry = NULL;
    
    if (fullName[0] == '_') {
      
      if (precursorMz < 0.001) {
	// no precursor m/z specified,
        unsigned int charge = 0;
        string::size_type slashPos = fullName.rfind('/');
        if (slashPos != string::npos && slashPos < fullName.length() - 1) {
          charge = atoi(fullName.substr(slashPos + 1).c_str());
        }
	// try to calculate from MW and charge
        if (charge > 0) {
	  precursorMz = mw / (double)charge;
	} else {
	  precursorMz = mw;
	}
      }

      entry = new SpectraSTLibEntry(fullName, precursorMz, comments, "Normal", NULL);
    
    } else {
      
      Peptide* pep = new Peptide(fullName, 0, mods);

      // check legality of peptide and mod strings -- will not insert later if illegal (parsing will continue anyway)
      if (pep->illegalMspModStr) {
	g_log->error("MSP IMPORT", "Illegal modification string : \"" + mods + "\". Entry skipped.");
      }
      if (pep->illegalPeptideStr) {
	g_log->error("MSP IMPORT", "Illegal peptide ID string: \"" + fullName + "\". Entry skipped.");
      }

      // create the entry - without the peak list yet
      entry = new SpectraSTLibEntry(pep, comments, "Normal", NULL);

      stringstream nttss;
      nttss << pep->NTT();
      entry->setOneComment("NTT", nttss.str());
      
      stringstream nmcss;
      nmcss << pep->NMC();
      entry->setOneComment("NMC", nmcss.str());

      stringstream naass;
      naass << pep->NAA();
      entry->setOneComment("NAA", naass.str());


      // NIST's format for the protein string is a mess. Instead of trying to parse it and make sense out of it,
      // we'll just take the first word as our protein and ignore the rest. To preserve the information however,
      // we will dump that messy protein string entirely into a new attribute called "NISTProtein"
	  
      // This has the unfortunate effect that all peptides will appear to map to a single protein, while in reality
      // it may be multi-mapping. But since NIST doesn't use nonredundant databases, it's impossible anyway to tell whether
      // a peptide maps to several different proteins, or to the same protein with different names.
      string proteinStr("");
      if (entry->getOneComment("Protein", proteinStr)) {
	
	entry->setOneComment("NISTProtein", proteinStr);

	string::size_type spacePos = 0;
	string firstProtein = nextToken(proteinStr, 0, spacePos, " \t\r\n", " ");
	entry->setOneComment("Protein", "1/" + firstProtein);

	/*

	int proCount = 0;
	string::size_type semiPos = 0;
	string pro("");
	string newProtein("");
	while (!((pro = nextToken(proteinStr, semiPos, semiPos, ";\t\r\n", "; ")).empty())) {
	  string::size_type spacePos = 0;
	  string shortPro = nextToken(pro, 0, spacePos, " \t\r\n", " ");
  //       if (shortPro.find("|sp|") != string::npos) {
	    newProtein += "/" + shortPro; 
	    proCount++;
  //       }
	}
	
	if (proCount > 0) {
	  stringstream newProteinss;
	  newProteinss << proCount << newProtein;
	  entry->setOneComment("Protein", newProteinss.str());
    //    } else {
    //      entry->setOneComment("Protein", "UNMAPPED_TO_SP");
	}

	*/
      }
      
      string tfRatioStr("");
      string::size_type tfRatioPos = 0;
      if ((tfRatioPos = comments.find("Tfratio=", 0)) != string::npos) {
	if (tfRatioPos + 8 > comments.length() - 1) {
	  g_log->error("MSP IMPORT", "Badly formatted .msp file. Library creation truncated at point of error.");
	  pc.done();
	  return;
	}
	tfRatioStr = nextToken(comments, tfRatioPos + 8, tfRatioPos, " \t\r\n");
	double tfRatio = atof(tfRatioStr.c_str());
	double prob = tfRatio / (1.0 + tfRatio);
	stringstream probss;
	probss.precision(4);
	probss << fixed << prob;
	entry->setOneComment("Prob", probss.str());
	
      } else {
	entry->setOneComment("Prob", "1.0000"); // somehow can't find Tfratio. Make it Prob=1.0...
      }
      
      entry->getPeakList()->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);    
      
    }
    
    if (!(m_params.setFragmentation.empty())) {
      entry->setFragType(m_params.setFragmentation);
    } else {
      string fragType("CID");
      entry->setFragType(fragType); // default to CID -- better than leaving it blank!
    }

    // line now starts with "Num peaks: "
    int numPeaks = atoi((nextToken(line, 10, pos, "\r\n")).c_str());
    
    while (nextLine(fin, line, "Name: ", "")) {
      if (line == "_EOF_") {
        if (passAllFilters(entry)) { 
          pc.increment();
          m_count++;
          			
          entry->annotatePeaks(true);
          m_lib->insertEntry(entry); 
            
        }
        
	delete entry;
 
        stringstream countss;
        countss << "Total of " << m_count << " spectra imported.";
        g_log->log("MSP IMPORT", countss.str());
        pc.done();
        return;
	
      }
      // here are the peaks		
      double mz = atof((nextToken(line, 0, pos)).c_str());
      // pos1 now stores the position of the first space after the mz
      float intensity = atof((nextToken(line, pos, pos)).c_str());
      // pos2 now stores the position of the first space after the intensity
      
      // annotation has quotes around it, remove them by adding the quote char to the skipover and delimiter strings passed into nextToken
      string annotation = nextToken(line, pos, pos, "\"\r\n", "\"\t");
      
      // break annotation into real annotation and info
      string info("");
      
      string::size_type spacePos = annotation.find_first_of(" \t", 0);
      string::size_type dummyPos = 0;
      if (spacePos != string::npos) {
        info = nextToken(annotation, spacePos + 1, dummyPos, "\r\n", " \t");
        annotation = annotation.substr(0, spacePos);
      }
      // annotation will get an empty string if there's no annotation
      entry->getPeakList()->insert(mz, intensity, annotation, info);
    }	
    
    if (passAllFilters(entry)) { 

      pc.increment();
      m_count++;	
      		
      entry->annotatePeaks(true);
      m_lib->insertEntry(entry); 
    }
    
    delete entry;
    
    
  }
  
  stringstream countss1;
  countss1 << "Total of " << m_count << " spectra imported.";
  g_log->log("MSP IMPORT", countss1.str());
  
  pc.done();
  
}





