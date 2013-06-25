#include "SpectraSTLib.hpp"
#include "SpectraSTLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"

#include "Peptide.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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

/* Class: SpectraSTLib
 *
 * The class that represents the library in memory. NOTE however that since the library
 * file .splib is huge, it is typically not read into memory in its entirety. Instead
 * SpectraSTLib only keeps the index and the fstream objects in memory, and only goes
 * to the .splib file to retrieve the entries when it is asked to
 * (by a call to retrieve() or cachedRetrieve()).
 *
 */

using namespace std;

extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog* g_log;

// Constructor for creation
SpectraSTLib::SpectraSTLib(vector<string>& impFileNames, SpectraSTCreateParams* createParams) :
    m_libFin(),
    m_libFout(),
    m_txtFout(),
    m_libFileName(),
    m_impFileNames(impFileNames),
    m_txtFileName(),
    m_mzIdxFileName(),
    m_pepIdxFileName(),
    m_newLibId(0),
    m_mzIndex(NULL),
    m_pepIndex(NULL),
    m_searchParams(NULL),
    m_createParams(createParams),
    m_count(0),
    m_noSptxt(false),
    m_mrmFout(NULL),
    m_mgfFout(NULL) {

    initializeLibCreateMode();

}

// Constructor for searching
SpectraSTLib::SpectraSTLib(string fullFileName, SpectraSTSearchParams* searchParams, bool loadPeptideIndex) :
    m_libFin(),
    m_libFout(),
    m_libFileName(fullFileName),
    m_impFileNames(),
    m_mzIdxFileName(),
    m_pepIdxFileName(),
    m_newLibId(0),
    m_mzIndex(NULL),
    m_pepIndex(NULL),
    m_searchParams(searchParams),
    m_createParams(NULL),
    m_count(0),
    m_noSptxt(false),
    m_mrmFout(NULL),
    m_mgfFout(NULL) {

//    initializeLibSearchMode(loadPeptideIndex);
    initializeDatabase();
    resetCache();
}

// Destructor
SpectraSTLib::~SpectraSTLib() {

    m_libFin.close();
    m_libFout.close();

    if (m_mzIndex) {
        delete m_mzIndex;
    }
    if (m_pepIndex) {
        delete m_pepIndex;
    }
    if (m_mrmFout) {
        delete m_mrmFout;
    }
    if (m_mgfFout) {
        delete m_mgfFout;
    }
    if (m_searchParams) {
        resetCache();
        shutdownDatabase();
    }

}

// initializeLibSearchMode - initializes the library for search mode. 
void SpectraSTLib::initializeLibSearchMode(bool loadPeptideIndex) {

    parseFileName(m_libFileName, m_libFileNameStruct);

    // check to make sure library file has right extension
    if (m_libFileNameStruct.ext != ".splib") {
        g_log->error("SEARCH", "Library specified is not a .splib file. No search is performed.");
        g_log->crash();
        return;
    }

    // open library file
    if (!myFileOpen(m_libFin, m_libFileName, true)) {
        g_log->error("SEARCH", "Cannot open SPLIB file \"" + m_libFileName + " for loading library. No search is performed.");
        g_log->crash();
        return;
    }

    // peek to see if it's a binary file or not
    bool binary = true;
    char firstChar = m_libFin.peek();
    if (firstChar == '#' || firstChar == 'N') {
        binary = false;
    }

    if (m_searchParams->databaseFile.empty()) {
        extractDatabaseFileFromPreamble(binary);
    }

    // assumes that the m/z index file has the same path and same base file name as the .splib file
    m_mzIdxFileName = m_libFileNameStruct.path + m_libFileNameStruct.name + ".spidx";

    // read the m/z index into memory for speedy lookup
    if (m_searchParams->indexCacheAll) {
        m_mzIndex = new SpectraSTMzLibIndex(m_mzIdxFileName, &m_libFin, SpectraSTMzLibIndex::CACHE_ALL, binary);
    } else {
        // *2 because it's +/- the tolerance. +4 to allow for some buffer.
        m_mzIndex = new SpectraSTMzLibIndex(m_mzIdxFileName, &m_libFin, m_searchParams->indexRetrievalMzTolerance * 2.0 + 4.0, binary);
    }

    // if asked, also read the peptide index into memory
    if (loadPeptideIndex) {
        m_pepIdxFileName = m_libFileNameStruct.path + m_libFileNameStruct.name + ".pepidx";
        m_pepIndex = new SpectraSTPeptideLibIndex(m_pepIdxFileName, &m_libFin, binary);
    }
}

// initializeLibSearchMode - initializes the library for create mode. 
void SpectraSTLib::initializeLibCreateMode() {

    // instantiate the appropriate importer (the extension of the file will
    // specify the importer) to read the file
    SpectraSTLibImporter* importer = SpectraSTLibImporter::createSpectraSTLibImporter(m_impFileNames, this, *m_createParams);
    if (!importer) {
        return;
    }

    // importer will generate the output library file name; get it here
    m_libFileName = importer->getOutputFileName();

    parseFileName(m_libFileName, m_libFileNameStruct);
    string pathPlusBaseName = m_libFileNameStruct.path + m_libFileNameStruct.name;

    // try to open the output .splib file for writing
    if (!myFileOpen(m_libFout, m_libFileName, true)) {
        g_log->error("CREATE", "Cannot open SPLIB file \"" + m_libFileName + "\" for writing library. No library creation performed.");
        return;
    }
    
    // if binary format is specified (which is default), also writes a text .sptxt file for human reading
    if (m_createParams->binaryFormat) {
        m_txtFileName = m_libFileNameStruct.path + m_libFileNameStruct.name + ".sptxt";
        if (!myFileOpen(m_txtFout, m_txtFileName)) {
            g_log->error("CREATE", "Cannot open SPTXT file \"" + m_libFileName + "\" for writing text library. Text library writing skipped.");
            m_noSptxt = true;  // somehow failed; flag so that we won't write the .sptxt later on.
        }
    }

    // place the indices in the same path as the .splib
    m_mzIdxFileName = pathPlusBaseName + ".spidx";
    m_pepIdxFileName = pathPlusBaseName + ".pepidx";

    // create an index object to prepare for entries to be inserted
    m_mzIndex = new SpectraSTMzLibIndex(m_mzIdxFileName);

    // create a peptide index object to prepare for entries to be inserted
    m_pepIndex = new SpectraSTPeptideLibIndex(m_pepIdxFileName);

    // make a directory to hold the dta's
    if (m_createParams->writeDtaFiles) {
        makeDir(pathPlusBaseName + "_dtas/");
    }

    if (!(m_createParams->printMRMTable.empty())) {
        string mrmFileName(pathPlusBaseName + ".mrm");
        m_mrmFout = new ofstream();
        if (!myFileOpen(*m_mrmFout, mrmFileName)) {
            g_log->error("CREATE", "Cannot open file \"" + mrmFileName + "\" for writing MRM transition list. Writing skipped.");
            delete (m_mrmFout);
            m_mrmFout = NULL;
        }
    }

    if (m_createParams->writeMgfFile) {
        string mgfFileName(pathPlusBaseName + ".mgf");
        m_mgfFout = new ofstream();
        if (!myFileOpen(*m_mgfFout, mgfFileName)) {
            g_log->error("CREATE", "Cannot open file \"" + mgfFileName + "\" for writing library in MGF format. Writing skilled.");
            delete (m_mgfFout);
            m_mgfFout = NULL;
        }
    }

    if (m_createParams->writePAIdent) {
        string PAIdentFileName(pathPlusBaseName + ".PAIdent");
        m_PAIdentFout = new ofstream();
        if (!myFileOpen(*m_PAIdentFout, PAIdentFileName)) {
            g_log->error("CREATE", "Cannot open file \"" + PAIdentFileName + "\" for writing library in PAIdent format. Writing skilled.");
            delete (m_PAIdentFout);
            m_PAIdentFout = NULL;
        }
    }

    // start importing: the importer will call SpectraSTLib::insertEntry
    // repeated as it parses through the file.
    importer->import();

    // by then, the index object should contain indexes to all the inserted entries. Now write
    // the entire index to the .spidx and .pepidx files for future use.
    m_mzIndex->writeToFile();
    m_pepIndex->writeToFile();

    // display done creation messages
    if (!g_quiet) {
        cout << endl;
        cout << "Library file ";
        if (m_createParams->binaryFormat) {
            cout << "(BINARY) ";
        }
        cout << "\"" << m_libFileName << "\" created." << endl;
        if (m_createParams->binaryFormat && !m_noSptxt) {
            cout << "Library file (TEXT) \"" << m_txtFileName << "\" created." << endl;
        }
        cout << "M/Z Index file \"" << m_mzIdxFileName << "\" created." << endl;
        cout << "Peptide Index file \"" << m_pepIdxFileName << "\" created." << endl;
        if (m_createParams->writeDtaFiles) {
            cout << "Dtas of library spectra created in directory \"" << pathPlusBaseName << "_dtas/\" ." << endl;
        }
        if (m_createParams->writeMgfFile && m_mgfFout) {
            cout << "MGF of library spectra \"" << pathPlusBaseName + ".mgf\" created." << endl;
        }
        if (m_createParams->writePAIdent && m_PAIdentFout) {
            cout << "PAIdent of library spectra \"" << pathPlusBaseName + ".PAIdent\" created." << endl;
        }
        if (!m_createParams->printMRMTable.empty() && m_mrmFout) {
            cout << "MRM Table file \"" << pathPlusBaseName + ".mrm" << "\" created." << endl;
        }

        cout << endl;

        m_pepIndex->printStats(cout);
        cout << endl;
    }

    delete (importer);

}

// insertEntry - insert a library entry. This is called by the SpectraSTLibImporter repeatedly
// as it parses through the imported file. Does a few things: (1) assign a unique libID to each entry,
// (2) writes the entry to the .splib file, while at the same time remembering the file offset, 
// (3) tells the index to store this file offset for this particular entry to future retrieval.
void SpectraSTLib::insertEntry(SpectraSTLibEntry* entry) {

    // check to make sure we're in the Create mode; only then we will try to insert entries!
    if (!m_createParams) {
        return;
    }

    // give this entry a LibID
    entry->setLibId(m_newLibId++);

    if (!(m_createParams->printMRMTable.empty())) {
        // to print MRM table, need accurate mass of fragments. Re-annotate and use theoretical mass
        entry->annotatePeaks(true, true);
    }
    
    
    // remember offset
    fstream::off_type offset = m_libFout.tellp();

    // update counts
    m_count++;

    // put in remark if asked
    if (!(m_createParams->remark.empty())) {
        entry->setOneComment("Remark", m_createParams->remark);
    }

    // write .dta files if asked
    if (m_createParams->writeDtaFiles) {
        // this basically just writes the library entry's
        // peak list out as a dta file, with a descriptive file name
        FileName fn;
        parseFileName(m_libFileName, fn);

        stringstream ss;
        ss << m_libFileNameStruct.path << m_libFileNameStruct.name << "_dtas/";
        ss << entry->getSafeName() << ".dta";

        /*
    // hijack
    ss << fn.name << '.';
    ss.width(5);
    ss.fill('0');
    ss << right << entry->getLibId() + 1;
    ss << '.';
    ss.width(5);
    ss.fill('0');
    ss << right << entry->getLibId() + 1;
    ss << '.' << entry->getCharge();
    ss << ".dta";
    // end hijack
    */

        entry->writeDtaFile(ss.str());
    }

    if (m_createParams->writeMgfFile && m_mgfFout) {
        entry->writeMgfFile(*m_mgfFout);
    }

    if (!(m_createParams->printMRMTable.empty()) && m_mrmFout) {
        entry->writeMRM(*m_mrmFout, m_createParams->printMRMTable);
    }

    if (m_createParams->writePAIdent && m_PAIdentFout) {
        entry->writePAIdent(*m_PAIdentFout, m_libFileNameStruct.name);
    }

    // writing to .splib file. Note that the actual entries are NOT kept in memory. They will be written as they are created, and
    // deleted immediately afterwards.
    if (m_createParams->binaryFormat) {

        // remember the file offset into the binary .splib file
        fstream::off_type binaryFileOffset = m_libFout.tellp();
        stringstream bfoss;
        bfoss << binaryFileOffset;

        // write the library in binary format to .splib
        entry->writeToBinaryFile(m_libFout);

        if (!m_noSptxt) {
            // note the binary file offset as a comment -- just as a convenience for people examining the .sptxt file to view
            // the spectrum using plotspectrast
            entry->setOneComment("BinaryFileOffset", bfoss.str());
            // write the library in text format to .sptxt
            entry->writeToFile(m_txtFout);
        }

    } else {
        // write to library in text format to .splib file
        entry->writeToFile(m_libFout);
    }

    // update m/z index object
    m_mzIndex->insertEntry(entry, offset);

    // update peptide index object
    m_pepIndex->insertEntry(entry, offset);
}


// retrieve - retrieves all library entries within a m/z tolerance of the target m/z, 
// and store them in the vector 'entries'. Basically calls SpectraSTLibIndex::retrieve
void SpectraSTLib::retrieve(vector<SpectraSTLibEntry*>& entries, double lowMz, double highMz) {

    // check to make sure we are in the Search mode
    if (!m_searchParams) {
        return;
    }

    retrieveSQL(entries, lowMz, highMz);
//    m_mzIndex->retrieve(entries, lowMz, highMz, true);
}

// writePreamble - writes some information about the library to the library file (.sptxt if binary library format is used, .splib otherwise)
void SpectraSTLib::writePreamble(vector<string>& lines) {

    // check to make sure we are in the Create mode
    if (!m_createParams) {
        return;
    }

    string fullLibFileName(m_libFileName);
    makeFullPath(fullLibFileName);

    FileName fn;
    parseFileName(fullLibFileName, fn);

    // hijack -- don't print detailed information
    // lines.clear();
    // END hijack

    if (m_createParams->binaryFormat) {
        fstream::off_type offset = 0;

        if (!m_noSptxt) {

            // writing to .sptxt file
            offset = m_txtFout.tellp();
            if (offset != 0) {
                // already something in the file, don't write preamble
                return;
            }
            // print header
            m_txtFout << "### " << fn.name + ".sptxt" << "  (Text version of " << fn.name + ".splib" << ")" << endl;
            m_txtFout << "### " << "SpectraST (version " << SPECTRAST_VERSION << '.' << SPECTRAST_SUB_VERSION << ", " << szTPPVersionInfo << ")" << endl;
            m_txtFout << "### " << endl;

            // print the rest of the preamble
            for (vector<string>::iterator i = lines.begin(); i != lines.end(); i++) {
                m_txtFout << "### " << (*i) << endl;
            }
            m_txtFout << "### ===" << endl;
        }

        // writing to .splib (binary) file
        offset = m_libFout.tellp();
        if (offset != 0) {
            return;
        }
        // print header
        int spectrastVersion = SPECTRAST_VERSION;
        int spectrastSubVersion = SPECTRAST_SUB_VERSION;
        m_libFout.write((char*)(&spectrastVersion), sizeof(int));
        m_libFout.write((char*)(&spectrastSubVersion), sizeof(int));
        m_libFout << fn.name + ".splib" << endl;
        // print the rest of the preamble
        unsigned int numLines = (unsigned int)lines.size();
        m_libFout.write((char*)(&numLines), sizeof(unsigned int));
        for (vector<string>::iterator j = lines.begin(); j != lines.end(); j++) {
            m_libFout << (*j) << endl;
        }


    } else {

        // writing to .splib file
        fstream::off_type offset = m_libFout.tellp();
        if (offset != 0) {
            return;
        }

        // print header
        m_libFout << "### " << fn.name + ".splib" << endl;
        m_libFout << "### " << "SpectraST (version " << SPECTRAST_VERSION << '.' << SPECTRAST_SUB_VERSION << ", " << szTPPVersionInfo << ")" << endl;
        m_libFout << "### " << endl;

        // print the rest of the preamble
        for (vector<string>::iterator i = lines.begin(); i != lines.end(); i++) {
            m_libFout << "### " << (*i) << endl;
        }
        m_libFout << "### ===" << endl;
    }

}

// extractDatabaseFileFromPreamble - try to guess what the sequence database should be from the library's preamble. This is 
// useful when the user did not specify the -sD option when searching. This routine will parse the searched library's preamble,
// figure out which sequence database is searched most often in the datasets used in building the library, and assume that 
// is the most appropriate sequence database to be fed to RefreshParser. 
void SpectraSTLib::extractDatabaseFileFromPreamble(bool binary) {

    map<string, int> allDatabaseFiles;

    if (!binary) {
        char firstChar = (char)(m_libFin.peek());
        if (firstChar != '#') return;

        string line;
        string::size_type pos = 0;
        while (nextLine(m_libFin, line, "### ===", "")) {
            string prefix = line.substr(0, 3);
            if (prefix != "###") {
                break;
            }
            string rest = nextToken(line, 3, pos, "\r\n", " \t");

            if (!rest.empty() && ((pos = rest.find(" against \"")) != string::npos)) {
                string databaseFile = nextToken(rest, pos + 10, pos, "\"\t\r\n", "\"\t\r\n");
                map<string, int>::iterator found = allDatabaseFiles.find(databaseFile);
                if (found != allDatabaseFiles.end()) {
                    found->second++;
                } else {
                    allDatabaseFiles[databaseFile] = 1;
                }

            }
        }


    } else {

        int spectrastVersion = 0;
        int spectrastSubVersion = 0;
        unsigned int numLines = 0;
        string line("");

        m_libFin.read((char*)(&spectrastVersion), sizeof(int));
        m_libFin.read((char*)(&spectrastSubVersion), sizeof(int));
        if (!nextLine(m_libFin, line)) {
            g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
            g_log->crash();
        }

        string fileName(line);

        m_libFin.read((char*)(&numLines), sizeof(unsigned int));

        for (unsigned int i = 0; i < numLines; i++) {

            if (!nextLine(m_libFin, line)) {
                g_log->error("GENERAL", "Corrupt .splib file from which to import entry.");
                g_log->crash();
            }

            string::size_type pos = 0;

            if (!line.empty() && ((pos = line.find(" against \"")) != string::npos)) {
                string databaseFile = nextToken(line, pos + 10, pos, "\"\t\r\n", "\"\t\r\n");
                map<string, int>::iterator found = allDatabaseFiles.find(databaseFile);
                if (found != allDatabaseFiles.end()) {
                    found->second++;
                } else {
                    allDatabaseFiles[databaseFile] = 1;
                }

            }

        }
    }

    int maxOcc = 0;
    string maxDatabaseFile("");
    for (map<string, int>::iterator d = allDatabaseFiles.begin(); d != allDatabaseFiles.end(); d++) {
        if (d->second > maxOcc) {
            maxOcc = d->second;
            maxDatabaseFile = d->first;
        }
    }
    if (!maxDatabaseFile.empty()) {
        m_searchParams->databaseFile = maxDatabaseFile;
    }

}

//================================================================================
// SQLite3 database method
//
// The workflow of SQLite C/C++ API has two modes:
//
// 1. EXEC: int sqlite3_exec() function. Execuate an SQL statement directly.
//
// 2. Prepare - (Bind values) - Step(s) - (Clear Bindings) [- Step(s) ] - Finalize
//
// In this implementation, the PREPARE method is put in constructor with DB initialization,
// and FINALIZE method is put in the destructor with DB shutting down.
//
// The CYCLE of [bind - step(s) - clear bindings] is execute for every miss in retrieveSQL().

void SpectraSTLib::initializeDatabase() {

    db = NULL;
    int rc = 0;

    sqlite3_initialize();

    sqlite3_open_v2(":memory:", &db, SQLITE_OPEN_READONLY, NULL); // in-memory mode

    string db_stmt = "ATTACH DATABASE '";
    db_stmt = db_stmt + m_libFileName + "' AS db";

    rc = sqlite3_exec(db, db_stmt.c_str(), NULL, NULL, NULL);

    if(rc != SQLITE_OK) {
        sqlite3_close(db);
        cerr << "Cannot open the given library." << endl;
        exit(1);
    }

    db_stmt = "PRAGMA cache_size = 128000"; /* change the default cache size (# of pages) of SQLite. 128000 PAGES = 128 MB */
    rc = sqlite3_exec(db, db_stmt.c_str(), NULL, NULL, NULL);

    stmt = NULL;

    sqlite3_prepare_v2(db, "SELECT * FROM db.Library WHERE Mz BETWEEN (:lowMz) AND (:highMz);",
                       -1, &stmt, NULL);

//    CREATE TABLE Library (PeptideName text, Mz real, Charge integer,
//    NumPeak integer, Mz1 real, Intensity1 real, Annotation1 text, ...);

    hit = 0;
    miss = 0;
}

void SpectraSTLib::retrieveSQL(vector<SpectraSTLibEntry *> &hits, double lowMz, double highMz)
{
    int idx_low = (int) (lowMz - MIN_MZ) / BLOCK_SIZE;
    int idx_high = (int) (highMz - MIN_MZ) / BLOCK_SIZE;

    for(int idx = idx_low; idx <= idx_high; ++idx)
    {
        if(m_cache.find(idx) != m_cache.end()) // update cache kick out queue
        {
            list<int>::iterator found = find(cache_queue.begin(), cache_queue.end(), idx);
            cache_queue.erase(found);
            cache_queue.push_back(idx);
            ++hit;
        }
        else // this query is a miss, go to sqlite and find it in database
        {
            cache_queue.push_back(idx);
            m_cache[idx];
            ++miss;

            if(m_cache.size() > CACHE_SIZE)
            {
                int kick_out = cache_queue.front();
                cache_queue.erase(cache_queue.begin());

                for(vector<Entry>::iterator iter = m_cache[kick_out].begin(); iter != m_cache[kick_out].end(); ++iter)
                    delete iter->second;

                m_cache.erase(kick_out);
            }

            int lowMz_idx = sqlite3_bind_parameter_index(stmt, ":lowMz");
            int highMz_idx = sqlite3_bind_parameter_index(stmt, ":highMz");

            sqlite3_bind_double(stmt, lowMz_idx, idx * BLOCK_SIZE + MIN_MZ);
            sqlite3_bind_double(stmt, highMz_idx, (idx + 1) * BLOCK_SIZE + MIN_MZ);

            while(sqlite3_step(stmt) == SQLITE_ROW)
            {
                const char* raw_name = (const char*) (sqlite3_column_text(stmt, 0));
                string name(raw_name);

                int parentCharge = sqlite3_column_int(stmt, 1);

                double parentMz = sqlite3_column_double(stmt, 2);

                unsigned int numPeaks = sqlite3_column_int(stmt, 3);

                SpectraSTPeakList* newPeaklist = new SpectraSTPeakList(parentMz,
                                                                       parentCharge,
                                                                       numPeaks);

                for(int i = 0; i != numPeaks; ++i)
                {
                    double mz = sqlite3_column_double(stmt, 3 * i + 4);
                    double intensity = sqlite3_column_double(stmt, 3 * i + 5);
                    const char* raw_anno = (const char*) sqlite3_column_text(stmt, 3 * i + 6);
                    string annotations(raw_anno);

                    newPeaklist->insertForSearch(mz, intensity, annotations);
                }

                SpectraSTLibEntry* newLibEntry = new SpectraSTLibEntry(name, parentMz, string(),
                                                                       string(), newPeaklist, string());

                int intMz = parentMz;
                m_cache[idx].push_back(Entry(intMz, newLibEntry));
            }

            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);
        }

        for(block::iterator iter = m_cache[idx].begin(); iter != m_cache[idx].end(); ++iter)
        {
            if(iter->first <= (int) (highMz) && iter->first >= (int) (lowMz))
                hits.push_back(iter->second);
        }
    }
}

void SpectraSTLib::resetCache()
{
    for(map<int, block>::iterator iter1 = m_cache.begin(); iter1 != m_cache.end(); ++iter1)
    {
        for(vector<Entry>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ++iter2)
            delete iter2->second;
    }

    m_cache.clear();
}

void SpectraSTLib::shutdownDatabase()
{
    sqlite3_finalize(stmt);

    sqlite3_close(db);
    sqlite3_shutdown();

    int block_size = BLOCK_SIZE;

    cout << "The num of block is " << m_cache.size() << endl;
    cout << "The block size is " << block_size << " (intMz)"<<endl;
    cout << "Total miss is " << miss << endl;
    cout << "Total hit is " << hit << endl;
    cout << "The hit rate is " << (double) hit / (hit + miss) << endl;
}
