###################################################################################
# Program       : Makefile                                                        #
# Author        : Henry Lam <hlam@systemsbiology.org>                      #
# Date          : 01/18/2003                                                      #
#                                                                                 #
# THE SOFTWARE IS PROVIDED BY THE INSTITUTE FOR SYSTEMS BIOLOGY (ISB)             #
# "AS IS" AND "WITH ALL FAULTS." ISB MAKES NO REPRESENTATIONS OR WARRANTI         #
# ES OF ANY KIND CONCERNING THE QUALITY, SAFETY OR SUITABILITY OF THE             #
# SOFTWARE, EITHER EXPRESSED OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IM-     #
# PLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,          #
# OR NON-INFRINGEMENT.                                                            #
#                                                                                 #
# ISB MAKES NO REPRESENTATIONS OR WARRANTIES AS TO THE TRUTH, ACCURACY            #
# OR COMPLETENESS OF ANY STATEMENTS, INFORMATION OR MATERIALS CONCERNING          #
#                                                                                 #
# THE SOFTWARE THAT IS CONTAINED IN ANY DOCUMENTATION INCLUDED WITH THE           #
# SOFTWARE OR ON AND WITHIN ANY OF THE WEBSITES OWNED AND OPERATED BY ISB         #
#                                                                                 #
# IN NO EVENT WILL ISB BE LIABLE FOR ANY INDIRECT, PUNITIVE, SPECIAL,             #
# INCIDENTAL OR CONSEQUENTIAL DAMAGES HOWEVER THEY MAY ARISE AND EVEN IF          #
# ISB HAVE BEEN PREVIOUSLY ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.            #
#                                                                                 #
###################################################################################

SYSLIBS= -lm -lz

LDFLAGS= $(SYSLIBS)

ARCH=linux_standalone

EXE=${ARCH}/spectrast #${ARCH}/plotspectrast ${ARCH}/plotspectrast.cgi ${ARCH}/Lib2HTML

# lfs support
LFSFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
WARNINGFLAGS = -Werror -Wformat -Wstrict-aliasing -Wno-deprecated -Wno-char-subscripts
IFLAGS=

ifeq (${LGPL},1)
CXX=g++ ${DEBUG} ${IFLAGS} ${WARNINGFLAGS} ${LFSFLAGS} ${OSFLAGS} -DSTANDALONE_LINUX -D__LGPL__
else
CXX=g++ ${DEBUG} ${IFLAGS} ${WARNINGFLAGS} ${LFSFLAGS} ${OSFLAGS} -DSTANDALONE_LINUX
endif

######################################################################
#
# do not edit past this point
#
######################################################################

MYCRAMP = ${ARCH}/SpectraST_cramp.o ${ARCH}/SpectraST_ramp.o ${ARCH}/SpectraST_base64.o ${ARCH}/SpectraST_util.o
#MYCRAMP = $(TPPLIB) 
ifeq (${LGPL},1)
# no kwset, no REFRESH function 
else
MYKWSET = ${ARCH}/SpectraST_kwset.o ${ARCH}/SpectraST_obstack.o
endif

OBJS= ${ARCH}/SpectraSTLib.o ${ARCH}/SpectraSTLibIndex.o \
	${ARCH}/SpectraSTLibEntry.o ${ARCH}/SpectraSTPeakList.o \
	${ARCH}/SpectraSTMzLibIndex.o ${ARCH}/SpectraSTCreateParams.o \
	${ARCH}/SpectraSTLibImporter.o ${ARCH}/SpectraSTMspLibImporter.o \
	${ARCH}/SpectraSTPepXMLLibImporter.o ${ARCH}/SpectraSTSpLibImporter.o \
	${ARCH}/SpectraSTXHunterLibImporter.o ${ARCH}/SpectraSTTsvLibImporter.o \
	${ARCH}/SpectraSTSearchTask.o \
	${ARCH}/SpectraSTDtaSearchTask.o ${ARCH}/SpectraSTMspSearchTask.o ${ARCH}/SpectraSTMgfSearchTask.o \
	${ARCH}/SpectraSTMzXMLSearchTask.o ${ARCH}/SpectraSTCandidate.o\
	${ARCH}/SpectraSTSearch.o ${ARCH}/SpectraSTReplicates.o \
	${ARCH}/SpectraSTSearchOutput.o ${ARCH}/SpectraSTTxtSearchOutput.o\
	${ARCH}/SpectraSTXlsSearchOutput.o ${ARCH}/SpectraSTPepXMLSearchOutput.o \
	${ARCH}/SpectraSTHtmlSearchOutput.o ${ARCH}/SpectraSTPeptideLibIndex.o ${ARCH}/SpectraSTSimScores.o \
	${ARCH}/SpectraSTSearchParams.o ${ARCH}/SpectraSTMain.o \
	${ARCH}/SpectraSTQuery.o ${ARCH}/SpectraSTFileList.o \
	${ARCH}/SpectraSTLog.o ${ARCH}/SpectraSTMs2LibImporter.o \
	${ARCH}/SpectraSTMzXMLLibImporter.o \
	${ARCH}/SpectraSTSearchTaskStats.o \
	${ARCH}/SpectraSTFastaFileHandler.o \
	${ARCH}/SpectraSTDenoiser.o \
	${ARCH}/FileUtils.o ${ARCH}/Peptide.o ${ARCH}/XMLWalker.o \
	${ARCH}/ProgressCount.o ${ARCH}/Predicate.o \
	${MYCRAMP} ${MYKWSET}


#
# rules
#
.SUFFIXES:	.o .cpp 

#
# note we have our own arch subdirs, needed
# to avoid filename conflicts in the OBJ_ARCH
# directory of the parent project
#

${ARCH}/%.o : %.cpp
	@ mkdir -p $(ARCH)
	$(CXX) $(OLD) $(DEBUG) $(IFLAGS) $(LFSFLAGS) ${WARNINGFLAGS} -O2 -c $< -o $@

#
# targets
#

all : ${EXE}

${ARCH}/spectrast : ${OBJS}
	$(CXX) $(DEBUG) -O2 $^ $(LDFLAGS) -o $@ 

${ARCH}/plotspectrast.cgi : ${ARCH}/plotspectrast_cgi.o ${ARCH}/SpectraSTLibEntry.o ${ARCH}/SpectraSTQuery.o ${ARCH}/SpectraSTPeakList.o ${ARCH}/SpectraSTDenoiser.o ${ARCH}/SpectraSTLog.o ${ARCH}/FileUtils.o ${ARCH}/Peptide.o ${MYCRAMP}
	$(CXX) -O2 $^ $(LDFLAGS) -o $@ 

${ARCH}/plotspectrast : ${ARCH}/plotspectrast.o ${ARCH}/SpectraSTLibEntry.o ${ARCH}/SpectraSTQuery.o ${ARCH}/SpectraSTPeakList.o ${ARCH}/SpectraSTDenoiser.o ${ARCH}/SpectraSTLog.o ${ARCH}/FileUtils.o ${ARCH}/Peptide.o ${MYCRAMP}
	$(CXX) -O2 $^ $(LDFLAGS) -o $@ 

${ARCH}/Lib2HTML : ${ARCH}/Lib2HTML.o ${ARCH}/SpectraSTPeptideLibIndex.o ${ARCH}/SpectraSTMzLibIndex.o ${ARCH}/SpectraSTLibIndex.o ${ARCH}/SpectraSTLibEntry.o ${ARCH}/SpectraSTPeakList.o ${ARCH}/SpectraSTDenoiser.o ${ARCH}/SpectraSTLog.o ${ARCH}/FileUtils.o ${ARCH}/Peptide.o ${MYCRAMP}
	$(CXX) -O2 $^ $(LDFLAGS) -o $@

${ARCH}/plotspectrast_cgi.o : plotspectrast.cpp 
	$(CXX) $(LFSFLAGS) $(IFLAGS) ${WARNINGFLAGS} -D RUN_AS_CGI -O2 -o ${ARCH}/plotspectrast_cgi.o -c plotspectrast.cpp 

clean:
	rm -rf ${ARCH}; rm -f $(EXE) core* *~

debug:
	DEBUG='-g'  make all

LGPL:
	LGPL='1' make -f Makefile_STANDALONE_LINUX clean all



#
# dependencies
#
${ARCH}/SpectraSTLib.o : SpectraSTLib.cpp  SpectraSTLib.hpp  SpectraSTLibIndex.hpp  SpectraSTPeptideLibIndex.hpp  SpectraSTLibImporter.hpp  FileUtils.hpp Peptide.hpp
${ARCH}/SpectraSTLibIndex.o : SpectraSTLibIndex.cpp  SpectraSTLibIndex.hpp  SpectraSTLibEntry.hpp  FileUtils.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTMzLibIndex.o : SpectraSTMzLibIndex.cpp  SpectraSTMzLibIndex.hpp  SpectraSTLibIndex.hpp SpectraSTLibEntry.hpp  FileUtils.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTLibEntry.o : SpectraSTLibEntry.cpp  SpectraSTLibEntry.hpp  SpectraSTPeakList.hpp  FileUtils.hpp Peptide.hpp  
${ARCH}/SpectraSTPeakList.o : SpectraSTPeakList.cpp  SpectraSTPeakList.hpp  SpectraSTSimScores.hpp  SpectraSTSearchParams.hpp  SpectraSTDenoiser.hpp FileUtils.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTLibImporter.o : SpectraSTLibImporter.cpp SpectraSTLibImporter.hpp  SpectraSTLib.hpp  SpectraSTMspLibImporter.hpp
${ARCH}/SpectraSTMspLibImporter.o : SpectraSTMspLibImporter.cpp SpectraSTMspLibImporter.hpp  SpectraSTLibImporter.hpp FileUtils.hpp Peptide.hpp
${ARCH}/SpectraSTPepXMLLibImporter.o : SpectraSTPepXMLLibImporter.cpp SpectraSTPepXMLLibImporter.hpp  SpectraSTLibImporter.hpp XMLWalker.hpp FileUtils.hpp Peptide.hpp
${ARCH}/SpectraSTSpLibImporter.o : SpectraSTSpLibImporter.cpp SpectraSTSpLibImporter.hpp  SpectraSTLibImporter.hpp XMLWalker.hpp  FileUtils.hpp Peptide.hpp ProgressCount.hpp SpectraSTLog.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTMs2LibImporter.o : SpectraSTMs2LibImporter.cpp SpectraSTMs2LibImporter.hpp  SpectraSTLibImporter.hpp FileUtils.hpp Peptide.hpp ProgressCount.hpp SpectraSTLog.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTXHunterLibImporter.o : SpectraSTXHunterLibImporter.cpp SpectraSTXHunterLibImporter.hpp  SpectraSTLibImporter.hpp FileUtils.hpp Peptide.hpp ProgressCount.hpp
${ARCH}/SpectraSTTsvLibImporter.o : SpectraSTTsvLibImporter.cpp SpectraSTTsvLibImporter.hpp  SpectraSTLibImporter.hpp FileUtils.hpp Peptide.hpp ProgressCount.hpp
${ARCH}/SpectraSTMzXMLLibImporter.o : SpectraSTMzXMLLibImporter.cpp SpectraSTMzXMLLibImporter.hpp  SpectraSTLibImporter.hpp SpectraSTReplicates.hpp FileUtils.hpp ProgressCount.hpp SpectraST_cramp.hpp SpectraST_ramp.h
${ARCH}/SpectraSTFastaFileHandler.o : SpectraSTFastaFileHandler.cpp SpectraSTFastaFileHandler.hpp SpectraSTLog.hpp FileUtils.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTSearchTask.o : SpectraSTSearchTask.cpp  SpectraSTSearchTask.hpp SpectraSTLib.hpp SpectraSTSearchOutput.hpp SpectraSTMzXMLSearchTask.hpp SpectraSTMspSearchTask.hpp SpectraSTDtaSearchTask.hpp SpectraSTMgfSearchTask.hpp SpectraSTSearchTaskStats.hpp
${ARCH}/SpectraSTSearchTaskStats.o : SpectraSTSearchTaskStats.cpp SpectraSTSearchTaskStats.hpp SpectraSTSearch.hpp SpectraSTMspSearchTask.hpp  SpectraSTMzXMLSearchTask.hpp FileUtils.hpp
${ARCH}/SpectraSTMspSearchTask.o : SpectraSTMspSearchTask.cpp SpectraSTMspSearchTask.hpp SpectraSTSearch.hpp SpectraSTPeakList.hpp FileUtils.hpp
${ARCH}/SpectraSTMzXMLSearchTask.o : SpectraSTMzXMLSearchTask.cpp SpectraSTMzXMLSearchTask.hpp SpectraSTSearch.hpp SpectraSTPeakList.hpp FileUtils.hpp SpectraST_cramp.hpp SpectraST_ramp.h
${ARCH}/SpectraSTDtaSearchTask.o : SpectraSTDtaSearchTask.cpp SpectraSTDtaSearchTask.hpp SpectraSTSearch.hpp SpectraSTPeakList.hpp FileUtils.hpp
${ARCH}/SpectraSTMgfSearchTask.o : SpectraSTMgfSearchTask.cpp SpectraSTMgfSearchTask.hpp SpectraSTSearch.hpp SpectraSTPeakList.hpp FileUtils.hpp
${ARCH}/SpectraSTCandidate.o : SpectraSTCandidate.cpp  SpectraSTCandidate.hpp  SpectraSTLibEntry.hpp  FileUtils.hpp
${ARCH}/SpectraSTQuery.o : SpectraSTQuery.cpp SpectraSTQuery.hpp SpectraSTPeakList.hpp
${ARCH}/SpectraSTFileList.o : SpectraSTFileList.cpp SpectraSTFileList.hpp SpectraSTCreateParams.hpp SpectraSTSearchParams.hpp SpectraSTLib.hpp SpectraSTSearchTask.hpp FileUtils.hpp
${ARCH}/SpectraSTSearch.o : SpectraSTSearch.cpp  SpectraSTSearch.hpp  SpectraSTLib.hpp SpectraSTCandidate.hpp  SpectraSTSearchOutput.hpp  SpectraSTLibEntry.hpp  FileUtils.hpp
${ARCH}/SpectraSTReplicates.o : SpectraSTReplicates.cpp  SpectraSTReplicates.hpp  SpectraSTLibEntry.hpp SpectraSTPeakList.hpp SpectraSTCreateParams.hpp SpectraSTDenoiser.hpp Peptide.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTSearchOutput.o : SpectraSTSearchOutput.cpp SpectraSTSearchOutput.hpp  SpectraSTSearchParams.hpp  SpectraSTSimScores.hpp  SpectraSTTxtSearchOutput.hpp  SpectraSTXlsSearchOutput.hpp  SpectraSTPepXMLSearchOutput.hpp  FileUtils.hpp
${ARCH}/SpectraSTTxtSearchOutput.o : SpectraSTTxtSearchOutput.cpp SpectraSTTxtSearchOutput.hpp  SpectraSTSearchOutput.hpp Peptide.hpp FileUtils.hpp
${ARCH}/SpectraSTXlsSearchOutput.o : SpectraSTXlsSearchOutput.cpp SpectraSTXlsSearchOutput.hpp  SpectraSTSearchOutput.hpp Peptide.hpp FileUtils.hpp
${ARCH}/SpectraSTPepXMLSearchOutput.o : SpectraSTPepXMLSearchOutput.cpp SpectraSTPepXMLSearchOutput.hpp  SpectraSTSearchOutput.hpp Peptide.hpp  SpectraSTConstants.hpp
${ARCH}/SpectraSTHtmlSearchOutput.o : SpectraSTHtmlSearchOutput.cpp SpectraSTHtmlSearchOutput.hpp  SpectraSTSearchOutput.hpp Peptide.hpp  SpectraSTConstants.hpp
${ARCH}/SpectraSTPeptideLibIndex.o : SpectraSTPeptideLibIndex.cpp SpectraSTPeptideLibIndex.hpp  SpectraSTLibIndex.hpp SpectraSTLibEntry.hpp  FileUtils.hpp  Peptide.hpp
${ARCH}/SpectraSTSimScores.o : SpectraSTSimScores.cpp SpectraSTSimScores.hpp
${ARCH}/SpectraSTSearchParams.o : SpectraSTSearchParams.cpp  SpectraSTSearchParams.hpp  FileUtils.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTCreateParams.o : SpectraSTCreateParams.cpp  SpectraSTCreateParams.hpp  FileUtils.hpp SpectraSTConstants.hpp
${ARCH}/SpectraSTLog.o : SpectraSTLog.cpp SpectraSTLog.hpp 
${ARCH}/SpectraSTDenoiser.o : SpectraSTDenoiser.cpp SpectraSTDenoiser.hpp SpectraSTPeakList.hpp
${ARCH}/SpectraSTMain.o : SpectraSTMain.cpp SpectraSTLib.hpp  SpectraSTLibEntry.hpp  SpectraSTLibIndex.hpp  SpectraSTSearchTask.hpp SpectraSTSearchParams.hpp SpectraSTLog.hpp SpectraSTConstants.hpp
${ARCH}/FileUtils.o : FileUtils.cpp FileUtils.hpp
${ARCH}/Peptide.o : Peptide.cpp Peptide.hpp
${ARCH}/XMLWalker.o : XMLWalker.cpp XMLWalker.hpp
${ARCH}/Predicate.o : Predicate.cpp Predicate.hpp
${ARCH}/ProgressCount.o : ProgressCount.cpp ProgressCount.hpp
${ARCH}/plotspectrast.o : plotspectrast.cpp SpectraSTLibEntry.hpp SpectraSTPeakList.hpp SpectraSTLog.hpp SpectraST_cramp.hpp SpectraST_constants.h
${ARCH}/plotspectrast_cgi.o : plotspectrast.cpp SpectraSTLibEntry.hpp SpectraSTPeakList.hpp SpectraSTLog.hpp SpectraST_cramp.hpp SpectraST_constants.h
${ARCH}/Lib2HTML.o : Lib2HTML.cpp SpectraSTLibEntry.hpp SpectraSTPeptideLibIndex.hpp SpectraSTLog.hpp Peptide.hpp FileUtils.hpp
${ARCH}/SpectraST_base64.o : SpectraST_base64.cpp SpectraST_base64.h
${ARCH}/SpectraST_ramp.o :  SpectraST_ramp.cpp SpectraST_ramp.h
${ARCH}/SpectraST_cramp.o : SpectraST_cramp.cpp SpectraST_cramp.hpp
${ARCH}/SpectraST_util.o : SpectraST_util.cpp SpectraST_util.h
