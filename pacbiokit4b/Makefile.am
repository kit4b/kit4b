bin_PROGRAMS = pacbiokit4b
pacbiokit4b_SOURCES= SQLiteSummaries.cpp SQLiteSummaries.h SSW.cpp SSW.h SWAlign.cpp SWAlign.h PBAssemb.cpp PBAssemb.h PBECContigs.cpp PBECContigs.h \
                     SeqStore.cpp SeqStore.h PBFilter.cpp PBFilter.h pacbiocommon.h PacBioUtility.cpp PacBioUtility.h pacbiokit4b.cpp pacbiokit4b.h \
                     PBErrCorrect.cpp PBErrCorrect.h MAConsensus.cpp MAConsensus.h AssembGraph.cpp AssembGraph.h \
                     MAFKMerDist.cpp MAFKMerDist.h PBSWService.cpp PBSWService.h BKSProvider.cpp BKSProvider.h BKSRequester.cpp BKSRequester.h BKScommon.h

# set the include path found by configure
AM_CPPFLAGS= $(all_includes)

# our library
LDADD = ../libkit4b/libkit4b.a ../libzlib/libzlib.a ../libBKPLPlot/libBKPLPlot.a

