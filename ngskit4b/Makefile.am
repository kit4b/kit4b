bin_PROGRAMS = ngskit4b
ngskit4b_SOURCES= KAlignerCL.cpp Blitz.cpp genDESeq.cpp kit4bax.cpp LocKMers.cpp mergeoverlaps.cpp Scaffolder.cpp SQLiteSummaries.cpp \
		KAligner.cpp csv2sqlite.cpp genhyperconserved.cpp ngskit4b.cpp MapLoci2Feat.cpp MergeReadPairs.cpp SimReads.cpp \
		SSRdiscovery.cpp AlignsBootstrap.cpp deNovoAssemb.cpp genkmarkers.cpp rnade.cpp maploci2features.cpp PEScaffold.cpp \
		SNPs2pgSNPs.cpp SNPs2pgSNPs.h ArtefactReduce.cpp fastaextract.cpp genmarkerseq.cpp kit4bdna.cpp MarkerKMers.cpp psl2sqlite.cpp \
		SQLiteDE.cpp AssembGraph.cpp FastaNxx.cpp genpseudogenome.cpp kmermarkers.cpp Markers.cpp ReadStats.cpp SQLiteMarkers.cpp \
		Assemble.cpp FilterSAMAlignments.cpp gensnpmarkers.cpp LocateROI.cpp MarkerSeq.cpp RemapLoci.cpp \
		Aligner.h KAligner.h AssembGraph.h FastaNxx.h kit4bdna.h LocKMers.h MarkerSeq.h PEScaffold.h Scaffolder.h SQLiteDE.h SQLiteSummaries.h \
		AlignsBootstrap.h Blitz.h FilterSAMAlignments.h ngskit4b.h MapLoci2Feat.h Markers.h ReadStats.h SimReads.h SQLiteMarkers.h \
		SSRdiscovery.h ArtefactReduce.h deNovoAssemb.h genhyperconserved.h LocateROI.h MarkerKMers.h MergeReadPairs.h RemapLoci.h \
		Benchmarker.cpp Benchmarker.h SQLitePSL.h genbioseq.cpp genbiobed.cpp gengoassoc.cpp gengoterms.cpp goassoc.cpp seghaplotypes.cpp seghaplotypes.h \
		CallHaplotypes.cpp CallHaplotypes.h pbautils.cpp pbautils.h CGenMLdatasets.cpp CGenMLdatasets.h CWIGutils.cpp CWIGutils.h xroiseqs.cpp rnaexpr.cpp rnaexpr.h \
		hammings.cpp fasta2bed.cpp GBSmapSNPs.cpp GBSmapSNPs.h LocHap2Bed.cpp LocHap2Bed.h pangenome.cpp pangenome.h repassemb.cpp repassemb.h sarscov2ml.cpp sarscov2ml.h

# set the include path found by configure
AM_CPPFLAGS= $(all_includes)

AM_CFLAGS= -g
# our library
LDADD = ../libkit4b/libkit4b.a ../libzlib/libzlib.a ../libBKPLPlot/libBKPLPlot.a

