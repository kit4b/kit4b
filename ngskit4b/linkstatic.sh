#!/bin/bash
# linking ngskit4b as a static image with very few library dependences - save grief with gcc verions not present on clusters etc
# ngskit4bs
g++  -static -g -O2   -o ngskit4bs KAlignerCL.o Blitz.o genDESeq.o kit4bax.o LocKMers.o mergeoverlaps.o Scaffolder.o SQLiteSummaries.o \
    KAligner.o csv2sqlite.o genhyperconserved.o ngskit4b.o MapLoci2Feat.o MergeReadPairs.o SimReads.o SSRdiscovery.o AlignsBootstrap.o \
    deNovoAssemb.o genkmarkers.o CallHaplotypes.o rnade.o maploci2features.o PEScaffold.o SNPs2pgSNPs.o ArtefactReduce.o fastaextract.o \
    genmarkerseq.o kit4bdna.o MarkerKMers.o psl2sqlite.o SQLiteDE.o AssembGraph.o FastaNxx.o genpseudogenome.o kmermarkers.o Markers.o \
    ReadStats.o SQLiteMarkers.o Assemble.o FilterSAMAlignments.o gensnpmarkers.o LocateROI.o MarkerSeq.o RemapLoci.o Benchmarker.o \
    genbioseq.o genbiobed.o gengoassoc.o gengoterms.o goassoc.o seghaplotypes.o hammings.o fasta2bed.o LocHap2Bed.o pangenome.o GBSmapSNPs.o \
    repassemb.o sarscov2ml.o haps2genotype.o \
    ../libkit4b/libkit4b.a ../libzlib/libzlib.a ../libBKPLPlot/libBKPLPlot.a -lrt -ldl -lpthread
cp ngskit4bs ~/bin

