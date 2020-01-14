# kit4b 
kit4b - the 'K-mer Informed Toolkit for Bioinformatics' - is an integrated toolkit of high performance bioinformatics subprocesses targeting the challenges of next generation sequencing analytics. 

This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in incompatibilty with 'BioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga' parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.
'kit4b' is being released under the Opensource Software License Agreement (GPLv3)

## Why YAT (Yet Another Toolkit)
Compared with other widely used bioinformatics toolkits providing alignment capabilities, kit4b provides substantial gains in both the proportion and quality of aligned sequence reads at competitive or increased computational efficiency. Unlike most other aligners, bkit and BioKanga utilises Hamming distances between putative alignments to the targeted genome assembly for any given read as the discrimative acceptance criteria rather than relying on sequencer generated quality scores.

Another primary differentiator for kit4b is that this toolkit can process billions of reads against targeted genomes containing 100 million contigs and totalling up to 100Gbp of sequence.

## Toolset Components
The kit4b toolset contains a number of subprocesses, each of which is targeting a specific bioinformatics analytics task. Primary subprocesses provide functionality for:
  simreads              Generate simulated NGS readsets
  ngsqc                 Process NGS reads and report quality scores with compositional distributions
  fasta2nxx             Generate N10..N90 over Fasta sequences
  filter                Filter NGS reads for sequencer errors and/or exact duplicates
  assemb                de Novo assemble filtered reads into contigs
  scaffold              Scaffold de Novo assembled contigs
  index                 Generate index over genome assembly or sequences
  kmarkers              NGS reads alignment-less K-mer derived marker sequences generation
  prekmarkers           NGS reads alignment-less prefix K-mer derived marker sequences generation
  pseudogenome          Concatenate sequences to create pseudo-genome assembly
  kalign                Align (under development) NGS reads to indexed genome assembly or sequences
  pescaffold            Scaffold assembly contigs using PE read alignments
  ssr                   Identify SSRs in multifasta sequences
  maploci               Map aligned reads loci to known features
  rnade                 RNA-seq differential expression analyser with optional Pearsons generation
  gendeseq              Generate tab delimited counts file for input to DESeq or EdgeR
  xfasta                Extract fasta sequences from multifasta file
  mergeoverlaps         Merge PE short insert overlap reads
  snpmarkers            SNP alignment derived marker sequences identification
  markerseqs            Generate marker sequences from SNP loci
  blitz                 Blat like local align genomic sequences
  remaploci             Remap alignment loci
  filtchrom             Filter SAM/BAM alignments by chromosome
  locateroi             Locate and report regions of interest
  alignsbs              Alignments bootstrapper
  ultras                Identify Utra/Hyper conserved elements
  psl2sqlite            Generate SQLite Blat alignment Database from Blat generated PSL alignments
  snpm2sqlite           Generate SQLite Marker Database from SNP markers
  snps2sqlite           Generate SQLite SNP Database from aligner identified SNPs
  de2sqlite             Generate SQLite DE Database from RNA-seq DE
  genbioseq             Generate bioseq pre-parsed sequence file
  goassoc               GO association inferencing
  gengoassoc            Generate biogoassoc pre-indexed GO associations
  gengoterms            Generate biogoterms pre-indexed GO terms
  hammings              Generate hamming distances for K-mer over sequences
  fasta2bed             Generate BED file from fasta containing sequence names and lengths

There are many other sperate processes additional to kit4b included in this release, all sharing the same common library modules which are utilised by kit4b. The most notable of these is the kit4bpacbio process which has subprocesses targeting denovo PacBio assemblies. Documentation on these additional processes will be provided at a later date when the primary source 'bkit' has matured.

## Build and installation
### Linux
To build on Ubuntu, clone this repository, run `autoreconf -f -i`, `configure` and `make`. The following example will install the biokanga toolkit to a `bin` directory underneath the user's home directory.

```
mkdir ~/bin
git clone https://github.com/kit4b/kit4b.git
cd kit4b
autoreconf -f -i
./configure --prefix=$HOME
make install
```

Alternatively, the binary built for the appropriate platform can be used directly.

### Windows
To build on Windows, the current version requires Visual Studio 2017 or 2019 **with build tools of at least v140**. 
1. Open the `kit4b.sln` file in Visual Studio. 
2. Under the Build menu, select Configuration Manager. 
3. For Active solution platform, select x64. 
4. The project can then be built. By default, executables will be copied into the `bin64` directory.


## Documentation
Documentation for the core functionality of kit4b and kit4pacbio is available under the `Docs` directory.

## Contributing
kit4b is maintained by Dr Stuart Stephen, formerly with the - BioKanga - Crop Bioinformatics and Data Science team at CSIRO based in Canberra, Australia. Although nominally now retired Dr Stephen is still actively engaged in bioinformatics related research. 

Contributions are most welcome. To contribute, follow these steps.

1. Fork kit4b into your own repository ([more information](https://help.github.com/articles/about-forks/))
2. Clone and enter the repository to your development machine
3. Checkout the `dev` branch
4. Make and checkout a new branch for your work (`git checkout -b great-new-feature`)
5. Make regular commits on your new branch
6. Push your branch back to your github repository (`git push origin great-new-feature`)
7. Create a pull request to the `dev` branch of the github kit4b repository ([more information](https://help.github.com/articles/creating-a-pull-request/))
8. If your contribution is related to a known existing issue, refer to the issue in the pull request comment


## Issues
Please report issues on the [github project](https://github.com/kit4b/kit4b/issues).

