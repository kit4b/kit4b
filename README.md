# kit4b 
**kit4b** - the 'K-mer Informed Toolkit for Bioinformatics' - is an integrated toolkit of high performance bioinformatics processes and subprocesses targeting the challenges of next generation sequencing analytics. 

This toolkit is a source base clone of **BioKanga** release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) containing significant functionality enhancements breaking backward compatibility with **BioKanga**.

**kit4b** is being released under the Opensource Software License Agreement (GPLv3) as is **BioKanga**

## Why YAT (Yet Another Toolkit)
Compared with other widely used bioinformatics toolkits, **kit4b** provides substantial gains in throughput and quality of processing with integrated functionality providing the majority of bioinformatics workflow processing stages without the need to incorporate third party processing into the workflow.

## Toolset Components
The **kit4b** toolset contains a number of processes and subprocesses, each of which is targeting a specific bioinformatics analytics task.

The primary process is **kitb** which has subprocesses providing functionality for:
  * simreads -            Generate simulated NGS readsets
  * ngsqc    -            Process NGS reads and report quality scores with compositional distributions
  * fasta2nxx             Generate N10..N90 over Fasta sequences
  * filter                Filter NGS reads for sequencer errors and/or exact duplicates
  * assemb                de Novo assemble filtered reads into contigs
  * scaffold              Scaffold de Novo assembled contigs
  * index                 Generate index over genome assembly or sequences
  * kmarkers              NGS reads alignment-less K-mer derived marker sequences generation
  * prekmarkers           NGS reads alignment-less prefix K-mer derived marker sequences generation
  * pseudogenome          Concatenate sequences to create pseudo-genome assembly
  * kalign                Align (under development) NGS reads to indexed genome assembly or sequences
  * pescaffold            Scaffold assembly contigs using PE read alignments
  * ssr                   Identify SSRs in multifasta sequences
  * maploci               Map aligned reads loci to known features
  * rnade                 RNA-seq differential expression analyser with optional Pearsons generation
  * gendeseq              Generate tab delimited counts file for input to DESeq or EdgeR
  * xfasta                Extract fasta sequences from multifasta file
  * mergeoverlaps         Merge PE short insert overlap reads
  * snpmarkers            SNP alignment derived marker sequences identification
  * markerseqs            Generate marker sequences from SNP loci
  * blitz                 Blat like local align genomic sequences
  * remaploci             Remap alignment loci
  * filtchrom             Filter SAM/BAM alignments by chromosome
  * locateroi             Locate and report regions of interest
  * alignsbs              Alignments bootstrapper
  * ultras                Identify Utra/Hyper conserved elements
  * psl2sqlite            Generate SQLite Blat alignment Database from Blat generated PSL alignments
  * snpm2sqlite           Generate SQLite Marker Database from SNP markers
  * snps2sqlite           Generate SQLite SNP Database from aligner identified SNPs
  * de2sqlite             Generate SQLite DE Database from RNA-seq DE
  * genbioseq             Generate bioseq pre-parsed sequence file
  * goassoc               GO association inferencing
  * gengoassoc            Generate biogoassoc pre-indexed GO associations
  * gengoterms            Generate biogoterms pre-indexed GO terms
  * hammings              Generate hamming distances for K-mer over sequences
  * fasta2bed             Generate BED file from fasta containing sequence names and lengths

There are many other sperate processes additional to **kit4b** included in this release, all sharing the same common library modules in **libkitb**. 

The most notable of these additional processes is **kit4bpacbio** which has subprocesses targeting PacBio long read error correction and denovo assemblies. 
Documentation on these additional processes will be provided at a later date when the primary process **bkit** has matured with complete documentation having been made available.

## Build and installation
### Linux
To build on Ubuntu, clone this repository, run `autoreconf -f -i`, `configure` and `make`. 

The following example will install the **kitb** toolkit processes to a `bin` directory underneath the user's home directory.

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
Incomplete documentation for the core functionality of **kit4b** is available under the `Docs` directory.

## Contributing
kit4b is maintained by Dr Stuart Stephen, formerly with the - **BioKanga** - Crop Bioinformatics and Data Science team at CSIRO based in Canberra, Australia. Although nominally now retired Dr Stephen is still actively engaged in bioinformatics related research. 

Contributions are most welcome. To contribute, follow these steps.

1. Fork **kit4b** into your own repository ([more information](https://help.github.com/articles/about-forks/))
2. Clone and enter the repository to your development machine
3. Checkout the `dev` branch
4. Make and checkout a new branch for your work (`git checkout -b great-new-feature`)
5. Make regular commits on your new branch
6. Push your branch back to your github repository (`git push origin great-new-feature`)
7. Create a pull request to the `dev` branch of the github kit4b repository ([more information](https://help.github.com/articles/creating-a-pull-request/))
8. If your contribution is related to a known existing issue, refer to the issue in the pull request comment


## Issues
Please report any issues to [github project](https://github.com/kit4b/kit4b/issues).
