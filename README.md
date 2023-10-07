# kit4b 1.98.0
**kit4b** - the 'K-mer Informed Toolkit for Bioinformatics' - is an integrated toolkit of high performance bioinformatics processes and subprocesses targeting the challenges of next generation sequencing analytics. 

## Why YAT (Yet Another Toolkit)
Compared with other widely used bioinformatics toolkits, **kit4b** provides substantial gains in throughput and quality of processing with integrated functionality providing the majority of bioinformatics workflow processing stages without the need to incorporate third party processing into the workflow.

## Toolset Components
The **kit4b** toolset contains a number of processes and subprocesses, each of which is targeting a specific bioinformatics analytics task. **kit4b** core functionality is written in C++, and can be natively compiled for execution in either a Windows and Linux hosting environment. 

The primary process is **ngskitb** targeting NGS workflows which has subprocesses providing functionality to:
```
  * Generate UCSC Personal Genome SNP format output ready for direct loading into the UCSC Genome Browser
  * Provide aligner benchmarking
  * Generate simulated NGS readsets
  * Process NGS reads and report quality scores with compositional distributions
  * Generate N10..N90 over Fasta sequences
  * Filter NGS reads for sequencer errors and/or exact duplicates
  * de Novo assemble NGS filtered reads into contigs
  * Scaffold de Novo assembled contigs
  * Generate suffix array index over genome assembly or sequences
  * NGS reads alignment-less K-mer derived marker sequences generation
  * NGS reads alignment-less prefix K-mer derived marker sequences generation
  * Concatenate sequences to create pseudo-genome assembly
  * Align NGS reads to indexed genome assembly or sequences
  * Scaffold assembly contigs using PE read alignments
  * Identify SSRs in multifasta sequences
  * Map aligned reads genomic loci to known features
  * RNA-seq differential expression analyser with optional Pearsons generation
  * Generate tab delimited counts file for input to DESeq or EdgeR
  * Extract fasta sequences from multifasta file
  * Merge PE short insert overlap reads
  * SNP alignment derived marker sequences identification
  * Generate marker sequences from SNP loci
  * Blat like local align genomic sequences
  * Remap alignment loci
  * Filter SAM/BAM alignments by chromosome
  * Locate and report regions of interest
  * Alignments bootstrapper
  * Identify Utra/Hyper conserved elements
  * Generate SQLite Blat alignment Database from Blat generated PSL alignments
  * Generate SQLite Marker Database from SNP markers
  * Generate SQLite SNP Database from aligner identified SNPs
  * Generate SQLite DE Database from RNA-seq DE
  * Generate bioseq pre-parsed sequence file
  * GO association inferencing
  * Generate biogoassoc pre-indexed GO associations
  * Generate biogoterms pre-indexed GO terms
  * Generate hamming distances for K-mers over targeted sequences
  * Generate BED file from multifasta file containing sequence names and lengths
```
There are many other processes additional to **kit4b** included in this release, all sharing the same common library modules in **libkitb**. 

The most notable of these additional processes is **pacbiokit4b** which has subprocesses targeting PacBio long read error correction and denovo assemblies. 
Documentation on these additional processes will be provided at a later date when the primary process **ngsbkit** has matured with complete documentation having been made available.

## Documentation
Limited documentation is located in the /kit4b/Docs directory

## Build and installation
### Linux
To build on Ubuntu firstly clone this repository. Dependent on your hosting environment it may be neccessary to firstly update the repository 'aclocal.m4' by running 'aclocal'.


The following example, for Ubuntu 18..22.04 LTS, will install the **kitb** toolkit processes to a user created `bin` directory in the user's home directory. It assumes that the Autotools toolchain ('autotools-dev', 'automake' and 'autoconf') have been installed and are accessible to the user. 

```
mkdir ~/bin
git clone https://github.com/kit4b/kit4b.git
cd kit4b
aclocal
autoreconf -f -i
./configure --prefix=$HOME
make install
```

Alternatively, the binary built for the appropriate platform can be used directly.

### Windows
To build on Windows, this release requires compiling with either Visual Studio 2019 or 2022. 
```
Open the `kit4b.sln` file in Visual Studio. 
Under the Build menu, select Configuration Manager. 
For Active solution platform, select x64. 
The project can then be built. By default, executables will be copied into the `bin64` directory.
```

## Documentation
Incomplete documentation for the core functionality of **kit4b** is available under the `Docs` directory.

## Contributing
kit4b is maintained by Dr Stuart Stephen, who developed **BioKanga** whilst with the Crop Bioinformatics and Data Science team at CSIRO based in Canberra, Australia. Although nominally now retired Dr Stephen is still actively engaged in bioinformatics related research. 

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

## Acknowlegement
This toolkit is a source base clone of **BioKanga** release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) which has been further developed to contain significant functionality enhancements, these enhancements have resulted in compromised backward compatibility with **BioKanga**.

**kit4b** is being released under the Opensource Software License Agreement (GPLv3), the same as **BioKanga** licencing.

