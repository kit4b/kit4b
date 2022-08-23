/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibility with 'BioKanga'.

Because of the potential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019, 2020
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Original 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cMaxExcludeChroms = 20;		// allow upto this many regexpr for specifying chroms to exclude
const int cMaxIncludeChroms = 20;		// allow upto this many regexpr for specifying chroms to include

class CFilterSAMAlignments
{
	CSAMfile *m_pInBAMfile;				// SAM/BAM input file
	CSAMfile *m_pOutBAMfile;			// SAM/BAM output file

	char m_szLine[cMaxReadLen  * 3];	// used to buffer input lines

	CUtility m_RegExprs;            // regular expression processing

	char m_szFiltChrom[_MAX_PATH];	// used to cache last chrom processed	
	bool m_bFiltChrom;				// and it's filtered status


	int
	SetChromFilters(int NumIncludeChroms,		 // number of chromosome regular expressions to include
						char **ppszIncludeChroms,	 // array of include chromosome regular expressions
						int NumExcludeChroms,		 // number of chromosome expressions to exclude
						char **ppszExcludeChroms);	 // array of exclude chromosome regular expressions

	bool					// true if chrom is accepted, false if chrom not accepted
		AcceptThisChromName(char* pszChrom);   // chromosome name

	char *TrimWhitespace(char *pTxt);	// trim whitespace
	bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
		FileReqWriteCompr(char *pszFile); // If last 3 chars of file name is ".gz" then this file is assumed to require compression

public:
	CFilterSAMAlignments();
	~CFilterSAMAlignments();

	void Reset(void);
	void Init(void);

	int								// number of alignments which were retained and written to output file after filtering was applied
		FilterSAMbyChrom(int NumIncludeChroms,		// number of retained chromosomes regular expressions
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to explicitly exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszInFile,			// input file containing alignments to be filtered
		char *pszOutFile);			// write filtered alignments to this output file
};

