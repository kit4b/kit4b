/*
This toolkit is a source base clone of 'BioKanga' release 4.4.2 (https://github.com/csiro-crop-informatics/biokanga) and contains
significant source code changes enabling new functionality and resulting process parameterisation changes. These changes have resulted in
incompatibilty with 'BioKanga'.

Because of the potentential for confusion by users unaware of functionality and process parameterisation changes then the modified source base
and resultant compiled executables have been renamed to 'kit4b' - K-mer Informed Toolkit for Bioinformatics.
The renaming will force users of the 'BioKanga' toolkit to examine scripting which is dependent on existing 'BioKanga'
parameterisations so as to make appropriate changes if wishing to utilise 'kit4b' parameterisations and functionality.

'kit4b' is being released under the Opensource Software License Agreement (GPLv3)
'kit4b' is Copyright (c) 2019
Please contact Dr Stuart Stephen < stuartjs@g3web.com > if you have any questions regarding 'kit4b'.

Orginal 'BioKanga' copyright notice has been retained and immediately follows this notice..
*/
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cAllocOutBuff = 0x100000;		// output buffering allocation size if BED or CSV format output

class CRemapLoci
{
	CSAMfile *m_pInBAMfile;				// if SAM/BAM input
	CSAMfile *m_pOutBAMfile;			// if SAM/BAM output

	CBEDfile *m_pMappingBED;			// BED containing remapping locii


	int m_hOutFile;						// output file handle for BED 
	int m_OutBuffIdx;					// current index into m_pszOutBuff at which to next write output formated for BED 
	int m_AllocOutBuff;				    // output buffer allocated to hold this many chars
	char *m_pszOutBuff;					// allocated output buffer

	char *TrimWhitespace(char *pTxt);	// trim whitespace
	bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
		FileReqWriteCompr(char *pszFile); // If last 3 chars of file name is ".gz" then this file is assumed to require compression

public:
	CRemapLoci();
	~CRemapLoci();

	void Reset(void);

	int
	RemapLocii(int PMode,					// processing mode
					 int FType,				// alignment file type
					char *pszInAlignFile,	// alignment file with loci to be remapped
					char *pszInBEDFile,     // BED file containing loci remapping
					char *pszRemappedFile);	// write remapped alignments to this file

	int
	RemapBEDLocii(char *pszInAlignFile,		// BED alignment file with loci to be remapped
				char *pszRemappedFile);		// write remapped alignments to this file

	int
	RemapSAMLocii(char *pszInAlignFile,		// SAM or BAM alignment file with loci to be remapped
				char *pszRemappedFile);		// write remapped alignments to this file

};

