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

#include "SQLiteSummaries.h"

typedef struct TAG_sSubProcess {
	const char *pszName;	// subprocess name
	const char *pszBriefDescr; // subprocess brief description
	const char *pszFullDescr; // subprocess full description
	int (* SubFunct)(int argc, char* argv[]);
} tsSubProcess;

extern const char *cpszProgVer;				// will be incremented with each release
extern CStopWatch gStopWatch;				// time keeper
extern CDiagnostics gDiagnostics;			// for writing diagnostics messages to log file
extern CSQLiteSummaries gSQLiteSummaries;	// for writing processing result summaries to SQLite database
extern int	gExperimentID;					// SQLite experiment identifier
extern int gProcessID;						// SQLite processor identifier
extern int	gProcessingID;					// SQLite processing identifier

extern char gszProcName[_MAX_FNAME];		// process name
extern tsSubProcess *gpszSubProcess;		// selected subprocess