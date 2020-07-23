#pragma once

const int cMaxLenPrefix = 10;					// prefixes must be no longer than this number of chars
const int cAllocPGBuffInSize = 0x1fffffff;		// allocating buffering for input file
const int cAllocPGBuffOutSize = cAllocPGBuffInSize;		// allocating buffering for output files
const int cAllocAsciiChkSize = 0x03fffff;		// check if first 4MB of input file is ascii - only ascii files are parsed

typedef enum TAG_eModePG {
	eMPGDefault = 0,		// default is to prefix fasta descriptors with originating genome identifier
	eMPGFilterPrefix,			// filter SAM alignments removing alignments to targets which are not prefixed
	eMPGPlaceHolder			// used to mark end of processing modes
	}eModePG;
