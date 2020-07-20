#pragma once

const int cAllocPGBuffInSize = 0x1fffffff;		// allocating buffering for input file
const int cAllocPGBuffOutSize = cAllocPGBuffInSize;		// allocating buffering for output files

typedef enum TAG_eModePG {
	eMPGDefault = 0,		// default is to prefix fasta descriptors with originating genome identifier
	eMPGFilterPrefix,			// filter SAM alignments removing alignments to targets which are not prefixed
	eMPGPlaceHolder			// used to mark end of processing modes
	}eModePG;
