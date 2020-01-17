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
#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CEndian::CEndian(void)
{
INT16 EndianTest = 0x4321;
m_bIsBigEndian = (*(char *)&EndianTest == 0x21) ? false : true;
}

CEndian::~CEndian(void)
{
}

// Swaps 16bit (2byte) endian order
UINT16 
CEndian::SwapUI16Endians(UINT16 Val)
{
return (((Val&0x00FF)<< 8)+((Val&0xFF00)>>8));
}

// Swaps 32bit (4byte) endian order
UINT32 
CEndian::SwapUI32Endians(UINT32 Val)
{
return (((Val&0x000000FF)<<24)+((Val&0x0000FF00)<<8)+
   ((Val&0x00FF0000)>>8)+((Val&0xFF000000)>>24));
}

// Swaps 64bit (8byte) endian order
UINT64 
CEndian::SwapUI64Endians(UINT64 Val)
{
UINT32 LoInt = SwapUI32Endians((UINT32)(Val >> 32));
UINT32 HiInt = SwapUI32Endians((UINT32)Val);
return(((UINT64)HiInt << 32) + (UINT64)LoInt);
}

