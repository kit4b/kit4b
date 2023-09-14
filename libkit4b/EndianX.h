#pragma once
// Note: this header file was originally named as 'Endian.h' but
// this naming convention caused conflicts when compiling with gcc 11.3.0 under
// Windows WSL - I suspect file name case folding - so have renamed to
// be EndianX.h
//
class CEndian
{
protected:
	bool m_bIsBigEndian;				// true if on a big-endian machine
public:
	CEndian(void);
	~CEndian(void);
	uint16_t SwapUI16Endians(uint16_t Val);	// Swaps 16bit (2byte) endian order
	uint32_t SwapUI32Endians(uint32_t Val);	// Swaps 32bit (4byte) endian order
	uint64_t SwapUI64Endians(uint64_t Val);	// Swaps 64bit (8byte) endian order
};
