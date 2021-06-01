#pragma once

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
