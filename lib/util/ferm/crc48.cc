/*! \file
 *  \brief 48-bit CRC
 */

#include "crc48.h"

namespace CRC48
{

/************************************************************************ 
  48 bit CRC routines. 
 
  These routines can be used to calculate a 48 bit CRC over an array of 
  characters. 
 
  Author: Gary Mussar    mussar@bnr.ca 
  Dec 1, 1991  These routines are public domain and may be freely use 
               in both academic projects and commercial products. 
 
               The author would appreciate receiving any improvements  
               (or bug fixes) 
**************************************************************************/ 
#include "crc48.h" 
 
/* 
   Define a number required items. 
 
   initializer  - use initialize the CRC accumulator at the start of frame 
   generator    - 48 bit CRC generating polynomial 
   complement   - use to modify CRC residual before transmission 
   CRCHIGHBIT   - bit position in base integer used when carrying a 
		  bit from on array element to the next 
   CRCBYTE0     - expression to extract byte 0 of CRC from structure 
   CRCBYTE1     - expression to extract byte 1 of CRC from structure 
      .               . 
      .               . 
*/ 
#ifdef CRCBASE8 
#	define CRCHIGHBIT            0x80 
	static CRC48_t initializer = {{0x33, 0xa5, 0x57, 0xdf, 0xf1, 0xec}}; 
	static CRC48_t generator   = {{0x28, 0x35, 0x09, 0x71, 0xdb, 0xea}}; 
	static CRC48_t complement  = {{0xcc, 0x5a, 0x57, 0xdf, 0x0e, 0x13}}; 
#	define CRCBYTE0 crc[0] 
#	define CRCBYTE1 crc[1] 
#	define CRCBYTE2 crc[2] 
#	define CRCBYTE3 crc[3] 
#	define CRCBYTE4 crc[4] 
#	define CRCBYTE5 crc[5] 
#endif 
#ifdef CRCBASE16 
#	define CRCHIGHBIT            0x8000 
	static CRC48_t initializer = {{0xa533, 0xdf57, 0xecf1}}; 
	static CRC48_t generator   = {{0x3528, 0x7109, 0xeadb}}; 
	static CRC48_t complement  = {{0x5acc, 0xdf57, 0x130e}}; 
#	define CRCBYTE0 crc[0] 
#	define CRCBYTE1 crc[0] >> 8 
#	define CRCBYTE2 crc[1] 
#	define CRCBYTE3 crc[1] >> 8 
#	define CRCBYTE4 crc[2] 
#	define CRCBYTE5 crc[2] >> 8 
#endif 
#ifdef CRCBASE32 
#	define CRCHIGHBIT            0x80000000L 
	static CRC48_t initializer = {{0xdf57a533L, 0xecf1L}}; 
	static CRC48_t generator   = {{0x71093528L, 0xeadbL}}; 
	static CRC48_t complement  = {{0xdf575accL, 0x130eL}}; 
#	define CRCBYTE0 crc[0] 
#	define CRCBYTE1 crc[0] >> 8 
#	define CRCBYTE2 crc[0] >> 16 
#	define CRCBYTE3 crc[0] >> 24 
#	define CRCBYTE4 crc[1] 
#	define CRCBYTE5 crc[1] >> 8 
#endif 
#ifdef CRCBASE64 
#	define CRCHIGHBIT            0x800000000000L 
	static CRC48_t initializer = {{0xecf1df57a533L}}; 
	static CRC48_t generator   = {{0xeadb71093528L}}; 
	static CRC48_t complement  = {{0x130edf575accL}}; 
#	define CRCBYTE0 crc[0] 
#	define CRCBYTE1 crc[0] >> 8 
#	define CRCBYTE2 crc[0] >> 16 
#	define CRCBYTE3 crc[0] >> 24 
#	define CRCBYTE4 crc[0] >> 32 
#	define CRCBYTE5 crc[0] >> 40 
#endif 
 
void initCRC48(CRC48_t& acc)
{ 
  /* Initialize accumulator */ 
  acc = initializer; 
} 
 
 
void calcCRC48(CRC48_t& acc, const void *dataPtr, int count)
{ 
  /* Perform CRC over every byte in array */ 
  for (unsigned char *p = (unsigned char *)dataPtr; count; count--) 
  { 
    /* Mathematically, it makes no difference whether all 
       the data bits of a byte are XORed at once and then 
       check CRC over the next 8 bits or whether the CRC is 
       done one bit at a time. */ 
 
    /* XOR whole byte into CRC accumulator */ 
    acc.crc[0] ^= *(p++); 
 
    /* Perform CRC check over next 8 bits */ 
    for(int bits=0; bits<8; bits++) 
    { 
      /* Remember value of least significant byte of CRC */ 
      int bit = acc.crc[0] & 0x01; 
      
      /* Shift whole CRC one bit to right */ 
      for (int i=0; i<CRCARRAYSIZE-1; i++) 
      { 
	acc.crc[i] >>= 1; 
	if (acc.crc[i+1] & 0x01) acc.crc[i] ^= CRCHIGHBIT; 
      } 
      acc.crc[CRCARRAYSIZE-1] >>= 1; 
 
      /* If least significant bit was set, we should subtract 
	 the generator polymonial */ 
      if (bit) 
      { 
	for(int i = 0; i<CRCARRAYSIZE; i++) 
	  acc.crc[i] ^= generator.crc[i]; 
      } 
    } 
  } 
} 
 
void getCRC48(const CRC48_t& acc, void *dataPtr, int count)
{ 
  CRC48_t xmit; 
 
  /* Sanity check. Who wants 0 or less bytes? */ 
  if (count < 1) return; 
 
  /* Copy CRC to temporary location and prepare for transmission by 
     XORing with complement */ 
  xmit = acc; 
  for(int i = 0; i < CRCARRAYSIZE; i++) 
    xmit.crc[i] ^= complement.crc[i]; 
 
  unsigned char *p = (unsigned char*)dataPtr; 
 
  /* Painfully extract each byte. This is uglier and less efficient than 
     it should be so that we can accomdate big endian and little endian 
     architectures as well as situations where the real size of the 
     integer used in the CRC array is larger than we have been told */ 
  p[0] = xmit.CRCBYTE0; 
  if (count > 1) 
    p[1] = xmit.CRCBYTE1; 
  if (count > 2) 
    p[2] = xmit.CRCBYTE2; 
  if (count > 3) 
    p[3] = xmit.CRCBYTE3; 
  if (count > 4) 
    p[4] = xmit.CRCBYTE4; 
  if (count > 5) 
    p[5] = xmit.CRCBYTE5; 
  
} 

} // namespace CRC
