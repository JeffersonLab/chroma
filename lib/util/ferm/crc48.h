// -*- C++ -*-
/*! \file
 *  \brief 48-bit CRC
 */

#ifndef  __crc48_h__
#define  __crc48_h__

namespace CRC48
{

/************************************************************************ 
  48 bit CRC routines. 
 
  These routines can be used to calculate a 48 bit CRC over an array of 
  characters. The calculated value can be used to generate data frames 
  that have correct CRC when received by systems which support either 
  CCITT-16 CRCs or CCITT-32 bit CRCs. 
 
  These routines use a structure to hold the 48 bit CRC. The structure 
  consists of an array of integers (or longs). The number of elements 
  required to hold the CRC is dependent on the base size of an integer 
  in the target system. The MACROs CRCBASE8, CRCBASE16, CRCBASE32 and 
  CRCBASE64 are used to indicate the appropriate size. Any base can be 
  used that is less than or equal to the size of integers (or longs), 
  however, the CRC calculations will be less efficient. 
 
  Author: Gary Mussar    mussar@bnr.ca 
  Dec 1, 1991  These routines are public domain and may be freely use 
               in both academic projects and commercial products. 
 
               The author would appreciate receiving any improvements  
               (or bug fixes) 
**************************************************************************/ 
 
/* Use 8 bit base */
#define CRCBASE8

/* 
   Make sure only one CRCBASE size is selected 
*/ 
#ifdef CRCBASE8 
#	undef CRCBASE16 
#	undef CRCBASE32 
#	undef CRCBASE64 
#else 
#	ifdef CRCBASE16 
#		undef CRCBASE32 
#		undef CRCBASE64 
#	else 
#		ifdef CRCBASE32 
#			undef CRCBASE64 
#		else 
#			ifndef CRCBASE64 
#				define CRCBASE32 
#			endif 
#		endif 
#	endif 
#endif 
 
/* 
   Assign appropriate variable types and the number of elements of the 
   array to hold the CRC. The variable type MUST be unsigned. 
*/ 
#ifdef CRCBASE8 
#	define CRCBASETYPE unsigned char 
#	define CRCARRAYSIZE 6 
#endif 
#ifdef CRCBASE16 
#	define CRCBASETYPE unsigned int 
#	define CRCARRAYSIZE 3 
#endif 
#ifdef CRCBASE32 
#	define CRCBASETYPE unsigned long 
#	define CRCARRAYSIZE 2 
#endif 
#ifdef CRCBASE64 
#	define CRCBASETYPE unsigned long 
#	define CRCARRAYSIZE 1 
#endif 
 
/* 
   Define the CRC structure and typedef for easier use. 
*/ 
struct CRC48_t
{ 
  CRCBASETYPE crc[CRCARRAYSIZE]; 
}; 
 
/* 
   Reset the CRC at the beginning of a frame. 
*/ 
void initCRC48(CRC48_t& acc); 
 
/* 
   Compute the CRC over an array of bytes. This routine can be called 
   multiple times as long as the data is transmitted in the same order. 
*/ 
void calcCRC48(CRC48_t& acc, const void *dataPtr, int count); 
 
/* 
   Retrieve part of all of the 48 bit CRC. For machines supporting 32 bit 
   CRCs, only the first two bytes or the 48 bit CRC are required. For 
   machines support 16 bit CRCs, the first 4 bytes are required. The data 
   is retrieved in the order required for transmission (ie. byte[0] is 
   transmitted first, followed by byte[1], etc.). 
*/ 
void getCRC48(const CRC48_t& acc, void *dataPtr, int count); 

} // namespace CRC48

#endif
