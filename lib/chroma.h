// $Id: chroma.h,v 3.0 2006-04-03 04:58:43 edwards Exp $
//
// Main chroma include file. If you include this, you include
// all prototypes
//
/*! \file
 * \brief Primary include file for CHROMA in application codes
 *
 * This is the only file needed for main programs. 
 * If you include this, you include all prototypes
 */

/*! \mainpage  CHROMA
 *
 * \section Description
 *
 * The Chroma package supports data-parallel programming constructs for
 * lattice field theory and in particular lattice QCD. It uses the
 * SciDAC QDP++ data-parallel programming (in C++) that presents a
 * single high-level code image to the user, but can generate highly
 * optimized code for many architectural systems including single node
 * workstations, multi-threaded SMP workstations (soon to come),
 * clusters of workstations via QMP, and classic vector computers.
 */

#ifndef CHROMA_INCLUDE
#define CHROMA_INCLUDE

#include "chromabase.h"
#include "chromainc.h"

#endif
