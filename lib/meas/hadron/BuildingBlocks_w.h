//#####################################################################################
//#####################################################################################
//
// BuildingBlocks.hh
//
//#####################################################################################
//#####################################################################################
//
// description:
//
// BuildingBlocks combines a gauge field, a forward propagator, and seqential
// (or backward) propagators to form the corresponding building blocks (three point
// functions).
//
// authors:
//
// Dru B. Renner, dru@mit.edu, 2002 - port of "MIT" code to qdp++
// Dru B. Renner, dru@mit.edu, July 2003 - analyze MILC configurations
//
//#####################################################################################
//#####################################################################################

#ifndef INCLUDE_BuildingBlocks_hh
#define INCLUDE_BuildingBlocks_hh

//#####################################################################################
//#####################################################################################

static const char* const CVSBuildingBlocks_hh =
  "$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/lib/meas/hadron/BuildingBlocks_w.h,v 1.1 2003-08-19 17:15:24 bjoo Exp $";

//#####################################################################################
//#####################################################################################

#include "qdp.h"
#include "chroma.h"

using namespace QDP;

//#####################################################################################
// Used to Set Requested Link Patterns
//#####################################################################################

typedef void (*BBLinkPattern)( bool &                          DoThisPattern,
                               bool &                          DoFurtherPatterns,
                               multi1d< unsigned short int > & LinkPattern );

//#####################################################################################
// Record CVS Info for BuildingBlocks.hh and BuildingBlocks.cc
//#####################################################################################

void CVSBuildingBlocks( TextWriter & Out );

//#####################################################################################
// Construct Building Blocks
//#####################################################################################

void BuildingBlocks( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
                     const SftMom &                        Phases,
	             const multi2d< const char* > &        BinaryDataFileNames );

//#####################################################################################
//#####################################################################################
//
// B is an array of several backward propagators.
//
// F is a forward propagator.
//
// U is the gauge field.
//
// MaxNLinks is the maximum number of links to include in operators.
//
// LinkPattern is a pointer to a function specifing which link patterns to include.
//
// BinaryDataFileNames is an array of file name patterns with an order corresponding
// to B.
//
// This function will fill many files with the operators corresponding to each
// element of B.
//
//#####################################################################################
//#####################################################################################

#endif

//#####################################################################################
//#####################################################################################
