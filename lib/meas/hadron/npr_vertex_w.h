
#ifndef INCLUDE_BuildingBlocks_hh
#define INCLUDE_BuildingBlocks_hh

namespace Chroma {
//###################################################################################//
// CVS Header                                                                        //
//###################################################################################//

static const char* const CVSBuildingBlocks_hh =
  "$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/lib/meas/hadron/npr_vertex_w.h,v 1.2 2006-10-30 22:32:01 edwards Exp $";

//###################################################################################//
// Record CVS Info for BuildingBlocks.hh and BuildingBlocks.cc                       //
//###################################################################################//

void CVSBuildingBlocks( TextWriter & Out );

//###################################################################################//
// Used to Set Requested Link Patterns                                               //
//###################################################################################//

typedef void (*BBLinkPattern)( bool &                          DoThisPattern,
                               bool &                          DoFurtherPatterns,
                               multi1d< unsigned short int > & LinkPattern );

//###################################################################################//
// Construct Building Blocks                                                         //
//###################################################################################//

void BuildingBlocks( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
                     const SftMom &                        Phases,
                     const SftMom &                        PhasesCanonical,
	             const multi2d< string > &             BinaryDataFileNames,
                     const signed short int                T1,
                     const signed short int                T2,
                     const signed short int                Tsrc,
                     const signed short int                Tsnk,
		     const std::string&                    SeqSourceType, 
		     const multi1d< int >&                 SnkMom, 
		     const signed short int                DecayDir,
		     const bool                            TimeReverse,
		     const bool                            Translate );

//###################################################################################//
// Arguments                                                                         //
//###################################################################################//
//                                                                                   //
// B is an array of several backward propagators.                                    //
//                                                                                   //
// F is a forward propagator.                                                        //
//                                                                                   //
// U is a gauge field.                                                               //
//                                                                                   //
// MaxNLinks is the maximum number of links to include in operators.                 //
//                                                                                   //
// LinkPattern is a pointer to a function specifing which link patterns to include.  //
//                                                                                   //
// BinaryDataFileNames is an array of file name patterns with an order corresponding //
// to B.                                                                             //
//                                                                                   //
// T1 is the earliest time slice for which building blocks will be calculated.       //
//                                                                                   //
// T2 is the latest time slice for which building blocks will be calculated.         //
//                                                                                   //
// This function will fill many files with the operators corresponding to each       //
// element of B.                                                                     //
//                                                                                   //
//###################################################################################//
//###################################################################################//

}  // end namespace Chroma

#endif

//###################################################################################//
//###################################################################################//
