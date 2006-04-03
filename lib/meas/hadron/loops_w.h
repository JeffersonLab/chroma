// -*- C++ -*-
// $Id: loops_w.h,v 3.0 2006-04-03 04:58:59 edwards Exp $
/*! \file
 *  \brief Quark loop via noise
 */


#ifndef  LOOPS_W_INC
#define  LOOPS_W_INC 

namespace Chroma {

//
//  Compute 
//             trace ( \Gamma M^-1 )
//  using noise sources.
//
//
void loops(const LatticeFermion &q_source,
	   const LatticeFermion &psi,
	   int length,
	   int t0,
	   XMLWriter& xml_gamma,
	   const string& xml_tag) ; 

}  // end namespace Chroma

#endif
