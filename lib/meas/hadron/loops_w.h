// -*- C++ -*-
// $Id: loops_w.h,v 1.1 2004-02-08 11:23:01 mcneile Exp $
/*! \file
 *  \brief Quark loop via noise
 */


//
//  Compute 
//             trace ( \Gamma M^-1 )
//  using noise sources.
//
//

#ifndef  LOOPS_W_INC
#define  LOOPS_W_INC 

void loops(const LatticeFermion &q_source,
	   const LatticeFermion &psi,
	   int length,
	   int t0,
	   XMLWriter& xml_gamma,
	   const string& xml_tag) ; 


#endif
