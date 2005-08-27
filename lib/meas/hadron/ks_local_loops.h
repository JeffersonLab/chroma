// -*- C++ -*-
// $Id: ks_local_loops.h,v 1.2 2005-08-27 11:27:47 mcneile Exp $
/*! \file
 *  \brief Wrapper routines for computing loops with staggeref fermions
 */        


//
//
//

#ifndef  KS_LOCAL_LOOPS_INC 
#define  KS_LOCAL_LOOPS_INC 

#include "enum_loops_s.h"

namespace Chroma {

void ks_local_loops(
		 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 const multi1d<LatticeColorMatrix> & u,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 int volume_source
		 ) ;



  void ks_local_loops(
		      Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		      LatticeStaggeredFermion & q_source,
		      LatticeStaggeredFermion & psi ,
		      multi1d<LatticeColorMatrix> & u,
		      XMLFileWriter & xml_out,
		      XMLReader & xml_in ,
		      int t_length,
		      Real Mass,
		      int Nsamp,
		      Real RsdCG,
		      int CFGNO,
                 int volume_source
		      ) ;





}  // end namespace Chroma

#endif
