// -*- C++ -*-
// $Id: ks_local_loops.h,v 1.3 2005-09-03 12:32:10 mcneile Exp $
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




void ks_fuzz_loops(
		 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 LatticeStaggeredFermion & psi_fuzz,
		 const multi1d<LatticeColorMatrix> & u,
		 const multi1d<LatticeColorMatrix> & u_smr,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 int volume_source,
		 int fuzz_width, 
		 int j_decay
		 )  ;



}  // end namespace Chroma

#endif
