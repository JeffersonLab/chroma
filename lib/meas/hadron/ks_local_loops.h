// -*- C++ -*-
// $Id: ks_local_loops.h,v 1.1 2005-06-27 15:54:21 mcneile Exp $
/*! \file
 *  \brief Wrapper routines for computing loops with staggeref fermions
 */        


//
//
//

#ifndef  KS_LOCAL_LOOPS_INC 
#define  KS_LOCAL_LOOPS_INC 

namespace Chroma {

  enum VolSrc {
    Z2NOISE ,
    GAUSSIAN
  };
  typedef   VolSrc  VolSrc_type ;


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
