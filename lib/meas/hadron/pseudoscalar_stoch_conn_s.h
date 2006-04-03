// $Id: pseudoscalar_stoch_conn_s.h,v 3.0 2006-04-03 04:59:00 edwards Exp $
#ifndef PSEUDOSCALAR_STOCH_CONN_H
#define PSEUDOSCALAR_STOCH_CONN_H

#include "meas/hadron/stoch_conn_s.h"

namespace Chroma {

  class stoch_conn_corr ; 

  class fourlink_pseudoscalar_stoch_conn  : public stoch_conn_corr
  {
  public :


    void compute(LatticeStaggeredFermion & q_source1,
		 LatticeStaggeredFermion & q_source2,
		 LatticeStaggeredFermion & psi1,
		 LatticeStaggeredFermion & psi2,
		 int isample);

    fourlink_pseudoscalar_stoch_conn(int t_len, int nsample,
			       const multi1d<LatticeColorMatrix> & uin,
			       Stag_shift_option type_of_shift_in)  
      : stoch_conn_corr(t_len,nsample,uin,type_of_shift_in)
      {

	outer_tag = "stoch_conn_g5_cross_1"  ; 
	inner_tag1 = "stoch_conn_corr1" ; 
	inner_tag2 = "stoch_conn_corr2" ; 
      }

    virtual ~fourlink_pseudoscalar_stoch_conn()
      {
      }


  protected:


  } ; 









}  // end namespace Chroma

#endif
