// $Id: pseudoscalar_stoch_conn_s.h,v 3.1 2007-05-14 13:40:39 egregory Exp $
#ifndef PSEUDOSCALAR_STOCH_CONN_H
#define PSEUDOSCALAR_STOCH_CONN_H

#include "meas/hadron/stoch_conn_s.h"

namespace Chroma {

  class stoch_conn_corr ; 

  class fourlink_pseudoscalar_stoch_conn  : public stoch_conn_corr
  {
  private: 
    bool Fsrc;
    bool Fsink;

  public :


    void compute(LatticeStaggeredFermion & q_source1,
		 LatticeStaggeredFermion & q_source2,
		 LatticeStaggeredFermion & psi1,
		 LatticeStaggeredFermion & psi2,
		 int isample);

    void compute(LatticeStaggeredFermion & q_source1,
		 LatticeStaggeredFermion & q_source2,
		 LatticeStaggeredFermion & psi1,
		 LatticeStaggeredFermion & psi2,
		 //		 bool fuzz_src, bool fuzz_sink,
		 const multi1d<LatticeColorMatrix> & u_smr,
		 int fuzz_width,
		 int j_decay,
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

    fourlink_pseudoscalar_stoch_conn(int t_len, int nsample,
			       const multi1d<LatticeColorMatrix> & uin,
				     Stag_shift_option type_of_shift_in, 
				     bool fuzz_src, bool fuzz_sink)  
      : stoch_conn_corr(t_len,nsample,uin,type_of_shift_in)
      {

	Fsrc=fuzz_src;
	Fsink=fuzz_sink;

	outer_tag = "stoch_conn_g5_cross_1"  ; 
	if((fuzz_src)&&(!fuzz_sink)){
	  outer_tag = "stoch_conn_g5_cross_1_fzsrc"  ;
	  printf("stoch_conn_g5_cross_1_fzsrc!!!!!!!\n");

	}
	if((!fuzz_src)&&(fuzz_sink)){
	  outer_tag = "stoch_conn_g5_cross_1_fzsink"  ;
	  printf("stoch_conn_g5_cross_1_fzsink!!!!!!!\n");
	}
	if((fuzz_src)&&(fuzz_sink)){
	  outer_tag = "stoch_conn_g5_cross_1_fzsrc_fzsink"  ;
	  printf("stoch_conn_g5_cross_1_fzsrc_fzsink!!!!!!!\n");
	}
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
