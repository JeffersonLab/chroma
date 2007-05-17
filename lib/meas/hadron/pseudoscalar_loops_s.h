// $Id: pseudoscalar_loops_s.h,v 3.2 2007-05-17 15:29:57 egregory Exp $
#ifndef PSEUDOSCALAR_LOOPS_S_H
#define PSEUDOSCALAR_LOOPS_S_H

#include "meas/hadron/loops_s.h"

namespace Chroma {

  class staggered_loops ; 

  class threelink_pseudoscalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    threelink_pseudoscalar_loop(int t_len, 
				int nsample, 
				const multi1d<LatticeColorMatrix> & uin,
				Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_gamma3gamma5_cross_one"  ; 
	inner_tag = "loop" ; 
      }


    ~threelink_pseudoscalar_loop()
      {
      }


  protected:


  } ; 



  class fourlink_pseudoscalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    fourlink_pseudoscalar_loop(int t_len, int nsample,
			       const multi1d<LatticeColorMatrix> & uin,
			       Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {

	outer_tag = "loop_gamma5_cross_one"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~fourlink_pseudoscalar_loop()
      {
      }


  protected:


  } ; 




  class fourlink_pseudoscalar_kilcup_loop  : public staggered_loops
  {
  public :

    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    void compute(LatticeStaggeredFermion & psi, int isample, Real mass) ; 

    fourlink_pseudoscalar_kilcup_loop(int t_len, int nsample,
				      const multi1d<LatticeColorMatrix> & uin, 
				      Stag_shift_option type_of_shift_in)
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_gamma5_cross_one_k"  ; 
	inner_tag = "loop" ; 
      }



    virtual ~fourlink_pseudoscalar_kilcup_loop()
      {
      }


  protected:


  } ; 




  class zerolink_pseudoscalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    zerolink_pseudoscalar_loop(int t_len, 
			       int nsample,
			       const multi1d<LatticeColorMatrix> & uin,
			       Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin, type_of_shift_in)
      {
	outer_tag = "loop_gamma5_cross_gamma5"  ; 
	inner_tag = "loop" ; 
      }

    ~zerolink_pseudoscalar_loop()
      {
      }


  protected:


  } ; 


  // Now fuzzed loops

  class threelink_pseudoscalar_loop_fuzz  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    threelink_pseudoscalar_loop_fuzz(int t_len, 
				int nsample, 
				const multi1d<LatticeColorMatrix> & uin,
				Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_gamma3gamma5_cross_one_fz"  ; 
	inner_tag = "loop" ; 
      }


    ~threelink_pseudoscalar_loop_fuzz()
      {
      }


  protected:


  } ; 



  class fourlink_pseudoscalar_loop_fuzz  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    fourlink_pseudoscalar_loop_fuzz(int t_len, int nsample,
			       const multi1d<LatticeColorMatrix> & uin,
			       Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {

	outer_tag = "loop_gamma5_cross_one_fz"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~fourlink_pseudoscalar_loop_fuzz()
      {
      }


  protected:


  } ; 




  class fourlink_pseudoscalar_kilcup_loop_fuzz  : public staggered_loops
  {
  public :

        void compute(LatticeStaggeredFermion & q_source, 
    		 LatticeStaggeredFermion & psi, int isample) ; 

    //    void compute(LatticeStaggeredFermion & psi, int isample, Real mass) ; 
    void compute(LatticeStaggeredFermion & psi_fuzz,
		 LatticeStaggeredFermion & psi,
		 int isample, Real mass);

    fourlink_pseudoscalar_kilcup_loop_fuzz(int t_len, int nsample,
				      const multi1d<LatticeColorMatrix> & uin, 
				      Stag_shift_option type_of_shift_in)
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_gamma5_cross_one_k_fz"  ; 
	inner_tag = "loop" ; 
      }



    virtual ~fourlink_pseudoscalar_kilcup_loop_fuzz()
      {
      }


  protected:


  } ; 







}  // end namespace Chroma

#endif
