// $Id: scalar_loops_s.h,v 3.1 2007/04/11 15:24:45 egregory Exp $
#ifndef SCALAR_LOOPS_S_H
#define SCALAR_LOOPS_S_H

#include "meas/hadron/loops_s.h"

namespace Chroma {

  class staggered_loops ; 

  class local_scalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    local_scalar_loop(int t_len, int nsample,
		      const multi1d<LatticeColorMatrix> & uin, 
		      Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_one_cross_one"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~local_scalar_loop()
      {
      }


  protected:

  
  } ; 

  //! Class for local (zero-link) \f$1\otimes1\f$ scalar loop using VKVR trick
  class local_scalar_kilcup_loop
    : public staggered_loops
  {
  public:
    //! Do the measurement (needs psi, not q_source!)
    void compute(LatticeStaggeredFermion& psi,
		 int isample, Real mass);
    
    //! Satisfy virtual function compute(LSF&, LSF&, int)
    void compute(LatticeStaggeredFermion& q_source,
		 LatticeStaggeredFermion& psi, 
		 int isample) {}

    //! Constructor (sets up staggered loop and XML tags)
    local_scalar_kilcup_loop(int t_len, int nsample,
			     const multi1d<LatticeColorMatrix>& uin,
			     Stag_shift_option type_of_shift_in) 
      : staggered_loops(t_len, nsample, uin, type_of_shift_in)
      {
	outer_tag = "loop_one_cross_one_k";
	inner_tag = "loop";
      }

      //! Virtual destructor
      virtual ~local_scalar_kilcup_loop() { }

  protected:
  };

  class non_local_scalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    non_local_scalar_loop(int t_len, int nsample,
			  const multi1d<LatticeColorMatrix> & uin,
			  Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_one_cross_gamma3"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~non_local_scalar_loop()
      {
      }


  protected:

  
  } ; 

  class fourlink_scalar_loop 
    : public staggered_loops
  {
  public:
    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 
  
    fourlink_scalar_loop(int t_len, int nsample,
			 const multi1d<LatticeColorMatrix> & uin,
			 Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_one_cross_gamma5"  ; 
	inner_tag = "loop" ; 
      }

      virtual ~fourlink_scalar_loop() {}
  protected:
  };

  class fourlink_scalar_kilcup_loop
    : public staggered_loops
  {
  public:
    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) {}
    void compute(LatticeStaggeredFermion& psi, int isample, Real Mass);
    
    fourlink_scalar_kilcup_loop(int t_len, int nsample,
				const multi1d<LatticeColorMatrix> & uin,
				Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_one_cross_gamma5_k"  ; 
	inner_tag = "loop" ; 
      }
      
      virtual ~fourlink_scalar_kilcup_loop() {}
  protected:
  };

  // fuzzed loops

  class local_scalar_loop_fuzz  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    local_scalar_loop_fuzz(int t_len, int nsample,
		      const multi1d<LatticeColorMatrix> & uin, 
		      Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_one_cross_one_fz"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~local_scalar_loop_fuzz()
      {
      }


  protected:


  } ; 

  //! Class for local (zero-link) \f$1\otimes1\f$ scalar loop, with VKVR and fuzzing
  class local_scalar_kilcup_loop_fuzz : public staggered_loops
  {
  public:
    //! Do the measurement (with VKVR)
    void compute(LatticeStaggeredFermion& psi_fuzz,
		 LatticeStaggeredFermion& psi,
		 int isample, Real mass);
    
    //! Empty compute (satisfies pure virtual compute in staggered_loops)
    void compute(LatticeStaggeredFermion& q_source,
		 LatticeStaggeredFermion& psi, int isample) { }

    //! Set up staggered loop, and set XML tags
    local_scalar_kilcup_loop_fuzz(int t_len, int nsample,
				  const multi1d<LatticeColorMatrix>& uin,
				  Stag_shift_option type_of_shift_in)
      : staggered_loops(t_len, nsample, uin, type_of_shift_in)
      {
	outer_tag = "loop_one_cross_one_k_fz";
	inner_tag = "loop";
      }

    //! Virtual destructor
    virtual ~local_scalar_kilcup_loop_fuzz() {}
  
  protected:
  };

  class non_local_scalar_loop_fuzz  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    non_local_scalar_loop_fuzz(int t_len, int nsample,
			  const multi1d<LatticeColorMatrix> & uin,
			  Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_one_cross_gamma3_fz"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~non_local_scalar_loop_fuzz()
      {
      }


  protected:


  } ; 


}  // end namespace Chroma

#endif
