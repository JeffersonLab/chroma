#ifndef BICGSTAB_KERNELS_SCALARSITE_H
#define BICGSTAB_KERNELS_SCALARSITE_H
#include "chromabase.h"
#include "chroma_config.h"

/* The funky kernels used by BiCGStab.
   These versions should always work... */

namespace Chroma { 

  namespace BiCGStabKernels { 
    void initScalarSiteKernels();
    void finishScalarSiteKernels();

    REAL64* getNormSpace();

    struct ord_xymz_normx_arg 
    {
      REAL64* x_ptr;
      REAL64* y_ptr;
      REAL64* z_ptr;
      REAL64* norm_ptr;
      int atom;
    };

#include "actions/ferm/invert/ord_xmyz_normx_kernel.h"


    template<>
      inline
      void xymz_normx(LatticeDiracFermionD& x, 
		      const LatticeDiracFermionD& y, 
		      const LatticeDiracFermionD& z, 
		      Double& x_norm,
		      const Subset& s)
    {
      REAL64 norm = 0;

      if ( s.hasOrderedRep() ) {
	
	REAL64* x_ptr = (REAL64*)&(x.elem(s.start()).elem(0).elem(0).real());
	REAL64* y_ptr = (REAL64*)&(y.elem(s.start()).elem(0).elem(0).real());
	REAL64* z_ptr = (REAL64*)&(z.elem(s.start()).elem(0).elem(0).real());
	REAL64* norms = getNormSpace();

	ord_xymz_normx_arg  arg={x_ptr,y_ptr,z_ptr,norms,4*3*2};
	int len = (s.end()-s.start()+1);
	
	dispatch_to_threads(len,arg,ord_xymz_normx_kernel);
	norm = norms[0];
	// Sum the norms...
	for(int i=1 ; i < qdpNumThreads(); i++) { 
	  norm += norms[i];
	}
      }
      else { 
	QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	QDP_abort(1);
      }

      QDPInternal::globalSum(norm);
      x_norm = Double(norm);

    }


    struct ord_yxpaymabz_arg
    {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL32* z_ptr;
      REAL32  a_re;
      REAL32  a_im;
      REAL32  b_re;
      REAL32  b_im;
      int atom;
    };

#include "actions/ferm/invert/ord_yxpaymabz_kernel.h"
    template<>
      inline
      void yxpaymabz(LatticeDiracFermionF& x, 
		     LatticeDiracFermionF& y, 
		     LatticeDiracFermionF& z, 
		     const ComplexF& a, const ComplexF& b, const Subset& s)
    {
      //QDPIO::cout << "Using optimized yxpaymabz" << endl;

      if( s.hasOrderedRep() ) { 
	REAL32* x_ptr = (REAL32*)&(x.elem(s.start()).elem(0).elem(0).real());
	REAL32* y_ptr = (REAL32*)&(y.elem(s.start()).elem(0).elem(0).real());
	REAL32* z_ptr = (REAL32*)&(z.elem(s.start()).elem(0).elem(0).real());
	REAL32 a_re = a.elem().elem().elem().real();
	REAL32 a_im = a.elem().elem().elem().imag();
	REAL32 b_re = b.elem().elem().elem().real();
	REAL32 b_im = b.elem().elem().elem().imag();
	ord_yxpaymabz_arg arg={x_ptr,y_ptr,z_ptr,a_re,a_im, b_re, b_im,4*3*2};

	int len = (s.end()-s.start()+1);
	dispatch_to_threads(len,arg,ord_yxpaymabz_kernel);
      }
      else {
	QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	QDP_abort(1);
      }
    }


    struct ord_norm2x_cdotxy_arg
    {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL64* norm_space;
      int atom;

    };

      
#include "actions/ferm/invert/ord_norm2x_cdotxy_kernel.h"

    template<>
      inline
      void norm2x_cdotxy(const LatticeDiracFermionF& x, 
			 const LatticeDiracFermionF& y, Double& norm2x, DComplex& cdotxy,  const Subset& s) 
      {

	REAL64 norm_array[3]={0,0,0};


	if ( s.hasOrderedRep() ) {
	  REAL32* x_ptr = (REAL32*)&(x.elem(s.start()).elem(0).elem(0).real());
	  REAL32* y_ptr = (REAL32*)&(y.elem(s.start()).elem(0).elem(0).real());
	  REAL64* norm_space = getNormSpace();
	  
	  int len = (s.end()-s.start()+1);
	  ord_norm2x_cdotxy_arg arg={x_ptr,y_ptr,norm_space,4*3*2};
	  dispatch_to_threads(len,arg, ord_norm2x_cdotxy_kernel);

	  for(int i=0; i < 3*qdpNumThreads(); i+=3) { 
	    norm_array[0] += norm_space[i];
	    norm_array[1] += norm_space[i+1];
	    norm_array[2] += norm_space[i+2];
	  }
	}
	else { 
	  QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	  QDP_abort(1);
	}
	
	QDPInternal::globalSumArray(norm_array,3);
	
	norm2x.elem().elem().elem().elem()  =norm_array[0];
	cdotxy.elem().elem().elem().real()= norm_array[1];
	cdotxy.elem().elem().elem().imag() = norm_array[2];
	
      }


    struct ord_xpaypbz_arg {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL32* z_ptr;
      REAL32 a_re;
      REAL32 a_im;
      REAL32 b_re;
      REAL32 b_im;
      int atom;
    };


#include "actions/ferm/invert/ord_xpaypbz_kernel.h"

    template<>
      inline
      void xpaypbz(LatticeDiracFermionF& x, 
		   LatticeDiracFermionF& y, 
		   LatticeDiracFermionF& z, 
		   ComplexF& a, 
		   ComplexF& b, 
		   const Subset& s)

    {

      if( s.hasOrderedRep() ) { 
	REAL32* x_ptr = (REAL32*)&(x.elem(s.start()).elem(0).elem(0).real());
	REAL32* y_ptr = (REAL32*)&(y.elem(s.start()).elem(0).elem(0).real());
	REAL32* z_ptr = (REAL32*)&(z.elem(s.start()).elem(0).elem(0).real());
	REAL32 a_re = a.elem().elem().elem().real();
	REAL32 a_im = a.elem().elem().elem().imag();
	REAL32 b_re = b.elem().elem().elem().real();
	REAL32 b_im = b.elem().elem().elem().imag();
	ord_xpaypbz_arg arg={x_ptr,y_ptr,z_ptr,a_re,a_im,b_re,b_im,4*3*2};

	int len=(s.end()-s.start()+1);
	dispatch_to_threads(len,arg, ord_xpaypbz_kernel);
      }
      else {
	QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	QDP_abort(1);
      }
      

    }



    struct ord_xmay_normx_cdotzx_arg {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL32* z_ptr;
      REAL32 a_re;
      REAL32 a_im;
      REAL64* norm_space;
      int atom;
    };


#include "actions/ferm/invert/ord_xmay_normx_cdotzx_kernel.h"


    template<> 
      inline
      void xmay_normx_cdotzx(LatticeDiracFermionF& x,
			     const LatticeDiracFermionF& y, 
			     const LatticeDiracFermionF& z, 
			     ComplexF& a, 
			     Double& normx, 
			     DComplex& cdotzx, const Subset& s) 
    {

      REAL64 norm_array[3]={0,0,0};

      if( s.hasOrderedRep() ) { 
	REAL32* x_ptr = (REAL32*)&(x.elem(s.start()).elem(0).elem(0).real());
	REAL32* y_ptr = (REAL32*)&(y.elem(s.start()).elem(0).elem(0).real());
	REAL32* z_ptr = (REAL32*)&(z.elem(s.start()).elem(0).elem(0).real());
	
	
	REAL32 a_re = a.elem().elem().elem().real();
	REAL32 a_im = a.elem().elem().elem().imag();
	REAL64* norm_space = getNormSpace();
	ord_xmay_normx_cdotzx_arg arg = {x_ptr,y_ptr,z_ptr,a_re,a_im, norm_space,4*3*2};


	int len =(s.end()-s.start()+1);
	dispatch_to_threads(len,arg,ord_xmay_normx_cdotzx_kernel);

	for(int i=0; i < 3*qdpNumThreads(); i+=3) { 
	  norm_array[0] += norm_space[i];
	  norm_array[1] += norm_space[i+1];
	  norm_array[2] += norm_space[i+2];
	}

      }
      else {
	QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	QDP_abort(1);
      }

      QDPInternal::globalSumArray(norm_array,3);
      normx.elem().elem().elem().elem()  =norm_array[0];
      cdotzx.elem().elem().elem().real()= norm_array[1];
      cdotzx.elem().elem().elem().imag() = norm_array[2];
      
      //      x[s] -= a*y;
      //normx = norm2(x,s);
      // cdotzx = innerProduct(z,x,s);
    }


    struct ord_cxmayf_arg {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL32 a_re;
      REAL32 a_im;
      int atom;
    };


#include "actions/ferm/invert/ord_cxmayf_kernel.h"
    template<>
      inline
      void cxmay(LatticeDiracFermionF& x, 
		 const LatticeDiracFermionF& y, 
		 const ComplexF& a, 
		 const Subset& s)

    {

      if( s.hasOrderedRep() ) { 
	REAL32* x_ptr = (REAL32*)&(x.elem(s.start()).elem(0).elem(0).real());
	REAL32* y_ptr = (REAL32*)&(y.elem(s.start()).elem(0).elem(0).real());
	REAL32 a_re = a.elem().elem().elem().real();
	REAL32 a_im = a.elem().elem().elem().imag();
	ord_cxmayf_arg arg={x_ptr,y_ptr,a_re,a_im,4*3*2};

	int len=(s.end()-s.start()+1);
	dispatch_to_threads(len,arg, ord_cxmayf_kernel);
      }
      else {
	QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	QDP_abort(1);
      }
      

    }


    /********************
     * IBICGSTAB KERNELS
     ******************** */


    template<typename F>
    struct ib_zvupdates_arg {
      F* r_ptr;
      F* z_ptr;
      F* v_ptr;
      F* u_ptr;
      F* q_ptr;
     
      F alpha_re;
      F alpha_im;

      F alpha_rat_beta_re;
      F alpha_rat_beta_im;
      
      F alpha_delta_re;
      F alpha_delta_im;
      
      F beta_re;
      F beta_im;

      F delta_re;
      F delta_im;
      int atom;
    };

#include "actions/ferm/invert/ord_ib_zvupdates_kernel.h"

#if 1
    template<>
      inline
    void 
      ibicgstab_zvupdates(const LatticeDiracFermionF3& r, 
			  LatticeDiracFermionF3& z, 
			  LatticeDiracFermionF3& v,
			  const LatticeDiracFermionF3& u, 
			  const LatticeDiracFermionF3& q,
			  const ComplexD& alpha,
			  const ComplexD& alpha_rat_beta,
			  const ComplexD& alpha_delta, 
			  const ComplexD& beta,
			  const ComplexD& delta,
			  const Subset& sub)
      {     
	if( sub.hasOrderedRep() ) { 
	  REAL32* r_ptr = (REAL32*)&(r.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* z_ptr = (REAL32*)&(z.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* v_ptr = (REAL32*)&(v.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* u_ptr = (REAL32*)&(u.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* q_ptr = (REAL32*)&(q.elem(sub.start()).elem(0).elem(0).real());
	
	  REAL32 a_re = (REAL32)alpha.elem().elem().elem().real();
	  REAL32 a_im = (REAL32)alpha.elem().elem().elem().imag();

	  REAL32 arb_re = (REAL32)alpha_rat_beta.elem().elem().elem().real();
	  REAL32 arb_im = (REAL32)alpha_rat_beta.elem().elem().elem().imag();

	  REAL32 ad_re = (REAL32)alpha_delta.elem().elem().elem().real();
	  REAL32 ad_im = (REAL32)alpha_delta.elem().elem().elem().imag();

	  REAL32 beta_re = (REAL32)beta.elem().elem().elem().real();
	  REAL32 beta_im = (REAL32)beta.elem().elem().elem().imag();

	  REAL32 delta_re = (REAL32)delta.elem().elem().elem().real();
	  REAL32 delta_im = (REAL32)delta.elem().elem().elem().imag();
	  
	  ib_zvupdates_arg<REAL32> arg={r_ptr,z_ptr,v_ptr,u_ptr,q_ptr, 
					a_re, a_im, arb_re, arb_im, ad_re, ad_im, 
					beta_re, beta_im, delta_re, delta_im,4*3*2};

	  int len=(sub.end()-sub.start()+1);
	  dispatch_to_threads(len,arg, ord_ib_zvupdates_kernel_real32);
      }
      else {
	QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	QDP_abort(1);
      }
      


    }
#endif    


    template<>
      inline
      void 
      ibicgstab_zvupdates(const LatticeDiracFermionD3& r, 
			  LatticeDiracFermionD3& z, 
			  LatticeDiracFermionD3 &v,
			  const LatticeDiracFermionD3& u, 
			  const LatticeDiracFermionD3& q,
			  const ComplexD& alpha,
			  const ComplexD& alpha_rat_beta,
			  const ComplexD& alpha_delta, 
			  const ComplexD& beta,
			  const ComplexD& delta,
			  const Subset& sub)
      {     
	if( sub.hasOrderedRep() ) { 
	  REAL64* r_ptr = (REAL64*)&(r.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* z_ptr = (REAL64*)&(z.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* v_ptr = (REAL64*)&(v.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* u_ptr = (REAL64*)&(u.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* q_ptr = (REAL64*)&(q.elem(sub.start()).elem(0).elem(0).real());
	
	  REAL64 a_re = alpha.elem().elem().elem().real();
	  REAL64 a_im = alpha.elem().elem().elem().imag();

	  REAL64 arb_re = alpha_rat_beta.elem().elem().elem().real();
	  REAL64 arb_im = alpha_rat_beta.elem().elem().elem().imag();

	  REAL64 ad_re = alpha_delta.elem().elem().elem().real();
	  REAL64 ad_im = alpha_delta.elem().elem().elem().imag();

	  REAL64 beta_re = beta.elem().elem().elem().real();
	  REAL64 beta_im = beta.elem().elem().elem().imag();

	  REAL64 delta_re = delta.elem().elem().elem().real();
	  REAL64 delta_im = delta.elem().elem().elem().imag();
	  
	  ib_zvupdates_arg<REAL64> arg={r_ptr,z_ptr,v_ptr,u_ptr,q_ptr, 
					a_re, a_im, arb_re, arb_im, ad_re, ad_im, 
					beta_re, beta_im, delta_re, delta_im,4*3*2};

	  int len=(sub.end()-sub.start()+1);
	  dispatch_to_threads(len,arg, ord_ib_zvupdates_kernel_real64);
	}
	else {
	  QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	  QDP_abort(1);
	}
      }




    template<typename F> 
      struct ib_rxupdate_arg { 
	F* s_ptr;
	F* t_ptr;
	F* z_ptr;
	F* r_ptr;
	F* x_ptr;

	F omega_re; // omega
	F omega_im; 
	int atom;
      };

#include "actions/ferm/invert/ord_ib_rxupdate_kernel.h"

    template<>
      inline
      void 
      ibicgstab_rxupdate(const ComplexD& omega,
			  const LatticeDiracFermionF3& s,
			  const LatticeDiracFermionF3& t,
			  const LatticeDiracFermionF3& z,
			  LatticeDiracFermionF3& r,
			  LatticeDiracFermionF3& x,
			  const Subset& sub)
      {     
	if( sub.hasOrderedRep() ) { 
	  REAL32* s_ptr = (REAL32*)&(s.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* t_ptr = (REAL32*)&(t.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* z_ptr = (REAL32*)&(z.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* r_ptr = (REAL32*)&(r.elem(sub.start()).elem(0).elem(0).real());
	  REAL32* x_ptr = (REAL32*)&(x.elem(sub.start()).elem(0).elem(0).real());
	
	  REAL32 omega_re = (REAL32)omega.elem().elem().elem().real();
	  REAL32 omega_im = (REAL32)omega.elem().elem().elem().imag();
	  
	  ib_rxupdate_arg<REAL32> arg={s_ptr,t_ptr,z_ptr,r_ptr,x_ptr, 
					omega_re, omega_im,4*3*2};

	  int len=(sub.end()-sub.start()+1);
	  dispatch_to_threads(len,arg, ord_ib_rxupdate_kernel_real32);
	}
	else {
	  QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	  QDP_abort(1);
	}
      }


    template<>
      inline
      void 
      ibicgstab_rxupdate(const ComplexD& omega,
			  const LatticeDiracFermionD3& s,
			  const LatticeDiracFermionD3& t,
			  const LatticeDiracFermionD3& z,
			  LatticeDiracFermionD3& r,
			  LatticeDiracFermionD3& x,
			  const Subset& sub)
      {     
	if( sub.hasOrderedRep() ) { 
	  REAL64* s_ptr = (REAL64*)&(s.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* t_ptr = (REAL64*)&(t.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* z_ptr = (REAL64*)&(z.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* r_ptr = (REAL64*)&(r.elem(sub.start()).elem(0).elem(0).real());
	  REAL64* x_ptr = (REAL64*)&(x.elem(sub.start()).elem(0).elem(0).real());
	
	  REAL64 omega_re = omega.elem().elem().elem().real();
				
	  REAL64 omega_im = omega.elem().elem().elem().imag();
	  
	  ib_rxupdate_arg<REAL64> arg={s_ptr,t_ptr,z_ptr,r_ptr,x_ptr, 
					omega_re, omega_im, 4*3*2};

	  int len=(sub.end()-sub.start()+1);
	  dispatch_to_threads(len,arg, ord_ib_rxupdate_kernel_real64);
	}
	else {
	  QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	  QDP_abort(1);
	}
      }


    template<typename F> 
      struct ib_stupdate_arg { 
	F a_r;   // Real + imaginary parts of alpha
	F a_i;   
	F* r;
	F* u;
	F* v;
	F* q;
	F* r0;
	F* f0;
	F* s;
	F* t;
	
	// Reduction results 
	REAL64* norm_space; // Space for 12 doubles
	int atom;

      };
#include "actions/ferm/invert/ord_ib_stupdates_reduces.h"


    template<>
      inline
      void ibicgstab_stupdates_reduces(const ComplexD& alpha,
				       const LatticeDiracFermionF3& r,
				       const LatticeDiracFermionF3& u,
				       const LatticeDiracFermionF3& v,
				       const LatticeDiracFermionF3& q,
				       const LatticeDiracFermionF3& r0,
				       const LatticeDiracFermionF3& f0,
				       LatticeDiracFermionF3& s,
				       LatticeDiracFermionF3& t,
				       ComplexD& phi,
				       ComplexD& pi,
				       ComplexD& gamma,
				       ComplexD& eta,
				       ComplexD& theta,
				       Double& kappa,
				       Double& rnorm,
				       const Subset& sub)
      
      {     

	REAL64 norm_array[12]={0,0,0,0,0,0,0,0,0,0,0,0};

	if( sub.hasOrderedRep() ) { 
	  ib_stupdate_arg<REAL32> arg = { 
	    (REAL32)alpha.elem().elem().elem().real(),
	    (REAL32)alpha.elem().elem().elem().imag(),
	    (REAL32*)&(r.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(u.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(v.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(q.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(r0.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(f0.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(s.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL32*)&(t.elem(sub.start()).elem(0).elem(0).real()),
	    getNormSpace(),
	    4*3*2};
	  
	  for(int i=0; i < 12*qdpNumThreads(); i++) { 
	    arg.norm_space[i] = (REAL64)0;
	  }


	  int len=(sub.end()-sub.start()+1);
	  dispatch_to_threads(len,arg,ord_ib_stupdates_kernel_real32);
	  for(int i=0; i < qdpNumThreads(); i++) { 
	    norm_array[0] += arg.norm_space[12*i];
	    norm_array[1] += arg.norm_space[12*i+1];
	    norm_array[2] += arg.norm_space[12*i+2];
	    norm_array[3] += arg.norm_space[12*i+3];
	    norm_array[4] += arg.norm_space[12*i+4];
	    norm_array[5] += arg.norm_space[12*i+5];
	    norm_array[6] += arg.norm_space[12*i+6];
	    norm_array[7] += arg.norm_space[12*i+7];
	    norm_array[8] += arg.norm_space[12*i+8];
	    norm_array[9] += arg.norm_space[12*i+9];
	    norm_array[10] += arg.norm_space[12*i+10];
	    norm_array[11] += arg.norm_space[12*i+11];
	  }
	  QDPInternal::globalSumArray(norm_array,12);
	 
	  phi.elem().elem().elem().real() = norm_array[0];
	  phi.elem().elem().elem().imag() = norm_array[1];

	  gamma.elem().elem().elem().real() = norm_array[2];
	  gamma.elem().elem().elem().imag() = norm_array[3];

	  pi.elem().elem().elem().real() = norm_array[4];
	  pi.elem().elem().elem().imag() = norm_array[5];

	  eta.elem().elem().elem().real() = norm_array[6];
	  eta.elem().elem().elem().imag() = norm_array[7];

	  theta.elem().elem().elem().real() = norm_array[8];
	  theta.elem().elem().elem().imag() = norm_array[9];

	  kappa.elem().elem().elem().elem() = norm_array[10];
	  rnorm.elem().elem().elem().elem() = norm_array[11];

	}
	else {
	  QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	  QDP_abort(1);
	}
      }


    template<>
      inline
      void ibicgstab_stupdates_reduces(const ComplexD& alpha,
				       const LatticeDiracFermionD3& r,
				       const LatticeDiracFermionD3& u,
				       const LatticeDiracFermionD3& v,
				       const LatticeDiracFermionD3& q,
				       const LatticeDiracFermionD3& r0,
				       const LatticeDiracFermionD3& f0,
				       LatticeDiracFermionD3& s,
				       LatticeDiracFermionD3& t,
				       ComplexD& phi,
				       ComplexD& pi,
				       ComplexD& gamma,
				       ComplexD& eta,
				       ComplexD& theta,
				       Double& kappa,
				       Double& rnorm,
				       const Subset& sub)
      
      {     

	REAL64 norm_array[12]={0,0,0,0,0,0,0,0,0,0,0,0};

	if( sub.hasOrderedRep() ) { 
	  
	  ib_stupdate_arg<REAL64> arg = { 
	    (REAL64)alpha.elem().elem().elem().real(),
	    (REAL64)alpha.elem().elem().elem().imag(),
	    (REAL64*)&(r.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(u.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(v.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(q.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(r0.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(f0.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(s.elem(sub.start()).elem(0).elem(0).real()),
	    (REAL64*)&(t.elem(sub.start()).elem(0).elem(0).real()),
	    getNormSpace(),
	    4*3*2};
	  
	  for(int i=0; i < 12*qdpNumThreads(); i++) { 
	    arg.norm_space[i] = (REAL64)0;
	  }

	  int len=(sub.end()-sub.start()+1);
	  dispatch_to_threads(len,arg,ord_ib_stupdates_kernel_real64);

	  
	  for(int i=0; i < qdpNumThreads(); i++) { 
	    norm_array[0] += arg.norm_space[12*i];
	    norm_array[1] += arg.norm_space[12*i+1];
	    norm_array[2] += arg.norm_space[12*i+2];
	    norm_array[3] += arg.norm_space[12*i+3];
	    norm_array[4] += arg.norm_space[12*i+4];
	    norm_array[5] += arg.norm_space[12*i+5];
	    norm_array[6] += arg.norm_space[12*i+6];
	    norm_array[7] += arg.norm_space[12*i+7];
	    norm_array[8] += arg.norm_space[12*i+8];
	    norm_array[9] += arg.norm_space[12*i+9];
	    norm_array[10] += arg.norm_space[12*i+10];
	    norm_array[11] += arg.norm_space[12*i+11];
	  }
	  QDPInternal::globalSumArray(norm_array,12);
	 
	  phi.elem().elem().elem().real() = norm_array[0];
	  phi.elem().elem().elem().imag() = norm_array[1];

	  gamma.elem().elem().elem().real() = norm_array[2];
	  gamma.elem().elem().elem().imag() = norm_array[3];

	  pi.elem().elem().elem().real() = norm_array[4];
	  pi.elem().elem().elem().imag() = norm_array[5];

	  eta.elem().elem().elem().real() = norm_array[6];
	  eta.elem().elem().elem().imag() = norm_array[7];

	  theta.elem().elem().elem().real() = norm_array[8];
	  theta.elem().elem().elem().imag() = norm_array[9];

	  kappa.elem().elem().elem().elem() = norm_array[10];
	  rnorm.elem().elem().elem().elem() = norm_array[11];

	}
	else {
	  QDPIO::cerr << "I only work for ordered subsets for now" << endl;
	  QDP_abort(1);
	}

      }





  }
}



#endif
