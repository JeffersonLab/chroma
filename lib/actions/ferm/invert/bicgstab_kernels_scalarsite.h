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
    };

#include "actions/ferm/invert/ord_xmyz_normx_kernel.h"

    template<>
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

	ord_xymz_normx_arg  arg={x_ptr,y_ptr,z_ptr,norms};
	int len = 4*3*2*(s.end()-s.start()+1);
	
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

      Internal::globalSum(norm);
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

    };

#include "actions/ferm/invert/ord_yxpaymabz_kernel.h"

    template<>
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
	ord_yxpaymabz_arg arg={x_ptr,y_ptr,z_ptr,a_re,a_im, b_re, b_im};

	int len = 4*3*2*(s.end()-s.start()+1);
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
     

    };

      
#include "actions/ferm/invert/ord_norm2x_cdotxy_kernel.h"

    template<>
      void norm2x_cdotxy(const LatticeDiracFermionF& x, 
			 const LatticeDiracFermionF& y, Double& norm2x, DComplex& cdotxy,  const Subset& s) 
      {

	REAL64 norm_array[3]={0,0,0};


	if ( s.hasOrderedRep() ) {
	  REAL32* x_ptr = (REAL32*)&(x.elem(s.start()).elem(0).elem(0).real());
	  REAL32* y_ptr = (REAL32*)&(y.elem(s.start()).elem(0).elem(0).real());
	  REAL64* norm_space = getNormSpace();
	  
	  int len = 4*3*2*(s.end()-s.start()+1);
	  ord_norm2x_cdotxy_arg arg={x_ptr,y_ptr,norm_space};
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
	
	Internal::globalSumArray(norm_array,3);
	
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
    };


#include "actions/ferm/invert/ord_xpaypbz_kernel.h"

    template<>
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
	ord_xpaypbz_arg arg={x_ptr,y_ptr,z_ptr,a_re,a_im,b_re,b_im};

	int len=4*3*2*(s.end()-s.start()+1);
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
    };

#include "actions/ferm/invert/ord_xmay_normx_cdotzx_kernel.h"
  
    template<> 
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
	ord_xmay_normx_cdotzx_arg arg = {x_ptr,y_ptr,z_ptr,a_re,a_im, norm_space};


	int len =4*3*2*(s.end()-s.start()+1);
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

      Internal::globalSumArray(norm_array,3);
      normx.elem().elem().elem().elem()  =norm_array[0];
      cdotzx.elem().elem().elem().real()= norm_array[1];
      cdotzx.elem().elem().elem().imag() = norm_array[2];
      
      //      x[s] -= a*y;
      //normx = norm2(x,s);
      // cdotzx = innerProduct(z,x,s);
    }


    

  }
}



#endif
