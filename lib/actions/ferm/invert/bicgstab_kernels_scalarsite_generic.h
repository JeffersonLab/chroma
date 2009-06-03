#ifndef BICGSTAB_KERNELS_SCALARSITE_GENERIC_H
#define BICGSTAB_KERNELS_SCALARSITE_GENERIC_H
#include "chromabase.h"


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

    void ord_xymz_normx_kernel(int lo, int hi, int my_id, ord_xymz_normx_arg* a)
    {
      REAL64* x_ptr;
      REAL64* y_ptr;
      REAL64* z_ptr;
      REAL64 norm=0;

      x_ptr = &(a->x_ptr[lo]);
      y_ptr = &(a->y_ptr[lo]);
      z_ptr = &(a->z_ptr[lo]);

      int len = hi-lo;
      for(int count = 0; count < len; count+=4) { 

	x_ptr[count] = y_ptr[count] - z_ptr[count];
	x_ptr[count+1] = y_ptr[count+1] - z_ptr[count+1];
	x_ptr[count+2] = y_ptr[count+2] - z_ptr[count+2];
	x_ptr[count+3] = y_ptr[count+3] - z_ptr[count+3];

	norm += x_ptr[count]*x_ptr[count];
	norm += x_ptr[count+1]*x_ptr[count+1];
	norm += x_ptr[count+2]*x_ptr[count+2];
	norm += x_ptr[count+3]*x_ptr[count+3];
      }
      
      a->norm_ptr[my_id] = norm;
    }

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

#if 1

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

    void ord_yxpaymabz_kernel(int lo, int hi, int my_id, ord_yxpaymabz_arg* a)
    {
      REAL32* x_ptr = &(a->x_ptr[lo]);
      REAL32* y_ptr = &(a->y_ptr[lo]);
      REAL32* z_ptr = &(a->z_ptr[lo]);

      REAL32 a_re = a->a_re;
      REAL32 a_im = a->a_im;
      REAL32 b_re = a->b_re;
      REAL32 b_im = a->b_im;

      int len = hi - lo;
      for(int count = 0; count < len; count+=4) { 
	REAL32 tmp_re, tmp_im;
	REAL32 tmp_re2, tmp_im2;
	tmp_re  = y_ptr[count] - b_re*z_ptr[count];
	tmp_re +=  b_im*z_ptr[count+1];
	tmp_im  = y_ptr[count+1] - b_re*z_ptr[count+1];
	tmp_im -=  b_im*z_ptr[count];
	

	tmp_re2  = y_ptr[count+2] - b_re*z_ptr[count+2];
	tmp_re2 +=  b_im*z_ptr[count+3];
	tmp_im2  = y_ptr[count+3] - b_re*z_ptr[count+3];
	tmp_im2 -=  b_im*z_ptr[count+2];
	
	y_ptr[count]   = x_ptr[count] + a_re*tmp_re ;
	y_ptr[count]  -=  a_im*tmp_im;
	y_ptr[count+1] = x_ptr[count+1] + a_re*tmp_im ;
	y_ptr[count+1]+= a_im*tmp_re; 
	y_ptr[count+2]   = x_ptr[count+2] + a_re*tmp_re2 ;
	y_ptr[count+2]  -=  a_im*tmp_im2;
	y_ptr[count+3] = x_ptr[count+3] + a_re*tmp_im2 ;
	y_ptr[count+3]+= a_im*tmp_re2; 
      }
    }

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
      // T tmp;
      // tmp[s] = y - b*z;
      // y[s] = x + a*tmp;
    }
#endif

#if 1

    struct ord_norm2x_cdotxy_arg
    {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL64* norm_space;
     

    };


    void ord_norm2x_cdotxy_kernel(int lo, int hi, int my_id, ord_norm2x_cdotxy_arg* a)
    {
      REAL32* x_ptr = &(a->x_ptr[lo]);
      REAL32* y_ptr = &(a->y_ptr[lo]);
      REAL64 norm_array[3] = {0,0,0};

      int len = hi-lo;
      for(int count = 0; count < len; count+=4) { 
	// norm
	norm_array[0] += x_ptr[count]*x_ptr[count];
	norm_array[0] += x_ptr[count+1]*x_ptr[count+1];
	norm_array[0] += x_ptr[count+2]*x_ptr[count+2];
	norm_array[0] += x_ptr[count+3]*x_ptr[count+3];
	
	// cdot re
	norm_array[1] += x_ptr[count]*y_ptr[count];
	norm_array[1] += x_ptr[count+1]*y_ptr[count+1];
	norm_array[1] += x_ptr[count+2]*y_ptr[count+2];
	norm_array[1] += x_ptr[count+3]*y_ptr[count+3];
	
	// cdot im
	norm_array[2] += x_ptr[count]*y_ptr[count+1];
	norm_array[2] -= x_ptr[count+1]*y_ptr[count];
	norm_array[2] += x_ptr[count+2]*y_ptr[count+3];
	norm_array[2] -= x_ptr[count+3]*y_ptr[count+2];
      }
      a->norm_space[3*my_id]=norm_array[0];
      a->norm_space[3*my_id+1]=norm_array[1];
      a->norm_space[3*my_id+2]=norm_array[2];
    }
      
      

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
	
	// norm2x = norm2(x,s);
	//cdotxy = innerProduct(x,y,s);
      }
 
#endif

#if 1

    struct ord_xpaypbz_arg {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL32* z_ptr;
      REAL32 a_re;
      REAL32 a_im;
      REAL32 b_re;
      REAL32 b_im;
    };

    void ord_xpaypbz_kernel(int lo, int hi, int my_id, ord_xpaypbz_arg* a)
    {
      REAL32* x_ptr = &(a->x_ptr[lo]);
      REAL32* y_ptr = &(a->y_ptr[lo]);
      REAL32* z_ptr = &(a->z_ptr[lo]);

      REAL32 a_re = a->a_re;
      REAL32 a_im = a->a_im;
      REAL32 b_re = a->b_re;
      REAL32 b_im = a->b_im;

      int len = hi - lo;
      for(int count = 0; count < len; count+=4) { 
	  REAL32 tmp_re, tmp_im;
	  REAL32 tmp_re2, tmp_im2;

	  tmp_re  = x_ptr[count] + a_re*y_ptr[count];
	  tmp_re -=  a_im*y_ptr[count+1];
	  tmp_im  = x_ptr[count+1] + a_re*y_ptr[count+1];
	  tmp_im +=  a_im*y_ptr[count];
	  


	  tmp_re2  = x_ptr[count+2] + a_re*y_ptr[count+2];
	  tmp_re2 -=  a_im*y_ptr[count+3];
	  tmp_im2  = x_ptr[count+3] + a_re*y_ptr[count+3];
	  tmp_im2 +=  a_im*y_ptr[count+2];


	  x_ptr[count]   = tmp_re + b_re*z_ptr[count] ;
	  x_ptr[count]  -=  b_im*z_ptr[count+1];
	  x_ptr[count+1] = tmp_im + b_re*z_ptr[count+1];
	  x_ptr[count+1]+= b_im *z_ptr[count]; 

	  
	  x_ptr[count+2]   = tmp_re2 + b_re*z_ptr[count+2] ;
	  x_ptr[count+2]  -=  b_im*z_ptr[count+3];
	  x_ptr[count+3] = tmp_im2 + b_re*z_ptr[count+3];
	  x_ptr[count+3]+= b_im *z_ptr[count+2]; 

      }
    }

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
      

    //      T tmp;
    // tmp[s] = x + a*y;
    //  x[s] = tmp + b*z;
    }

#endif

#if 1    

    struct ord_xmay_normx_cdotzx_arg {
      REAL32* x_ptr;
      REAL32* y_ptr;
      REAL32* z_ptr;
      REAL32 a_re;
      REAL32 a_im;
      REAL64* norm_space;
    };

    void ord_xmay_normx_cdotzx_kernel(int lo, int hi, int my_id, ord_xmay_normx_cdotzx_arg *a) 
    {

      REAL32* x_ptr=&(a->x_ptr[lo]);
      REAL32* y_ptr=&(a->y_ptr[lo]);
      REAL32* z_ptr=&(a->z_ptr[lo]);
      REAL32 a_re = a->a_re;
      REAL32 a_im = a->a_im;
      REAL64 norm_array[3]={0,0,0};

      int len = hi-lo;
      for(int count = 0; count < len; count+=4) { 

	x_ptr[count] -= a_re*y_ptr[count];
	x_ptr[count] += a_im*y_ptr[count+1];
	
	x_ptr[count+1] -= a_im*y_ptr[count];
	x_ptr[count+1] -= a_re*y_ptr[count+1];
	
	x_ptr[count+2] -= a_re*y_ptr[count+2];
	x_ptr[count+2] += a_im*y_ptr[count+3];
	
	x_ptr[count+3] -= a_im*y_ptr[count+2];
	x_ptr[count+3] -= a_re*y_ptr[count+3];
	
	norm_array[0] += x_ptr[count]*x_ptr[count];
	norm_array[0] += x_ptr[count+1]*x_ptr[count+1];
	norm_array[0] += x_ptr[count+2]*x_ptr[count+2];
	norm_array[0] += x_ptr[count+3]*x_ptr[count+3];
	
	norm_array[1] += z_ptr[count]*x_ptr[count];
	norm_array[1] += z_ptr[count+1]*x_ptr[count+1];
	norm_array[1] += z_ptr[count+2]*x_ptr[count+2];
	norm_array[1] += z_ptr[count+3]*x_ptr[count+3];
	
	norm_array[2] += z_ptr[count]*x_ptr[count+1];
	norm_array[2] -= z_ptr[count+1]*x_ptr[count];
	norm_array[2] += z_ptr[count+2]*x_ptr[count+3];
	norm_array[2] -= z_ptr[count+3]*x_ptr[count+2];
	
      }
      a->norm_space[3*my_id]=norm_array[0];
      a->norm_space[3*my_id+1]=norm_array[1];
      a->norm_space[3*my_id+2]=norm_array[2];
    }      

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


    
#endif
  }
}



#endif
