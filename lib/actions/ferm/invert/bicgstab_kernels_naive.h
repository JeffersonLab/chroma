#ifndef BICGSTAB_KERNELS_NAIVE_H
#define BICGSTAB_KERNELS_NAIVE_H
#include "chromabase.h"


/* The funky kernels used by BiCGStab.
   These versions should always work... */

namespace Chroma { 

  namespace BiCGStabKernels { 

   
    template<typename T>
      void xymz_normx(T& x, const T& y, const T& z, Double& x_norm,
			 const Subset& s)
    {
      x[s] = y-z;
      x_norm = norm2(x,s);
    }

    template<typename T, typename C>
      void yxpaymabz(T& x, T&y, T&z, const C& a, const C& b, const Subset& s)
    {
      T tmp;
      tmp[s] = y - b*z;
      y[s] = x + a*tmp;
    }

    template<typename T>
      void norm2x_cdotxy(const T&x, const T&y, Double& norm2x, DComplex& cdotxy,  const Subset& s) 
      {
	norm2x = norm2(x,s);
	cdotxy = innerProduct(x,y,s);
      }
    
    template<typename T, typename C>
      void xpaypbz(T& x, T& y, T&z, C& a, C& b, const Subset& s)
    {
      T tmp;
      tmp[s] = x + a*y;
      x[s] = tmp + b*z;
    }
    
    template<typename T, typename C> 
      void xmay_normx_cdotzx(T& x,const T& y, const T& z, C& a, Double& normx, DComplex& cdotzx, const Subset& s) 
    {
      x[s] -= a*y;
      normx = norm2(x,s);
      cdotzx = innerProduct(z,x,s);
    }

  }
    


}



#endif
