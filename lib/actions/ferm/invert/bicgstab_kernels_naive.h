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


    template<typename T, typename C>
      void cxmay(T& x, const T& y, const C& a, const Subset& s)
    {
      x[s] -= a*y;
    } 



    template<typename T, typename C>
      inline
      void ibicgstab_zvupdates(const T& r, T& z, T &v,
			       const T& u, const T& q,
			       const C& alpha,
			       const C& alpha_rat_beta,
			       const C& alpha_delta, 
			       const C& beta,
			       const C& delta,
			       const Subset& s)
    {
      T tmp;
      tmp[s] = alpha_rat_beta*z;
      z[s] = tmp  + alpha*r ;
      z[s] -= alpha_delta*v;
      
      tmp = v;
      v[s] = u+beta*tmp;
      v[s]-= delta*q;
    }

    template<typename T, typename C, typename F>
      void ibicgstab_stupdates_reduces(const C& alpha,
				       const T& r,
				       const T& u,
				       const T& v,
				       const T& q,
				       const T& r0,
				       const T& f0,
				       T& s,
				       T& t,
				       C& phi,
				       C& pi,
				       C& gamma,
				       C& eta,
				       C& theta,
				       F& kappa,
				       F& rnorm,
				       const Subset& sub)
    {
      s[sub] = r - alpha*v;
      t[sub] = u - alpha*q;

      phi = innerProduct(r0,s,sub);
      gamma = innerProduct(f0,s,sub);
      pi = innerProduct(r0,q,sub);
      eta = innerProduct(f0,t,sub);
      theta = innerProduct(t,s,sub);
      kappa=norm2(t,sub);
      rnorm = norm2(r,sub);
    }
	
    template<typename T, typename C>
      void ibicgstab_rxupdate(const C& omega,
			      const T& s,
			      const T& t,
			      const T& z,
			      T& r,
			      T& x,
			      const Subset& sub)
    {
      r[sub] = s - omega*t;
      x[sub] += omega*s;
      x[sub] += z;
    }
		
  }
    


}



#endif
