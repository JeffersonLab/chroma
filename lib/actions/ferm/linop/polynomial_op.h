// -*- C++ -*-
// $Id: polynomial_op.h,v 3.5 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Polynomial filter for 4D Dirac operators. It creates an approximation
 *    to 1/Q^2 in the range \f$[\mu, \lambda_{max}]\f$ with \f$Q = \gamma5 M\f$
 *    where M is a dirac operator of some kind (EO preconditioned is accepted).
 *    Can handle any 4D operator but not 5D since gamma_5 hermiticity is more
 *    involved there. 
 *
 *   Initial version Feb. 6, 2006 (Kostas Orginos)
 */

#ifndef _LPOLY_H
#define _LPOLY_H

#include "handle.h"
#include "linearop.h"
#include "polylinop.h"

using namespace QDP::Hints;

namespace Chroma 
{
  //! Polynomial operator
  /*! \ingroup linop */
  template<typename T, typename P, typename Q>
  class lpoly: public PolyLinearOperator<T,P,Q>
  {
  private:
    Handle< DiffLinearOperator<T,P,Q> > M;   // this is the operator

    Double LowerBound ;
    Double UpperBound ;

    int degree ;

    // The polynomium is: c_Zero * Prod_i ( Q^2 - root[i])*inv_c[i] 

    multi1d<DComplex> root ;
    multi1d<DComplex> inv_c ;

    Real c_Zero ; // the zeroth order approximation to the inverse  ie. a constant


    int bitrevers(int n,int maxBits){
      // int bits[maxBits] ;
      // int br_bits[maxBits] ;
      multi1d<int> bits(maxBits);
      multi1d<int> br_bits(maxBits);

      int br(0);
      for(int i(0);i<maxBits;i++){
	bits[i] = n % 2 ;
	n /= 2;
	int br_i = maxBits-i-1 ;
	br_bits[br_i] = bits[i] ;
	br = br + (br_bits[br_i]<<br_i) ;
      }
      return br ;
    }

    void  GetRoots(multi1d<DComplex>& r,     multi1d<DComplex>& ic){

      Double eps = LowerBound/UpperBound ;
      // 2PI/(N+1) ;
      Double www = 6.28318530717958647696/(degree+1.0) ;
      
      // complex conjugate pairs are at i and N-1-i
      for(int i(0);i<degree;i++){
	//cout<<"i: "<<i<<endl ;
	Real ii = i+1.0 ;
	r[i] = UpperBound*cmplx(0.5*(1.0+eps)*(1.0-cos(www*ii)) , -sqrt(eps)*sin(www*ii));
	ic[i] = 1.0/(UpperBound*0.5*(1.0+eps) - r[i]) ;
      }
      c_Zero = 2.0/(UpperBound+LowerBound) ;
      // The c_Zero constant is computed here
      // The polynomium is: c_Zero * Prod_i ( Q^2 - root[i])*inv_c[i] 
    }


  public:
    // constructor
    // need to modify the contructor to pass down the roots
    // this is doing my own stupid ordering...
    lpoly(Handle< DiffLinearOperator<T,P,Q> > m_,
	  int degree_,
	  const Real& lower_bound_,
	  const Real& upper_bound_, 
	  int ord) : M(m_), LowerBound(lower_bound_),  UpperBound(upper_bound_), degree(degree_)
      { //This sets up Chebyshev Polynomials. But we should have a general
	// class that allows us to use different types of polynomials
     
	//Qsq = new lmdagm<T>(*M) ;

	//QDPIO::cout<<"lpoly constructor\n" ;

	if(degree_%2 !=0 ){
	  QDPIO::cout<<"lpoly: Polynomium degree must be even.\n" ;
	  degree_++ ;
	  degree++  ;
	  QDPIO::cout<<"lpoly: Using degree:" << degree<<endl ;       
	}
	//UpperBound = upper_bound_ ;
	root.resize(degree_);
	inv_c.resize(degree_);
	
	// get the natural order roots and constants
	multi1d<DComplex> r(degree);
	multi1d<DComplex> ic(degree);
	GetRoots(r,ic);

	// Arrange the roots in complex conjugate pairs     
	int j(0) ;
	// complex conjugate pairs are at i and N-1-i
	for(int k(1);k<=ord;k++){
	  //QDPIO::cout<<"k: "<<k<<endl ;
	  for(int i(k-1);i<degree/2;i+=ord){
	    //QDPIO::cout<<"i: "<<i<<endl ;
	    root[j] = r[i] ;
	    inv_c[j] = ic[i] ;
	    root [degree-1-j] = conj(root[j]    ) ;
	    inv_c[degree-1-j] = conj(inv_c[j]) ;
	    j++ ;
	  }
	}
	// The polynomium is: c_Zero * Prod_i ( Q^2 - root[i])*inv_c[i] 
     
      }

    // This is doing the Jensen bit reversal ordering
    lpoly(Handle< DiffLinearOperator<T,P,Q> > m_,
	  int degree_,
	  const Real& lower_bound_,
	  const Real& upper_bound_) : M(m_), LowerBound(lower_bound_),  UpperBound(upper_bound_), degree(degree_){
      //This sets up Chebyshev Polynomials. But we should have a general
      // class that allows us to use different types of polynomials
      
      //QDPIO::cout<<"tlpoly constructor\n" ;
      //QDPIO::cout<<"tlpoly: degree: " << degree<<endl ;  

      int maxBits = 0 ;
      while(degree_%2==0)
	{
	  degree_/=2 ;
	  maxBits ++ ;
	}
      
      if(degree_ !=1 ){
	QDPIO::cout<<"tlpoly: Polynomium degree must be power of two.\n" ;
	while(degree_!=1){
	  degree_ = degree + 1 ;
	  degree = degree_ ;
	  QDPIO::cout<<"tlpoly: Using degree: " << degree<<endl ;  
	  maxBits=0 ;
	  while(degree_%2==0){
	    degree_/=2 ;
	    maxBits++ ;
	  }
	}
	
	QDPIO::cout<<"tlpoly: Using degree: " << degree<<endl ;        
      }
      
      //QDPIO::cout<<"tlpoly: maxBits: " << maxBits<<endl ;        
      //UpperBound = upper_bound_ ;
      root.resize(degree);
      inv_c.resize(degree);
      
      // get the natural order roots and constants
      multi1d<DComplex> r(degree);
      multi1d<DComplex> ic(degree);
      GetRoots(r,ic);
      
      // complex conjugate pairs are at i and N-1-i
      for(int i(0);i<degree;i++){
	//cout<<"i: "<<i<<endl ;
	int ii = bitrevers(i,maxBits) ;
	root[i] = r[ii] ;
	inv_c[i] = ic[ii] ;
      }
      // The polynomium is: c_Zero * Prod_i ( Q^2 - root[i])*inv_c[i] 
    }

    ~lpoly(){}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return M->getFermBC();}
    
    //! Subset comes from underlying operator
    inline const Subset& subset() const {return M->subset();}

    //! Returns  the roots
    multi1d<DComplex> Roots() const {
      return root ;
    }

    //! Returns  the roots
    multi1d<DComplex> MonomialNorm() const {
      return inv_c ;
    }


    //! Apply the underlying (Mdagger M) operator
    void Qsq(T& y, const T& x) const {
      T t ;   moveToFastMemoryHint(t);
      (*M)(t,x,PLUS) ;
      (*M)(y,t,MINUS);
    }

    void applyChebInv(T& x, const T& b) const 
      {
	//Chi has the initial guess as it comes in
	//b is the right hand side
	Real a = 2/(UpperBound-LowerBound) ; 
	Real d = (UpperBound+LowerBound)/(UpperBound-LowerBound) ; 
	T QsqX ;
	//*Qsq(QsqX,x,PLUS) ;
	Qsq(QsqX,x) ;
	T r ;
	r = b - QsqX ;
	T y(x) ;
	int m(1) ;
	Real sigma_m(d) ;
	Real sigma_m_minus_one(1.0) ;
	Real sigma_m_plus_one ;
     
	x = x + a/d*r ;
	//*Qsq(QsqX,x,PLUS) ;
	Qsq(QsqX,x) ;
	r = b - QsqX ;
	T p ;
	while (m<degree+1){
	  //QDPIO::cout<<"iter: "<<m<<" 1/sigma: "<<1.0/sigma_m<<endl ;
	  sigma_m_plus_one = 2.0*d*sigma_m - sigma_m_minus_one ;
	  p = 2.0*sigma_m/sigma_m_plus_one *(d*x+a*r) - 
	    sigma_m_minus_one/sigma_m_plus_one*y ;
	  y = x; x= p ;
	  //*Qsq(QsqX,x,PLUS) ;
	  Qsq(QsqX,x) ;
	  r = b - QsqX ;

	  m++;
	  sigma_m_minus_one = sigma_m ;
	  sigma_m = sigma_m_plus_one ;
	}
	QDPIO::cout<<"applyChebInv: 1/sigma("<<m<<"): "<<1.0/sigma_m<<endl ;

      }

    //! Here is your apply function
    // use the Golub algorithm here
    void operator()(T& chi, const T& b, PlusMinus isign) const 
      {
	chi = zero ; // zero the initial guess. This way P(Qsq)*b is produced     
	applyChebInv(chi,b) ;
      }
   
    //! Apply the A or A_dagger piece: P(Qsq) = A_dagger(Qsq) * A(Qsq) 
    // use the root representation here
    void applyA(T& chi,const T& b, PlusMinus isign) const{

      int start, end ;
      chi = b ;
      switch (isign){
      case PLUS:
	start = 0 ;
	end = degree/2 ;
	break ;
      case MINUS:
	start = degree/2  ;
	end   = degree ;
	break ;
      }

      //QDPIO::cout<<start<<" "<<end<<endl ;

     
      Real c0 = sqrt(c_Zero) ;

      T tt = c0*b ;

      for(int i(start);i<end;i++){
	//QDPIO::cout<<" root, norm: "<<root[i]<<" "<<inv_c[i]<<endl ;
	//*Qsq(chi, tt, PLUS);
	Qsq(chi, tt);
	chi -= root[i] * tt;
	chi *= inv_c[i] ;
	tt = chi ;
      }

    }

    //! Apply the derivative of the linop
    void deriv(P& ds_u, const T& chi, const T& psi, PlusMinus isign) const
      {
	// There's a problem here. The interface wants to do something like
	//    \f$\chi^\dag * \dot(Qsq) \psi\f$
	// the derivs all expect to do the derivative on the Qsq leaving
	// open some color indices. The chi^dag and psi multiplications will
	// trace over the spin indices.
	// However, there is no guarantee in the interface that  chi = psi,
	// so you **MUST** apply all the parts of the monomial each to psi
	// and to chi. Sorry, current limitation. In the future, we could
	// add another deriv function that takes only one fermion and you
	// could optimize here.

	// So do something like
	ds_u.resize(Nd);
	ds_u = zero;
	P  ds_tmp;

	multi1d<T> chi_products(degree+1);
	multi1d<T> psi_products(degree+1);

	multi1d<T> M_chi_products(degree+1);
	multi1d<T> M_psi_products(degree+1);

	// Build up list of products of roots of poly against psi and chi each
	// play some trick trying to pad boundaries


	// Qsq is M^dag M 
	chi_products[degree] = sqrt(c_Zero)*chi;
	for(int n(degree-1); n >= 0; --n){
	  Qsq(chi_products[n], chi_products[n+1]);
	  chi_products[n] -= conj(root [n]) * chi_products[n+1];
	  chi_products[n] *= conj(inv_c[n]) ;
	}
 
	psi_products[0] = sqrt(c_Zero)*psi;  
	for(int n(0); n < degree; ++n){
	  Qsq(psi_products[n+1], psi_products[n]);
	  psi_products[n+1] -= root[n] * psi_products[n];
	  psi_products[n+1] *= inv_c[n] ;
	}

	for(int n(0); n < degree+1; ++n)
	{
	  (*M)(M_chi_products[n],chi_products[n],PLUS);
	  (*M)(M_psi_products[n],psi_products[n],PLUS);
	}

	// Now do force computation
	// be more careful than me about bouds
	for(int n(0); n < degree; ++n)
	{
	  M->deriv(ds_tmp, M_chi_products[n+1], psi_products[n], PLUS);
	  for(int d(0);d<Nd;d++)
	    ds_tmp[d] *= inv_c[n] ; 
	  ds_u += ds_tmp;   // build up force
       
	  M->deriv(ds_tmp, chi_products[n+1], M_psi_products[n], MINUS);
	  for(int d(0);d<Nd;d++)
	    ds_tmp[d] *= inv_c[n] ; 
	  ds_u += ds_tmp;   // build up force
	}

	//Check The force:
	/**
	   QDPIO::cout<<"CheckProducts: ";
	   QDPIO::cout<< innerProduct(psi_products[0],psi_products[degree])<<endl ;
	   for(int n(0); n < degree+1; ++n){     
	   QDPIO::cout<<"CheckProducts: "<<n<< ": " ;
	   QDPIO::cout<< innerProduct(chi_products[n],psi_products[n])<<endl ;
     
	   }
	**/
      }

  };

}// End Namespace Chroma

#endif 
