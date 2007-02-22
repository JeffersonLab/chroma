// -*- C++ -*-
// $Id: polprec_op.h,v 3.2 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Maybe usefull for polynomial preconditioning. It contructs an
 *     operator that is \f$QP(Q^2)Q\f$, where \f$Q = \gamma_5 M\f$, with M any 4D dirac
 *     operator (including EO preconditioned variants). P(Q^2) is could be any 
 *     function of Q^2.
 *
 *  Initial version Feb. 6, 2006 (Kostas Orginos)
 */

#ifndef _LPOLPREC_H
#define _LPOLPREC_H

#include "handle.h"
#include "linearop.h"

using namespace QDP::Hints;

namespace Chroma 
{
  //! Polynomial preconditioner
  /*! \ingroup linop */
  template<typename T, typename P, typename Q>
  class PolyPrec: public DiffLinearOperator<T,P,Q>
  {
  private:
    Handle< DiffLinearOperator<T,P,Q> > M   ;   // this is the operator
    Handle< DiffLinearOperator<T,P,Q> > Pol ;   // this is the preconditioner
    
  public:
    // constructor
    PolyPrec(Handle< DiffLinearOperator<T,P,Q> > m_,
	     Handle< DiffLinearOperator<T,P,Q> > p_) : M(m_), Pol(p_) {    }
   
   
    ~PolyPrec(){}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return M->getFermBC();}
    
    //! Subset comes from underlying operator
    // Hopefully the subsets of the preconditioner and the matrix
    // are the same...
    inline const Subset& subset() const {return M->subset();}

    inline void g5M(T& out, const T& in,  enum PlusMinus isign) const 
      {
	const int G5=Ns*Ns-1;
      
	T  tmp;
	const Subset& sub = M->subset();
      
	// [ Gamma(5) D ]^{dag} = Gamma(5) D
	// using D = gamma_5 D^{dag} gamma_5
      
	(*M)(tmp, in, PLUS);
	out[sub] = Gamma(G5)*tmp;
      }


    //! The operator is Q P(Q^2) Q
    // Q is g5M 
    // If I ever do this for DWF the g5M definition has to change....
    void operator()(T& chi, const T& psi, PlusMinus isign) const
      {
	T tmp ;
      
	g5M(chi,psi,PLUS) ;
	(*Pol)(tmp,chi,PLUS) ;
	g5M(chi,tmp,PLUS) ;

      }
    

    //! Apply the derivative of the linop
    void deriv(P& ds_u, const T& chi, const T& psi, PlusMinus isign) const 
      {
	const int G5=Ns*Ns-1;
	// So do something like
	ds_u.resize(Nd);
	ds_u = zero;
	P  ds_tmp;
      
      
	const Subset& sub = M->subset();

	T g5Mchi,g5Mpsi ;
	T g5Pg5Mchi, g5Pg5Mpsi ;
	{
	  T tt ;
	  g5M(g5Mchi,chi,PLUS);
	  (*Pol)(tt,g5Mchi,PLUS) ;
	  g5Pg5Mchi[sub] = Gamma(G5)*tt ;

	  g5M(g5Mpsi,psi,PLUS) ;
	  (*Pol)(tt,g5Mpsi,PLUS) ;
	  g5Pg5Mpsi[sub] = Gamma(G5)*tt ;
	}

	M->deriv(ds_tmp, g5Pg5Mchi, psi, PLUS);
	ds_u += ds_tmp;   // build up force
	Pol->deriv(ds_tmp, g5Mchi, g5Mpsi, PLUS);
	ds_u += ds_tmp;   // build up force
	M->deriv(ds_tmp, chi, g5Pg5Mpsi, MINUS);
	ds_u += ds_tmp;   // build up force
      }

  };

}// End Namespace Chroma

#endif 
