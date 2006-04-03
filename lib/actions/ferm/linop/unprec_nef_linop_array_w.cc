// $Id: unprec_nef_linop_array_w.cc,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine
  /*! \ingroup linop
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param b5_           NEF parameter (Read)
   * \param c5_           NEF parameter (Read)
   * \param m_q_          quark mass    (Read)
   */
  void UnprecNEFDWLinOpArray::create(Handle< FermState<T,P,Q> > state,
				     const Real& WilsonMass_, 
				     const multi1d<Real>& b5_, const multi1d<Real>& c5_, 
				     const Real& m_q_, int N5_)
  {
    WilsonMass = WilsonMass_;
    b5 = b5_ ;
    c5 = c5_ ;
    m_q = m_q_ ;
    N5  = N5_ ;

    fb5.resize(N5);
    fc5.resize(N5);
    Real ff = Nd - WilsonMass;
    for(int s=0; s < N5; ++s)
    {
      fb5[s] = b5[s]*ff + 1;
      fc5[s] = c5[s]*ff - 1;
    }

    D.create(state);
    fbc = state->getFermBC();
  }



  //! Apply unpreconditioned domain-wall fermion linear operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecNEFDWLinOpArray::operator() (multi1d<LatticeFermion>& chi, 
					  const multi1d<LatticeFermion>& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    if( chi.size() != N5 )  chi.resize(N5);

    //
    //  Chi   =  D' Psi
    //
    LatticeFermion  tmp,c_tmp,cmb ;
    moveToFastMemoryHint(tmp);
    moveToFastMemoryHint(c_tmp);
    moveToFastMemoryHint(cmb);

    switch (isign)
    {
    case PLUS:
    {
      for(int n=0; n < N5; ++n)
      {
	if (n == 0)
	  //c_tmp=0.5*(psi[1]-
	  //   m_q*psi[N5-1] - GammaConst<Ns,Ns*Ns-1>()*(m_q*psi[N5-1]+psi[1])) ;
	  c_tmp = chiralProjectMinus(psi[1]) -m_q*chiralProjectPlus(psi[N5-1]);
	else if (n == N5-1)
	  //c_tmp = 0.5*(psi[N5-2]  - 
	  //   m_q*psi[0] + GammaConst<Ns,Ns*Ns-1>()*(psi[N5-2] + m_q*psi[0])) ;
	  c_tmp = chiralProjectPlus(psi[N5-2]) - m_q*chiralProjectMinus(psi[0]) ;
	else
	  //c_tmp=0.5*(psi[n+1]+
	  //	   psi[n-1]+GammaConst<Ns,Ns*Ns-1>()*(psi[n-1]-psi[n+1]));
	  c_tmp = chiralProjectPlus(psi[n-1]) + chiralProjectMinus(psi[n+1]);
      
	cmb = b5[n]*psi[n] + c5[n]*c_tmp;
	D(tmp, cmb, isign);
	chi[n]  = fb5[n]*psi[n] + fc5[n]*c_tmp;
	chi[n] -= 0.5*tmp;
      }
    }
    break;

    case MINUS:
    {
      multi1d<LatticeFermion>  Dpsi(N5), c5Dpsi(N5);
      for(int n=0; n < N5; ++n)
      {
	D(Dpsi[n], psi[n], isign);

	Real cc5 = -0.5*c5[n];
	c5Dpsi[n] = fc5[n]*psi[n] + cc5*Dpsi[n];
      }

      for(int n=0; n < N5; ++n)
      {
 	if (n == 0){
	  c_tmp = chiralProjectPlus(c5Dpsi[1]) - m_q*chiralProjectMinus(c5Dpsi[N5-1]);
	}
	else if (n == N5-1)
	  c_tmp = chiralProjectMinus(c5Dpsi[N5-2]) - m_q*chiralProjectPlus(c5Dpsi[0]);
	else
	  c_tmp = chiralProjectMinus(c5Dpsi[n-1]) + chiralProjectPlus(c5Dpsi[n+1]);
      
	Real bb5 = -0.5*b5[n];
	chi[n]  = fb5[n]*psi[n] + bb5*Dpsi[n];
	chi[n] += c_tmp;
      }
    }
    break;
    }

    getFermBC().modifyF(chi);

    END_CODE();
  }

  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void UnprecNEFDWLinOpArray::Dminus(LatticeFermion& chi,
				     const LatticeFermion& psi,
				     enum PlusMinus isign,
				     int s5) const
  {
    LatticeFermion tt ;
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);

    chi = fc5[s5]*psi + Real(-0.5*c5[s5])*tt;
    chi *= -1;
  }


  //! Derivative
  void 
  UnprecNEFDWLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign) const
  {
    START_CODE();

    // THIS IS A VERY INEFFICIENT IMPLEMENTATION - JUST WANT IT FOR DEBUGGING

    ds_u.resize(Nd);
    ds_u = zero;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);

    LatticeFermion  tmp,c_tmp,cmb;

    switch (isign)
    {
    case PLUS:
    {
      for(int n=0; n < N5; ++n)
      {
	if (n == 0)
	{
	  tmp = b5[0]*psi[0];
	  D.deriv(ds_tmp, chi[0], tmp, isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectMinus(psi[1]);
	  tmp *= c5[0];
	  D.deriv(ds_tmp, chi[0], tmp, isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectPlus(psi[N5-1]);
	  tmp *= -m_q*c5[0];
	  D.deriv(ds_tmp, chi[0], tmp, isign);
	  ds_u += ds_tmp;
	}
	else if (n == N5-1)
	{
	  tmp = b5[n]*psi[n];
	  D.deriv(ds_tmp, chi[n], tmp, isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectPlus(psi[n-1]);
	  tmp *= c5[n];
	  D.deriv(ds_tmp, chi[n], tmp, isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectMinus(psi[0]);
	  tmp *= -m_q*c5[n];
	  D.deriv(ds_tmp, chi[n], tmp, isign);
	  ds_u += ds_tmp;
	}
	else
	{
	  tmp = b5[n]*psi[n];
	  D.deriv(ds_tmp, chi[n], tmp, isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectPlus(psi[n-1]);
	  tmp *= c5[n];
	  D.deriv(ds_tmp, chi[n], tmp, isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectMinus(psi[n+1]);
	  tmp *= c5[n];
	  D.deriv(ds_tmp, chi[n], tmp, isign);
	  ds_u += ds_tmp;
	}
      }
    }
    break;

    case MINUS:
    {
      for(int n=0; n < N5; ++n)
      {
	if (n == 0)
	{
	  tmp = b5[0]*chi[0];
	  D.deriv(ds_tmp, tmp, psi[0], isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectPlus(chi[0]);
	  tmp *= c5[1];
	  D.deriv(ds_tmp, tmp, psi[1], isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectMinus(chi[0]);
	  tmp *= -m_q*c5[N5-1];
	  D.deriv(ds_tmp, tmp, psi[N5-1], isign);
	  ds_u += ds_tmp;
	}
	else if (n == N5-1)
	{
	  tmp = b5[n]*chi[n];
	  D.deriv(ds_tmp, tmp, psi[n], isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectMinus(chi[n]);
	  tmp *= c5[n-1];
	  D.deriv(ds_tmp, tmp, psi[n-1], isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectPlus(chi[n]);
	  tmp *= -m_q*c5[0];
	  D.deriv(ds_tmp, tmp, psi[0], isign);
	  ds_u += ds_tmp;
	}
	else
	{
	  tmp = b5[n]*chi[n];
	  D.deriv(ds_tmp, tmp, psi[n], isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectMinus(chi[n]);
	  tmp *= c5[n-1];
	  D.deriv(ds_tmp, tmp, psi[n-1], isign);
	  ds_u += ds_tmp;

	  tmp = chiralProjectPlus(chi[n]);
	  tmp *= c5[n+1];
	  D.deriv(ds_tmp, tmp, psi[n+1], isign);
	  ds_u += ds_tmp;
	}
      }
    }
    break;
    }

    for(int m=0; m < Nd; ++m)
      ds_u[m] *= Real(-0.5);

    getFermBC().zero(ds_u);

    END_CODE();
  }

}; // End Namespace Chroma

