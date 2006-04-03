// $Id: unprec_ovdwf_linop_array_w.cc,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param m_q_          quark mass    (Read)
   */
  void 
  UnprecOvDWLinOpArray::create(Handle< FermState<T,P,Q> > state,
			       const Real& WilsonMass_, const Real& m_q_, int N5_)
  {
    WilsonMass = WilsonMass_;
    m_q = m_q_;
    a5  = 1.0;
    N5  = N5_;

    D.create(state);
//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
  }


  //! Apply the operator onto a source vector
  /*!
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  UnprecOvDWLinOpArray::operator() (multi1d<LatticeFermion>& chi,
				    const multi1d<LatticeFermion>& psi, 
				    enum PlusMinus isign) const
  {
    START_CODE();

    if( chi.size() != N5 )  chi.resize(N5);

    //
    //  Chi   =  D' Psi
    //
    Real fact1 = a5*(Nd - WilsonMass);
    Real fact2 = -0.5*a5;

    if (isign == PLUS)
    {
      LatticeFermion  tmp1, tmp2;
      moveToFastMemoryHint(tmp1);
      moveToFastMemoryHint(tmp2);

      for(int n=0; n < N5; ++n)
      {
	if (n == 0)
	{
	  tmp1   = psi[n] - m_q*chiralProjectPlus(psi[N5-1]) + chiralProjectMinus(psi[1]);
	  D(tmp2, tmp1, isign);
	  chi[n] = fact1*tmp1 + fact2*tmp2 + psi[n] 
	    + m_q*chiralProjectPlus(psi[N5-1]) - chiralProjectMinus(psi[1]);
	}
	else if (n == N5-1)
	{
	  tmp1   = psi[n] + chiralProjectPlus(psi[N5-2]) - m_q*chiralProjectMinus(psi[0]);
	  D(tmp2, tmp1, isign);
	  chi[n] = fact1*tmp1 + fact2*tmp2 + psi[n] 
	    - chiralProjectPlus(psi[N5-2]) + m_q*chiralProjectMinus(psi[0]);
	}
	else
	{
	  tmp1   = psi[n] + chiralProjectPlus(psi[n-1]) + chiralProjectMinus(psi[n+1]);
	  D(tmp2, tmp1, isign);
	  chi[n] = fact1*tmp1 + fact2*tmp2 + psi[n] 
	    - chiralProjectPlus(psi[n-1]) - chiralProjectMinus(psi[n+1]);
	}
      }          
    }
    else   // isign = MINUS   case
    {
      multi1d<LatticeFermion>  tmp(N5);   // should be more clever and reduce temporaries
      LatticeFermion  tmp1;
      moveToFastMemoryHint(tmp);
      moveToFastMemoryHint(tmp1);

      for(int n=0; n < N5; ++n)
      {
	D(tmp1, psi[n], isign);
	tmp[n] = fact1*psi[n] + fact2*tmp1;
      }

      for(int n=0; n < N5; ++n)
      {
	if (n == 0)
	{
	  chi[0] = tmp[0] + psi[0] + chiralProjectPlus(tmp[1]) - chiralProjectPlus(psi[1])
	    - m_q*(chiralProjectMinus(tmp[N5-1]) - chiralProjectMinus(psi[N5-1]));
	}
	else if (n == N5-1)
	{
	  chi[n] = tmp[n] + psi[n] + chiralProjectMinus(tmp[N5-2]) - chiralProjectMinus(psi[N5-2])
	    - m_q*(chiralProjectPlus(tmp[0]) - chiralProjectPlus(psi[0]));
	}
	else
	{
	  chi[n] = tmp[n] + psi[n] + chiralProjectMinus(tmp[n-1]) - chiralProjectMinus(psi[n-1])
	    + chiralProjectPlus(tmp[n+1]) - chiralProjectPlus(psi[n+1]);
	}
      }          
    }

    getFermBC().modifyF(chi);

    END_CODE();
  }

  //! Apply the Dminus operator on a lattice fermion. See my notes ;-)
  void 
  UnprecOvDWLinOpArray::Dminus(LatticeFermion& chi,
			       const LatticeFermion& psi,
			       enum PlusMinus isign,
			       int s5) const
  {
    LatticeFermion tt ; moveToFastMemoryHint(tt);
    D.apply(tt,psi,isign,0);
    D.apply(tt,psi,isign,1);
    chi = (1.0 - (Nd-WilsonMass))*psi +0.5*tt ; //really -(-.5)D
  }
  

} // End Namespace Chroma

