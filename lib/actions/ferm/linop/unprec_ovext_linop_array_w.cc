/* $Id: unprec_ovext_linop_array_w.cc,v 1.8 2005-03-02 16:27:15 bjoo Exp $
/*! \file
*  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
*/

#include "chromabase.h"
#include "actions/ferm/linop/unprec_ovext_linop_array_w.h"

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
  UnprecOvExtLinOpArray::create(const multi1d<LatticeColorMatrix>& u_, 
				const Real& WilsonMass_, const Real& m_q_, int N5_)
  {
    WilsonMass = WilsonMass_;
    m_q = m_q_;
    a5  = 1.0;
    N5  = N5_;

    if ((N5 & 1) == 0)
    {
      QDPIO::cerr << "UnprecOvExtLinOpArray: require odd N5" << endl;
      QDP_abort(1);
    }

    NN5 = (N5-1) >> 1;

    D.create(u_, -WilsonMass_);

    fact1 = -0.5*(1 + m_q);
    fact2 = sqrt((1 - m_q)/(2*NN5));
    cc.resize(NN5);
    ss.resize(NN5);

    // Use polar form for coeff. for now - can generalize
    for(int n=1; n <= NN5; ++n)
    {
      Real theta = twopi*(n - 0.5)/(4*NN5);
      cc[n-1] = pow(cos(theta),2);
      ss[n-1] = sin(theta);
    }
  }


  //! Apply the operator onto a source vector
  /*!
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  UnprecOvExtLinOpArray::operator() (multi1d<LatticeFermion>& chi,
				     const multi1d<LatticeFermion>& psi, 
				     enum PlusMinus isign) const
  {
    START_CODE();

    // chi.resize(N5);
    int G5 = Ns*Ns - 1;

    // Run through all the pseudofermion fields
    LatticeFermion  tmp;

    for(int n=0; n < N5; ++n)
    {
      if (n == 0)
      {
	chi[0] = fact1 * (Gamma(G5) * psi[0]);
	for(int s=1; s < N5; s+=2)
	  chi[0] += fact2 * psi[s];
      }
      else if ((n & 1) == 1)
      {
	int nn = (n-1) >> 1;
	D(tmp, psi[n], PLUS);
	chi[n] = fact2*psi[0] + (a5*cc[nn])*(Gamma(G5)*tmp) + ss[nn]*psi[n+1];
      }
      else
      {
	int nn = (n >> 1) - 1;
	D(tmp, psi[n], PLUS);
	chi[n] = ss[nn]*psi[n-1] - a5*(Gamma(G5)*tmp);
      }
    }

    END_CODE();
  }


  //! Derivative
  void 
  UnprecOvExtLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);
    ds_u = zero;

    int G5 = Ns*Ns - 1;

    // Run through all the pseudofermion fields
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    LatticeFermion  tmp;

    for(int n=0; n < N5; ++n)
    {
      if ((n & 1) == 1)
      {
	int nn = (n-1) >> 1;
	tmp = (a5*cc[nn])*(Gamma(G5)*chi[n]);
	D.deriv(ds_tmp, tmp, psi[n], PLUS);
      }
      else
      {
	int nn = (n >> 1) - 1;
	tmp = - a5*(Gamma(G5)*chi[n]);
	D.deriv(ds_tmp, tmp, psi[n], PLUS);
      }

      ds_u += ds_tmp;
    }

    END_CODE();
  }
}; // End Namespace Chroma

