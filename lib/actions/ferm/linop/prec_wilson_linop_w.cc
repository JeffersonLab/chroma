// $Id: prec_wilson_linop_w.cc,v 1.3 2004-01-30 04:22:20 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"

//! Creation routine
/*!
 * \param u_ 	  gauge field     	       (Read)
 * \param Mass_   fermion kappa   	       (Read)
 */
void EvenOddPrecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				    const Real& Mass_)
{
  Mass = Mass_;
  D.create(u_);

  fact = Nd + Mass;
  invfact = 1/fact;
}


//! Creation routine with Anisotropy
/*!
 * \param u_ 	  gauge field     	       (Read)
 * \param Mass_   fermion kappa   	       (Read)
 * \param aniso   anisotropy struct   	       (Read)
 */
void EvenOddPrecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				    const Real& Mass_,
				    const AnisoParam_t& aniso)
{
  Mass = Mass_;

  multi1d<LatticeColorMatrix> u = u_;
  Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
  
  if (aniso.anisoP)
  {
    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      if (mu != aniso.t_dir)
	u[mu] *= ff;
    }
  }
  D.create(u);

  fact = 1 + (Nd-1)*ff + Mass;
  invfact = 1/fact;
}


//! Apply even-odd linop component
/*!
 * The operator acts on the entire even sublattice
 *
 * \param chi 	  Pseudofermion field     	       (Write)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void 
EvenOddPrecWilsonLinOp::evenOddLinOp(LatticeFermion& chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
{
  START_CODE("EvenOddPrecWilsonLinOp::evenOddLinOp");

  Real mhalf = -0.5;

  D.apply(chi, psi, isign, 0);
  chi[rb[0]] *= mhalf;
  
  END_CODE("EvenOddPrecWilsonLinOp::evenOddLinOp");
}

//! Apply odd-even linop component
/*!
 * The operator acts on the entire odd sublattice
 *
 * \param chi 	  Pseudofermion field     	       (Write)
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 */
void 
EvenOddPrecWilsonLinOp::oddEvenLinOp(LatticeFermion& chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
{
  START_CODE("EvenOddPrecWilsonLinOp::oddEvenLinOp");

  Real mhalf = -0.5;

  D.apply(chi, psi, isign, 1);
  chi[rb[1]] *= mhalf;
  
  END_CODE("EvenOddPrecWilsonLinOp::oddEvenLinOp");
}

