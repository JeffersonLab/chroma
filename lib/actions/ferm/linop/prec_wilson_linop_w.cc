// $Id: prec_wilson_linop_w.cc,v 1.7 2004-07-28 02:38:02 edwards Exp $
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
  START_CODE();

  Real mhalf = -0.5;

  D.apply(chi, psi, isign, 0);
  chi[rb[0]] *= mhalf;
  
  END_CODE();
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
  START_CODE();

  Real mhalf = -0.5;

  D.apply(chi, psi, isign, 1);
  chi[rb[1]] *= mhalf;
  
  END_CODE();
}

void EvenOddPrecWilsonLinOp::operator()(LatticeFermion & chi, 
					 const LatticeFermion& psi, 
					 enum PlusMinus isign) const
{
  LatticeFermion tmp1, tmp2, tmp3;  // if an array is used here, 
                                    // the space is not reserved
  

  Real mquarterinvfact = -0.25*invfact;

  // tmp1[0] = D_eo psi[1]
  D.apply(tmp1, psi, isign, 0);

  // tmp2[1] = D_oe tmp1[0]
  D.apply(tmp2, tmp1, isign, 1);

  // Now we have tmp2[1] = D_oe D_eo psi[1]

  // now scale tmp2[1] with (-1/4)/fact = (-1/4)*(1/(Nd + m))
  // with a vscale -- using tmp2 on both sides should be OK, but
  // just to be safe use tmp3
  tmp3[rb[1]] = mquarterinvfact*tmp2;

  // now tmp3[1] should be = (-1/4)*(1/(Nd + m) D_oe D_eo psi[1]

  // Now get chi[1] = fact*psi[1] + tmp3[1]
  // with a vaxpy3 
  // chi[1] = (Nd + m) - (1/4)*(1/(Nd + m)) D_oe D_eo psi[1]
  //
  // in this order, this last job could be replaced with a 
  // vaxpy3_norm if we wanted the || M psi ||^2 over the subset.
  chi[rb[1]] = fact*psi + tmp3;
}

 
