// $Id: prec_asq_mdagm_s.cc,v 1.2 2003-12-11 17:11:17 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
#include "prec_asq_mdagm_s.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param _u_fat         fat7 links                        (Read)
 * \param _u_triple    triple links                        (Read)
 * \param _Mass        fermion mass   	                   (Read)
 */
void PrecAsqtadMdagM::create(const multi1d<LatticeColorMatrix>& _u_fat, const multi1d<LatticeColorMatrix>& _u_triple, const Real& _Mass)
{
  Mass = _Mass;
  u_fat = _u_fat;
  u_triple = _u_triple;
  D.create(u_fat,u_triple);
}


//! Apply Asqtad staggered fermion linear operator
/*!
 * \ingroup linop
 *
 * The operator acts on the checkerboarded lattice
 *
 * \param psi 	  Pseudofermion field     	       (Read)
 * \param isign   Flag ( PLUS | MINUS )   	       (Read)
 * \param cb      Checkerboard of OUTPUT VECTOR        (Read)
 */
void PrecAsqtadMdagM::operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const
{

  START_CODE("AsqtadMdagM");

  Real mass_sq = Mass*Mass;
  LatticeFermion tmp1, tmp2;

  //
  //  Chi     =  4m**2 Psi     -  D'  D'      Psi
  //     E                E        EO  OE   


  D.apply(tmp1, psi, isign, 1);
  D.apply(tmp2, tmp1, isign, 0);


  chi[rb[0]] = 4*mass_sq*psi - tmp2;
  
  END_CODE("AsqtadMdagM");

//  return chi;
}
