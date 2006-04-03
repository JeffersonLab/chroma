// $Id: asqtad_mdagm_s.cc,v 3.0 2006-04-03 04:58:49 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */
// NEW $Id: asqtad_linop.cc 2003/11/13 steve

#include "chromabase.h"
#include "actions/ferm/linop/asqtad_mdagm_s.h"

using namespace QDP::Hints;
namespace Chroma 
{ 
  //! Creation routine
  /*! \ingroup fermact
   *
   * \param _u_fat         fat7 links                        (Read)
   * \param _u_triple    triple links                        (Read)
   * \param _Mass        fermion mass   	                   (Read)
   */
  void AsqtadMdagM::create(Handle<AsqtadConnectStateBase> state,
			   const Real& Mass_)
  {
    Mass = Mass_;

    // We wire them into D so we don't need to keep them
    // u_fat = _u_fat;
    // u_triple = _u_triple;
    D.create(state);
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
  void AsqtadMdagM::operator()(LatticeStaggeredFermion& chi, 
			       const LatticeStaggeredFermion& psi, 
			       enum PlusMinus isign) const
  {
    START_CODE();

    Real mass_sq = Mass*Mass;
    LatticeStaggeredFermion tmp1, tmp2;
    moveToFastMemoryHint(tmp1);
    moveToFastMemoryHint(tmp2);

    tmp1 = tmp2 = zero;

    //
    //  Chi     =  4m**2 Psi     -  D'  D'      Psi
    //     E                E        EO  OE   

    D.apply(tmp1, psi, isign, 1);
    D.apply(tmp2, tmp1, isign, 0);

    chi[rb[0]] = 4*mass_sq*psi - tmp2;
  
    END_CODE();
  }

} // End Namespace Chroma

