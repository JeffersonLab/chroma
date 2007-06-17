// $Id: eoprec_clover_extfield_linop_w.cc,v 3.3 2007-06-17 02:25:16 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion linear operator in an external field
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/extfield_fermstate_w.h"
#include "actions/ferm/linop/eoprec_clover_extfield_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecCloverExtFieldLinOp::create(Handle< FermState<T,P,Q> > fs, 
					      const CloverFermActParams& param_)
  {
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    Handle< ExtFieldFermState<T,P,Q> > efs = fs.cast< ExtFieldFermState<T,P,Q> >();

    param = param_;

    clov.create(efs->getOriginalState(), param);
 
    invclov.create(efs->getOriginalState(),param, clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part

    D.create(efs->getU1State(), param.anisoParam);

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
  }

  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecCloverExtFieldLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					      enum PlusMinus isign) const
  {
    clov.apply(chi, psi, isign, 1);
  }


  //! Apply the the even-even block onto a source vector
  void 
  EvenOddPrecCloverExtFieldLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
						enum PlusMinus isign) const
  {
    // Nuke for testing
    clov.apply(chi, psi, isign, 0);
  }

  //! Apply the inverse of the even-even block onto a source vector
  void 
  EvenOddPrecCloverExtFieldLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
						   enum PlusMinus isign) const
  {
    invclov.apply(chi, psi, isign, 0);
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
  EvenOddPrecCloverExtFieldLinOp::evenOddLinOp(LatticeFermion& chi, 
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
  EvenOddPrecCloverExtFieldLinOp::oddEvenLinOp(LatticeFermion& chi, 
					       const LatticeFermion& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;
  
    END_CODE();
  }


  //! Apply even-odd preconditioned Clover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecCloverExtFieldLinOp::operator()(LatticeFermion & chi, 
						  const LatticeFermion& psi, 
						  enum PlusMinus isign) const
  {

    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);
    invclov.apply(tmp2, tmp1, isign, 0);
    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);


    chi[rb[1]] += mquarter*tmp1;
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecCloverExtFieldLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecCloverExtFieldLinOp::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }
} // End Namespace Chroma
