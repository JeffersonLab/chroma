/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_clover_back_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine with Anisotropy
  /*!
   * \param fs 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecCloverBackLinOp::create(Handle< FermState<T,P,Q> > fs,
				 const CloverBackFermActParams& param_)
  {
    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    CloverFermActParams* foo= new CloverFermActParams;
    *foo = param_ ;

    UnprecCloverLinOp::create(fs,*foo);
   
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
  }



  //! Apply unpreconditioned Clover fermion linear operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field          (Write)
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecCloverBackLinOp::operator()(LatticeFermion & chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
  {
       
    UnprecCloverLinOp::operator()(chi,psi,isign) ;
    LatticeFermion tmp; moveToFastMemoryHint(tmp);
    chi +=  lambda*(Gamma(gamma)*chi) ;
    
  }

  //! Return flops performed by the operator()
  unsigned long UnprecCloverBackLinOp::nFlops() const
  {
    unsigned long site_flops= UnprecCloverBackLinOp::nFlops();
    // add the site flops for the vector addition 
    site_flops += Ns*Nc*5*Layout::sitesOnNode(); 
    return site_flops ;
  }

} // End Namespace Chroma
