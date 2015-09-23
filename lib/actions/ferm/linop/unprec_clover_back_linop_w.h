// -*- C++ -*-
/*! \file
 *  \brief Unpreconditioned Clover fermion linear operator
 */

#ifndef __unprec_clover_back_linop_w_h__
#define __unprec_clover_back_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "actions/ferm/fermacts/clover_back_fermact_params_w.h"



namespace Chroma 
{ 
  //! Unpreconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  
  class UnprecCloverBackLinOp : public  UnprecCloverLinOp
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecCloverBackLinOp() {}

    //! Full constructor
    UnprecCloverBackLinOp(Handle< FermState<T,P,Q> > fs,
			  const CloverBackFermActParams& param_);

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
                const CloverBackFermActParams& param_);

    //! Apply the operator onto a source std::vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;
    
    unsigned long nFlops() const ;
    

    //! Destructor is automatic
    ~UnprecCloverBackLinOp() {}

  private:
    int gamma ;
    Complex lambda ;
  };

} // End Namespace Chroma


#endif
