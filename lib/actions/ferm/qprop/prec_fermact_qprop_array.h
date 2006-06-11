// -*- C++ -*-
// $Id: prec_fermact_qprop_array.h,v 3.1 2006-06-11 06:30:32 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#ifndef __prec_fermact_qprop_array_h__
#define __prec_fermact_qprop_array_h__

#include "chromabase.h"
#include "fermact.h"
#include "io/aniso_io.h"

namespace Chroma 
{
  //! Propagator of a generic even-odd preconditioned 5D fermion linear operator
  /*! \ingroup qprop
   *
   * This routine is actually generic to all even-odd 5D preconditioned fermions
   */
  template<typename T, typename P, typename Q>
  class PrecFermAct5DQprop : public SystemSolverArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    PrecFermAct5DQprop(Handle< EvenOddPrecLinearOperatorArray<T,P,Q > > A_,
		       const InvertParam_t& invParam_) : A(A_), invParam(invParam_) 
      {}

    //! Another constructor for compatibility
    /*!
     * This constructor exists so this class can be compatible with
     * optimized inverter constructors
     *
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    PrecFermAct5DQprop(Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > A_,
		       Handle< FermState<T,P,Q> > state_,  // throw away
		       const Real& OverMass_,    // throw away
		       const Real& Mass_,        // throw away
		       const AnisoParam_t& anisoParam_, // throw away
 		       const InvertParam_t& invParam_) : A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~PrecFermAct5DQprop() {}

    //! Expected length of array index
    int size() const {return A->size();}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const;

  private:
    Handle< EvenOddPrecLinearOperatorArray<T,P,Q> > A;
    const InvertParam_t invParam;
  };


}; // namespace Chroma


#endif
