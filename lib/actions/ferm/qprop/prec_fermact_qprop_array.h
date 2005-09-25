// $Id: prec_fermact_qprop_array.h,v 2.0 2005-09-25 21:04:30 edwards Exp $
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
  template<typename T, typename P>
  class PrecFermAct5DQprop : public SystemSolver< multi1d<T> >
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    PrecFermAct5DQprop(Handle< const EvenOddPrecLinearOperator< multi1d<T>, P > > A_,
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
    PrecFermAct5DQprop(Handle< const EvenOddPrecLinearOperator< multi1d<T>, P > > A_,
		       Handle<const ConnectState> state_,  // throw away
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
    int operator() (multi1d<T>& psi, const multi1d<T>& chi) const;

  private:
    Handle< const EvenOddPrecLinearOperator< multi1d<T>, P > > A;
    const InvertParam_t invParam;
  };


}; // namespace Chroma


#endif
