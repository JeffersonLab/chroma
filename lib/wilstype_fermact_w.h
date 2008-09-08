// -*- C++ -*-
// $Id: wilstype_fermact_w.h,v 3.4 2008-09-08 16:00:19 bjoo Exp $

/*! @file
 * @brief Wilson-like fermion actions
 */

#ifndef __wilstype_fermact_w_h__
#define __wilstype_fermact_w_h__

#include "fermact.h"
#include "eoprec_linop.h"

namespace Chroma
{

  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class WilsonTypeFermAct : public DiffFermAct4D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~WilsonTypeFermAct() {}

    //! Produce a linear operator M^dag.M for this action
    /*! Default implementation */
    virtual DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const
      {
	return new DiffMdagMLinOp<T,P,Q>(this->linOp(state));
      }

    //! Produce a hermitian version of the linear operator
    virtual LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    /*! Default implementation */
    virtual LinOpSystemSolver<T>* invLinOp(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    /*! Default implementation */
    virtual MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
					   const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (M+shift)*psi=chi 
    /*! Default implementation */
    virtual LinOpMultiSystemSolver<T>* mInvLinOp(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    /*! Default implementation */
    virtual MdagMMultiSystemSolver<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const;

    /*! Default implementation */
    virtual MdagMMultiSystemSolverAccumulate<T>* mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
						 const GroupXML_t& invParam) const;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol         quark propagator ( Write )
     * \param q_src         source ( Read )
     * \param xml_out       diagnostic output ( Modify )
     * \param state         gauge connection state ( Read )
     * \param invParam      inverter parameters ( Read )
     * \param quarkSpinType compute only a non-relativistic prop ( Read )
     * \param ncg_had       number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam,
			   QuarkSpinType quarkSpinType,
			   int& ncg_had) const;
  };


  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class WilsonTypeFermAct5D : public DiffFermAct5D<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~WilsonTypeFermAct5D() {}

    //! Produce a linear operator M^dag.M for this action
    /*! Default implementation */
    virtual DiffLinearOperatorArray<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const
      {
	return new DiffMdagMLinOpArray<T,P,Q>(this->linOp(state));
      }

    //! Produce a hermitian version of the linear operator
    virtual LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const = 0;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    /*! Default implementation provided */
    virtual LinOpSystemSolverArray<T>* invLinOp(Handle< FermState<T,P,Q> > state,
						const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    /*! Default implementation provided */
    virtual MdagMSystemSolverArray<T>* invMdagM(Handle< FermState<T,P,Q> > state,
						const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve PV*psi=chi 
    /*! 
     * Default implementation provided
     *
     * Do we need this critter? 
     */
    virtual LinOpSystemSolverArray<T>* invLinOpPV(Handle< FermState<T,P,Q> > state,
						  const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve PV^dag*PV*psi=chi 
    /*! Default implementation provided */
    virtual MdagMSystemSolverArray<T>* invMdagMPV(Handle< FermState<T,P,Q> > state,
						  const GroupXML_t& invParam) const;
 
    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    /*! Default implementation provided */
    virtual MdagMMultiSystemSolverArray<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
						      const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    /*! Default implementation provided */
    virtual MdagMMultiSystemSolverAccumulateArray<T>* mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
						      const GroupXML_t& invParam) const;


    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
    /*! Default implementation provided */
    virtual MdagMMultiSystemSolverArray<T>* mInvMdagMPV(Handle< FermState<T,P,Q> > state,
						      const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (PV^dag*PV+shift)*psi=chi 
    /*! Default implementation provided */
    virtual MdagMMultiSystemSolverAccumulateArray<T>* mInvMdagMPVAcc(Handle< FermState<T,P,Q> > state,
							const GroupXML_t& invParam) const;

    //! Produce an unpreconditioned linear operator projecting 5D to 4D (the inverse of qprop below)
    virtual LinearOperator<T>* linOp4D(Handle< FermState<T,P,Q> > state,
				       const Real& m_q,
				       const GroupXML_t& invParam) const = 0;

    //! Produce a  DeltaLs = 1-epsilon^2(H) operator
    virtual LinearOperator<T>* DeltaLs(Handle<  FermState<T,P,Q> > state,
				       const GroupXML_t& invParam) const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol         quark propagator ( Write )
     * \param q_src         source ( Read )
     * \param xml_out       diagnostic output ( Modify )
     * \param state         gauge connection state ( Read )
     * \param invParam      inverter parameters ( Read )
     * \param quarkSpinType compute only a non-relativistic prop ( Read )
     * \param ncg_had       number of solver iterations ( Write )
     */
    virtual void quarkProp(typename PropTypeTraits<T>::Type_t& q_sol,
			   XMLWriter& xml_out,
			   const typename PropTypeTraits<T>::Type_t& q_src,
			   Handle< FermState<T,P,Q> > state,
			   const GroupXML_t& invParam,
			   QuarkSpinType quarkSpinType,
			   int& ncg_had) const;
  };







}

#endif
