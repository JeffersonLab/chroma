#ifndef CENTRAL_TPREC_FERMACT_W_H
#define CENTRAL_TPREC_FERMACT_W_H

#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "wilstype_fermact_w.h"
#include "central_tprec_linop.h"

namespace Chroma { 

  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class CentralTimePrecFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~CentralTimePrecFermAct() {}

    virtual CentralTimePrecLinearOperator<T,P,Q>* linOp( Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    /*! Default implementation */
    virtual DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const
      {
	return new DiffMdagMLinOp<T,P,Q>(this->linOp(state));
      }

  };

  //-------------------------------------------------------------------------------------------
  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */
  template<typename T, typename P, typename Q>
  class Central2TimePrecFermAct : public WilsonTypeFermAct<T,P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~Central2TimePrecFermAct() {}

    virtual Central2TimePrecLinearOperator<T,P,Q>* linOp( Handle< FermState<T,P,Q> > state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    /*! Default implementation */
    virtual DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const
      {
	return new DiffMdagMLinOp<T,P,Q>(this->linOp(state));
      }

  };




}

#endif
#endif
#endif


#endif

