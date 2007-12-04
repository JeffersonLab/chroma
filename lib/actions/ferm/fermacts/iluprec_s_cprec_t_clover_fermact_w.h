#ifndef ILUPREC_S_CPREC_T_CLOVER_FERMACT_W_H
#define ILUPREC_S_CPREC_T_CLOVER_FERMACT_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/aniso_io.h"
#include "iluprec_s_cprec_t_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma { 
  namespace ILUPrecSpaceCentralPrecTimeCloverFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  class ILUPrecSpaceCentralPrecTimeCloverFermAct :
    public ILUPrecSpaceCentralPrecTimeWilsonTypeFermAct< 
    LatticeFermion, 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> >
  {
  public:
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;

    //! Destructor is automagic
    ~ILUPrecSpaceCentralPrecTimeCloverFermAct(){}


    //! Construct from Params
    ILUPrecSpaceCentralPrecTimeCloverFermAct( Handle< CreateFermState<T,P,Q> > cfs_,
      const CloverFermActParams& param_ ) : cfs(cfs_), param(param_) {}

    // Copy Constructor
    ILUPrecSpaceCentralPrecTimeCloverFermAct( const ILUPrecSpaceCentralPrecTimeCloverFermAct&  a) :
      cfs(a.cfs), param(a.param) {}

    //! Assignment 
    ILUPrecSpaceCentralPrecTimeCloverFermAct& operator=(const ILUPrecSpaceCentralPrecTimeCloverFermAct& a) 
    {
      cfs = a.cfs; 
      param = a.param; 
      return *this;
    }
    
    //! Produce a linear Operator for this action
    ILUPrecSpaceCentralPrecTimeLinearOperator<T,P,Q>* linOp( Handle< FermState<T,P,Q> > state) const;
   
    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

  protected:
    //! Return the fermion create state object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}
    
  private:
    ILUPrecSpaceCentralPrecTimeCloverFermAct(){}   

  private:

    Handle< CreateFermState<T,P,Q> >  cfs;
    CloverFermActParams param;
  };



}


#endif

#endif
#endif
#endif
