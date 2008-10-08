#ifndef ILU2PREC_S_CPREC_T_CLOVER_FERMACT_W_H
#define ILU2PREC_S_CPREC_T_CLOVER_FERMACT_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/aniso_io.h"
#include "iluprec_s_cprec_t_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma { 
  namespace ILU2PrecSpaceCentralPrecTimeCloverFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  class ILU2PrecSpaceCentralPrecTimeCloverFermAct :
    public ILU2PrecSpaceCentralPrecTimeWilsonTypeFermAct< 
    LatticeFermion, 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> >
  {
  public:
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;

    //! Destructor is automagic
    ~ILU2PrecSpaceCentralPrecTimeCloverFermAct(){}


    //! Construct from Params
    ILU2PrecSpaceCentralPrecTimeCloverFermAct( Handle< CreateFermState<T,P,Q> > cfs_,
      const CloverFermActParams& param_ ) : cfs(cfs_), param(param_) {}

    // Copy Constructor
    ILU2PrecSpaceCentralPrecTimeCloverFermAct( const ILU2PrecSpaceCentralPrecTimeCloverFermAct&  a) :
      cfs(a.cfs), param(a.param) {}

    //! Assignment 
    ILU2PrecSpaceCentralPrecTimeCloverFermAct& operator=(const ILU2PrecSpaceCentralPrecTimeCloverFermAct& a) 
    {
      cfs = a.cfs; 
      param = a.param; 
      return *this;
    }
    
    //! Produce a linear Operator for this action
    ILU2PrecSpaceCentralPrecTimeLinearOperator<T,P,Q>* linOp( Handle< FermState<T,P,Q> > state) const;
   
    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

  protected:
    //! Return the fermion create state object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}
    
  private:
    ILU2PrecSpaceCentralPrecTimeCloverFermAct(){}   

  private:

    Handle< CreateFermState<T,P,Q> >  cfs;
    CloverFermActParams param;
  };



}


#endif

#endif
#endif
#endif
