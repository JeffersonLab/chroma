#ifndef UNPREC_S_CPREC_T_WILSON_FERMACT_W_H
#define UNPREC_S_CPREC_T_WILSON_FERMACT_W_H

#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "io/aniso_io.h"
#include "unprec_s_cprec_t_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma { 
  namespace UnprecSpaceCentralPrecTimeWilsonFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  class UnprecSpaceCentralPrecTimeWilsonFermAct :
    public UnprecSpaceCentralPrecTimeWilsonTypeFermAct< 
    LatticeFermion, 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> >
  {
  public:
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;

    //! Destructor is automagic
    ~UnprecSpaceCentralPrecTimeWilsonFermAct(){}

    //! General Constructor
    UnprecSpaceCentralPrecTimeWilsonFermAct( Handle< CreateFermState<T,P,Q> > cfs_,
			 const Real& Mass_,
			 const AnisoParam_t& aniso_) : cfs(cfs_) {
      param.Mass = Mass_;
      param.anisoParam = aniso_;
    }

    //! Construct from Params
    UnprecSpaceCentralPrecTimeWilsonFermAct( Handle< CreateFermState<T,P,Q> > cfs_,
      const WilsonFermActParams& param_ ) : cfs(cfs_), param(param_) {}

    // Copy Constructor
    UnprecSpaceCentralPrecTimeWilsonFermAct( const UnprecSpaceCentralPrecTimeWilsonFermAct&  a) :
      cfs(a.cfs), param(a.param) {}

    //! Assignment 
    UnprecSpaceCentralPrecTimeWilsonFermAct& operator=(const UnprecSpaceCentralPrecTimeWilsonFermAct& a) 
    {
      cfs = a.cfs; 
      param = a.param; 
      return *this;
    }
    
    //! Produce a linear Operator for this action
    UnprecSpaceCentralPrecTimeLinearOperator<T,P,Q>* linOp( Handle< FermState<T,P,Q> > state) const;
   
    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

  protected:
    //! Return the fermion create state object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}
    
  private:
    UnprecSpaceCentralPrecTimeWilsonFermAct(){}   

  private:

    Handle< CreateFermState<T,P,Q> >  cfs;
    WilsonFermActParams param;
  };



}

#endif 
#endif
#endif

#endif
