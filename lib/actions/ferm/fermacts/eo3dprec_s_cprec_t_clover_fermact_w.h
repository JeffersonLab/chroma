#ifndef EO3DPREC_S_CPREC_T_CLOVER_FERMACT_W_H
#define EO3DPREC_S_CPREC_T_CLOVER_FERMACT_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/aniso_io.h"
#include "eo3dprec_s_cprec_t_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma { 
  namespace EO3DPrecSpaceCentralPrecTimeCloverFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  class EO3DPrecSpaceCentralPrecTimeCloverFermAct :
    public EO3DPrecSpaceCentralPrecTimeConstDetWilsonTypeFermAct< 
    LatticeFermion, 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> >
  {
  public:
    typedef LatticeFermion T;
    typedef multi1d<LatticeColorMatrix> P;
    typedef multi1d<LatticeColorMatrix> Q;

    //! Destructor is automagic
    ~EO3DPrecSpaceCentralPrecTimeCloverFermAct(){}

    //! General Constructor
    EO3DPrecSpaceCentralPrecTimeCloverFermAct( Handle< CreateFermState<T,P,Q> > cfs_,
			 const Real& Mass_,
			 const AnisoParam_t& aniso_) : cfs(cfs_) {
      param.Mass = Mass_;
      param.anisoParam = aniso_;
    }

    //! Construct from Params
    EO3DPrecSpaceCentralPrecTimeCloverFermAct( Handle< CreateFermState<T,P,Q> > cfs_,
      const CloverFermActParams& param_ ) : cfs(cfs_), param(param_) {}

    // Copy Constructor
    EO3DPrecSpaceCentralPrecTimeCloverFermAct( const EO3DPrecSpaceCentralPrecTimeCloverFermAct&  a) :
      cfs(a.cfs), param(a.param) {}

    //! Assignment 
    EO3DPrecSpaceCentralPrecTimeCloverFermAct& operator=(const EO3DPrecSpaceCentralPrecTimeCloverFermAct& a) 
    {
      cfs = a.cfs; 
      param = a.param; 
      return *this;
    }
    
    //! Produce a linear Operator for this action
    EO3DPrecSpaceCentralPrecTimeLinearOperator<T,P,Q>* linOp( Handle< FermState<T,P,Q> > state) const;
   
    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

  protected:
    //! Return the fermion create state object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}
    
  private:
    EO3DPrecSpaceCentralPrecTimeCloverFermAct(){}   

  private:

    Handle< CreateFermState<T,P,Q> >  cfs;
    CloverFermActParams param;
  };



}

#endif
#endif
#endif

#endif
