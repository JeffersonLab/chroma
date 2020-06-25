/*
 * projector_clover_mg_proto.h
 *
 *  Created on: Jun 8, 2020
 *      Author: eloy (so you know who to blame)
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_PROJECTOR_CLOVER_MG_PROTO_QPHIX_H_
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_PROJECTOR_CLOVER_MG_PROTO_QPHIX_H_

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/mg_proto/mgproto_solver_params.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "lattice/qphix/qphix_mgdeflation.h"
#include <stdexcept>
#include <memory>

using namespace QDP;

namespace Chroma {

//! Registration and other yuckies
  namespace ProjectorMGProtoQPhiXCloverEnv
  {
    //! Register the projector
    bool registerAll();


  }

  class ProjectorMGProtoQPhiXClover : public Projector<LatticeFermion>
  {
  public:
	  using T = LatticeFermion;
	  using Q = multi1d<LatticeColorMatrix>;
	  using Ts = const std::vector<std::shared_ptr<T>>;
	  using const_Ts = const std::vector<std::shared_ptr<const T>>;

	  ProjectorMGProtoQPhiXClover(Handle< LinearOperator<T> > A_,
			  Handle< FermState<T,Q,Q> > state_,
			  const MGProtoMGDeflationParams& deflation_params_);

	  ~ProjectorMGProtoQPhiXClover();

    	  //! Apply the oblique projector A*V*inv(U^H*A*V)*U^H
    	  /*! 
    	   * Returns A*V*inv(U^H*A*V)*U^H*chi = psi
    	   */
    	  void AVUObliqueProjector(Ts& psi, const_Ts& chi) const override;

    	  //! Apply the oblique projector V*inv(U^H*A*V)*U^H*A
    	  /*! 
    	   * Returns V*inv(U^H*A*V)*U^H*A*chi = psi
    	   */
    	  void VUAObliqueProjector(Ts& psi, const_Ts& chi) const override;

    	  //! Rank of the projector, which is the rank of U and V also
    	  unsigned int rank() const override;

    	  //! Return U[i]
    	  void U(unsigned int i, T& psi) const override;

     	  //! Return V[i]
    	  void V(unsigned int i, T& psi) const override;

    	  //! Return U[i]^H*A*V[i]
    	  void lambda(unsigned int i, DComplex& lambda) const override;

	  //! Return the subset on which the operator acts
	  const Subset& subset() const override { return _subset; }

  private:
	  void apply(Ts& psi, const_Ts& chi, bool do_VVA=true) const;
	  const Subset _subset;
	  const Handle< LinearOperator< T > > A;
	  const std::shared_ptr<MG::MGDeflation> deflation;
	  std::vector<std::complex<double>> diag_inv_VAV;
  };

};




#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_PROJECTOR_CLOVER_MG_PROTO_H_ */
