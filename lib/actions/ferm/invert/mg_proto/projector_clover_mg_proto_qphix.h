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

	  ProjectorMGProtoQPhiXClover(Handle< LinearOperator<T> > A_,
			  Handle< FermState<T,Q,Q> > state_,
			  const MGProtoMGDeflationParams& deflation_params_);

	  ~ProjectorMGProtoQPhiXClover();

    	  //! Apply the orthonormal projector
    	  /*! 
    	   * Returns   V*V^H*chi = psi at some accuracy.
    	   */
    	  void orthonormalProjector(T& psi, const T& chi) const;

    	  //! Apply the oblique projector A*V*inv(V^H*A*V)*V^H
    	  /*! 
    	   * Returns A*V*inv(V^H*A*V)*V^H*chi = psi
    	   */
    	  void AVVObliueProjector(T& psi, const T& chi) const;

    	  //! Apply the oblique projector V*inv(V^H*A*V)*V^H*A
    	  /*! 
    	   * Returns V*inv(V^H*A*V)*V^H*A*chi = psi
    	   */
    	  void VVAObliueProjector(T& psi, const T& chi) const;

    	  //! Rank of the projector, which is the rank of V also
    	  unsigned int rank() const;

    	  //! Return v_i
    	  void V(unsigned int i, T& psi) const;

    	  //! Return v_i^H*A*V_i
    	  void lambda(unsigned int i, DComplex& lambda) const;

	  //! Return the subset on which the operator acts
	  const Subset& subset() const { return _subset; }

  private:
	  void apply(T& psi, const T& chi, bool do_VVA=true) const;
	  const Subset _subset;
	  const Handle< LinearOperator< T > > A;
	  const std::shared_ptr<MG::MGDeflation> deflation;
	  std::vector<std::complex<double>> diag_inv_VAV;
  };

};




#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_PROJECTOR_CLOVER_MG_PROTO_H_ */
