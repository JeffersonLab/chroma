/*
 * projector_clover_mg_proto.h
 *
 *  Created on: Jun 8, 2020
 *      Author: eloy (so you know who to blame)
 */

#include "chromabase.h"
#include "actions/ferm/invert/mg_proto/projector_clover_mg_proto_qphix.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/mg_proto/mg_proto_qphix_helpers.h"

#include "lattice/solver.h"
#include "lattice/fgmres_common.h"
#include "lattice/qphix/invfgmres_qphix.h"
#include "lattice/qphix/qphix_qdp_utils.h"
#include "lattice/qphix/qphix_clover_linear_operator.h"
#include "actions/ferm/invert/mg_solver_exception.h"

#include <memory>
#include <stdexcept>
using namespace QDP;

namespace Chroma
{
  namespace ProjectorMGProtoQPhiXCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("MG_PROTO_QPHIX_CLOVER_PROJECTOR");

      //! Local registration flag
      bool registered = false;
    }



    // Double precision
    Projector<LatticeFermion>* createProjector(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new ProjectorMGProtoQPhiXClover(A,state,MGProtoMGDeflationParams(xml_in, path));
    }

    //! Register all the factories
     bool registerAll()
     {
       bool success = true;
       if (! registered)
       {

    	   success &= Chroma::TheLinOpFermProjectorFactory::Instance().registerObject(name, createProjector);
    	   registered = true;
       }
       return success;
     }
  };

  // Save me typing, by exposing this file level from here
  using T = ProjectorMGProtoQPhiXClover::T;
  using Q = ProjectorMGProtoQPhiXClover::Q;
  using Ts = ProjectorMGProtoQPhiXClover::Ts;
  using const_Ts = ProjectorMGProtoQPhiXClover::const_Ts;

  // Constructor
  ProjectorMGProtoQPhiXClover::ProjectorMGProtoQPhiXClover(Handle< LinearOperator<T> > A_,
		  Handle< FermState<T,Q,Q> > state_,
		  const MGProtoMGDeflationParams& param_) :  A(A_), deflation(MGProtoHelpersQPhiX::createMGDeflation(param_, state_->getLinks())), _subset(A_->subset())
  {
	  diag_inv_VAV = deflation->GetDiagInvCoarse();
  }

  // Destructor
  ProjectorMGProtoQPhiXClover::~ProjectorMGProtoQPhiXClover(){}

  //! Apply the oblique projector A*V*inv(U^H*A*V)*U^H
  /*! 
   * Returns A*V*inv(U^H*A*V)*U^H*chi = psi
   */
  void ProjectorMGProtoQPhiXClover::AVUObliqueProjector(Ts& psi, const_Ts& chi) const {
	  apply(psi, chi, false);
  }

  //! Apply the oblique projector V*inv(U^H*A*V)*U^H*A
  /*! 
   * Returns V*inv(U^H*A*V)*U^H*A*chi = psi
   */
  void ProjectorMGProtoQPhiXClover::VUAObliqueProjector(Ts& psi, const_Ts& chi) const {
	  apply(psi, chi, true);
  }

  //! Rank of the projector, which is the rank of V also
  unsigned int ProjectorMGProtoQPhiXClover::rank() const {
	  return deflation->GetRank();
  }

  //! Return U[i]
  void ProjectorMGProtoQPhiXClover::U(unsigned int i, T& psi) const {
	  const LatticeInfo& info = deflation->GetInfo();
	  QPhiXSpinor qphix_out(info);
	  deflation->V(i, qphix_out);
	  QPhiXSpinorToQDPSpinor(qphix_out,0,psi);
  }

  //! Return V[i]
  void ProjectorMGProtoQPhiXClover::V(unsigned int i, T& psi) const {
	  const LatticeInfo& info = deflation->GetInfo();
	  QPhiXSpinor qphix_out(info);
	  deflation->g5V(i, qphix_out);
	  QPhiXSpinorToQDPSpinor(qphix_out,0,psi);
  }

  //! Return U_i^H*A*V_i
  void ProjectorMGProtoQPhiXClover::lambda(unsigned int i, DComplex& lambda) const {
	  assert(i < diag_inv_VAV.size());
	  std::complex<double> vav = 1./diag_inv_VAV[i];
	  lambda.elem().elem().elem() = RComplex<double>(std::real(vav), std::imag(vav));
  }


  void
  ProjectorMGProtoQPhiXClover::apply(Ts& psi, const_Ts& chi, bool do_VUA) const
  {
	  assert(psi.size() == chi.size());
	  int ncols = psi.size();

	  QDPIO::cout << "Jolly Greetings from Multigridland for deflation" << std::endl;
	  StopWatch swatch;
	  StopWatch swatch2;

	  swatch.reset();
	  swatch.start();

	  const LatticeInfo& info = deflation->GetInfo();
	  QPhiXSpinor qphix_in(info, ncols);
	  QPhiXSpinor qphix_out(info, ncols);

	  for (int col=0; col<ncols; ++col) QDPSpinorToQPhiXSpinor(*chi[col],qphix_in,col);
	  ZeroVec(qphix_out);

	  swatch2.reset();
	  swatch2.start();
	  if (do_VUA) {
	  	  deflation->VVA(qphix_out, qphix_in);
	  } else {
	  	  deflation->AVV(qphix_out, qphix_in);
	  }
	  swatch2.stop();

	  for (int col=0; col<ncols; ++col) {
	    *psi[col] = zero;
	    QPhiXSpinorToQDPSpinor(qphix_out,col,*psi[col]);
	  }


	  swatch.stop();
	  QDPIO::cout << "MG_PROTO_CLOVER_PROJECTOR_TIME: call_time = "<< swatch2.getTimeInSeconds() << " sec.  total_time=" << swatch.getTimeInSeconds() << " sec." << std::endl;
  }
};

