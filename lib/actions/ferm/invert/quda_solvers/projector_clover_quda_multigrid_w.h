// -*- C++ -*-
/*! \file
 *  \QUDA MULTIGRID Clover solver.
 */

#ifndef __projector_quda_multigrid_clover_h__
#define __projector_quda_multigrid_clover_h__

#include "chroma_config.h"

#ifdef BUILD_MUGIQ
#  include "mugiq.h"
#  include <quda.h>

#  include "actions/ferm/fermbcs/simple_fermbc.h"
#  include "actions/ferm/fermstates/periodic_fermstate.h"
#  include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#  include "actions/ferm/linop/clover_term_w.h"
#  include "eigsolve_mugiq.h"
#  include "handle.h"
#  include "io/aniso_io.h"
#  include "linearop.h"
#  include "meas/gfix/temporal_gauge.h"
#  include "quda_mg_utils.h"
#  include "state.h"
#  include "syssolver.h"
#  include "util/gauge/reunit.h"
#  include <sstream>
#  include <string>
#  ifdef QDP_IS_QDPJIT
#    include "actions/ferm/invert/quda_solvers/qdpjit_memory_wrapper.h"
#  endif

//#include <util_quda.h>

namespace Chroma
{

  //! Richardson system solver namespace
  namespace ProjectorMugiqMULTIGRIDCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  }

  //! Project onto a Clover Fermion System using the QUDA inverter
  /*! \ingroup invert
   *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */

  class ProjectorMugiqMULTIGRIDClover : public Projector<LatticeFermion>,
					public QUDAMGUtils::AbstractQUDAClover
  {
  public:
    using T = LatticeFermion;
    using Q = multi1d<LatticeColorMatrix>;
    using Ts = const std::vector<std::shared_ptr<T>>;
    using const_Ts = const std::vector<std::shared_ptr<const T>>;


    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param deflation_params  deflation parameters ( Read )
     */
    ProjectorMugiqMULTIGRIDClover(Handle<LinearOperator<T>> A_, Handle<FermState<T, Q, Q>> state_,
				  const MugiqMGDeflationCloverParams& deflation_params_);

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
    const Subset& subset() const override
    {
      return A->subset();
    }

  private:
    void apply(Ts& psi, const_Ts& chi, bool do_VVA = true) const;
    void createMg(const MugiqMGDeflationCloverParams& param);
    void getVector(unsigned int i, T& psi, bool g5 = false) const;

    const Handle<LinearOperator<T>> A;
    QudaEigParam qudaEigParam;
    std::shared_ptr<mugiq::MugiqEigParam> mugiqEigParam;
    std::shared_ptr<mugiq::Eigsolve_Mugiq> mugiqEigsolver;
    std::shared_ptr<mugiq::MG_Mugiq> mugiqMg;
  };

} // End namespace

#endif // BUILD_MUGIQ
#endif
