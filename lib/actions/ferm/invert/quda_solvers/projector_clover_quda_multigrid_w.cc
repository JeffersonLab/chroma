/*! \file
 *  \QUDA MULTIGRID Clover solver.
 */
// comment
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/invert/quda_solvers/projector_clover_quda_multigrid_w.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "io/aniso_io.h"

#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "handle.h"
#include "meas/glue/mesplq.h"
// QUDA Headers
#include <memory>
#include <stdexcept>
#include <quda.h>
// #include <util_quda.h>
#include "actions/ferm/invert/quda_solvers/quda_mg_utils.h"

// Quda function from lib/interface_quda.cpp
extern quda::cudaGaugeField* checkGauge(QudaInvertParam* param);

namespace Chroma
{
  namespace ProjectorMugiqMULTIGRIDCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QUDA_MULTIGRID_CLOVER_PROJECTOR");

      //! Local registration flag
      bool registered = false;
    }

    // Save me typing, by exposing this file level from here
    using T = ProjectorMugiqMULTIGRIDClover::T;
    using Q = ProjectorMugiqMULTIGRIDClover::Q;
    using Ts = ProjectorMugiqMULTIGRIDClover::Ts;
    using const_Ts = ProjectorMugiqMULTIGRIDClover::const_Ts;

    Projector<LatticeFermion>* createProjector(
      XMLReader& xml_in, const std::string& path,
      Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>
	state,
      Handle<LinearOperator<LatticeFermion>> A)
    {
      return new ProjectorMugiqMULTIGRIDClover(A, state,
					       MugiqMGDeflationCloverParams(xml_in, path));
    }

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
	success &=
	  Chroma::TheLinOpFermProjectorFactory::Instance().registerObject(name, createProjector);
	registered = true;
      }
      return success;
    }
  }


  // Constructor
  ProjectorMugiqMULTIGRIDClover::ProjectorMugiqMULTIGRIDClover(
    Handle<LinearOperator<T>> A, Handle<FermState<T, Q, Q>> state,
    const MugiqMGDeflationCloverParams& param)
    : QUDAMGUtils::AbstractQUDAClover(state, param), A(A)
  {
    createMg(param);
  }

  quda::TimeProfile profileEigensolveMuGiq("computeEvecsMuGiq");

  void ProjectorMugiqMULTIGRIDClover::createMg(const MugiqMGDeflationCloverParams& param)
  {

    QudaEigParam eig_param = newQudaEigParam();
    QudaInvertParam eig_inv_param = newQudaInvertParam();
    eig_inv_param = quda_inv_param;
    eig_inv_param.input_location = QUDA_CUDA_FIELD_LOCATION;
    eig_inv_param.output_location = QUDA_CUDA_FIELD_LOCATION;
    eig_inv_param.cuda_prec = QUDA_SINGLE_PRECISION;
    
    eig_param.invert_param = &eig_inv_param;

    eig_param.eig_type = QUDA_EIG_TR_LANCZOS;
    eig_param.spectrum = QUDA_SPECTRUM_SR_EIG;

    eig_param.nConv = param.EigenSolverMaxRank;

    eig_param.nEv = param.EigenSolverMaxRank;
    eig_param.nKr = param.EigenSolverMaxRestartSize;
    eig_param.tol =
      param.RsdTarget.elem().elem().elem().elem().elem(); // for the love of some god...
    eig_param.batched_rotate = 0;
    eig_param.require_convergence = QUDA_BOOLEAN_YES;
    eig_param.check_interval = 10;
    eig_param.max_restarts = 100;
    eig_param.cuda_prec_ritz = QUDA_DOUBLE_PRECISION;

    eig_param.use_norm_op = QUDA_BOOLEAN_YES;
    eig_param.use_dagger =  QUDA_BOOLEAN_NO;
    eig_param.compute_svd = QUDA_BOOLEAN_YES;

    eig_param.use_poly_acc = QUDA_BOOLEAN_NO;
    eig_param.poly_deg = 0;
    eig_param.a_min = 0;
    eig_param.a_max = 1;
    eig_param.arpack_check = QUDA_BOOLEAN_NO;

    mugiqMg =
      std::make_shared<mugiq::MG_Mugiq>(&subspace_pointers->mg_param, profileEigensolveMuGiq);
    //- Create the eigensolver environment
    mugiqEigParam = std::make_shared<mugiq::MugiqEigParam>(&eig_param);
    mugiqEigsolver = std::make_shared<mugiq::Eigsolve_Mugiq>(mugiqEigParam.get(), mugiqMg.get(),
							     &profileEigensolveMuGiq);

    mugiqEigsolver->printInfo();

    //- Compute eigenvectors and (local) eigenvalues
    mugiqEigsolver->computeEvecs();
    mugiqEigsolver->computeEvals();
    //mugiqEigsolver->printEvals();
  }

  //! Apply the oblique projector A*V*inv(U^H*A*V)*U^H
  /*! 
   * Returns A*V*inv(U^H*A*V)*U^H*chi = psi
   */
  void ProjectorMugiqMULTIGRIDClover::AVUObliqueProjector(Ts& psi, const_Ts& chi) const
  {
    apply(psi, chi, false);
  }

  //! Apply the oblique projector V*inv(U^H*A*V)*U^H*A
  /*! 
   * Returns V*inv(U^H*A*V)*U^H*A*chi = psi
   */
  void ProjectorMugiqMULTIGRIDClover::VUAObliqueProjector(Ts& psi, const_Ts& chi) const
  {
    apply(psi, chi, true);
  }

  //! Rank of the projector, which is the rank of V also
  unsigned int ProjectorMugiqMULTIGRIDClover::rank() const
  {
    return mugiqEigParam->nEv;
  }

  //! Return U[i]
  void ProjectorMugiqMULTIGRIDClover::U(unsigned int i, T& psi) const
  {
    getVector(i, psi, true);
  }

  //! Return V[i]
  void ProjectorMugiqMULTIGRIDClover::V(unsigned int i, T& psi) const
  {
    getVector(i, psi, false);
  }

  //! Return U_i^H*A*V_i
  void ProjectorMugiqMULTIGRIDClover::lambda(unsigned int i, DComplex& lambda) const
  {
    assert(i < rank());
    std::complex<double> vav = (*mugiqEigsolver->getEvals())[i];
    lambda.elem().elem().elem() = RComplex<double>(std::real(vav), std::imag(vav));
  }

  void ProjectorMugiqMULTIGRIDClover::apply(Ts& psi, const_Ts& chi, bool do_VUA) const
  {
    assert(psi.size() == chi.size());
    int ncols = psi.size();

    if (invParam.asymmetricP)
      throw std::runtime_error("Unsupported asymmetricP");
    if (!do_VUA)
      throw std::runtime_error("Mugiq isn't supporting AVU projector");

    StopWatch swatch;
    StopWatch swatch2;

    swatch.reset();
    swatch.start();

    T mod_chi;

    /// Copy from function invertQuda in quda (lib/interface_quda.cpp)

    // check the gauge fields have been created???
    const int* X = checkGauge(const_cast<QudaInvertParam*>(&quda_inv_param))->X();

    bool pc_solution = (quda_inv_param.solution_type == QUDA_MATPC_SOLUTION) ||
		       (quda_inv_param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);

    for (int i = 0; i < ncols; ++i)
    {
      // Copy source into mod_chi, and zero the off-parity
      mod_chi[rb[0]] = zero;
      mod_chi[rb[1]] = *chi[i];

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
      void* spinorIn = (void*)&(mod_chi.elem(rb[1].start()).elem(0).elem(0).real());
      void* spinorOut = (void*)&(psi[i]->elem(rb[1].start()).elem(0).elem(0).real());
#else
      void* spinorIn;
      void* spinorOut;
      GetMemoryPtr2(spinorIn, spinorOut, mod_chi.getId(), psi[i]->getId())
#endif

      // wrap CPU host side pointers
      quda::ColorSpinorParam cpuParam(spinorIn, *const_cast<QudaInvertParam*>(&quda_inv_param), X,
				      pc_solution, quda_inv_param.input_location);
      quda::ColorSpinorField* h_in = quda::ColorSpinorField::Create(cpuParam);

      cpuParam.v = spinorOut;
      cpuParam.location = quda_inv_param.output_location;
      quda::ColorSpinorField* h_out = quda::ColorSpinorField::Create(cpuParam);

      mugiqEigsolver->projectVector(*h_out, *h_in);

      delete h_in;
      delete h_out;
    }

    swatch.stop();
    QDPIO::cout << "MG_PROTO_CLOVER_PROJECTOR_TIME: call_time = " << swatch2.getTimeInSeconds()
		<< " sec.  total_time=" << swatch.getTimeInSeconds() << " sec." << std::endl;
  }

  void ProjectorMugiqMULTIGRIDClover::getVector(unsigned int i, T& psi, bool g5) const
  {
    /// Copy from function invertQuda in quda (lib/interface_quda.cpp)

    // check the gauge fields have been created???
    const int* X = checkGauge(const_cast<QudaInvertParam*>(&quda_inv_param))->X();

    bool pc_solution = (quda_inv_param.solution_type == QUDA_MATPC_SOLUTION) ||
		       (quda_inv_param.solution_type == QUDA_MATPCDAG_MATPC_SOLUTION);

#ifndef BUILD_QUDA_DEVIFACE_SPINOR
    void* spinorOut = (void*)&(psi.elem(rb[1].start()).elem(0).elem(0).real());
#else
    void* spinorOut;
    GetMemoryPtr1(spinorOut, psi.getId());
#endif

    // wrap CPU host side pointers
    quda::ColorSpinorParam cpuParam(spinorOut, *const_cast<QudaInvertParam*>(&quda_inv_param), X,
				    pc_solution, quda_inv_param.output_location);
    quda::ColorSpinorField* h_out = quda::ColorSpinorField::Create(cpuParam);

    mugiqEigsolver->getVector(i, *h_out);

    if (g5)
      quda::gamma5(*h_out, *h_out);

    delete h_out;
  }
}
