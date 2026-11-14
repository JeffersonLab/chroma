// -*- C++ -*-
/*! \file                                                                    
 * \brief Multigrid prototype next
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_MGPROTON__
#define __INCLUDE_MGPROTON__

#include "chromabase.h"
#include "util/ferm/superb_contractions.h"
#include "util/ferm/superb_options.h"

#include <set>

#ifdef BUILD_SB
namespace Chroma
{
  namespace SB
  {
    /// Either a Chroma solver or a superb solver
    struct ChimeraSolver {
      /// Action type
      using Action = Handle<
	FermionAction<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>;

      /// Action
      Action S;

      /// State type
      using State =
	Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>;

      /// State
      State state;

      /// Chroma solver (optional)
      Handle<SystemSolver<LatticeFermion>> PP;

      /// Operator on scxyztX (optional)
      Operator<Nd + 7, Complex> op;

      /// Constructor
      /// \param fermAction: XML for the fermion action
      /// \param invParam: XML for the quark propagator
      /// \param u: gauge fields
      ChimeraSolver(const GroupXML_t& fermAction, const GroupXML_t& invParam,
		    const multi1d<LatticeColorMatrix>& u);
    };

    /// Callback function for doInversion
    /// Arguments of the callback:
    /// \param tensor: output tensor with order cSxyztXns,
    ///        (c: color, S: sink spin, xyzt: lattice dims, X: even-odd, n: colorvec index, s: source spin)
    /// \param first_s: first spin source (s) returned
    /// \param first_n: first colorvec index (n) returned

    template <typename COMPLEX = Complex>
    using doInversionOnTSliceFn =
      std::function<void(Tensor<Nd + 5, COMPLEX>, int first_s, int first_n)>;

    template <typename COMPLEX_CHI>
    void doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
		     const Tensor<2, SB::Complex>& spin_sources, int max_rhs,
		     const doInversionOnTSliceFn<SB::Complex>& call);
    template <typename COMPLEX_CHI>
    Tensor<Nd + 5, SB::Complex>
    doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
		int first_tslice_out, int n_tslice_out, const std::vector<int>& spin_sources,
		int max_rhs, const std::string& order_out = "cSxyztXns");

    /// Multiple spin-color lattice fields
    using MultipleLatticeFermions = std::vector<std::shared_ptr<LatticeFermion>>;
    using ConstMultipleLatticeFermions = std::vector<std::shared_ptr<const LatticeFermion>>;

    void doInversion(const ChimeraSolver& sol, MultipleLatticeFermions& psis,
		     const ConstMultipleLatticeFermions& chis, int max_rhs = 0);
    Operator<Nd + 7, Complex> getOperator(const ChimeraSolver& sol, int max_rhs = 1);

    template <typename COMPLEX>
    EigensolverFun<Nd + 7, COMPLEX> getInexactEigensolverGD(Operator<Nd + 7, COMPLEX> op,
							    const Options& ops);

    /// Function that returns V[from..from+size-1]
    /// \param from: first index to return
    /// \param size: number of vectors to return starting from `from`
    /// \param label: label indicating the columns

    template <std::size_t NOp, typename COMPLEX>
    using VectorFun = std::function<Tensor<NOp + 1, COMPLEX>(unsigned int, unsigned int, char)>;

    /// Returns the pieces of a projector of the form:
    ///   V*inv(U'*op*V)*U'*op, where U'*op*V = diag(lambda)
    ///
    /// It can work as a projector on the right of inv(op), that is P*inv(op), where
    ///   P = V*inv(U'*op*V)*U'*op
    /// or it can work as a projector on the left of inv(op), that is inv(op)*Q, where
    ///   Q = op*V*inv(U'*op*V)*U'
    ///
    /// Then, the trace of either P*inv(op) or inv(op)*Q is
    ///   tr V*inv(U'*op*V)*U' = \sum_i u_i'*vi/lambda_i

    template <std::size_t NOp, typename COMPLEX>
    struct Projector {
      /// Function that applies V * inv(U'*op*V) * U'
      Operator<NOp, COMPLEX> V_inv_Ut;

      /// Function that returns the i-th left base of the projector
      VectorFun<NOp, COMPLEX> V;

      /// Function that returns the i-th right base of the projector
      VectorFun<NOp, COMPLEX> U;

      /// Inner products with the operator, lambda_i = U[i]'*op*V[i]
      std::vector<COMPLEX> lambdas;

      /// Function that applies the operator
      Operator<NOp, COMPLEX> op;
    };

    /// Chroma or mgproton projector, that is, P*P*x = P*x

    struct ChimeraProjector {
      /// Action type
      using Action = Handle<
	FermionAction<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>;

      /// Action
      Action S;

      /// State type
      using State =
	Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>;

      /// State
      State state;

      /// Chroma projector (optional)
      Handle<Chroma::Projector<LatticeFermion>> chroma_proj;

      /// Operator on scxyztX (optional)
      SB::Projector<Nd + 7, Complex> op;

      /// Constructor
      /// \param fermAction: XML for the fermion action
      /// \param projParam: XML for the projector
      /// \param u: gauge fields
      ChimeraProjector(const GroupXML_t& fermAction, const GroupXML_t& projParam,
		    const multi1d<LatticeColorMatrix>& u);
    };

    void doVUAObliqueProjector(const ChimeraProjector& proj, MultipleLatticeFermions& psis,
			       const ConstMultipleLatticeFermions& chis, int max_rhs = 0);
    unsigned int getProjectorRank(const ChimeraProjector& proj);
    void getV(const ChimeraProjector& proj, unsigned int from, MultipleLatticeFermions& psis);
    void getU(const ChimeraProjector& proj, unsigned int from, MultipleLatticeFermions& psis);
    DComplex getLambda(const ChimeraProjector& proj, unsigned int index);
  }
}
#endif // BUILD_SB

#endif // __INCLUDE_MGPROTON__
