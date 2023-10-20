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

    /// Multiple spin-color lattice fields
    template <typename COMPLEX_CHI>
    Tensor<Nd + 5, SB::Complex>
    doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
		int first_tslice_out, int n_tslice_out, const std::vector<int>& spin_sources,
		int max_rhs, const std::string& order_out = "cSxyztXns");

    using MultipleLatticeFermions = std::vector<std::shared_ptr<LatticeFermion>>;
    using ConstMultipleLatticeFermions = std::vector<std::shared_ptr<const LatticeFermion>>;

    void doInversion(const ChimeraSolver& sol, MultipleLatticeFermions& psis,
		     const ConstMultipleLatticeFermions& chis, int max_rhs = 0);
    Operator<Nd + 7, Complex> getOperator(const ChimeraSolver& sol, int max_rhs = 1);

    template <typename COMPLEX>
    EigensolverFun<Nd + 7, COMPLEX> getInexactEigensolverGD(Operator<Nd + 7, COMPLEX> op,
							    const Options& ops);
  }
}
#endif // BUILD_SB

#endif // __INCLUDE_MGPROTON__
