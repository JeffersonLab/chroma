// -*- C++ -*-
/*! \file                                                                    
 * \brief Multigrid prototype next
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#include "util/ferm/mgproton.h"
#include "chromabase.h"
#include "util/ferm/superb_contractions.h"
#include "util/ferm/superb_options.h"

#include <set>

#ifdef BUILD_SB
namespace Chroma
{
  namespace SB
  {

    /// Verbosity level for solvers
    enum Verbosity {
      NoOutput = 0,    ///< Print nothing
      JustSummary = 1, ///< Print summary at the end of the solver execution
      Detailed = 2,    ///< Print progress
      VeryDetailed = 3 ///< Print whatever
    };

    /// Expected solver domain and image
    enum SolverSpace {
      FullSpace,     ///< input and output vector with full space support
      OnlyEvensSpace ///< input and output vector restricted to even sites
    };

    /// Return the map from string to Verbosity values
    inline const std::map<std::string, Verbosity>& getVerbosityMap()
    {
      static const std::map<std::string, Verbosity> m{{"NoOutput", NoOutput},
						      {"Summary", JustSummary},
						      {"Detailed", Detailed},
						      {"VeryDetailed", VeryDetailed},
						      {"false", NoOutput}};
      return m;
    }

    /// Return the map from string to Verbosity values
    inline const std::map<std::string, ColOrdering>& getColOrderingMap()
    {
      static const std::map<std::string, ColOrdering> m{{"row", RowMajor}, {"column", ColumnMajor}};
      return m;
    }

    namespace detail
    {
      /// Check that the common dimensions have the same size
      /// \param V: tensor to check
      /// \param W: other tensor to check
      /// \param check_dims: dimension labels to check (or none for all of them)
      /// \param not_check_dims: dimension labels to not check

      template <std::size_t NV, std::size_t NW, typename TV, typename TW>
      void check_compatible(Tensor<NV, TV> V, Tensor<NW, TW> W,
			    Maybe<std::string> check_dims = none,
			    Maybe<std::string> not_check_dims = none)
      {
	auto wdims = W.kvdim();
	for (auto it : V.kvdim())
	  if ((!check_dims.hasSome() || check_dims.getSome().find(it.first) != std::string::npos) &&
	      (!not_check_dims.hasSome() ||
	       not_check_dims.getSome().find(it.first) != std::string::npos) &&
	      wdims.count(it.first) > 0 && wdims.at(it.first) != it.second)
	    throw std::runtime_error("check_compatible: some label does not match");
      }

      /// Returns max ||I - C||_F, which is a heuristic of the orthogonality level of C=V^\dagger*V
      /// \param C: matrix to test
      /// \param order_t: dimension labels that do not participate in the orthogonalization
      /// \param order_cols: dimension labels that are the columns of the matrices V and W
      /// \param m_to_rows: map from columns to the alternative label for columns

      template <std::size_t Ncols, std::size_t N, typename COMPLEX>
      double ortho_level(Tensor<N, COMPLEX> C, const std::string& order_t,
			 const std::string& order_cols, const remap& m_to_rows)
      {
	// Check Nrows
	if (order_cols.size() != Ncols)
	  throw std::runtime_error("ortho_level: invalid template argument `Ncols`");

	// Check that the matrix is square
	std::string order_rows = detail::update_order(order_cols, m_to_rows);
	auto dim = C.kvdim();
	std::size_t m = volume(dim, order_rows);
	if (m != volume(dim, order_cols))
	  throw std::runtime_error("ortho_level: the input tensor is not square");

	// Compute ||I - C||_F^2
	C = C.clone();
	identity<N, COMPLEX>(dim, m_to_rows, OnMaster).scale(-1).addTo(C);
	auto fnorm2 = contract<N - Ncols * 2>(C.conj(), C, order_rows + order_cols, OnHost,
					      OnEveryoneReplicated);

	// Return the largest square root
	double tol = 0;
	fnorm2.foreachWithCPUFun(
	  [&](const COMPLEX& t) { tol = std::max(tol, (double)std::sqrt(std::real(t))); });
	return tol;
      }
    }

    /// Whether to check orthogonalization level, used by `ortho`

    enum class CheckOrthoLevel { dontCheckOrthoLevel, doCheckOrthoLevel };

    /// Orthonormalize W against V, W_out <- (I-V*V')*W_in*R, with R such that W_out'*W_out = I,
    /// for each combination of values of the order_t dimensions, or once if it is empty.
    /// \param V: orthogonal matrix to orthogonalize against
    /// \param W: matrix to orthogonalize
    /// \param order_t: dimension labels that do not participate in the orthogonalization
    /// \param order_rows: dimension labels that are the rows of the matrices V and W
    /// \param order_cols: dimension labels that are the columns of the matrices V and W
    /// \param max_its: maximum number of iterations
    /// \param checkOrthoLevel: whether to stop earlier if the orthogonality level is good enough
    /// \param verb: verbosity level
    /// \param prefix: prefix for logging
    ///
    /// NOTE: ortho is mostly used to normalize vectors, that is, tensors whose columns are
    /// aren't distributed. The function that checks the orthogonality level assumes that,
    /// and it's very inefficient when columns are distributed.

    template <std::size_t Nrows, std::size_t Ncols, std::size_t NV, typename COMPLEX,
	      std::size_t NW>
    void ortho(Tensor<NV, COMPLEX> V, Tensor<NW, COMPLEX> W, const std::string& order_t,
	       const std::string& order_rows, const std::string& order_cols,
	       unsigned int max_its = 4,
	       CheckOrthoLevel checkOrthoLevel = CheckOrthoLevel::doCheckOrthoLevel,
	       Verbosity verb = NoOutput, const std::string& prefix = "")
    {
      // Check Nrows
      if (order_rows.size() != Nrows)
	throw std::runtime_error("ortho: invalid template argument `Nrows`");
      if (order_cols.size() != Ncols)
	throw std::runtime_error("ortho: invalid template argument `Ncols`");

      // Quick shortcut
      if (!W)
	return;

      // Check that V and W are compatible, excepting the column dimensions
      if (V)
	detail::check_compatible(V, W, none, order_cols);

      Tracker _t(std::string("ortho ") + prefix);

      // Find the ordering
      std::string Wcorder = detail::remove_dimensions(W.order, order_t + order_rows);
      std::string Vcorder =
	V ? detail::remove_dimensions(V.order, order_t + order_rows) : std::string{};

      // Create an alternative view of W with different labels for the column
      remap Wac = detail::getNewLabels(Wcorder, (V ? V.order : std::string{}) + W.order);
      std::string Wacorder = detail::update_order(Wcorder, Wac);
      auto Wa = W.rename_dims(Wac);

      constexpr std::size_t Nt = NW - Nrows - Ncols;
      double l = 0;
      unsigned int i = 0;
      for (; i <= max_its;)
      {
	// W = W - V*(V'*W)
	if (V)
	  contract(V.scale(-1), contract<NV + NW - Nrows * 2>(V.conj(), W, order_rows), Vcorder,
		   AddTo, W);

	// Compute Wa'*W and the orthogonality level of the basis
	Tensor<Nt + 2 * Ncols, COMPLEX> C;
	if (checkOrthoLevel == CheckOrthoLevel::doCheckOrthoLevel)
	{
	  C = contract<Nt + 2 * Ncols>(Wa.conj(), W, order_rows, none, OnEveryoneReplicated);
	  l = detail::ortho_level<Ncols>(C, order_t, Wcorder, Wac);

	  if (verb >= Detailed)
	    QDPIO::cout << prefix << " ortho #its: " << i << " |I-V'*V|_F: " << detail::tostr(l)
			<< std::endl;

	  // If ||I-C|| < 1, then the basis has no linear dependencies; the conditioning of the basis is
	  // also related to that norm, but we choose an arbitrary close but smaller value than one.
	  if (l < 0.7)
	    break;
	}
	else
	{
	  C = contract<Nt + 2 * Ncols>(Wa.conj(), W, order_rows, none);

	  if (verb >= Detailed)
	    QDPIO::cout << prefix << " ortho #its: " << i << " (no checking ortho level)"
			<< std::endl;
	}

	if (i >= max_its)
	  throw std::runtime_error("ortho: failing in orthonormalizing the basis");

	// W = W/chol(Wa'*W), where Wa'*W has dimensions (rows,cols)=(Wacorder,Wcorder)
	cholInv<NW, Nt + 2 * Ncols, NW, COMPLEX>(std::move(C), Wacorder, Wcorder, Wa, Wacorder,
						 CopyTo, W);

	++i;

	// Stop by number of iterations when don't checking orthogonality level
	if (checkOrthoLevel == CheckOrthoLevel::dontCheckOrthoLevel && i >= max_its)
	  break;
      }

      if (verb >= JustSummary)
      {
	QDPIO::cout << prefix << " ortho summary rank: " << detail::volume(W.kvdim(), Wcorder)
		    << " #its: " << i;
	if (checkOrthoLevel == CheckOrthoLevel::doCheckOrthoLevel)
	  QDPIO::cout << " |I-V'*V|_F: " << detail::tostr(l);
	QDPIO::cout << std::endl;
      }
    }

    /// Orthonormalize W, W_out <- W_in*R, with R such that W_out'*W_out = I,
    /// for each combination of values of the order_t dimensions, or once if it is empty.
    /// \param W: matrix to orthogonalize
    /// \param order_t: dimension labels that do not participate in the orthogonalization
    /// \param order_rows: dimension labels that are the rows of the matrices V and W
    /// \param order_cols: dimension labels that are the columns of the matrices V and W
    /// \param max_its: maximum number of iterations
    /// \param checkOrthoLevel: whether to stop earlier if the orthogonality level is good enough
    /// \param verb: verbosity level
    /// \param prefix: prefix for logging
    ///
    /// NOTE: ortho is mostly used to normalize vectors, that is, tensors whose columns are
    /// aren't distributed. The function that checks the orthogonality level assumes that,
    /// and it's very inefficient when columns are distributed.

    template <std::size_t Nrows, std::size_t Ncols, std::size_t NW, typename COMPLEX>
    void ortho(Tensor<NW, COMPLEX> W, const std::string& order_t, const std::string& order_rows,
	       const std::string& order_cols, unsigned int max_its = 4,
	       CheckOrthoLevel checkOrthoLevel = CheckOrthoLevel::doCheckOrthoLevel,
	       Verbosity verb = NoOutput, const std::string& prefix = "")
    {
      ortho<Nrows, Ncols, NW, COMPLEX, NW>(Tensor<NW, COMPLEX>{}, W, order_t, order_rows,
					   order_cols, max_its, checkOrthoLevel, verb, prefix);
    }

    /// Solve iteratively op * y = x using FGMRES
    /// \param op: problem matrix
    /// \param prec: left preconditioner
    /// \param x: input right-hand-sides
    /// \param y: the solution vectors
    /// \param max_basis_size: maximum rank of the search subspace basis per t
    /// \param tol: maximum tolerance
    /// \param max_its: maximum number of iterations
    /// \param error_if_not_converged: throw an error if the tolerance was not satisfied
    /// \param ortho_each_its: orthogonalize every this number of iterations
    /// \param max_residual_updates: recompute residual vector every this number of restarts
    /// \param passing_initial_guess: whether `y` contains a solution guess
    /// \param verb: verbosity level
    /// \param prefix: prefix printed before every line
    ///
    /// For FGMRES we find an approximation of op^{-1}*b in a space Z by minimizing
    /// || op*Z*x - b ||_2:
    ///
    /// argmin_x || op*Z*x - b ||_2 = (U'*op*Z)^{-1} * (U'*b), where U = op*Z
    ///
    /// In FGMRES, it just happens that Z_0 = prec * r and Z_i = prec * U_{i-1} if i>0.

    template <std::size_t NOp, typename COMPLEX>
    void fgmres(const Operator<NOp, COMPLEX>& op, Operator<NOp, COMPLEX> prec,
		const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX>& y,
		unsigned int max_basis_size, double tol, unsigned int max_its = 0,
		bool error_if_not_converged = true, unsigned int ortho_each_its = 0,
		unsigned int max_residual_updates = 0, bool passing_initial_guess = false,
		Verbosity verb = NoOutput, std::string prefix = "")
    {
      detail::log(1, prefix + " starting fgmres");

      // TODO: add optimizations for multiple operators
      if (op.order_t.size() > 0)
	throw std::runtime_error("Not implemented");

      // Check options
      if (max_its == 0 && tol <= 0)
	throw std::runtime_error("fgmres: please give a stopping criterion, either a tolerance or "
				 "a maximum number of iterations");
      if (max_basis_size == 0)
	max_basis_size = 5;
      if (max_its == 0)
	max_its = std::numeric_limits<unsigned int>::max();
      if (ortho_each_its == 0)
	ortho_each_its = (std::is_same<COMPLEX, double>::value ||
			  std::is_same<COMPLEX, std::complex<double>>::value)
			   ? 8
			   : 4;
      if (max_residual_updates == 0)
	max_residual_updates = (std::is_same<COMPLEX, double>::value ||
				std::is_same<COMPLEX, std::complex<double>>::value)
				 ? 4
				 : 2;

      // Check that the operator and the preconditioner are compatible with the input and output vectors
      if (!op.d.is_compatible(x) || !op.d.is_compatible(y) || (prec && !prec.d.is_compatible(op.d)))
	throw std::runtime_error("Either the input or the output vector isn't compatible with the "
				 "operator or the preconditioner");

      // Get an unused label for the search subspace columns
      char Vc = detail::get_free_label(x.order);
      char Vac = detail::get_free_label(x.order + std::string(1, Vc));
      std::string order_cols = detail::remove_dimensions(x.order, op.i.order);
      std::string order_rows = detail::remove_dimensions(op.d.order, op.order_t);
      std::size_t num_cols = x.volume(order_cols);
      if (num_cols == 0)
	return;

      // Counting op and prec applications
      unsigned int nops = 0, nprecs = 0;

      // Compute residual, r = op * y - x
      Tensor<NOp + 1, COMPLEX> r;
      if (passing_initial_guess)
      {
	r = op(y);
	nops += num_cols;
      }
      else
      {
	r = op.template make_compatible_img<NOp + 1>(order_cols, x.kvdim());
	r.set_zero();
	y.set_zero();
      }
      x.scale(-1).addTo(r);
      auto normr0 = norm<1>(r, op.order_t + order_cols); // type std::vector<real of T>
      if (max(normr0) == 0)
	return;

      // Allocate the search subspace U (onto the left singular space), and
      // Z (onto the right singular space)
      auto U = r.template make_compatible<NOp + 2>(std::string({'%', Vc}), '%', "",
						   {{Vc, max_basis_size}});

      // Allocate the search subspace Z (onto the right singular space)
      auto Z = U.make_compatible();

      // Extend r with Vc
      auto r_Vc = r.append_dimension(Vc);

      // Do the iterations
      auto normr = normr0.clone();	 ///< residual norms
      unsigned int it = 0;		 ///< iteration number
      double max_tol = HUGE_VAL;	 ///< maximum residual norm
      unsigned int ires = 0;		 ///< vector index for the last restart starting
      unsigned int residual_updates = 0; ///< number of residual updates
      for (it = 0; it < max_its;)
      {
	unsigned int expansion_size =
	  std::min(std::min(max_basis_size - ires, ortho_each_its), max_its - it);

	// Expand the search subspace from residual
	if (prec)
	{
	  for (unsigned int i = 0; i < expansion_size; ++i)
	  {
	    // Z(:,ires+i) = prec * U(:,ires+i-1)
	    prec(i == 0 ? r_Vc : U.kvslice_from_size({{Vc, ires + i - 1}}, {{Vc, 1}}),
		 Z.kvslice_from_size({{Vc, ires + i}}, {{Vc, 1}}));
	    nprecs += num_cols;

	    // U(:,ires+i) = op * Z(:,ires+i)
	    op(Z.kvslice_from_size({{Vc, ires + i}}, {{Vc, 1}}),
	       U.kvslice_from_size({{Vc, ires + i}}, {{Vc, 1}}));
	    nops += num_cols;

	    ++it;
	  }
	}
	else
	{
	  // U(:,ires+i) = op ^ i * r for i=1..expansion_size-1
	  op(r_Vc, U.kvslice_from_size({{Vc, ires}}, {{Vc, expansion_size}}), Vc);
	  r_Vc.copyTo(Z.kvslice_from_size({{Vc, ires}}, {{Vc, 1}}));
	  U.kvslice_from_size({{Vc, ires}}, {{Vc, expansion_size - 1}})
	    .copyTo(Z.kvslice_from_size({{Vc, ires + 1}}, {{Vc, expansion_size - 1}}));
	  nops += num_cols * expansion_size;
	  it += expansion_size;
	}

	// Orthogonalize U and put it into W: W = orth(U(:,2:end))
	// NOTE: for small max_basis_size and assuming r is far from a left singular vector,
	//       a light or none orthogonalization should be enough
	unsigned int basis_size = ires + expansion_size;
	auto Up = U.kvslice_from_size({}, {{Vc, basis_size}});
	auto Uo = Up.rename_dims({{Vc, Vac}});
	Uo = Uo.clone(); /// TODO: Avoid this
	ortho<NOp, 1>(Uo, op.order_t + order_cols, order_rows, std::string(1, Vac), 4,
		      CheckOrthoLevel::doCheckOrthoLevel, verb, prefix);

	// Restrict to Uo: [x_rt H_rt] = Uo'*U = Uo'*[r Up]
	auto x_rt = contract<2>(Uo.conj(), r, order_rows);
	auto H_rt = contract<3>(Uo.conj(), Up, order_rows);

	// Solve the projected problem: y_rt = (Uo'*U(:2:end))\(Uo'*r);
	auto y_rt = solve<1, 2>(H_rt, std::string(1, Vac), std::string(1, Vc),
				x_rt.rename_dims({{Vac, Vc}}), std::string(1, Vc))
		      .rename_dims({{Vac, Vc}})
		      .make_sure(none, none, OnEveryoneReplicated);

	// Update solution: y += -Z*y_rt
	contract(Z.kvslice_from_size({}, {{Vc, basis_size}}).scale(-1), y_rt, std::string(1, Vc),
		 AddTo, y);

	// Compute residual
	if (residual_updates < max_residual_updates)
	{
	  // Update residual by saving a matvec: r += -U*y_rt
	  contract(Up.scale(-1), y_rt, std::string(1, Vc), AddTo, r);
	  residual_updates++;
	}
	else
	{
	  op(y, r); // r = op(y)
	  x.scale(-1).addTo(r);
	  nops += num_cols;
	  residual_updates = 0;
	}

	// Compute the norm
	auto normr = norm<1>(r, op.order_t + order_cols);

	// Check final residual
	if (superbblas::getDebugLevel() > 0)
	{
	  auto rd = r.make_compatible();
	  op(y, rd); // rd = op(y)
	  nops += num_cols;
	  x.scale(-1).addTo(rd);
	  r.scale(-1).addTo(rd);
	  auto normrd = norm<1>(rd, op.order_t + order_cols);
	  double max_tol_d = max(div(normrd, normr));
	  QDPIO::cout << prefix
		      << " MGPROTON FGMRES error in residual vector: " << detail::tostr(max_tol_d)
		      << std::endl;
	}

	// Get the worse tolerance
	max_tol = max(div(normr, normr0));

	// Report iteration
	if (verb >= Detailed)
	  QDPIO::cout << prefix << " MGPROTON FGMRES iteration #its.: " << it
		      << " max rel. residual: " << detail::tostr(max_tol, 2) << std::endl;

	// Stop if the residual tolerance is satisfied
	if (max_tol <= tol)
	  break;

	// Start a new expansion after the current one if there's space left; otherwise do a full restart
	ires += expansion_size;
	if (ires >= max_basis_size)
	  ires = 0;
      }

      // Check final residual
      if (error_if_not_converged)
      {
	op(y, r); // r = op(y)
	nops += num_cols;
	x.scale(-1).addTo(r);
	auto normr = norm<1>(r, op.order_t + order_cols);
	max_tol = max(div(normr, normr0));
	if (tol > 0 && max_tol > tol)
	  throw std::runtime_error("fgmres didn't converged and you ask for checking the error");
      }

      // Report iteration
      if (verb >= JustSummary)
	QDPIO::cout << prefix << " MGPROTON FGMRES summary #its.: " << it
		    << " max rel. residual: " << detail::tostr(max_tol, 2) << " matvecs: " << nops
		    << " precs: " << nprecs << std::endl;
    }

    /// Solve iteratively op * y = x using BICGSTAB
    /// \param op: problem matrix
    /// \param prec: left preconditioner
    /// \param x: input right-hand-sides
    /// \param y: the solution vectors
    /// \param tol: maximum tolerance
    /// \param max_its: maximum number of iterations
    /// \param error_if_not_converged: throw an error if the tolerance was not satisfied
    /// \param max_residual_updates: recompute residual vector every this number of restarts
    /// \param passing_initial_guess: whether `y` contains a solution guess
    /// \param verb: verbosity level
    /// \param prefix: prefix printed before every line
    ///
    /// Method:
    /// r0 = p = r = b - A * x  // compute initial residual
    /// rho = r0' * r
    /// for i=1,2,..
    ///   Kp = prec * p
    ///   AKp = A * Kp
    ///   alpha = rho / (r0' * AKp)
    ///   s = r - AKp * alpha
    ///   Ks = prec * s
    ///   AKs = A * Ks
    ///   omega = (AKs' * Ks) / (AKs' * AKs)
    ///   x = x + Kp * alpha + Ks * omega
    ///   r = s - AKs * omega
    ///   exit if |r| < |b|*tol
    ///   prev_rho = rho
    ///   rho = r0' * r
    ///   beta = (rho/prev_rho)*(alpha/omega)
    ///   p = r + beta*(p - omega*AKp)
    /// end

    template <std::size_t NOp, typename COMPLEX>
    void bicgstab(const Operator<NOp, COMPLEX>& op, Operator<NOp, COMPLEX> prec,
		  const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX>& y, double tol,
		  unsigned int max_its = 0, bool error_if_not_converged = true,
		  unsigned int max_residual_updates = 0, bool passing_initial_guess = false,
		  Verbosity verb = NoOutput, std::string prefix = "")
    {
      detail::log(1, prefix + " starting bicgstab");

      // TODO: add optimizations for multiple operators
      if (op.order_t.size() > 0)
	throw std::runtime_error("Not implemented");

      // Check options
      if (max_its == 0 && tol <= 0)
	throw std::runtime_error(
	  "bicgstab: please give a stopping criterion, either a tolerance or "
	  "a maximum number of iterations");
      if (max_its == 0)
	max_its = std::numeric_limits<unsigned int>::max();
      if (max_residual_updates == 0)
	max_residual_updates = (std::is_same<COMPLEX, double>::value ||
				std::is_same<COMPLEX, std::complex<double>>::value)
				 ? 100
				 : 100;

      // Check that the operator and the preconditioner are compatible with the input and output vectors
      if (!op.d.is_compatible(x) || !op.d.is_compatible(y) || (prec && !prec.d.is_compatible(op.d)))
	throw std::runtime_error("Either the input or the output vector isn't compatible with the "
				 "operator or the preconditioner");

      // Get an unused label for the search subspace columns
      std::string order_cols = detail::remove_dimensions(x.order, op.i.order);
      std::string order_rows = detail::remove_dimensions(op.d.order, op.order_t);
      std::size_t num_cols = x.volume(order_cols);
      if (num_cols == 0)
	return;

      // Counting op and prec applications
      unsigned int nops = 0, nprecs = 0;

      // Compute residual, r = x - op * y
      Tensor<NOp + 1, COMPLEX> r;
      if (passing_initial_guess)
      {
	r = op(y).scale(-1);
	nops += num_cols;
      }
      else
      {
	r = op.template make_compatible_img<NOp + 1>(order_cols, x.kvdim());
	r.set_zero();
	y.set_zero();
      }
      x.addTo(r);
      auto normr0 = norm<1>(r, op.order_t + order_cols); // type std::vector<real of T>
      if (max(normr0) == 0)
	return;

      // Choose an arbitrary vector that isn't orthogonal to r
      auto r0 = op.template make_compatible_img<NOp + 1>(order_cols, x.kvdim());
      r.copyTo(r0);

      // Do the iterations
      auto normr = normr0.clone();	 ///< residual norms
      unsigned int it = 0;		 ///< iteration number
      double max_tol = HUGE_VAL;	 ///< maximum residual norm
      unsigned int residual_updates = 0; ///< number of residual updates
      auto p = r0.clone();		 ///< p0 = r0
      auto rho = contract<1>(r, r0.conj(), order_rows); // rho = r0' * r
      auto Kp = prec ? r0.make_compatible() : p;
      auto AKp = r0.make_compatible();
      auto Ks = prec ? r0.make_compatible() : r;
      auto AKs = r0.make_compatible();
      for (it = 0; it < max_its;)
      {
	// Kp = prec * p
	if (prec)
	{
	  prec(p, Kp);
	  nprecs += num_cols;
	}

	// AKp = A * Kp
	op(Kp, AKp);
	nops += num_cols;

	// alpha = rho / (r0' * AKp)
	auto alpha = div(rho, contract<1>(AKp, r0.conj(), order_rows));

	// s = r - alpha * AKp
	auto s = r;
	contract<NOp + 1>(AKp, alpha.scale(-1), "", AddTo, s);

	// Ks = prec * s
	if (prec)
	{
	  prec(s, Ks);
	  nprecs += num_cols;
	}

	// AKs = A * Ks
	op(Ks, AKs);
	nops += num_cols;

	// omega = (AKs' * Ks) / (AKs' * AKs)
	auto omega =
	  div(contract<1>(AKs.conj(), Ks, order_rows), contract<1>(AKs, AKs.conj(), order_rows));

	// y = y + alpha * Kp + omega * Ks
	contract<NOp + 1>(Kp, alpha, "", AddTo, y);
	contract<NOp + 1>(Ks, omega, "", AddTo, y);

	// r = s - omega * AKs (NOTE: s == r)
	contract<NOp + 1>(AKs, omega.scale(-1), "", AddTo, r);

	// prev_rho = rho, rho = r0' * r 
	auto prev_rho = rho;
	contract<1>(r, r0.conj(), order_rows, CopyTo, rho);

	// beta = rho / prev_rho * (alpha / omega)
	auto rho_ratio = div(rho, prev_rho);
	auto beta = mult(rho_ratio, div(alpha, omega));

	// beta_omega = beta * omega = rho / prev_rho * alpha
	auto beta_omega = mult(rho_ratio, alpha);

	// p = r + beta * p - beta_omega * AKp
	auto aux = AKs;
	contract<NOp + 1>(p, beta, "", CopyTo, aux);
	aux.copyTo(p);
	contract<NOp + 1>(AKp, beta_omega.scale(-1), "", AddTo, p);
	r.addTo(p);

	// Recompute residual vector if needed
	if (residual_updates < max_residual_updates)
	{
	  residual_updates++;
	}
	else
	{
	  op(y.scale(-1), r); // r = op(-y)
	  x.addTo(r);
	  r.copyTo(r0);
	  r.copyTo(p);
	  contract<1>(r, r0.conj(), order_rows, CopyTo, rho);
	  nops += num_cols;
	  residual_updates = 0;
	}

	// Compute the norm
	auto normr = norm<1>(r, op.order_t + order_cols);

	// Show error in the residual
	if (superbblas::getDebugLevel() > 0)
	{
	  auto rd = r.make_compatible();
	  op(y, rd); // rd = op(y)
	  nops += num_cols;
	  x.scale(-1).addTo(rd);
	  r.addTo(rd);
	  auto normrd = norm<1>(rd, op.order_t + order_cols);
	  double max_tol_d = max(div(normrd, normr));
	  QDPIO::cout << prefix
		      << " MGPROTON BICGSTAB error in residual vector: " << detail::tostr(max_tol_d)
		      << std::endl;
	}

	// Get the worse tolerance
	max_tol = max(div(normr, normr0));

	// Report iteration
	if (verb >= Detailed)
	  QDPIO::cout << prefix << " MGPROTON BICGSTAB iteration #its.: " << it
		      << " max rel. residual: " << detail::tostr(max_tol, 2) << std::endl;

	// Increase iterator counter
	++it;

	// Stop if the residual tolerance is satisfied
	if (max_tol <= tol)
	  break;
      }

      // Check final residual
      if (error_if_not_converged)
      {
	op(y, r); // r = op(y)
	nops += num_cols;
	x.scale(-1).addTo(r);
	auto normr = norm<1>(r, op.order_t + order_cols);
	max_tol = max(div(normr, normr0));
	if (tol > 0 && max_tol > tol)
	  throw std::runtime_error("bicgstab didn't converged and you ask for checking the error");
      }

      // Report iteration
      if (verb >= JustSummary)
	QDPIO::cout << prefix << " MGPROTON BICGSTAB summary #its.: " << it
		    << " max rel. residual: " << detail::tostr(max_tol, 2) << " matvecs: " << nops
		    << " precs: " << nprecs << std::endl;
    }

    /// Solve iteratively op * y = x using Minimal Residual (valid for positive definite linear systems)
    /// \param op: problem matrix
    /// \param prec: left preconditioner
    /// \param x: input right-hand-sides
    /// \param y: the solution vectors
    /// \param tol: maximum tolerance
    /// \param max_its: maximum number of iterations
    /// \param error_if_not_converged: throw an error if the tolerance was not satisfied
    /// \param max_residual_updates: recompute residual vector every this number of restarts
    /// \param passing_initial_guess: whether `y` contains a solution guess
    /// \param verb: verbosity level
    /// \param prefix: prefix printed before every line
    ///
    /// Method:
    /// r = b - A * x  // compute initial residual
    /// while |r| > |b|*tol
    ///   kr = prec * r
    ///   p = A * kr
    ///   alpha = (p' * r) / (p' * p)
    ///   x = x + kr * alpha
    ///   r = r - p * alpha
    /// end

    template <std::size_t NOp, typename COMPLEX>
    void mr(const Operator<NOp, COMPLEX>& op, const Operator<NOp, COMPLEX>& prec,
	    const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX>& y, double tol,
	    unsigned int max_its = 0, bool error_if_not_converged = true,
	    unsigned int max_residual_updates = 0, bool passing_initial_guess = false,
	    Verbosity verb = NoOutput, std::string prefix = "")
    {
      detail::log(1, prefix + " starting minimal residual");

      // TODO: add optimizations for multiple operators
      if (op.order_t.size() > 0)
	throw std::runtime_error("Not implemented");

      // Check options
      if (max_its == 0 && tol <= 0)
	throw std::runtime_error("mr: please give a stopping criterion, either a tolerance or "
				 "a maximum number of iterations");
      if (max_its == 0)
	max_its = std::numeric_limits<unsigned int>::max();
      if (max_residual_updates == 0)
	max_residual_updates = (std::is_same<COMPLEX, double>::value ||
				std::is_same<COMPLEX, std::complex<double>>::value)
				 ? 100
				 : 100;

      // Check that the operator is compatible with the input and output vectors
      if (!op.d.is_compatible(x) || !op.d.is_compatible(y))
	throw std::runtime_error(
	  "mr: Either the input or the output vector isn't compatible with the "
	  "operator");

      // Get an unused label for the search subspace columns
      std::string order_cols = detail::remove_dimensions(x.order, op.i.order);
      std::string order_rows = detail::remove_dimensions(op.d.order, op.order_t);
      std::size_t num_cols = x.volume(order_cols);
      if (num_cols == 0)
	return;

      // Counting op applications
      unsigned int nops = 0, nprecs = 0;

      // Compute residual, r = x - op * y
      Tensor<NOp + 1, COMPLEX> r;
      if (passing_initial_guess)
      {
	r = op(y).scale(-1);
	nops += num_cols;
      }
      else
      {
	r = op.template make_compatible_img<NOp + 1>(order_cols, x.kvdim());
	r.set_zero();
	y.set_zero();
      }
      x.addTo(r);
      auto normr0 = norm<1>(r, op.order_t + order_cols); // type std::vector<real of T>
      if (max(normr0) == 0)
	return;

      // Do the iterations
      auto normr = normr0.clone();		///< residual norms
      unsigned int it = 0;			///< iteration number
      double max_tol = HUGE_VAL;		///< maximum residual norm
      unsigned int residual_updates = 0;	///< number of residual updates
      auto p = r.make_compatible();		///< p will hold A * prec * r
      auto kr = prec ? r.make_compatible() : r; ///< p will hold prec * r
      for (it = 0; it < max_its;)
      {
	// kr = prec * r
	if (prec)
	{
	  prec(r, kr);
	  nprecs += num_cols;
	}

	// p = A * kr
	op(kr, p);
	nops += num_cols;

	// alpha = (p' * r) / (p' * p)
	auto alpha =
	  div(contract<1>(r, p.conj(), order_rows), contract<1>(p, p.conj(), order_rows));

	// y = y + alpha * kr
	contract<NOp + 1>(kr, alpha, "", AddTo, y);

	// r = r - alpha * p
	contract<NOp + 1>(p.scale(-1), alpha, "", AddTo, r);

	// Update residual if needed
	if (residual_updates < max_residual_updates)
	{
	  residual_updates++;
	}
	else
	{
	  op(y.scale(-1), r); // r = op(-y)
	  x.addTo(r);
	  nops += num_cols;
	  residual_updates = 0;
	}

	// Compute the norm
	auto normr = norm<1>(r, op.order_t + order_cols);

	// Show residual error
	if (superbblas::getDebugLevel() > 0)
	{
	  auto rd = r.make_compatible();
	  op(y, rd); // rd = op(y)
	  nops += num_cols;
	  x.scale(-1).addTo(rd);
	  r.addTo(rd);
	  auto max_norm = max(div(norm<1>(rd, op.order_t + order_cols), normr));
	  QDPIO::cout << prefix
		      << " MGPROTON MR error in residual vector: " << detail::tostr(max_norm)
		      << std::endl;
	}

	// Get the worse tolerance
	max_tol = max(div(normr, normr0));

	// Report iteration
	if (verb >= Detailed)
	  QDPIO::cout << prefix << " MGPROTON MR iteration #its.: " << it
		      << " max rel. residual: " << detail::tostr(max_tol, 2) << std::endl;

	// Increase iterator counter
	++it;

	// Stop if the residual tolerance is satisfied
	if (max_tol <= tol)
	  break;
      }

      // Check final residual
      if (error_if_not_converged)
      {
	op(y, r); // r = op(y)
	nops += num_cols;
	x.scale(-1).addTo(r);
	auto normr = norm<1>(r, op.order_t + order_cols);
	max_tol = max(div(normr, normr0));
	if (tol > 0 && max_tol > tol)
	  throw std::runtime_error("mr didn't converged and you ask for checking the error");
      }

      // Report iteration
      if (verb >= JustSummary)
	QDPIO::cout << prefix << " MGPROTON MR summary #its.: " << it
		    << " max rel. residual: " << detail::tostr(max_tol, 2) << " matvecs: " << nops
		    << " precs: " << nprecs << std::endl;
    }

    /// Solve iteratively op * y = x using Generalized Conjugate Residual
    /// \param op: problem matrix
    /// \param x: input right-hand-sides
    /// \param y: the solution vectors
    /// \param tol: maximum tolerance
    /// \param max_its: maximum number of iterations
    /// \param error_if_not_converged: throw an error if the tolerance was not satisfied
    /// \param max_residual_updates: recompute residual vector every this number of restarts
    /// \param passing_initial_guess: whether `y` contains a solution guess
    /// \param verb: verbosity level
    /// \param prefix: prefix printed before every line
    ///
    /// Method:
    /// r = b - A * x  // compute initial residual
    /// P = AP = []
    /// while |r| > |b|*tol
    ///   kr = prec * r
    ///   Akr = A * kr
    ///   beta = (AP' * AP)\(AP' * Akr)
    ///   p = (I - P/(P'A'AP)*P'A'A) * prec * r = kr - P * beta
    ///   Ap = A * p = Akr - AP * beta
    ///   alpha = (Ap' * r) / (Ap' * Ap)
    ///   x = x + p * alpha
    ///   r = r - Ap * alpha
    ///   P = [P p]
    ///   AP = [AP Ap]
    /// end

    template <std::size_t NOp, typename COMPLEX>
    void gcr(const Operator<NOp, COMPLEX>& op, const Operator<NOp, COMPLEX>& prec,
	     const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX>& y,
	     unsigned int max_basis_size, double tol, unsigned int max_its = 0,
	     bool error_if_not_converged = true, unsigned int max_residual_updates = 0,
	     bool passing_initial_guess = false, Verbosity verb = NoOutput, std::string prefix = "")
    {
      detail::log(1, prefix + " starting minimal residual");

      // TODO: add optimizations for multiple operators
      if (op.order_t.size() > 0)
	throw std::runtime_error("Not implemented");

      // Check options
      if (max_its == 0 && tol <= 0)
	throw std::runtime_error("mr: please give a stopping criterion, either a tolerance or "
				 "a maximum number of iterations");
      if (max_its == 0)
	max_its = std::numeric_limits<unsigned int>::max();
      if (max_residual_updates == 0)
	max_residual_updates = (std::is_same<COMPLEX, double>::value ||
				std::is_same<COMPLEX, std::complex<double>>::value)
				 ? 100
				 : 100;

      // Check that the operator is compatible with the input and output vectors
      if (!op.d.is_compatible(x) || !op.d.is_compatible(y))
	throw std::runtime_error(
	  "mr: Either the input or the output vector isn't compatible with the "
	  "operator");

      // Get an unused label for the search subspace columns
      char Vc = detail::get_free_label(x.order);
      std::string order_cols = detail::remove_dimensions(x.order, op.i.order);
      std::string order_rows = detail::remove_dimensions(op.d.order, op.order_t);
      std::size_t num_cols = x.volume(order_cols);
      if (num_cols == 0)
	return;

      // Counting op applications
      unsigned int nops = 0, nprecs = 0;

      // Compute residual, r = x - op * y
      Tensor<NOp + 1, COMPLEX> r;
      if (passing_initial_guess)
      {
	r = op(y).scale(-1);
	nops += num_cols;
      }
      else
      {
	r = op.template make_compatible_img<NOp + 1>(order_cols, x.kvdim());
	r.set_zero();
	y.set_zero();
      }
      x.addTo(r);
      auto normr0 = norm<1>(r, op.order_t + order_cols); // type std::vector<real of T>
      if (max(normr0) == 0)
	return;

      // Allocate the search subspace P (onto the left singular space), and
      // AP (= op * P)
      auto P = r.template make_compatible<NOp + 2>(std::string({'%', Vc}), '%', "",
						   {{Vc, max_basis_size}});
      auto AP = P.make_compatible();

      // Do the iterations
      auto normr = normr0.clone();		///< residual norms
      unsigned int it = 0;			///< iteration number
      double max_tol = HUGE_VAL;		///< maximum residual norm
      unsigned int residual_updates = 0;	///< number of residual updates
      auto p = r.make_compatible();		///< p will hold the next column in P
      auto Ap = r.make_compatible();		///< Ap will hold the next column in AP
      auto kr = prec ? r.make_compatible() : r; ///< p will hold prec * r
      auto Akr = r.make_compatible();		///< Akr will hold A * kr
      unsigned int active_P = 0;		///< number of active columns in P and AP
      auto AP_norm2 = P.template like_this<2>(order_cols + std::string{Vc}, {}, none,
					      detail::compatible_replicated_distribution(P.dist));
      for (it = 0; it < max_its;)
      {
	// kr = prec * r
	if (prec)
	{
	  prec(r, kr);
	  nprecs += num_cols;
	}

	// Akr = A * kr
	op(kr, Akr);
	nops += num_cols;

	if (active_P > 0)
	{
	  auto P_act = P.kvslice_from_size({}, {{Vc, active_P}});
	  auto AP_act = AP.kvslice_from_size({}, {{Vc, active_P}});
	  auto AP_norm2_act = AP_norm2.kvslice_from_size({}, {{Vc, active_P}});

	  // beta = (AP' * Akr) / AP_norm2
	  auto beta = div(contract<2>(Akr, AP_act.conj(), order_rows), AP_norm2_act);

	  // p = kr - P * beta
	  kr.copyTo(p);
	  contract<NOp + 1>(P_act, beta.scale(-1), std::string{Vc}, AddTo, p);

	  // Ap = Akr - AP * beta
	  Akr.copyTo(Ap);
	  contract<NOp + 1>(AP_act, beta.scale(-1), std::string{Vc}, AddTo, Ap);
	}
	else
	{
	  kr.copyTo(p);
	  Akr.copyTo(Ap);
	}

	// alpha = (Ap' * r) / (Ap' * Ap)
	auto alpha =
	  div(contract<1>(r, Ap.conj(), order_rows), contract<1>(Ap, Ap.conj(), order_rows));

	// y = y + alpha * p
	contract<NOp + 1>(p, alpha, "", AddTo, y);

	// r = r - alpha * Ap
	contract<NOp + 1>(Ap.scale(-1), alpha, "", AddTo, r);

	// Update residual if needed
	if (residual_updates < max_residual_updates)
	{
	  residual_updates++;
	}
	else
	{
	  op(y.scale(-1), r); // r = op(-y)
	  x.addTo(r);
	  nops += num_cols;
	  active_P = 0;
	  residual_updates = 0;
	}

	// Compute the norm
	auto normr = norm<1>(r, op.order_t + order_cols);

	// Check final residual
	if (superbblas::getDebugLevel() > 0)
	{
	  auto rd = r.make_compatible();
	  op(y, rd); // rd = op(y)
	  nops += num_cols;
	  x.scale(-1).addTo(rd);
	  r.addTo(rd);
	  auto max_norm = max(div(norm<1>(rd, op.order_t + order_cols), normr));
	  QDPIO::cout << prefix
		      << " MGPROTON GCR error in residual vector: " << detail::tostr(max_norm)
		      << std::endl;
	}

	// Get the worse tolerance
	max_tol = max(div(normr, normr0));

	// Report iteration
	if (verb >= Detailed)
	  QDPIO::cout << prefix << " MGPROTON GCR iteration #its.: " << it
		      << " max rel. residual: " << detail::tostr(max_tol, 2) << std::endl;

	// Increase iterator counter
	++it;

	// Stop if the residual tolerance is satisfied
	if (max_tol <= tol)
	  break;

	// Update P, AP, and AP_norm2
	if (active_P >= max_basis_size)
	  active_P = 1;
	else
	  ++active_P;
	auto P_new = P.kvslice_from_size({{Vc, active_P - 1}}, {{Vc, 1}});
	auto AP_new = AP.kvslice_from_size({{Vc, active_P - 1}}, {{Vc, 1}});
	p.copyTo(P_new);
	Ap.copyTo(AP_new);
	contract<2>(AP_new.conj(), AP_new, order_rows, CopyTo,
		    AP_norm2.kvslice_from_size({{Vc, active_P - 1}}, {{Vc, 1}}));
      }

      // Check final residual
      if (error_if_not_converged)
      {
	op(y, r); // r = op(y)
	nops += num_cols;
	x.scale(-1).addTo(r);
	auto normr = norm<1>(r, op.order_t + order_cols);
	max_tol = max(div(normr, normr0));
	if (tol > 0 && max_tol > tol)
	  throw std::runtime_error("gcr didn't converged and you ask for checking the error");
      }

      // Report iteration
      if (verb >= JustSummary)
	QDPIO::cout << prefix << " MGPROTON GCR summary #its.: " << it
		    << " max rel. residual: " << detail::tostr(max_tol, 2) << " matvecs: " << nops
		    << " precs: " << nprecs << std::endl;
    }

    template <std::size_t NOp, typename COMPLEX>
    Operator<NOp, COMPLEX> getSolver(const Operator<NOp, COMPLEX>& op, const Options& ops,
				     const Operator<NOp, COMPLEX>& prec = Operator<NOp, COMPLEX>(),
				     SolverSpace solverSpace = FullSpace);

    template <std::size_t NOp, typename COMPLEX>
    Projector<NOp, COMPLEX> getProjector(const Operator<NOp, COMPLEX>& op, const Options& ops);

    namespace detail
    {
      /// Apply the given function to all the given columns in x
      /// \param x: tensor to apply to fun
      /// \param y: output tensor
      /// \param max_rhs: maximum number of columns to apply `fun` at once
      /// \param fun: function to apply
      /// \param d: column dimension

      template <std::size_t N, typename T, typename Func>
      inline void foreachInChuncks(const Tensor<N, T>& x, Tensor<N, T>& y, unsigned int max_rhs,
				   const Func& fun, char d = 'n')
      {
	unsigned int n = x.kvdim().at(d);
	if (max_rhs == 0)
	  max_rhs = n;
	for (unsigned int i = 0, step = std::min(max_rhs, n); i < n;
	     i += step, step = std::min(max_rhs, n - i))
	  fun(x.kvslice_from_size({{d, i}}, {{d, step}}),
	      y.kvslice_from_size({{d, i}}, {{d, step}}));
      }

      /// Returns a FGMRES solver
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getFGMRESSolver(Operator<NOp, COMPLEX> op, const Options& ops,
					     Operator<NOp, COMPLEX> prec_)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup fgmres ") + prefix);

	// Get preconditioner
	Operator<NOp, COMPLEX> prec;
	Maybe<const Options&> precOps = getOptionsMaybe(ops, "prec");
	if (precOps && prec_)
	  throw std::runtime_error(
	    "getFGMRESSolver: invalid `prec` tag: the solver got a preconditioner already");
	if (precOps)
	  prec = getSolver(op, precOps.getSome());
	else if (prec_)
	  prec = prec_;

	// Get the remainder options
	unsigned int max_basis_size = getOption<unsigned int>(ops, "max_basis_size", 0);
	double tol = getOption<double>(ops, "tol", 0.0);
	unsigned int max_its = getOption<unsigned int>(ops, "max_its", 0);
	if (max_its == 0 && tol <= 0)
	  ops.throw_error("set either `tol` or `max_its`");
	bool error_if_not_converged = getOption<bool>(ops, "error_if_not_converged", true);
	unsigned int ortho_each_its = getOption<unsigned int>(ops, "ortho_each_its", 0);
	unsigned int max_residual_updates = getOption<unsigned int>(ops, "max_residual_updates", 0);
	unsigned int max_simultaneous_rhs = getOption<unsigned int>(ops, "max_simultaneous_rhs", 0);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  Tracker _t(std::string("fgmres ") + prefix);
		  _t.arity = x.kvdim().at('n');
		  foreachInChuncks(
		    x, y, max_simultaneous_rhs,
		    [=](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		      fgmres(op, prec, x, y, max_basis_size, tol, max_its, error_if_not_converged,
			     ortho_each_its, max_residual_updates, false /* no init guess */, verb,
			     prefix);
		    },
		    'n');
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.imgLayout,
		op.domLayout,
		DenseOperator(),
		op.preferred_col_ordering,
		false /* no Kronecker blocking */};
      }

      /// Returns a BICGSTAB solver
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getBicgstabSolver(Operator<NOp, COMPLEX> op, const Options& ops,
					       Operator<NOp, COMPLEX> prec_)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup bicgstab ") + prefix);

	// Get preconditioner
	Operator<NOp, COMPLEX> prec;
	Maybe<const Options&> precOps = getOptionsMaybe(ops, "prec");
	if (precOps && prec_)
	  throw std::runtime_error(
	    "getBicgstabSolver: invalid `prec` tag: the solver got a preconditioner already");
	if (precOps)
	  prec = getSolver(op, precOps.getSome());
	else if (prec_)
	  prec = prec_;

	// Get the remainder options
	double tol = getOption<double>(ops, "tol", 0.0);
	unsigned int max_its = getOption<unsigned int>(ops, "max_its", 0);
	if (max_its == 0 && tol <= 0)
	  ops.throw_error("set either `tol` or `max_its`");
	bool error_if_not_converged = getOption<bool>(ops, "error_if_not_converged", true);
	unsigned int max_residual_updates = getOption<unsigned int>(ops, "max_residual_updates", 0);
	unsigned int max_simultaneous_rhs = getOption<unsigned int>(ops, "max_simultaneous_rhs", 0);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  Tracker _t(std::string("bicgstab ") + prefix);
		  _t.arity = x.kvdim().at('n');
		  foreachInChuncks(
		    x, y, max_simultaneous_rhs,
		    [=](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		      bicgstab(op, prec, x, y, tol, max_its, error_if_not_converged,
			       max_residual_updates, false /* no init guess */, verb, prefix);
		    },
		    'n');
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.imgLayout,
		op.domLayout,
		DenseOperator(),
		op.preferred_col_ordering,
		false /* no Kronecker blocking */};
      }

      /// Returns a MR solver
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getMRSolver(Operator<NOp, COMPLEX> op, const Options& ops,
					 Operator<NOp, COMPLEX> prec_)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup mg ") + prefix);

	// Get preconditioner
	Operator<NOp, COMPLEX> prec;
	Maybe<const Options&> precOps = getOptionsMaybe(ops, "prec");
	if (precOps && prec_)
	  throw std::runtime_error(
	    "getMRSolver: invalid `prec` tag: the solver got a preconditioner already");
	if (precOps)
	  prec = getSolver(op, precOps.getSome());
	else if (prec_)
	  prec = prec_;

	// Get the remainder options
	double tol = getOption<double>(ops, "tol", 0.0);
	unsigned int max_its = getOption<unsigned int>(ops, "max_its", 0);
	if (max_its == 0 && tol <= 0)
	  ops.throw_error("set either `tol` or `max_its`");
	bool error_if_not_converged = getOption<bool>(ops, "error_if_not_converged", true);
	unsigned int max_residual_updates = getOption<unsigned int>(ops, "max_residual_updates", 0);
	unsigned int max_simultaneous_rhs = getOption<unsigned int>(ops, "max_simultaneous_rhs", 0);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  Tracker _t(std::string("mr ") + prefix);
		  _t.arity = x.kvdim().at('n');
		  foreachInChuncks(
		    x, y, max_simultaneous_rhs,
		    [=](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		      mr(op, prec, x, y, tol, max_its, error_if_not_converged, max_residual_updates,
			 false /* no init guess */, verb, prefix);
		    },
		    'n');
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.imgLayout,
		op.domLayout,
		DenseOperator(),
		op.preferred_col_ordering,
		false /* no Kronecker blocking */};
      }

      /// Returns a GCR solver
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getGCRSolver(Operator<NOp, COMPLEX> op, const Options& ops,
					  Operator<NOp, COMPLEX> prec_)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup gcr ") + prefix);

	// Get preconditioner
	Operator<NOp, COMPLEX> prec;
	Maybe<const Options&> precOps = getOptionsMaybe(ops, "prec");
	if (precOps && prec_)
	  throw std::runtime_error(
	    "getMRSolver: invalid `prec` tag: the solver got a preconditioner already");
	if (precOps)
	  prec = getSolver(op, precOps.getSome());
	else if (prec_)
	  prec = prec_;

	// Get the remainder options
	unsigned int max_basis_size = getOption<unsigned int>(ops, "max_basis_size", 3);
	double tol = getOption<double>(ops, "tol", 0.0);
	unsigned int max_its = getOption<unsigned int>(ops, "max_its", 0);
	if (max_its == 0 && tol <= 0)
	  ops.throw_error("set either `tol` or `max_its`");
	bool error_if_not_converged = getOption<bool>(ops, "error_if_not_converged", true);
	unsigned int max_residual_updates = getOption<unsigned int>(ops, "max_residual_updates", 0);
	unsigned int max_simultaneous_rhs = getOption<unsigned int>(ops, "max_simultaneous_rhs", 0);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  Tracker _t(std::string("gcr ") + prefix);
		  _t.arity = x.kvdim().at('n');
		  foreachInChuncks(
		    x, y, max_simultaneous_rhs,
		    [=](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		      gcr(op, prec, x, y, max_basis_size, tol, max_its, error_if_not_converged,
			  max_residual_updates, false /* no init guess */, verb, prefix);
		    },
		    'n');
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.imgLayout,
		op.domLayout,
		DenseOperator(),
		op.preferred_col_ordering,
		false /* no Kronecker blocking */};
      }

      /// Return the block diagonal of an operator
      /// \param op: operator to extract the block diagonal
      /// \param block_labels: labels that compose the blocks (and will be the rows)
      /// \param m: map from the operator block labels to the new labels for the columns
      /// \return: tensor with the block labels as rows and the renamed labels as columns

      template <std::size_t N, std::size_t NOp, typename COMPLEX>
      Tensor<N, COMPLEX>
      getBlockDiag(const Operator<NOp, COMPLEX>& op, const std::string& block_labels,
		   const remap& m,
		   BlockingAsSparseDimensions blockingAsSparseDimensions = ConsiderBlockingSparse)
      {
	// Get a sparse tensor representation of the operator
	remap rd;
	SpTensor<NOp, NOp, COMPLEX> sop;
	int power = 0; // clone only the block diagonal
	ColOrdering coBlk = RowMajor;
	const std::string prefix = "block diag";
	if (blockingAsSparseDimensions == ConsiderBlockingSparse)
	{
	  auto t = cloneOperatorToSpTensor(op, power, coBlk, false, prefix);
	  sop = t.first;
	  rd = t.second;
	}
	else
	{
	  auto t = cloneUnblockedOperatorToSpTensor(op, power, coBlk, false, prefix);
	  sop = t.first;
	  rd = t.second;
	}

	// Create tensor to return
	const std::string cols = update_order(block_labels, m);
	std::string order =
	  block_labels + cols + detail::remove_dimensions(op.d.order, block_labels);
	auto dims = op.d.kvdim();
	for (const auto& it : m)
	  dims[it.second] = dims[it.first];
	Tensor<N, COMPLEX> r = op.d.template make_compatible<N>(order, dims);

	// Copy the blocks
	remap m_cols = getNewLabels(cols, sop.data.order);
	remap m_blk;
	for (const auto& it : m)
	  m_blk[rd.at(it.first)] = m_cols.at(it.second);
	sop.data.rename_dims(m_blk).copyTo(r.rename_dims(m_cols));

	// Return tensor
	return r;
      }

      enum class SpinSplitting {
	None,	   // one spin output
	Chirality, // two spin output
	Full	   // four spin output (homoiconic coarse operator)
      };

      enum class NullVectorsNormalization {
	None,	 // don't normalize them
	Full,	 // normalize the whole vectors
	Blocking // normalize each block
      };

      /// Returns the prolongator constructed from the given operator
      /// \param op: operator to make the prolongator of
      /// \param num_null_vecs: total number of random vectors used to construct the prolongator
      /// \param max_num_null_vecs: maximum number of vectors to consider at once
      /// \param null_vectors_normalization: normalization to use before analysing frequencies
      /// \param mg_blocking: lattice dimension blocking for the prolongator
      /// \param layout_blocking: lattice dimension blocking for layout
      /// \param spin_splitting: treatment of the spin
      /// \param null_vecs_ops: solver options for the computing the null vecs
      /// \param solver_space: target space

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX>
      getMGProlongator(const Operator<NOp, COMPLEX>& op, unsigned int num_null_vecs,
		       unsigned int max_num_null_vecs, unsigned int num_colors,
		       NullVectorsNormalization null_vectors_normalization,
		       const std::map<char, unsigned int>& mg_blocking,
		       const std::map<char, unsigned int>& layout_blocking,
		       SpinSplitting spin_splitting, const Options& null_vecs_ops,
		       SolverSpace solverSpace)
      {
	Tracker _t("setup mg prolongator");
	detail::log(1, "starting getMGProlongator");

	// Check that the blocking divide the lattice sizes. For an operator with support only
	// on the even sites, make sure that the blocking on x is divisible by two
	auto opdims = op.d.kvdim();
	int X = op.imgLayout == EvensOnlyLayout ? 2 : opdims.at('X');
	bool x_blocking_divide_X = (mg_blocking.at('x') % X == 0);
	if (!x_blocking_divide_X && op.imgLayout == EvensOnlyLayout)
	  throw std::runtime_error(
	    "When using even-odd preconditioning, the blocking on x should be divisible by 2");
	for (const auto it : getNatLatticeDims(opdims, op.imgLayout))
	{
	  if (it.second % mg_blocking.at(it.first) != 0)
	    throw std::runtime_error("The operator dimensions are not divisible by the blocking");
	}

	// Check that there is enough spins to do splitting
	int ns = opdims.at('s');
	if (ns != 1 && ns != 2 && ns != Ns)
	  throw std::runtime_error("Error in getMGProlongator: Unsupported spin number");
	if (ns == 1)
	  spin_splitting = SpinSplitting::None;

	// Ignore max_num_null_vecs if it is smaller than num_colors
	max_num_null_vecs = std::max(max_num_null_vecs, num_colors + 1);

	const Operator<NOp, COMPLEX> null_solver =
	  getSolver(op, getOptions(null_vecs_ops, "solver"));
	auto eigensolver_ops = getOptionsMaybe(null_vecs_ops, "eigensolver");
	if (eigensolver_ops && max_num_null_vecs < num_null_vecs)
	  eigensolver_ops.getSome().throw_error(
	    "unsupported `eigensolver` when `max_num_null_vecs` is smaller than `num_null_vecs`");

	auto get_null_vecs_fn = [&](unsigned int num_null_vecs) {
	  Tensor<NOp + 1, COMPLEX> nv;
	  if (!eigensolver_ops)
	  {
	    // Create the random initial guesses to be used in solving Ax=0
	    auto b = op.template make_compatible_img<NOp + 1>("n", {{'n', num_null_vecs}});
	    if (solverSpace == FullSpace)
	    {
	      nrand(b);
	    }
	    else
	    {
	      b.set_zero();
	      nrand(b.kvslice_from_size({}, {{'X', 1}}));
	    }

	    // Solve Ax=0 with the random initial guesses
	    nv = null_solver(op(b));
	    b.scale(-1).addTo(nv);
	    b.release();
	  }
	  else
	  {
	    // Compute the eigenpairs of inv(op) * g5, the right singular vectors of op
	    auto eigensolver = SB::getInexactEigensolverGD(null_solver, eigensolver_ops.getSome());
	    double tol = getOption<double>(null_vecs_ops, "tol", 0.1);
	    auto values_vectors = eigensolver(num_null_vecs, tol);
	    nv = std::get<1>(values_vectors);
	  }
	  return nv;
	};

	Operator<NOp, COMPLEX> V;
	auto opdims_nat = opdims;
	if (opdims_nat.at('X') != 1)
	{
	  opdims_nat['x'] *= opdims_nat['X'];
	  opdims_nat['X'] = 1;
	}
	if (spin_splitting != SpinSplitting::Full)
	{
	  Tensor<NOp + 1 + 5, COMPLEX> nv_blk;
	  std::map<std::string, std::string> m_blk, m_blk_rev, m_blk_nv;
	  for (unsigned int nv_created = 0; nv_created < num_null_vecs;)
	  {
	    unsigned int current_nv = std::min(
	      num_null_vecs - nv_created, max_num_null_vecs - (nv_created == 0 ? 0 : num_colors));

	    // Get the null vectors
	    auto nv = get_null_vecs_fn(current_nv);

	    // Do chirality splitting nv2 = [ nv * gpos, nv * gneg ]
	    auto nv2 = nv;
	    if (spin_splitting == SpinSplitting::Chirality)
	    {
	      nv2 = nv.like_this(none, {{'n', current_nv * 2}});
	      auto g5 = getGamma5<COMPLEX>(ns, OnHost, nv.dist), g5pos = g5.cloneOn(OnHost),
		   g5neg = g5.cloneOn(OnHost);
	      for (int i = 0; i < ns; ++i) // make diagonal entries of gpos all positive or zero
		g5pos.set({{i, i}}, g5.get({{i, i}}) + COMPLEX{1});
	      for (int i = 0; i < ns; ++i) // make diagonal entries of gneg all negative or zero
		g5neg.set({{i, i}}, g5.get({{i, i}}) - COMPLEX{1});
	      nv2.kvslice_from_size({}, {{'n', current_nv}})
		.contract(g5pos, {{'i', 's'}}, NotConjugate, nv, {{'s', 'j'}}, NotConjugate);
	      nv2.kvslice_from_size({{'n', current_nv}}, {{'n', current_nv}})
		.contract(g5neg, {{'i', 's'}}, NotConjugate, nv, {{'s', 'j'}}, NotConjugate);
	    }
	    nv.release();

	    // Do the blocking, which encompasses the following transformations:
	    //  X0x -> WX0x, 1y -> Y1y, 2z -> Z2z, 3t -> T3t,
	    // where output X is a singlet dimension, and W,Y,Z, and T have size mg_blocking,
	    // and output's 0,1,2, and 3 have size layout_blocking, and the output's x,y,z, and t
	    // have the remaining

	    m_blk_nv = {{"X0x", "WX0x"}, {"1y", "Y1y"}, {"2z", "Z2z"}, {"3t", "T3t"},
			{"n", "cs"},	 {"c", "C"},	{"s", "S"}};
	    m_blk = {{"X0x", "WX0x"}, {"1y", "Y1y"}, {"2z", "Z2z"},
		     {"3t", "T3t"},   {"c", "C"},    {"s", "S"}};
	    m_blk_rev = {{"WX0x", "X0x"}, {"Y1y", "1y"}, {"Z2z", "2z"},
			 {"T3t", "3t"},	  {"C", "c"},	 {"S", "s"}};

	    if (!x_blocking_divide_X)
	      nv2 = toNaturalOrdering(nv2);
	    auto nv_blk0 = nv2.template reshape_dimensions<NOp + 1 + 5>(
	      m_blk_nv,
	      {{'X', 1}, // we don't do even-odd layout on the coarse operator space
	       {'W', mg_blocking.at('x') / (op.imgLayout == EvensOnlyLayout ? 2 : 1)},
	       {'Y', mg_blocking.at('y')},
	       {'Z', mg_blocking.at('z')},
	       {'T', mg_blocking.at('t')},
	       {'0', layout_blocking.at('x')},
	       {'1', layout_blocking.at('y')},
	       {'2', layout_blocking.at('z')},
	       {'3', layout_blocking.at('t')},
	       {'c', current_nv}},
	      true);
	    nv2.release();

	    // User even-odd ordering for nv_blk0
	    if (nv_blk0.kvdim().at('0') != 1)
	      throw std::runtime_error("getMGProlongator: unsupported blocking on the x direction");

	    // Cut down the extra null vectors
	    if (num_null_vecs > num_colors)
	    {
	      // Pre-normalize the vectors if the user asks
	      if (null_vectors_normalization == NullVectorsNormalization::Full)
	      {
		nv_blk0 =
		  vecnorm<NOp + 1 + 5 - 1>(nv_blk0, detail::union_dimensions(nv_blk0.order, "", "c"));
	      }
	      else if (null_vectors_normalization == NullVectorsNormalization::Blocking)
	      {
		nv_blk0 = vecnorm<6>(nv_blk0, "WYZTSC");
	      }

	      // Extend nv_blk with nv_blk0
	      nv_blk = concat(nv_blk, nv_blk0, 'c');
	      nv_blk0.release();

	      // Compute the left singular vectors
	      // NOTE: use the label '^' as the singular triple index
	      auto svd_result = svd<NOp + 1 + 5, NOp, NOp + 1>(nv_blk, "WYZTSC", "c", '^');
	      const auto& svd_u = std::get<0>(svd_result); // get the left singular vectors
	      auto svd_s =
		std::get<1>(svd_result).template cast<COMPLEX>(); // get the singular values
	      std::get<2>(svd_result).release(); // release the right singular vectors

	      if (nv_created + current_nv >= num_null_vecs)
	      {
		// If this is the last iteration, take the most common `num_colors` directions
		nv_blk = svd_u.rename_dims({{'^', 'c'}})
			   .kvslice_from_size({}, {{'c', num_colors}})
			   .clone();
	      }
	      else
	      {
		// Otherwise, take the first `num_colors`+1 most frequent directions and
		// scale them with the frequency
		nv_blk =
		  contract<NOp + 1 + 5>(svd_u.kvslice_from_size({}, {{'^', num_colors + 1}}),
					svd_s.kvslice_from_size({}, {{'^', num_colors + 1}}), "")
		    .rename_dims({{'^', 'c'}});
	      }
	    }
	    else
	    {
	      nv_blk = nv_blk0;
	    }

	    // Update loop var
	    nv_created += current_nv;
	  }

	  // Do the orthogonalization on each block and chirality
	  // NOTE: don't check the orthogonality level, it will replicate the inner products on all processes
	  ortho<6, 1>(nv_blk, "X0123xyzts", "WYZTSC", "c", 4, CheckOrthoLevel::dontCheckOrthoLevel,
		      JustSummary, "prolongator");

	  // Return the operator
	  auto nv_blk_eo_dim = nv_blk.kvdim();
	  if (X == 2 && nv_blk_eo_dim.at('x') % 2 == 0 &&
	      (nv_blk_eo_dim.at('y') == 1 || nv_blk_eo_dim.at('y') % 2 == 0) &&
	      (nv_blk_eo_dim.at('z') == 1 || nv_blk_eo_dim.at('z') % 2 == 0) &&
	      (nv_blk_eo_dim.at('t') == 1 || nv_blk_eo_dim.at('t') % 2 == 0))
	  {
	    nv_blk_eo_dim['x'] /= 2;
	    nv_blk_eo_dim['X'] = 2;
	  }
	  Tensor<NOp, COMPLEX> d = op.d.like_this(none, nv_blk_eo_dim).make_eg(), i = op.i;
	  V = Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto y0 = contract<NOp + 1 + 4>(nv_blk, toNaturalOrdering(x), "cs")
			  .template reshape_dimensions<NOp + 1>(
			    m_blk_rev, x_blocking_divide_X ? opdims : opdims_nat, true);
	      if (!x_blocking_divide_X)
		toEvenOddOrdering(y0).copyTo(y);
	      else
		y0.copyTo(y);
	    },
	    d,
	    i,
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = x_blocking_divide_X ? x : toNaturalOrdering(x);
	      auto x_blk = x0.template reshape_dimensions<NOp + 1 + 4>(m_blk, nv_blk.kvdim(), true);
	      if (nv_blk_eo_dim.at('X') == 1)
		contract<NOp + 1>(nv_blk.conj(), x_blk, "WYZTSC", CopyTo, y);
	      else
		toEvenOddOrdering(contract<NOp + 1>(nv_blk.conj(), x_blk, "WYZTSC")).copyTo(y);
	    },
	    op.order_t,
	    nv_blk_eo_dim.at('X') == 2 ? XEvenOddLayout : NaturalLayout,
	    op.imgLayout,
	    getNeighborsAfterBlocking(mg_blocking, op.d.kvdim(), op.neighbors, op.imgLayout),
	    op.preferred_col_ordering,
	    false /* no Kronecker format */};
	}
	else
	{
	  Tensor<NOp + 4, COMPLEX> nv_blk;
	  std::map<std::string, std::string> m_blk, m_blk_rev, m_blk_nv;
	  for (unsigned int nv_created = 0; nv_created < num_null_vecs;)
	  {
	    unsigned int current_nv = std::min(
	      num_null_vecs - nv_created, max_num_null_vecs - (nv_created == 0 ? 0 : num_colors));

	    // Get the null vectors
	    auto nv = get_null_vecs_fn(current_nv);

	    // Do the blocking, which encompasses the following transformations:
	    //  X0x -> WX0x, 1y -> Y1y, 2z -> Z2z, 3t -> T3t,
	    // where output X is a singlet dimension, and W,Y,Z, and T have size mg_blocking,
	    // and output's 0,1,2, and 3 have size layout_blocking, and the output's x,y,z, and t
	    // have the remaining.
	    // Also, enforce the values to be the same across different spins so that the coarse operator links
	    // are also the tensor product of a spin matrix and a color matrix

	    m_blk_nv = {{"X0x", "WX0x"}, {"1y", "Y1y"}, {"2z", "Z2z"},
			{"3t", "T3t"},	 {"ns", "c"},	{"c", "C"}};
	    m_blk = {{"X0x", "WX0x"}, {"1y", "Y1y"}, {"2z", "Z2z"}, {"3t", "T3t"}, {"c", "C"}};
	    m_blk_rev = {{"WX0x", "X0x"}, {"Y1y", "1y"}, {"Z2z", "2z"}, {"T3t", "3t"}, {"C", "c"}};

	    if (!x_blocking_divide_X)
	      nv = toNaturalOrdering(nv);
	    auto nv_blk0 = nv.template reshape_dimensions<NOp + 4>(
	      m_blk_nv,
	      {{'X', 1}, // we don't do even-odd layout on the coarse operator space
	       {'W', mg_blocking.at('x') / (op.imgLayout == EvensOnlyLayout ? 2 : 1)},
	       {'Y', mg_blocking.at('y')},
	       {'Z', mg_blocking.at('z')},
	       {'T', mg_blocking.at('t')},
	       {'0', layout_blocking.at('x')},
	       {'1', layout_blocking.at('y')},
	       {'2', layout_blocking.at('z')},
	       {'3', layout_blocking.at('t')},
	       {'c', current_nv * ns}},
	      true);
	    nv.release();

	    // User even-odd ordering for nv_blk0
	    if (nv_blk0.kvdim().at('0') != 1)
	      throw std::runtime_error("getMGProlongator: unsupported blocking on the x direction");

	    // Cut down the extra null vectors
	    if (num_null_vecs > num_colors)
	    {
	      // Pre-normalize the vectors if the user asks
	      if (null_vectors_normalization == NullVectorsNormalization::Full)
	      {
		nv_blk0 =
		  vecnorm<NOp + 4 -1>(nv_blk0, detail::union_dimensions(nv_blk0.order, "", "c"));
	      }
	      else if (null_vectors_normalization == NullVectorsNormalization::Blocking)
	      {
		nv_blk0 = vecnorm<5>(nv_blk0, "WYZTC");
	      }

	      // Extend nv_blk with nv_blk0
	      nv_blk = concat(nv_blk, nv_blk0, 'c');
	      nv_blk0.release();

	      // Compute the left singular vectors
	      // NOTE: use the label '^' as the singular triple index
	      auto svd_result = svd<NOp + 4, NOp - 1, NOp>(nv_blk, "WYZTC", "c", '^');
	      const auto& svd_u = std::get<0>(svd_result); // get the left singular vectors
	      auto svd_s =
		std::get<1>(svd_result).template cast<COMPLEX>(); // get the singular values
	      std::get<2>(svd_result).release(); // release the right singular vectors

	      if (nv_created + current_nv >= num_null_vecs)
	      {
		// If this is the last iteration, take the most common `num_colors` directions
		nv_blk = svd_u.rename_dims({{'^', 'c'}})
			   .kvslice_from_size({}, {{'c', num_colors}})
			   .clone();
	      }
	      else
	      {
		// Otherwise, take the first `num_colors`+1 most frequent directions
		nv_blk = contract<NOp + 4>(svd_u.kvslice_from_size({}, {{'^', num_colors + 1}}),
					   svd_s.kvslice_from_size({}, {{'^', num_colors + 1}}), "")
			   .rename_dims({{'^', 'c'}});
	      }
	    }
	    else
	    {
	      nv_blk = nv_blk0;
	    }

	    // Update loop var
	    nv_created += current_nv;
	  }

	  // Do the orthogonalization on each block and chirality
	  // NOTE: don't check the orthogonality level, it will replicate the inner products on all processes
	  ortho<5, 1>(nv_blk, "X0123xyzt", "WYZTC", "c", 4, CheckOrthoLevel::dontCheckOrthoLevel,
		      JustSummary, "prolongator");

	  // Return the operator
	  auto nv_blk_eo_dim = nv_blk.kvdim();
	  if (X == 2 && nv_blk_eo_dim.at('x') % 2 == 0 &&
	      (nv_blk_eo_dim.at('y') == 1 || nv_blk_eo_dim.at('y') % 2 == 0) &&
	      (nv_blk_eo_dim.at('z') == 1 || nv_blk_eo_dim.at('z') % 2 == 0) &&
	      (nv_blk_eo_dim.at('t') == 1 || nv_blk_eo_dim.at('t') % 2 == 0))
	  {
	    nv_blk_eo_dim['x'] /= 2;
	    nv_blk_eo_dim['X'] = 2;
	  }
	  Tensor<NOp, COMPLEX> d = op.d.like_this(none, nv_blk_eo_dim), i = op.i;
	  V = Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto y0 = contract<NOp + 1 + 4>(nv_blk, toNaturalOrdering(x), "c")
			  .template reshape_dimensions<NOp + 1>(
			    m_blk_rev, x_blocking_divide_X ? opdims : opdims_nat, true);
	      if (!x_blocking_divide_X)
		toEvenOddOrdering(y0).copyTo(y);
	      else
		y0.copyTo(y);
	    },
	    d,
	    i,
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = x_blocking_divide_X ? x : toNaturalOrdering(x);
	      auto x_blk = x0.template reshape_dimensions<NOp + 1 + 4>(m_blk, nv_blk.kvdim(), true);
	      if (nv_blk_eo_dim.at('X') == 1)
		contract<NOp + 1>(nv_blk.conj(), x_blk, "WYZTC", CopyTo, y);
	      else
		toEvenOddOrdering(contract<NOp + 1>(nv_blk.conj(), x_blk, "WYZTC")).copyTo(y);
	    },
	    op.order_t,
	    nv_blk_eo_dim.at('X') == 2 ? XEvenOddLayout : NaturalLayout,
	    op.imgLayout,
	    getNeighborsAfterBlocking(mg_blocking, op.d.kvdim(), op.neighbors, op.imgLayout),
	    op.preferred_col_ordering,
	    true /* using Kronecker format */};
	}

	// Test that the transpose of the prolongator is consistent
	{
	  auto xd = V.template make_compatible_dom<NOp + 2>("nk", {{'n', 4}, {'k', 1}});
	  auto xi = V.template make_compatible_img<NOp + 1>("m", {{'m', 4}});
	  nrand(xd);
	  nrand(xi);
	  auto xid = contract<3>(xi.conj(), V(xd), V.i.order).make_sure(none, OnHost, OnMaster);
	  auto xdi =
	    contract<3>(xd.conj(), V.tconj()(xi), V.d.order).make_sure(none, OnHost, OnMaster);
	  xid.conj().scale(-1).addTo(xdi);
	  double err = std::fabs(norm<1>(xdi, "k").get(Coor<1>{0}));
	  double rel = std::fabs(norm<1>(xid, "k").get(Coor<1>{0}));
	  if (err > rel * 1e-4)
	    throw std::runtime_error("The prolongator has an inconsistent transpose");
	}

	return V;
      }

      /// Return a list of destroy callbacks after setting a solver

      inline std::vector<DestroyFun>& getOperatorsCacheDestroyList()
      {
	static std::vector<DestroyFun> list;
	return list;
      }

      /// Call the destroy callbacks set up in `getOperatorsCacheDestroyList`
      inline void cleanOperatorsCache()
      {
	for (const auto& f : getOperatorsCacheDestroyList())
	  f();
      }

      /// Tuple storing the operator even-odd and odd-even parts and the block diagonal and its inverse

      template <std::size_t NOp, typename COMPLEX>
      using EvenOddOperatorParts = std::tuple<Operator<NOp, COMPLEX>,	///< op_eo
					      Operator<NOp, COMPLEX>,	///< op_oe
					      Tensor<NOp + 2, COMPLEX>, ///< opDiagE
					      Tensor<NOp + 2, COMPLEX>, ///< opDiagO
					      Tensor<NOp + 2, COMPLEX>, ///< opInvDiagE
					      Tensor<NOp + 2, COMPLEX>	///< opInvDiagO
					      >;

      /// Return the cache of block diagonals for even-odd operators generated by getEvenOddPrec
      /// NOTE: this one is destroyed by `cleanOperatorsCache`

      template <std::size_t NOp, typename COMPLEX>
      std::map<void*, EvenOddOperatorParts<NOp, COMPLEX>>& getEvenOddOperatorsPartsCache()
      {
	static std::map<void*, EvenOddOperatorParts<NOp, COMPLEX>> m = []() {
	  getOperatorsCacheDestroyList().push_back(
	    []() { getEvenOddOperatorsPartsCache<NOp, COMPLEX>().clear(); });
	  return std::map<void*, EvenOddOperatorParts<NOp, COMPLEX>>{};
	}();
	return m;
      }

      /// Return a cache for the prolongators
      /// NOTE: this one is destroyed by `cleanOperatorsCache`

      template <std::size_t NOp, typename COMPLEX>
      static std::map<std::string, Operator<NOp, COMPLEX>>& getProlongatorCache()
      {
	static std::map<std::string, Operator<NOp, COMPLEX>> m = [] {
	  getDestroyList().push_back([] { getProlongatorCache<NOp, COMPLEX>().clear(); });
	  return std::map<std::string, Operator<NOp, COMPLEX>>{};
	}();
	return m;
      }

      /// Return the prolongator and the coarse operator and spin splitting
      template <std::size_t NOp, typename COMPLEX>
      std::tuple<Operator<NOp, COMPLEX>, Operator<NOp, COMPLEX>, SpinSplitting>
      getProlongatorAndCoarse(Operator<NOp, COMPLEX> op, const Options& ops,
			      SolverSpace solverSpace)
      {
	Tracker _t("setup mg prolongator and coarse");

	// Get prolongator, V
	unsigned int num_null_vecs = getOption<unsigned int>(ops, "num_null_vecs");
	unsigned int max_num_null_vecs =
	  getOption<unsigned int>(ops, "max_num_null_vecs", num_null_vecs);
	unsigned int num_colors = getOption<unsigned int>(ops, "num_colors", num_null_vecs);
	std::vector<unsigned int> mg_blocking_v =
	  getOption<std::vector<unsigned int>>(ops, "blocking");
	if (mg_blocking_v.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getMGPrec: the blocking should be a vector with four elements");
	std::map<char, unsigned int> mg_blocking{{'x', mg_blocking_v[0]},
						 {'y', mg_blocking_v[1]},
						 {'z', mg_blocking_v[2]},
						 {'t', mg_blocking_v[3]}};
	std::vector<unsigned int> layout_blocking_v = getOption<std::vector<unsigned int>>(
	  ops, "coarse_layout_blocking", std::vector<unsigned int>{{1, 1, 1, 1}});
	if (layout_blocking_v.size() != Nd)
	  ops.getValue("coarse_layout_blocking")
	    .throw_error("getMGPrec: the blocking should be a vector with four elements");
	std::map<char, unsigned int> layout_blocking{{'x', layout_blocking_v[0]},
						     {'y', layout_blocking_v[1]},
						     {'z', layout_blocking_v[2]},
						     {'t', layout_blocking_v[3]}};
	static const std::map<std::string, SpinSplitting> m_spin_splitting{
	  {"none", SpinSplitting::None},
	  {"chirality_splitting", SpinSplitting::Chirality},
	  {"full", SpinSplitting::Full}};
	SpinSplitting spin_splitting = getOption<SpinSplitting>(
	  ops, "spin_splitting", m_spin_splitting, SpinSplitting::Chirality);

	static const std::map<std::string, NullVectorsNormalization> m_null_vecs_normalization{
	  {"none", NullVectorsNormalization::None},
	  {"full", NullVectorsNormalization::Full},
	  {"blocking", NullVectorsNormalization::Blocking}};
	NullVectorsNormalization null_vecs_normalization = getOption<NullVectorsNormalization>(
	  ops, "null_vecs_norm", m_null_vecs_normalization, NullVectorsNormalization::None);

	// Grab the prolongator from cache if the user name it
	Operator<NOp, COMPLEX> V;
	std::string prolongator_id = getOption<std::string>(ops, "prolongator_id", "");
	if (prolongator_id.size() == 0 ||
	    getProlongatorCache<NOp, COMPLEX>().count(prolongator_id) == 0)
	{
	  V = getMGProlongator(op, num_null_vecs, max_num_null_vecs, num_colors,
			       null_vecs_normalization, mg_blocking, layout_blocking,
			       spin_splitting, getOptions(ops, "null_vecs"), solverSpace);
	  if (prolongator_id.size() > 0)
	    getProlongatorCache<NOp, COMPLEX>()[prolongator_id] = V;
	}
	else
	{
	  QDPIO::cout << "Found prolongator for id " << prolongator_id << std::endl;
	  V = getProlongatorCache<NOp, COMPLEX>().at(prolongator_id);
	}

	// Compute the coarse operator, either V' * op * V or V' * op * g5 * V
	unsigned int create_coarse_max_rhs =
	  getOption<unsigned int>(ops, "create_coarse_max_rhs", 0);
	ColOrdering co = getOption<ColOrdering>(ops, "operator_ordering", getColOrderingMap(),
						op.preferred_col_ordering);
	ColOrdering co_blk =
	  getOption<ColOrdering>(ops, "operator_block_ordering", getColOrderingMap(), RowMajor);
	int ns = op.d.kvdim().at('s');
	auto g5 = getGamma5<COMPLEX>(ns, op.d.getDev(), op.d.dist);
	const Operator<NOp, COMPLEX> op_c = cloneOperator(
	  Operator<NOp, COMPLEX>{
	    [&](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      foreachInChuncks(x, y, create_coarse_max_rhs,
			       [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
				 if (spin_splitting != SpinSplitting::None || ns == 1)
				 {
				   V.tconj()(op(V(x)), y);
				 }
				 else
				 {
				   V.tconj()(
				     contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), op(V(x)), "s")
				       .rename_dims({{'i', 's'}}),
				     y);
				 }
			       });
	    },
	    V.d, V.d, nullptr, op.order_t, V.domLayout, V.domLayout, V.neighbors,
	    op.preferred_col_ordering, V.is_kronecker()},
	  co, co_blk, ConsiderBlockingSparse, "coarse");

	return {V, op_c, spin_splitting};
      }

      /// Returns a MG preconditioner.
      ///
      /// It returns an approximation of Op^{-1} = Op^{-1}*Q + Op^{-1}(I-Q), where Q is a projector
      /// on the left singular space of Op. (MG can be derived also from P*Op^{-1} + (I-P)*Op^{-1},
      /// where P is on the right singular space, producing pre-smoothers instead.)
      /// The approximation is constructed using an oblique projector and doing the inversions
      /// approximately:
      ///   1) Q = Op*V*(W'*Op*V)^{-1}*W', where W and V are left and right singular spaces.
      ///   2) [Op^{-1}*Q + Op^{-1}(I-Q)]*x \approx
      ///                       V*solver(W'*Op*V, W'*x) + solver(Op, x - Op*V*solver(W'*Op*V, W'*x))
      /// Note that if Op is \gamma_5-Hermitian, and W=\gamma_5*V, and \gamma_5 commutes with V
      /// (\gamma_5*V = V*\gamma_5^coarse), then the projector reduces to Q=Op*V(V'*Op*V)^{-1}*V',
      /// and, this coarse operator is easier to solve than V'\gamma_5*Op*V.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getMGPrec(Operator<NOp, COMPLEX> op, const Options& ops,
				       Operator<NOp, COMPLEX> prec_, SolverSpace solverSpace)
      {
	if (prec_)
	  throw std::runtime_error("getMGPrec: unsupported input preconditioner");

	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup mg solver ") + prefix);

	// Get prolongator and coarse operator
	auto prolongator_coarse_spin_splitting = getProlongatorAndCoarse(op, ops, solverSpace);
	Operator<NOp, COMPLEX> V = std::get<0>(prolongator_coarse_spin_splitting);
	Operator<NOp, COMPLEX> op_c = std::get<1>(prolongator_coarse_spin_splitting);
	SpinSplitting spin_splitting = std::get<2>(prolongator_coarse_spin_splitting);

	// Get the solver for the projector
	const Operator<NOp, COMPLEX> coarseSolver =
	  getSolver(op_c, getOptions(ops, "solver_coarse"));

	// Get the solver for the smoother
	const Operator<NOp, COMPLEX> opSolver = getSolver(op, getOptions(ops, "solver_smoother"));

	OperatorFun<NOp, COMPLEX> solver;
	int ns = op.d.kvdim().at('s');
	auto g5 = getGamma5<COMPLEX>(ns, op.d.getDev(), op.d.dist);
	if (getEvenOddOperatorsPartsCache<NOp, COMPLEX>().count(op.id.get()) == 0)
	{
	  solver = [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("mg solver ") + prefix);

	    // x0 = g5 * x if !do_chirality_splitting && ns > 1
	    auto x0 = x;
	    if (spin_splitting == SpinSplitting::None && ns > 1)
	    {
	      x0 =
		contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), x, "s").rename_dims({{'i', 's'}});
	    }

	    // y0 = V*solver(V'*Op*V, V'x0)
	    auto y0 = V(coarseSolver(V.tconj()(x0)));

	    // x1 = x - op*y0
	    auto x1 = op(y0.scale(-1));
	    x.addTo(x1);

	    // y = y0 + solver(Op, x1)
	    opSolver(std::move(x1), y);
	    y0.addTo(y);
	  };
	}
	else
	{
	  // Get the block diagonal of the operator with rows cs and columns CS
	  remap m_sc{{'s', 'S'}, {'c', 'C'}};
	  auto t = getEvenOddOperatorsPartsCache<NOp, COMPLEX>().at(op.id.get());
	  Operator<NOp, COMPLEX> op_eo = std::get<0>(t);
	  Operator<NOp, COMPLEX> op_oe = std::get<1>(t);
	  Tensor<NOp + 2, COMPLEX> opDiagE = std::get<2>(t);
	  Tensor<NOp + 2, COMPLEX> opDiagO = std::get<3>(t);
	  solver = [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("mg solver ") + prefix);

	    // x0 = g5 * x if !do_chirality_splitting && ns > 1
	    auto x0 = x;
	    if (spin_splitting == SpinSplitting::None && ns > 1)
	    {
	      x0 =
		contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), x, "s").rename_dims({{'i', 's'}});
	    }

	    // y0 = V*solver(V'*Op*V, V'x0)
	    auto y0 = V(coarseSolver(V.tconj()(x0)));

	    // y1 = x - op*y0
	    auto y0m = y0.scale(-1);
	    auto y0me = y0m.kvslice_from_size({{'X', 0}}, {{'X', 1}}).clone();
	    auto y0mo = y0m.kvslice_from_size({{'X', 1}}, {{'X', 1}}).clone();
	    auto y1 = x.clone();
	    auto y1e = y1.kvslice_from_size({{'X', 0}}, {{'X', 1}});
	    auto y1o = y1.kvslice_from_size({{'X', 1}}, {{'X', 1}});
	    contract(opDiagE, y0me.rename_dims(m_sc), "CS", AddTo, y1e);
	    contract(opDiagO, y0mo.rename_dims(m_sc), "CS", AddTo, y1o);
	    op_oe(y0me).addTo(y1o);
	    op_eo(y0mo).addTo(y1e);

	    // y = y0 + solver(Op, y1)
	    opSolver(std::move(y1), y);
	    y0.addTo(y);
	  };
	}

	// Return the solver
	return {solver,
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.imgLayout,
		op.domLayout,
		DenseOperator(),
		op.preferred_col_ordering,
		false /* no Kronecker blocking */};
      }

      /// Returns a general deflation preconditioner
      ///
      /// It returns an approximation of Op^{-1} = Op^{-1}*Q + Op^{-1}(I-Q), where Q is a projector
      /// on the left singular space of Op. 
      /// The approximation is constructed using an oblique projector and doing the inversions
      /// approximately:
      ///   1) Q = Op*V*(W'*Op*V)^{-1}*W', where W and V are left and right singular spaces.
      ///   2) [Op^{-1}*Q + Op^{-1}(I-Q)]*x \approx
      ///                       V*solver(W'*Op*V, W'*x) + solver(Op, x - Op*V*solver(W'*Op*V, W'*x))
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getProjPrec(Operator<NOp, COMPLEX> op, const Options& ops,
					 Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getProjPrec: unsupported input preconditioner");

	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup proj ") + prefix);

	// Getting a projector
	auto proj = getProjector(op, getOptions(ops, "proj"));

	// Get the solver for the smoother
	const Operator<NOp, COMPLEX> opSolver = getSolver(op, getOptions(ops, "solver_smoother"));

	OperatorFun<NOp, COMPLEX> solver;
	if (getEvenOddOperatorsPartsCache<NOp, COMPLEX>().count(op.id.get()) == 0)
	{
	  solver = [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("proj solver ") + prefix);

	    // y0 = V*inv(U'*Op*V) U' * x
	    auto y0 = proj.V_inv_Ut(x);

	    // y1 = (I - Q)*x = x - op*y0
	    auto y1 = op(y0.scale(-1));
	    x.addTo(y1);

	    // y = y0 + solver(Op, y1)
	    opSolver(std::move(y1), y);
	    y0.addTo(y);
	  };
	}
	else
	{
	  // Get the block diagonal of the operator with rows cs and columns CS
	  remap m_sc{{'s', 'S'}, {'c', 'C'}};
	  auto t = getEvenOddOperatorsPartsCache<NOp, COMPLEX>().at(op.id.get());
	  Operator<NOp, COMPLEX> op_eo = std::get<0>(t);
	  Operator<NOp, COMPLEX> op_oe = std::get<1>(t);
	  Tensor<NOp + 2, COMPLEX> opDiagE = std::get<2>(t);
	  Tensor<NOp + 2, COMPLEX> opDiagO = std::get<3>(t);
	  solver = [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("proj solver ") + prefix);

	    // y0 = V*inv(U'*Op*V) U' * x
	    auto y0 = proj.V_inv_Ut(x);

	    // y1 = x - op*y0
	    auto y0m = y0.scale(-1);
	    auto y0me = y0m.kvslice_from_size({{'X', 0}}, {{'X', 1}}).clone();
	    auto y0mo = y0m.kvslice_from_size({{'X', 1}}, {{'X', 1}}).clone();
	    auto y1 = x.clone();
	    auto y1e = y1.kvslice_from_size({{'X', 0}}, {{'X', 1}});
	    auto y1o = y1.kvslice_from_size({{'X', 1}}, {{'X', 1}});
	    contract(opDiagE, y0me.rename_dims(m_sc), "CS", AddTo, y1e);
	    contract(opDiagO, y0mo.rename_dims(m_sc), "CS", AddTo, y1o);
	    op_oe(y0me).addTo(y1o);
	    op_eo(y0mo).addTo(y1e);

	    // y = y0 + solver(Op, y1)
	    opSolver(std::move(y1), y);
	    y0.addTo(y);
	  };
	}

	// Return the solver
	return {solver,
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.imgLayout,
		op.domLayout,
		DenseOperator(),
		op.preferred_col_ordering,
		false /* no Kronecker blocking */};
      }

      /// Returns an even-odd preconditioner.
      ///
      /// It returns an approximation of Op^{-1} by splitting the rows and columns into the two
      /// (red-black, even-odd) and doing:
      ///
      ///  Op^{-1} = R^{-1} * A^{-1} * L^{-1} ->
      ///
      /// [ Op_ee Op_eo ]^{-1} = [       I            0 ] * [ Op_ee-Op_eo*Op_oo^{-1}*Op_oe   0   ]^{-1} *
      /// [ Op_oe Op_oo ]        [ -Op_oo^{-1}*Op_oe  I ]   [              0               Op_oo ]
      ///
      ///                        * [ I  -Op_eo*Op_oo^{-1} ]
      ///                          [ 0         I          ]
      ///
      /// The matrix Op_oo^{-1} is block diagonal and is computed directly, while
      /// (Op_ee-Op_eo*Op_oo^{-1}*Op_oe)^{-1} is approximated by an iterative solver. Note that
      /// the residual norm while solving A_ee=Op_ee-Op_eo*Op_oo^{-1}*Op_oe is the same as the original
      /// residual norm. To prove that, notice that the global residual will be zero for the odds,
      /// and that the global residual for a solution x'_e of A_ee is
      ///   r = L*A*R * R^{-1}*[x'_e x'_o] - b = [r_e; 0];
      /// therefore ||L^{-1}r|| = ||r|| because of the form of L, making
      ///   ||r|| = || A_e*x'_e - (b_e-Op_eo*Op_oo^{-1}*b_o) ||.
      ///
      /// Notice that A_ee^{-1} is Op^{-1}_ee. So a way for preconditioning
      /// the solution of A_ee is with a preconditioner on Op restricted to the even sites. The
      /// solver does that when options for `prec_ee' are given.
      ///
      /// Also, A_ee can be additionally preconditioned by the left with Op_ee. Note that left
      /// preconditioning will not change the original residual norm. This option is activated
      /// when use_Aee_prec is true.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getEvenOddPrec(Operator<NOp, COMPLEX> op, const Options& ops,
					    Operator<NOp, COMPLEX> prec_, SolverSpace solverSpace)
      {
	auto dims = op.d.kvdim();
	if (dims.count('X') == 0 || dims.at('X') != 2 || op.imgLayout != XEvenOddLayout)
	  ops.throw_error(
	    "getEvenOddPrec: only supported on explicitly colored operators with two colors");

	if (getFurthestNeighborDistance(op) != 1)
	  throw std::runtime_error(
	    "getEvenOddPrec: not supported for operators with other than distance-1 neighbors");

	if (prec_)
	  throw std::runtime_error("getEvenOddPrec: unsupported input preconditioner");

	bool use_Aee_prec = getOption<bool>(ops, "use_Aee_prec", false);
	std::string prefix = getOption<std::string>(ops, "prefix", "");

	Tracker _t(std::string("setup eo solver ") + prefix);

	// Get the block diagonal of the operator with rows cs and columns CS
	remap m_sc{{'s', 'S'}, {'c', 'C'}};
	Tensor<NOp + 2, COMPLEX> opDiagE, opDiagO, opInvDiagE, opInvDiagO;
	Operator<NOp, COMPLEX> op_eo, op_oe;
	if (getEvenOddOperatorsPartsCache<NOp, COMPLEX>().count(op.id.get()) == 0)
	{
	  auto opDiag = getBlockDiag<NOp + 2>(op, "cs", m_sc); // rows cs, cols CS
	  opDiagE = opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}).clone();
	  opDiagO = opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}).clone();
	  auto opInvDiag = inv(opDiag, "cs", "CS");
	  opInvDiagE = opInvDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}).clone();
	  opInvDiagO = opInvDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}).clone();
	  op_eo = op.kvslice_from_size({{'X', 1}}, {{'X', 1}}, {{'X', 0}}, {{'X', 1}});
	  op_oe = op.kvslice_from_size({{'X', 0}}, {{'X', 1}}, {{'X', 1}}, {{'X', 1}});
	  getEvenOddOperatorsPartsCache<NOp, COMPLEX>()[op.id.get()] = {
	    op_eo, op_oe, opDiagE, opDiagO, opInvDiagE, opInvDiagO};
	}
	else
	{
	  auto t = getEvenOddOperatorsPartsCache<NOp, COMPLEX>().at(op.id.get());
	  op_eo = std::get<0>(t);
	  op_oe = std::get<1>(t);
	  opDiagE = std::get<2>(t);
	  opDiagO = std::get<3>(t);
	  opInvDiagE = std::get<4>(t);
	  opInvDiagO = std::get<5>(t);
	}

	// Get an explicit form for A:=Op_ee-Op_eo*Op_oo^{-1}*Op_oe
	unsigned int create_operator_max_rhs =
	  getOption<unsigned int>(ops, "create_operator_max_rhs", 0);
	ColOrdering co = getOption<ColOrdering>(ops, "operator_ordering", getColOrderingMap(),
						op.preferred_col_ordering);
	ColOrdering co_blk =
	  getOption<ColOrdering>(ops, "operator_block_ordering", getColOrderingMap(), RowMajor);
	unsigned int max_dist_neighbors_opA = 2;

	Operator<NOp, COMPLEX> opA{
	  use_Aee_prec ?
		       // Do opA = I - Op_eo * Op_oo^{-1} * Op_oe * Op_ee^{-1} if use_Aee_prec
	    OperatorFun<NOp, COMPLEX>(
	      [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		foreachInChuncks(
		  x, y, create_operator_max_rhs,
		  [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		    Tracker _t(std::string("eo matvec Aee_prec ") + prefix);

		    // y = x
		    x.copyTo(y);

		    // y0 = Op_ee^{-1} * x
		    auto y0 = contract<NOp + 1>(opInvDiagE, x.rename_dims(m_sc), "CS");

		    // y1 = Op_oe * y0
		    auto y1 = op_oe(std::move(y0));

		    // y2 = Op_oo^{-1} * y1
		    auto y2 = contract<NOp + 1>(opInvDiagO, std::move(y1).rename_dims(m_sc), "CS");

		    // y += -Op_eo * y2
		    op_eo(std::move(y2)).scale(-1).addTo(y);
		  });
	      })
		       :
		       // Otherwise, do opA = Op_ee - Op_eo * Op_oo^{-1} * Op_oe
	    OperatorFun<NOp, COMPLEX>(
	      [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		foreachInChuncks(x, y, create_operator_max_rhs,
				 [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
				   Tracker _t(std::string("eo matvec ") + prefix);

				   // y = Op_ee * x
				   contract(opDiagE, x.rename_dims(m_sc), "CS", CopyTo, y);

				   // y1 = Op_oe * x
				   auto y1 = op_oe(x);

				   // y2 = Op_oo^{-1} * y1
				   auto y2 = contract<NOp + 1>(
				     opInvDiagO, std::move(y1).rename_dims(m_sc), "CS");

				   // y += -Op_eo * y2
				   op_eo(std::move(y2)).scale(-1).addTo(y);
				 });
	      }),
	  op.d.kvslice_from_size({{'X', 0}}, {{'X', 1}}),
	  op.i.kvslice_from_size({{'X', 0}}, {{'X', 1}}),
	  nullptr,
	  op.order_t,
	  EvensOnlyLayout,
	  EvensOnlyLayout,
	  getNeighbors(op.i.kvdim(), max_dist_neighbors_opA, EvensOnlyLayout),
	  op.preferred_col_ordering,
	  false /* no Kronecker blocking */};

	// Get preconditioner
	auto precOps = getOptionsMaybe(ops, "prec_ee");
	Operator<NOp, COMPLEX> prec_ee;
	if (precOps)
	{
	  auto prec_ee0 = getSolver(op, precOps.getSome(), {}, OnlyEvensSpace)
			    .kvslice_from_size({}, {{'X', 1}}, {}, {{'X', 1}});
	  if (use_Aee_prec)
	  {
	    prec_ee = Operator<NOp, COMPLEX>{
	      [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		Tracker _t(std::string("eo prec_ee ") + prefix);

		// y = Op_ee * prec_ee0(x)
		contract(opDiagE, prec_ee0(x).rename_dims(m_sc), "CS", CopyTo, y);
	      },
	      prec_ee0.d, prec_ee0.i, nullptr, prec_ee0};
	  }
	  else
	  {
	    prec_ee = prec_ee0;
	  }
	}

	// Get solver on opA
	const Operator<NOp, COMPLEX> solver = getSolver(opA, getOptions(ops, "solver"), prec_ee);

	// Create the solver
	Operator<NOp, COMPLEX> rop{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("eo solver ") + prefix);

	    // be = x_e - Op_eo*Op_oo^{-1}*x_o
	    Tensor<NOp + 1, COMPLEX> be;
	    if (solverSpace == FullSpace)
	    {
	      be = solver.template make_compatible_img<NOp + 1>("n", {{'n', x.kvdim().at('n')}});
	      x.kvslice_from_size({{'X', 0}}, {{'X', 1}}).copyTo(be);
	      op_eo(contract<NOp + 1>(opInvDiagO,
				      x.kvslice_from_size({{'X', 1}}, {{'X', 1}}).rename_dims(m_sc),
				      "CS"))
		.scale(-1)
		.addTo(be);
	    }
	    else
	    {
	      be = x.kvslice_from_size({{'X', 0}}, {{'X', 1}});
	    }

	    // Solve opA * y_e = be
	    if (solverSpace != FullSpace)
	      y.set_zero();
	    auto ye = y.kvslice_from_size({{'X', 0}}, {{'X', 1}});
	    solver(be, ye);

	    // Do y_e = Op_ee^{-1} y_e if use_Aee_prec
	    if (use_Aee_prec)
	    {
	      contract<NOp + 1>(opInvDiagE, ye.rename_dims(m_sc), "CS").copyTo(ye);
	    }

	    // y_o = Op_oo^{-1}*(-Op_oe*y_e + x_o)
	    if (solverSpace == FullSpace)
	    {
	      auto yo0 = be;
	      x.kvslice_from_size({{'X', 1}}, {{'X', 1}}).copyTo(yo0);
	      op_oe(ye).scale(-1).addTo(yo0);
	      contract<NOp + 1>(opInvDiagO,
				yo0.rename_dims(m_sc), "CS", CopyTo,
				y.kvslice_from_size({{'X', 1}}, {{'X', 1}}));
	    }
	  },
	  op.i,
	  op.d,
	  nullptr,
	  op.order_t,
	  op.imgLayout,
	  op.domLayout,
	  DenseOperator(),
	  op.preferred_col_ordering,
	  false /* no Kronecker blocking */};

	// Do a test
	if (superbblas::getDebugLevel() > 0)
	{
	  auto x = op.template make_compatible_img<NOp + 1>("n", {{'n', 2}});
	  urand(x, -1, 1);
	  auto y = op(rop(x));
	  x.scale(-1).addTo(y);
	  auto normx = norm<1>(x, "n");
	  auto normdiff = norm<1>(y, "n");
	  double max_err = 0;
	  for (int i = 0, vol = normdiff.volume(); i < vol; ++i)
	    max_err = std::max(max_err, (double)normdiff.get({{i}}) / normx.get({{i}}));
	  QDPIO::cout << " eo prec error: " << detail::tostr(max_err) << std::endl;
	}

	return rop;
      }

      /// Returns Schur complement preconditioner splitting over the spin components
      ///
      /// It returns an approximation of Op^{-1} by splitting the rows and columns into the two
      /// (even-odd spin components) and doing:
      ///
      ///  Op^{-1} = R^{-1} * A^{-1} * L^{-1} ->
      ///
      /// [ Op_ee Op_eo ]^{-1} = [       I            0 ] * [ Op_ee-Op_eo*Op_oo^{-1}*Op_oe   0   ]^{-1} *
      /// [ Op_oe Op_oo ]        [ -Op_oo^{-1}*Op_oe  I ]   [              0               Op_oo ]
      ///
      ///                        * [ I  -Op_eo*Op_oo^{-1} ]
      ///                          [ 0         I          ]
      ///
      /// The matrix Op_oo^{-1} is block diagonal and is computed directly, while
      /// (Op_ee-Op_eo*Op_oo^{-1}*Op_oe)^{-1} is approximated by an iterative solver. Note that
      /// the residual norm while solving A_ee=Op_ee-Op_eo*Op_oo^{-1}*Op_oe is the same as the original
      /// residual norm. To prove that, notice that the global residual will be zero for the odds,
      /// and that the global residual for a solution x'_e of A_ee is
      ///   r = L*A*R * R^{-1}*[x'_e x'_o] - b = [r_e; 0];
      /// therefore ||L^{-1}r|| = ||r|| because of the form of L, making
      ///   ||r|| = || A_e*x'_e - (b_e-Op_eo*Op_oo^{-1}*b_o) ||.
      ///
      /// Notice that A_ee^{-1} is Op^{-1}_ee. So a way for preconditioning
      /// the solution of A_ee is with a preconditioner on Op restricted to the even sites. The
      /// solver does that when options for `prec_ee' are given.
      ///
      /// Also, A_ee can be additionally preconditioned by the left with Op_ee. Note that left
      /// preconditioning will not change the original residual norm. This option is activated
      /// when use_Aee_prec is true.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getSpinEvenOddPrec(Operator<NOp, COMPLEX> op, const Options& ops,
						Operator<NOp, COMPLEX> prec_)
      {
	auto dims = op.d.kvdim();
	if (dims.count('s') == 0 || dims.at('s') % 2 != 0)
	  ops.throw_error(
	    "getSpinEvenOddPrec: only supported on operators with even spin components");

	if (prec_)
	  throw std::runtime_error("getSpinEvenOddPrec: unsupported input preconditioner");

	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup eo spin solver ") + prefix);

	// Partition the operator: s -> Ss, where S is the spin oddity
	char S = detail::get_free_label(op.d.order);
	std::string Ss = std::string({S, 's'});
	int ns = dims.at('s');
	auto eg_d = op.d.like_this(none, {{'s', ns / 2}}).make_eg();
	auto eg_i = op.i.like_this(none, {{'s', ns / 2}}).make_eg();
	Operator<NOp, COMPLEX> op_ee{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    auto x0 = op.template make_compatible_img<NOp + 1>("n", {{'n', x.kvdim().at('n')}});
	    x0.set_zero();
	    x.copyTo(x0.template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
		       .kvslice_from_size({{S, 0}}, {{S, 1}}));
	    op(x0)
	      .template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
	      .kvslice_from_size({{S, 0}}, {{S, 1}})
	      .copyTo(y);
	  },
	  eg_d, eg_i, nullptr, op};
	Operator<NOp, COMPLEX> op_eo{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    auto x0 = op.template make_compatible_img<NOp + 1>("n", {{'n', x.kvdim().at('n')}});
	    x0.set_zero();
	    x.copyTo(x0.template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
		       .kvslice_from_size({{S, 1}}, {{S, 1}}));
	    op(x0)
	      .template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
	      .kvslice_from_size({{S, 0}}, {{S, 1}})
	      .copyTo(y);
	  },
	  eg_d, eg_i, nullptr, op};

	Operator<NOp, COMPLEX> op_oe{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    auto x0 = op.template make_compatible_img<NOp + 1>("n", {{'n', x.kvdim().at('n')}});
	    x0.set_zero();
	    x.copyTo(x0.template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
		       .kvslice_from_size({{S, 0}}, {{S, 1}}));
	    op(x0)
	      .template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
	      .kvslice_from_size({{S, 1}}, {{S, 1}})
	      .copyTo(y);
	  },
	  eg_d, eg_i, nullptr, op};

	Operator<NOp, COMPLEX> op_oo{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    auto x0 = op.template make_compatible_img<NOp + 1>("n", {{'n', x.kvdim().at('n')}});
	    x0.set_zero();
	    x.copyTo(x0.template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
		       .kvslice_from_size({{S, 1}}, {{S, 1}}));
	    op(x0)
	      .template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}})
	      .kvslice_from_size({{S, 1}}, {{S, 1}})
	      .copyTo(y);
	  },
	  eg_d, eg_i, nullptr, op};


	// Get solver on op_oo
	const auto solver_oo = getSolver(op_oo, getOptions(ops, "solver_oo"));

	// Get an explicit form for A = Op_ee-Op_eo*Op_oo^{-1}*Op_oe
	unsigned int max_dist_neighbors_opA = 2;
	Operator<NOp, COMPLEX> opA{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("spin eo matvec ") + prefix);

	    // y = Op_ee * x
	    op_ee(x, y);

	    // y1 = Op_oe * x
	    auto y1 = op_oe(x);

	    // y2 = Op_oo^{-1} * y1
	    auto y2 = solver_oo(y1);

	    // y += -Op_eo * y2
	    op_eo(std::move(y2)).scale(-1).addTo(y);
	  },
	  eg_d,
	  eg_i,
	  nullptr,
	  op.order_t,
	  op.domLayout,
	  op.imgLayout,
	  getNeighbors(op.i.kvdim(), max_dist_neighbors_opA, op.imgLayout),
	  op.preferred_col_ordering,
	  false /* no Kronecker blocking */};

	// Get solver on opA
	const auto solver = getSolver(opA, getOptions(ops, "solver_A"));

	// Create the solver
	Operator<NOp, COMPLEX> rop{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("spin eo solver ") + prefix);

	    auto x_eo = x.template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}});
	    auto y_eo = y.template split_dimension<NOp + 2>('s', Ss, {{S, 2}, {'s', ns / 2}});

	    // be = x_e - Op_eo*Op_oo^{-1}*x_o
	    auto x_e = x_eo.kvslice_from_size({{S, 0}}, {{S, 1}});
	    auto x_o = x_eo.kvslice_from_size({{S, 1}}, {{S, 1}});
	    auto b_e = op_eo(solver_oo(x_o).scale(-1));
	    x_e.addTo(b_e);

	    // Solve opA * y_e = be
	    auto y_e = y_eo.kvslice_from_size({{S, 0}}, {{S, 1}});
	    auto y_o = y_eo.kvslice_from_size({{S, 1}}, {{S, 1}});
	    solver(b_e, y_e);

	    // y_o = Op_oo^{-1}*(-Op_oe*y_e + x_o)
	    auto yo0 = b_e;
	    x_o.copyTo(yo0);
	    op_oe(y_e).scale(-1).addTo(yo0);
	    solver_oo(yo0, y_o);
	  },
	  op.i,
	  op.d,
	  nullptr,
	  op.order_t,
	  op.imgLayout,
	  op.domLayout,
	  DenseOperator(),
	  op.preferred_col_ordering,
	  false /* no Kronecker blocking */};

	// Do a test
	if (superbblas::getDebugLevel() > 0)
	{
	  auto x = op.template make_compatible_img<NOp + 1>("n", {{'n', 2}});
	  urand(x, -1, 1);
	  auto y = op(rop(x));
	  x.scale(-1).addTo(y);
	  auto normx = norm<1>(x, "n");
	  auto normdiff = norm<1>(y, "n");
	  double max_err = 0;
	  for (int i = 0, vol = normdiff.volume(); i < vol; ++i)
	    max_err = std::max(max_err, (double)normdiff.get({{i}}) / normx.get({{i}}));
	  QDPIO::cout << " spin eo prec error: " << detail::tostr(max_err) << std::endl;
	}

	return rop;
      }

      enum LandSeaDomain { Land, Sea, All };

      /// Returns a generic two domains preconditioner.
      ///
      /// It returns an approximation of Op^{-1} by splitting the rows and columns into the two
      /// domains (islands, domains without holes, and sea, the domain which connect all the domains) and doing:
      ///
      /// [ Op_ss Op_si ]^{-1} = [       I            0 ] * [ Op_ss-Op_si*Op_ii^{-1}*Op_is   0   ]^{-1} *
      /// [ Op_is Op_ii ]        [ -Op_ii^{-1}*Op_is  I ]   [              0               Op_ii ]
      ///
      ///                        * [ I  -Op_si*Op_ii^{-1} ]
      ///                          [ 0         I          ]
      ///
      /// The matrix Op_ii^{-1} is block diagonal, while
      /// (Op_ss-Op_si*Op_ii^{-1}*Op_is)^{-1} should be close to Op_ss^{-1}. Note that
      /// the residual norm while solving A_ss=Op_ss-Op_si*Op_ii^{-1}*Op_is is the same as the original
      /// residual norm if the linear systems with Op_ii are solved exactly. See a prove in the comment of `getEvenOddPrec`.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getHierarchicalPrec(Operator<NOp, COMPLEX> op, const Options& ops,
						 Operator<NOp, COMPLEX> prec_,
						 LandSeaDomain land_sea_domain = All,
						 std::array<unsigned int, 4> given_divisions = {},
						 std::array<unsigned int, 4> given_land = {})
      {
	auto dims = op.d.kvdim();

	if (prec_)
	  throw std::runtime_error("getHierarchicalPrec: unsupported input preconditioner");
	if (!op.sp)
	  throw std::runtime_error("getHierarchicalPrec: unsupported implicit operators");

	std::string prefix = getOption<std::string>(ops, "prefix", "");
	auto divisions = land_sea_domain == All
			   ? getOption<std::array<unsigned int, 4>>(ops, "divisions")
			   : given_divisions;
	auto land =
	  land_sea_domain == All ? getOption<std::array<unsigned int, 4>>(ops, "land") : given_land;
	bool use_red_black = getOption<bool>(ops, "use_red_black", false);

	// Check divisions and land
	if (divisions == std::array<unsigned int, 4>{{}})
	  ops.getValue("divisions").throw_error("the four values cannot be zero");
	std::array<int, 4> latt_size{dims.at('x'), dims.at('y'), dims.at('z'), dims.at('t')};
	for (std::size_t i = 0; i < 4; ++i)
	  if (divisions[i] > latt_size[i])
	    divisions[i] = latt_size[i];
	for (std::size_t i = 0; i < 4; ++i)
	{
	  if (divisions[i] == 0)
	    continue;
	  if (land[i] == 0 || land[i] >= latt_size[i] / divisions[i])
	    ops.getValue("land").throw_error(
	      "values cannot be zero nor the total size of the divisision");
	}

	// Function
	// The i-th division starts at T / divisions * i + std::min(i, T % divisions).
	const auto part_ith_starts = [&](int i, int d) {
	  return latt_size[d] / divisions[d] * i + std::min(i, latt_size[d] % (int)divisions[d]);
	};
	// So all coordinates before division T % divisions have T / divisions + 1 elements,
	// and the others have T / divisions elements.
	std::array<int, 4> first_coor_within_small_partition, partition_size_large_partition,
	  partition_size_small_partition;
	for (int d = 0; d < 4; ++d)
	{
	  if (divisions[d] == 0)
	    continue;
	  first_coor_within_small_partition[d] = part_ith_starts(latt_size[d] % divisions[d], d);
	  partition_size_large_partition[d] = latt_size[d] / divisions[d] + 1;
	  partition_size_small_partition[d] = latt_size[d] / divisions[d];
	}

	const auto is_true_island = [&](const std::array<int, 4>& coor) {
	  // The coordinate is island if for all components is island
	  int num_sea = 0;
	  for (int d = 0; d < 4; ++d)
	  {
	    if (divisions[d] == 0)
	      continue;
	    int c = coor[d], disp, land_size;
	    if (c < first_coor_within_small_partition[d])
	    {
	      disp = c - part_ith_starts(c / partition_size_large_partition[d], d);
	      land_size = land[d];
	    }
	    else
	    {
	      disp = c - part_ith_starts(latt_size[d] % divisions[d] +
					   (c - first_coor_within_small_partition[d]) /
					     partition_size_small_partition[d],
					 d);
	      land_size = land[d];
	    }
	    if (disp >= land_size)
	    {
	      if (!use_red_black)
		return false;
	      else
		++num_sea;
	    }
	  }
	  if (!use_red_black)
	    return true;
	  else
	    return num_sea % 2 == 1;
	};
	// Border of the island
	const auto is_coast = [&](const std::array<int, 4>& coor) {
	  // The coordinate is land border if for all components is island and for some component is border
	  bool r = false;
	  for (int d = 0; d < 4; ++d)
	  {
	    if (divisions[d] == 0)
	      continue;
	    int c = coor[d], disp, land_size;
	    if (c < first_coor_within_small_partition[d])
	    {
	      disp = c - part_ith_starts(c / partition_size_large_partition[d], d);
	      land_size = land[d];
	    }
	    else
	    {
	      disp = c - part_ith_starts(latt_size[d] % divisions[d] +
					   (c - first_coor_within_small_partition[d]) /
					     partition_size_small_partition[d],
					 d);
	      land_size = land[d];
	    }
	    if (disp >= land_size)
	      return false;
	    if (disp == 0 || disp == land_size - 1)
	      r = true;
	  }
	  return r;
	};

	// Border of the sea
	const auto is_littoral = [&](const std::array<int, 4>& coor) {
	  // The coordinate is sea border if for all components is sea and for some component is border
	  bool r = false;
	  for (int d = 0; d < 4; ++d)
	  {
	    if (divisions[d] == 0)
	      continue;
	    int c = coor[d], disp, land_size, size;
	    if (c < first_coor_within_small_partition[d])
	    {
	      disp = c - part_ith_starts(c / partition_size_large_partition[d], d);
	      land_size = land[d];
	      size = partition_size_large_partition[d];
	    }
	    else
	    {
	      disp = c - part_ith_starts(latt_size[d] % divisions[d] +
					   (c - first_coor_within_small_partition[d]) /
					     partition_size_small_partition[d],
					 d);
	      land_size = land[d];
	      size = partition_size_small_partition[d];
	    }
	    if (disp < land_size)
	      return false;
	    if (disp == land_size || disp == size - 1)
	      r = true;
	  }
	  return r;
	};

	std::array<char, 4> latt_labels{'x', 'y', 'z', 't'};
	std::array<int, 4> idx_dom, idx_img;
	for (std::size_t i = 0; i < 4; ++i)
	{
	  idx_dom[i] =
	    std::find(op.d.order.begin(), op.d.order.end(), latt_labels[i]) - op.d.order.begin();
	  idx_img[i] =
	    std::find(op.i.order.begin(), op.i.order.end(), latt_labels[i]) - op.i.order.begin();
	}
	if (idx_dom != idx_img)
	  throw std::runtime_error("getHierarchicalPrec: unsupported operator, it should have the "
				   "same domain and image orderings");

	const auto is_island = [&](const Coor<NOp>& op_coor) {
	  std::array<int, 4> latt_coor;
	  for (std::size_t i = 0; i < 4; ++i)
	    latt_coor[i] = op_coor[idx_dom[i]];
	  switch (land_sea_domain)
	  {
	  case All: return is_true_island(latt_coor);
	  case Land: return is_coast(latt_coor);
	  case Sea: return is_littoral(latt_coor);
	  }
	  return false; // avoid warning
	};

	// return whether the coordinate is image part: island, domain part: island
	const auto is_island_island = [&](const Coor<NOp>& dom_coor, const Coor<NOp>& img_coor) {
	  return is_island(dom_coor) & is_island(img_coor);
	};
	// return whether the coordinate is image part: island, domain part: sea
	const auto is_island_sea = [&](const Coor<NOp>& dom_coor, const Coor<NOp>& img_coor) {
	  return !is_island(dom_coor) & is_island(img_coor);
	};
	// return whether the coordinate is image part: sea, domain part: island
	const auto is_sea_island = [&](const Coor<NOp>& dom_coor, const Coor<NOp>& img_coor) {
	  return is_island(dom_coor) & !is_island(img_coor);
	};
	// return whether the coordinate is image part: sea, domain part: sea
	const auto is_sea_sea = [&](const Coor<NOp>& dom_coor, const Coor<NOp>& img_coor) {
	  return !is_island(dom_coor) & !is_island(img_coor);
	};

	// Get the block diagonal of the operator with rows cs and columns CS
	Operator<NOp, COMPLEX> op_ss = op.kvslice_from_size(is_sea_sea);
	Operator<NOp, COMPLEX> op_si = op.kvslice_from_size(is_sea_island);
	Operator<NOp, COMPLEX> op_is = op.kvslice_from_size(is_island_sea);
	Operator<NOp, COMPLEX> op_ii = op.kvslice_from_size(is_island_island);
	Operator<NOp, COMPLEX> id_ii = op.get_identiy().kvslice_from_size(is_island_island);

	// Get solvers for the island and the sea operators
	const auto& solver_ssii_ops = getOptions(ops, "sea_land_solver");
	const Operator<NOp, COMPLEX> solver_ss = getSolver(op_ss, solver_ssii_ops);
	const Operator<NOp, COMPLEX> solver_ii = getSolver(op_ii, solver_ssii_ops);

	// Get solvers for inside opA
	const auto solver_opA_ops = getOptionsMaybe(ops, "schur_solver");
	const auto solver_ssii_schur_ops =
	  solver_opA_ops ? getOptionsMaybe(ops, "sea_land_schur_solver") : none;
	bool is_solver_ssii_schur_hie =
	  solver_ssii_schur_ops && land_sea_domain == All
	    ? getOption<std::string>(solver_ssii_schur_ops.getSome(), "type") == std::string("hie")
	    : false;
	const Operator<NOp, COMPLEX> solver_ss_schur =
	  solver_ssii_schur_ops
	    ? (is_solver_ssii_schur_hie
		 ? getHierarchicalPrec(op_ss, solver_ssii_schur_ops.getSome(),
				       Operator<NOp, COMPLEX>{}, Sea, divisions, land)
		 : getSolver(op_ss, solver_ssii_schur_ops.getSome()))
	    : solver_ss;
	const Operator<NOp, COMPLEX> solver_ii_schur =
	  solver_ssii_schur_ops
	    ? (is_solver_ssii_schur_hie
		 ? getHierarchicalPrec(op_ii, solver_ssii_schur_ops.getSome(),
				       Operator<NOp, COMPLEX>{}, Land, divisions, land)
		 : getSolver(op_ii, solver_ssii_schur_ops.getSome()))
	    : solver_ii;

	unsigned int create_operator_max_rhs =
	  getOption<unsigned int>(ops, "create_operator_max_rhs", 0);
	ColOrdering co = getOption<ColOrdering>(ops, "operator_ordering", getColOrderingMap(),
						op.preferred_col_ordering);
	ColOrdering co_blk =
	  getOption<ColOrdering>(ops, "operator_block_ordering", getColOrderingMap(), RowMajor);
	unsigned int max_dist_neighbors_opA = 2;

	// Do opA = I - Op_si * Op_ii^{-1} * Op_is * Op_ss^{-1}
	Operator<NOp, COMPLEX> opA{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    foreachInChuncks(x, y, create_operator_max_rhs,
			     [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
			       Tracker _t(std::string("hierarchical matvec ") + prefix);

			       // y = x
			       x.copyTo(y);

			       // y -= Op_si * Op_ii^{-1} * Op_is * Op_ss^{-1} * x
			       op_si(solver_ii_schur(op_is(solver_ss_schur(x)))).scale(-1).addTo(y);
			     });
	  },
	  op_ss.d, op_ss.i, nullptr, op_ss};
	const Operator<NOp, COMPLEX> opA_solver =
	  solver_opA_ops ? getSolver(opA, solver_opA_ops.getSome()) : Operator<NOp, COMPLEX>{};

	// Create the solver
	Operator<NOp, COMPLEX> rop{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Tracker _t(std::string("hierarchical solver ") + prefix);

	    // b_s = x_s - Op_si*Op_ii^{-1}*x_i = x - x_i - Op_si*Op_ii^{-1}*x_i
	    auto x_i = id_ii(x);
	    Tensor<NOp + 1, COMPLEX> b_s = op_si(solver_ii(x_i.scale(-1)));
	    x.addTo(b_s);
	    x_i.scale(-1).addTo(b_s);

	    // y_s = opA^{-1} * b_s
	    const auto y_s = opA_solver ? opA_solver(b_s) : b_s;

	    // y = Op_ss^{-1} * y_s
	    solver_ss(y_s, y);

	    // y += Op_ii^{-1}*(-Op_is*y + x_i)
	    op_is(y).scale(-1).addTo(x_i);
	    solver_ii(x_i, y_s);
	    y_s.addTo(y);
	  },
	  op.i,
	  op.d,
	  nullptr,
	  op.order_t,
	  op.imgLayout,
	  op.domLayout,
	  DenseOperator(),
	  op.preferred_col_ordering,
	  false /* no Kronecker blocking */};

	// Do a test
	if (superbblas::getDebugLevel() > 0)
	{
	  auto x = op.template make_compatible_img<NOp + 1>("n", {{'n', 2}});
	  urand(x, -1, 1);
	  auto y = op(rop(x));
	  x.scale(-1).addTo(y);
	  auto normx = norm<1>(x, "n");
	  auto normdiff = norm<1>(y, "n");
	  double max_err = 0;
	  for (int i = 0, vol = normdiff.volume(); i < vol; ++i)
	    max_err = std::max(max_err, (double)normdiff.get({{i}}) / normx.get({{i}}));
	  QDPIO::cout << " hierarchical prec error: " << detail::tostr(max_err) << std::endl;
	}

	return rop;
      }

      /// Returns an approximation of the inverse of the domains local to the processes
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation
      /// \param prec_: preconditioner
      ///
      /// By default (`with_correction` being false), the solver returns D^{-1}*x,
      /// where D is the local part of the operator. If asking for extra correction,
      /// `with_correciton` being true, then  the solver returns instead:
      ///
      ///    A^{-1} = (A - D + D)^{-1} = (I - (I - D^{-1}*A))^{-1}*D^{-1}
      ///           \approx (2*I - D^{-1}*A)*D^{-1}

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getDomainDecompositionPrec(Operator<NOp, COMPLEX> op,
							const Options& ops,
							Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getDomainDecompositionPrec: unsupported input preconditioner");

	// Get options
	bool with_correction = getOption<bool>(ops, "with_correction", false);
	std::string prefix = getOption<std::string>(ops, "prefix", "");

	Tracker _t(std::string("setup dd solver ") + prefix);

	// Get the local operator and make the solver
	auto local_op = op.getGlocal();
	const Operator<NOp, COMPLEX> solver = getSolver(local_op, getOptions(ops, "solver"));

	// Return the solver
	if (!with_correction)
	  return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		    // y = D^{-1}*x
		    solver(x.getGlocal(), y.getGlocal());
		  },
		  op.i,
		  op.d,
		  nullptr,
		  op.order_t,
		  op.imgLayout,
		  op.domLayout,
		  DenseOperator(),
		  op.preferred_col_ordering,
		  false /* no Kronecker blocking */};
	else
	  return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		    // y0 = D^{-1}*x
		    auto y0 = y.make_compatible();
		    solver(x.getGlocal(), y0.getGlocal());
		    // y1 = A*y0
		    auto y1 = op(y0);
		    // y = -D^{-1}*y1
		    solver(y1.scale(-1).getGlocal(), y.getGlocal());
		    // y += 2*y0
		    y0.scale(2).addTo(y);
		  },
		  op.i,
		  op.d,
		  nullptr,
		  op.order_t,
		  op.imgLayout,
		  op.domLayout,
		  DenseOperator(),
		  op.preferred_col_ordering,
		  false /* no Kronecker blocking */};
      }

      /// Returns the inverse of the block diagonal
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getBlockJacobi(Operator<NOp, COMPLEX> op, const Options& ops,
					    Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getBlockJacobi: unsupported input preconditioner");

	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup bj solver ") + prefix);

	// Get the blocking
	auto dim = op.d.kvdim();
	std::vector<unsigned int> default_blocking{
	  {op.imgLayout == EvensOnlyLayout ? 2u : 1u, 1u, 1u, 1u}};
	std::vector<unsigned int> blocking =
	  getOption<std::vector<unsigned int>>(ops, "blocking", default_blocking);
	if (blocking.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getBlocking: the blocking should be a vector with four elements");
	std::map<char, unsigned int> mblk{
	  {'x', blocking[0]}, {'y', blocking[1]}, {'z', blocking[2]}, {'t', blocking[3]}};

	// Shortcut for default blocking
	if (blocking == default_blocking)
	{
	  // Get the block diagonal of the operator with rows cs and columns CS
	  const std::string blk_rows = "cs"; // order of the block of rows to invert
	  remap m_blk = getNewLabels(blk_rows, op.d.order + op.i.order); // column labels
	  const std::string blk_cols =
	    update_order(blk_rows, m_blk); // order of the block of columns to invert
	  Tensor<NOp + 2, COMPLEX> opDiag = getBlockDiag<NOp + 2>(op, blk_rows, m_blk);

	  // Return the solver
	  return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		    // y = Op_diag^{-1} * x
		    solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		      opDiag, blk_rows, blk_cols, x.rename_dims(m_blk), blk_cols, CopyTo, y);
		  },
		  op.d, op.i, nullptr, op};
	}

	// Check that the blocking divide the lattice sizes. For an operator with support only
	// on the even sites, make sure that the blocking on x is divisible by two
	int X = op.imgLayout == EvensOnlyLayout ? 2 : dim.at('X');
	bool x_blocking_divide_X = (mblk.at('x') % X == 0);
	if (!x_blocking_divide_X && op.imgLayout == EvensOnlyLayout)
	  ops.getValue("blocking")
	    .throw_error(
	      "When using even-odd preconditioning, the blocking on x should be divisible by 2");
	if (!x_blocking_divide_X && blocking[0] != default_blocking[0])
	  ops.getValue("blocking")
	    .throw_error(
	      "unsupported blocking which is neither one nor divisible by 2 on the x direction");
	for (const auto it : getNatLatticeDims(dim, op.imgLayout))
	{
	  if (it.second % mblk.at(it.first) != 0)
	    ops.getValue("blocking")
	      .throw_error("The operator dimensions are not divisible by the blocking");
	}

	// Shortcut when no blocking on x direction
	if (default_blocking[0] == blocking[0])
	{
	  // Get the block diagonal of the operator with rows cs and columns CS
	  const std::string blk_rows = "cs0123"; // order of the block of rows to invert
	  remap m_blk = getNewLabels(blk_rows, op.d.order + op.i.order); // column labels
	  const std::string blk_cols =
	    update_order(blk_rows, m_blk); // order of the block of columns to invert

	  std::map<char, int> blk{
	    {'1', (int)blocking[1]}, {'2', (int)blocking[2]}, {'3', (int)blocking[3]}};
	  std::map<char, unsigned int> blk_u{
	    {'x', 1u}, {'y', blocking[1]}, {'z', blocking[2]}, {'t', blocking[3]}};
	  auto new_d = op.d
			 .template reshape_dimensions<NOp>(
			   {{"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true)
			 .make_eg();
	  auto new_i = op.i
			 .template reshape_dimensions<NOp>(
			   {{"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true)
			 .make_eg();
	  auto blk_op = Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = x.template reshape_dimensions<NOp + 1>(
		{{"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, dim, true);
	      auto y0 = op(std::move(x0));
	      y0.template reshape_dimensions<NOp + 1>({{"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}},
						      blk, true)
		.copyTo(y);
	    },
	    new_d,
	    new_i,
	    nullptr,
	    op.order_t,
	    op.domLayout == EvensOnlyLayout ? NaturalLayout : op.domLayout,
	    op.imgLayout == EvensOnlyLayout ? NaturalLayout : op.imgLayout,
	    getNeighborsAfterBlocking(blk_u, op.d.kvdim(), op.neighbors, op.imgLayout),
	    op.preferred_col_ordering,
	    op.is_kronecker()};

	  Tensor<NOp + 6, COMPLEX> opDiag =
	    getBlockDiag<NOp + 6>(blk_op, blk_rows, m_blk, ConsiderBlockingDense);

	  // Return the solver
	  return Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = x.template reshape_dimensions<NOp + 1>(
		{{"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true);

	      // y = Op_diag^{-1} * x
	      auto y0 = solve<6, NOp + 1, NOp + 6, NOp + 1, COMPLEX>(
		opDiag, blk_rows, blk_cols, std::move(x0).rename_dims(m_blk), blk_cols);

	      y0.template reshape_dimensions<NOp + 1>({{"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}},
						      dim, true)
		.copyTo(y);
	    },
	    op.d, op.i, nullptr, op};
	}
	else
	{
	  // Get the block diagonal of the operator with rows cs and columns CS
	  const std::string blk_rows = "csX0123"; // order of the block of rows to invert
	  remap m_blk = getNewLabels(blk_rows, op.d.order + op.i.order); // column labels
	  const std::string blk_cols =
	    update_order(blk_rows, m_blk); // order of the block of columns to invert

	  std::map<char, int> blk{{'0', (int)(blocking[1] / X)},
				  {'X', 1},
				  {'1', (int)blocking[1]},
				  {'2', (int)blocking[2]},
				  {'3', (int)blocking[3]}};
	  std::map<char, unsigned int> blk_u{
	    {'x', blocking[1] / X}, {'y', blocking[1]}, {'z', blocking[2]}, {'t', blocking[3]}};

	  auto new_d = op.d
			 .template reshape_dimensions<NOp>(
			   {{"X0x", "X0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true)
			 .make_eg();
	  auto new_i = op.i
			 .template reshape_dimensions<NOp>(
			   {{"X0x", "X0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true)
			 .make_eg();
	  auto blk_op = Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = x.template reshape_dimensions<NOp + 1>(
		{{"X0x", "X0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, dim, true);
	      auto y0 = op(std::move(x0));
	      y0.template reshape_dimensions<NOp + 1>(
		  {{"X0x", "X0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true)
		.copyTo(y);
	    },
	    new_d,
	    new_i,
	    nullptr,
	    op.order_t,
	    NaturalLayout,
	    NaturalLayout,
	    getNeighborsAfterBlocking(blk_u, op.d.kvdim(), op.neighbors, op.imgLayout),
	    op.preferred_col_ordering,
	    op.is_kronecker()};

	  Tensor<NOp + 7, COMPLEX> opDiag =
	    getBlockDiag<NOp + 7>(blk_op, blk_rows, m_blk, ConsiderBlockingDense);

	  // Return the solver
	  return Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = x.template reshape_dimensions<NOp + 1>(
		{{"X0x", "X0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blk, true);

	      // y = Op_diag^{-1} * x
	      auto y0 = solve<7, NOp + 1, NOp + 7, NOp + 1, COMPLEX>(
		opDiag, blk_rows, blk_cols, std::move(x0).rename_dims(m_blk), blk_cols);

	      y0.template reshape_dimensions<NOp + 1>(
		  {{"X0x", "X0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, dim, true)
		.copyTo(y);
	    },
	    op.d, op.i, nullptr, op};
	}
      }

      /// Returns a blocking, which should enhanced the performance of the sparse-dense tensor contraction.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getBlocking(Operator<NOp, COMPLEX> op, const Options& ops,
					 Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getBlocking: unsupported input preconditioner");

	auto dim = op.d.kvdim();

	std::vector<unsigned int> default_blocking{
	  {dim.count('0') == 1 ? (unsigned int)dim.at('0') : 1u,
	   dim.count('1') == 1 ? (unsigned int)dim.at('1') : 1u,
	   dim.count('2') == 1 ? (unsigned int)dim.at('2') : 1u,
	   dim.count('3') == 1 ? (unsigned int)dim.at('3') : 1u}};
	std::vector<unsigned int> blocking =
	  getOption<std::vector<unsigned int>>(ops, "blocking", default_blocking);
	if (blocking.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getBlocking: the blocking should be a vector with four elements");
	std::map<char, unsigned int> mblk{
	  {'x', blocking[0]}, {'y', blocking[1]}, {'z', blocking[2]}, {'t', blocking[3]}};

	// Check that the blocking divide the lattice sizes. For an operator with support only
	// on the even sites, make sure that the blocking on x is divisible by two
	auto opdims = op.d.kvdim();
	int X = op.imgLayout == EvensOnlyLayout ? 2 : opdims.at('X');
	bool x_blocking_divide_X = (mblk.at('x') % X == 0);
	if (!x_blocking_divide_X && op.imgLayout == EvensOnlyLayout)
	  ops.getValue("blocking")
	    .throw_error(
	      "When using even-odd preconditioning, the blocking on x should be divisible by 2");
	for (const auto it : getNatLatticeDims(opdims, op.imgLayout))
	{
	  if (it.second % mblk.at(it.first) != 0)
	    ops.getValue("blocking")
	      .throw_error("The operator dimensions are not divisible by the blocking");
	}

	// Construct the blocked operator
	ColOrdering co = getOption<ColOrdering>(ops, "operator_ordering", getColOrderingMap(),
						op.preferred_col_ordering);
	ColOrdering co_blk =
	  getOption<ColOrdering>(ops, "operator_block_ordering", getColOrderingMap(), RowMajor);
	unsigned int power = getOption<unsigned int>(ops, "power", 1);
	bool make_explicit = getOption<bool>(ops, "make_explicit", true);

	auto blkd = op.d
		      .template reshape_dimensions<NOp>(
			{{"0x", "0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}},
			{{'0', std::max(X, (int)mblk.at('x')) / X},
			 {'1', mblk.at('y')},
			 {'2', mblk.at('z')},
			 {'3', mblk.at('t')}},
			true)
		      .reorder("%0123Xxyzt", '%');

	const auto blkdim = blkd.kvdim();
	const Operator<NOp, COMPLEX> sop = cloneOperator(
	  Operator<NOp, COMPLEX>{
	    [&](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      op(x.template reshape_dimensions<NOp + 1>(
		   {{"0x", "0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, dim, true))
		.template reshape_dimensions<NOp + 1>(
		  {{"0x", "0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blkdim, true)
		.copyTo(y);
	    },
	    blkd, blkd, nullptr, op},
	  getFurthestNeighborDistance(op) * power, co, co_blk, ConsiderBlockingSparse, "blocking");

	auto solverOps = getOptionsMaybe(ops, "solver");
	const Operator<NOp, COMPLEX> solver =
	  solverOps.hasSome() ? getSolver(sop, solverOps.getSome()) : sop;

	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  solver(x.template reshape_dimensions<NOp + 1>(
			   {{"0x", "0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, blkdim, true))
		    .template reshape_dimensions<NOp + 1>(
		      {{"0x", "0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}}, dim, true)
		    .copyTo(y);
		},
		solverOps.hasSome() ? op.i : op.d, solverOps.hasSome() ? op.d : op.i, nullptr, sop};
      }
    }

    /// Returns an inexact Generalized Davidson on op*g5, that is the left SV of op
    /// NOTE: this is an eigensolver not a linear solver
    ///
    /// \param op: operator to make the inverse of
    /// \param ops: options to select the solver and the null-vectors creation

    template <typename COMPLEX>
    EigensolverFun<Nd + 7, COMPLEX> getInexactEigensolverGD(Operator<Nd + 7, COMPLEX> op,
							    const Options& ops)
    {
      using namespace detail;

      Tracker _t(std::string("setup gd solver "));

      // Get eigensolver properties
      unsigned int max_basis_size = getOption<unsigned int>(ops, "max_basis_size", 0);
      unsigned int max_block_size = getOption<unsigned int>(ops, "max_block_size", 1);
      Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);
      auto op_double = op.template cast<ComplexD>();

      // Return the solver
      return [=](int numEvals,
		 double tol) -> std::tuple<std::vector<double>, Tensor<Nd + 8, COMPLEX>> {
	DeviceHost primme_dev = OnHost;
#  if defined(SUPERBBLAS_USE_GPU)
	primme_dev = OnDefaultDevice;
#  endif
#  if defined(SUPERBBLAS_USE_HIP)
	hipblasHandle_t gpublas_handle;
	ns_getColorvecs::gpuBlasCheck(hipblasCreate(&gpublas_handle));
#  endif

	// Create an auxiliary struct for the PRIMME's matvec
	// NOTE: Please keep 'n' as the slowest index; the rows of vectors taken by PRIMME's matvec has dimensions 'cxyztX',
	// and 'n' is the dimension for the columns.
	ns_getColorvecs::OperatorAux<Nd + 7> opaux{op_double, primme_dev};

	// Make a bigger structure holding
	primme_params primme;
	primme_initialize(&primme);

	// Primme solver setup
	primme.numEvals = numEvals;
	primme.printLevel =
	  (verb == NoOutput ? 0 : verb == JustSummary ? 1 : verb == Detailed ? 3 : 5);
	primme.n = op_double.d.volume();
	primme.eps = tol;
	primme.target = primme_largest_abs;
	double zeros = 0;
	primme.targetShifts = &zeros;
	primme.numTargetShifts = 1;

	// Set parallel settings
	primme.nLocal = op_double.d.getLocal().volume();
	primme.numProcs = QDP::Layout::numNodes();
	primme.procID = QDP::Layout::nodeNumber();
	primme.globalSumReal = ns_getColorvecs::primmeGlobalSum;

	// No preconditioner for my matrix
	primme.matrixMatvec = ns_getColorvecs::primmeMatvecFermion;
	primme.matrix = &opaux;

	// Set block size
	primme.maxBasisSize = max_basis_size;
	primme.maxBlockSize = max_block_size;
	primme.ldOPs = primme.nLocal;

	// Should set lots of defaults
	if (primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme) < 0)
	{
	  QDPIO::cerr << __func__ << ": invalid preset method\n";
	  QDP_abort(1);
	}

	// Print primme options
	if (verb >= VeryDetailed && Layout::nodeNumber() == 0)
	  primme_display_params(primme);

	// Allocate space for converged Ritz values and residual norms
	std::vector<double> evals(primme.numEvals);
	std::vector<double> rnorms(primme.numEvals);
	Tensor<Nd + 8, ComplexD> evecs = op_double.d.template make_compatible<Nd + 8, ComplexD>(
	  "%n", '%', "", {{'n', (int)numEvals}}, primme_dev);
#  if defined(SUPERBBLAS_USE_GPU)
#    if defined(SUPERBBLAS_USE_CUDA)
	superbblas::detail::GpuBlasHandle gpublas_handle =
	  superbblas::detail::getGpuBlasHandle(evecs.ctx().toGpu(0));
	// Make sure cublas handle operates on legacy stream for primme
	ns_getColorvecs::gpuBlasCheck(cublasSetStream(gpublas_handle, 0));

#    endif
	primme.queue = &gpublas_handle;
#  endif

	// Call primme
#  if defined(SUPERBBLAS_USE_GPU)
	int ret = cublas_zprimme(evals.data(), evecs.data(), rnorms.data(), &primme);
#  else
	int ret = zprimme(evals.data(), evecs.data(), rnorms.data(), &primme);
#  endif

	if (verb != NoOutput)
	{
	  QDPIO::cout << "Eigenpairs converged: " << primme.initSize << std::endl;
	  QDPIO::cout << "Tolerance : " << primme.aNorm * primme.eps << std::endl;
	  QDPIO::cout << "Iterations: " << (int)primme.stats.numOuterIterations << std::endl;
	  QDPIO::cout << "Restarts  : " << (int)primme.stats.numRestarts << std::endl;
	  QDPIO::cout << "Matvecs   : " << (int)primme.stats.numMatvecs << std::endl;
	  QDPIO::cout << "Preconds  : " << (int)primme.stats.numPreconds << std::endl;
	  QDPIO::cout << "T. ortho  : " << primme.stats.timeOrtho << std::endl;
	  QDPIO::cout << "T. matvec : " << primme.stats.timeMatvec << std::endl;
	  QDPIO::cout << "Total time: " << primme.stats.elapsedTime << std::endl;
	}

	if (ret != 0)
	{
	  QDPIO::cerr << "Error: primme returned with nonzero exit status\n";
	  QDP_abort(1);
	}

	// Check the residuals, |op*v-lambda*v|_2<=|op|*tol
	if (evals.size() > 0)
	{
	  auto g5 = getGamma5<ComplexD>(op.d.kvdim().at('s'));
	  auto r =
	    op(contract<Nd + 8, 2, Nd + 8, ComplexD>(g5.rename_dims({{'j', 's'}}), evecs, "s")
		 .rename_dims({{'i', 's'}}));
	  std::vector<ComplexD> evals_cmpl(evals.begin(), evals.end());
	  contract<Nd + 8, Nd + 8, 1, ComplexD>(
	    evecs, asTensorView(evals_cmpl).rename_dims({{'i', 'n'}}).scale(-1), "", AddTo, r);
	  auto rnorm = norm<1>(r, "n");
	  for (int i = 0, vol = rnorm.volume(); i < vol; ++i)
	  {
	    if (rnorm.get({{i}}) > primme.stats.estimateLargestSVal * primme.eps * 10)
	    {
	      QDPIO::cerr << "Warning (or error): primme returned eigenpairs with too much error; "
			     "residual norm "
			  << i << " : " << rnorm.get({{i}}) << ", expected less than "
			  << primme.stats.estimateLargestSVal * primme.eps * 10 << std::endl;
	    }
	  }
	}

	// Cleanup
	primme_free(&primme);

#  if defined(SUPERBBLAS_USE_HIP)
	ns_getColorvecs::gpuBlasCheck(hipblasDestroy(gpublas_handle));
#  endif

	// Return
	return std::make_tuple(evals, evecs.template cast<COMPLEX>());
      };
    }

    template EigensolverFun<Nd + 7, SB::ComplexD>
    getInexactEigensolverGD(Operator<Nd + 7, SB::ComplexD> op, const Options& ops);
    template EigensolverFun<Nd + 7, SB::ComplexF>
    getInexactEigensolverGD(Operator<Nd + 7, SB::ComplexF> op, const Options& ops);

    namespace detail
    {
      /// Returns an inexact Generalized Davidson.
      /// NOTE: this is an eigensolver not a linear solver
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getInexactGD(Operator<NOp, COMPLEX> op, const Options& ops,
					  Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getInexactGD: unsupported input preconditioner");

	const Operator<NOp, COMPLEX> solver = getSolver(op, getOptions(ops, "solver"));
	double tol = getOption<double>(ops, "tol");
	auto eigensolver = getInexactEigensolverGD(solver, ops);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  auto values_vectors = eigensolver(x.kvdim().at('n'), tol);
		  std::get<1>(values_vectors).copyTo(y);
		},
		solver.i,
		solver.d,
		nullptr,
		solver.order_t,
		solver.imgLayout,
		solver.domLayout,
		DenseOperator(),
		solver.preferred_col_ordering,
		false /* no Kronecker form */};
      }

      /// Returns a shifted operator
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getShiftedOp(Operator<NOp, COMPLEX> op, const Options& ops,
					  Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getShiftedOp: unsupported input preconditioner");

	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup shift ") + prefix);

	// Get the remainder options
	double shift_I = getOption<double>(ops, "shift_I", 0.0);
	double shift_ig5 = getOption<double>(ops, "shift_ig5", 0.0);

	// Produced the shifted operator
	auto shifted_op = op;
	if (shift_I != 0 || shift_ig5 != 0)
	{
	  int ns = op.d.kvdim().at('s');
	  auto g5 = getGamma5<COMPLEX>(ns, op.d.getDev(), op.d.dist);
	  shifted_op = {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
			  // y = op * x
			  op(x, y);

			  // y +=  x * shift_I
			  if (shift_I != 0)
			    x.scale(shift_I).addTo(y);

			  // y += shift_ig5 * i * g5
			  if (shift_ig5 != 0)
			  {
			    COMPLEX ishift =
			      static_cast<COMPLEX>(std::complex<double>{0.0, shift_ig5});
			    if (ns == 1)
			    {
			      x.scale(ishift).addTo(y);
			    }
			    else
			    {
			      contract<NOp + 1>(g5.rename_dims({{'j', 's'}}).scale(ishift), x, "s",
						AddTo, y, {{'s', 'i'}});
			    }
			  }
			},
			op.d, op.i, nullptr, op};
	}

	// Get the solver
	Maybe<const Options&> solverOps = getOptionsMaybe(ops, "solver");
	if (solverOps)
	  return getSolver(shifted_op, solverOps.getSome());
	return shifted_op;
      }

      /// Returns the conjugate transpose of an operator
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getDagger(Operator<NOp, COMPLEX> op, const Options& ops,
				       Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getDagger: unsupported input preconditioner");

	int ns = op.d.kvdim().at('s');
	auto g5 = getGamma5<COMPLEX>(ns, op.d.getDev(), op.d.dist);

	// Return the solver
	return {
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    if (ns == 1)
	    {
	      x.copyTo(y);
	    }
	    else
	    {
	      // y = g5 * op * g5 * x
	      auto y0 = op(
		contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), x, "s").rename_dims({{'i', 's'}}));
	      contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), y0, "s", CopyTo,
				y.rename_dims({{'s', 'i'}}));
	    }
	  },
	  op.d, op.i, nullptr, op};
      }

      /// Returns an operator that applies \gamma_5
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getG5(Operator<NOp, COMPLEX> op, const Options& ops,
				   Operator<NOp, COMPLEX> prec_)
      {
	if (prec_)
	  throw std::runtime_error("getG5: unsupported input preconditioner");

	int ns = op.d.kvdim().at('s');
	auto g5 = getGamma5<COMPLEX>(ns, op.d.getDev(), op.d.dist);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  if (ns == 1)
		  {
		    x.copyTo(y);
		  }
		  else
		  {
		    // y = g5 * x
		    contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), x, "s", CopyTo,
				      y.rename_dims({{'s', 'i'}}));
		  }
		},
		op.d,
		op.i,
		nullptr,
		op.order_t,
		op.domLayout,
		op.imgLayout,
		getNeighbors(op.d.kvdim(), 0, op.domLayout),
		op.preferred_col_ordering,
		op.is_kronecker()};
      }

      /// Returns a solver with possible different precision than the operator's
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getCasting(Operator<NOp, COMPLEX> op, const Options& ops,
					Operator<NOp, COMPLEX> prec_, SolverSpace solverSpace)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup casting ") + prefix);

	// Get the current precision and the requested by the user
	enum Precision { Single, Double, Default };
	static const std::map<std::string, Precision> precisionMap{
	  {"default", Default}, {"single", Single}, {"float", Single}, {"double", Double}};
	Precision defaultPrecision = std::is_same<COMPLEX, ComplexD>::value ? Double : Single;
	Precision requestedPrecision =
	  getOption<Precision>(ops, "precision", precisionMap, Default);
	if (requestedPrecision == Default)
	  requestedPrecision = defaultPrecision;

	// Get the solver options
	const Options& solverOps = getOptions(ops, "solver");

	if (requestedPrecision == Double)
	{
	  return getSolver(op.template cast<ComplexD>(), solverOps, prec_.template cast<ComplexD>(),
			   solverSpace)
	    .template cast<COMPLEX>();
	}
	else
	{
	  return getSolver(op.template cast<ComplexF>(), solverOps, prec_.template cast<ComplexF>(),
			   solverSpace)
	    .template cast<COMPLEX>();
	}
      }
    }

    /// Returns an operator that approximate the inverse of a given operator
    /// \param op: operator to make the inverse of
    /// \param ops: options to select the solver from `solvers` and influence the solver construction

    template <std::size_t NOp, typename COMPLEX>
    Operator<NOp, COMPLEX> getSolver(const Operator<NOp, COMPLEX>& op, const Options& ops,
				     const Operator<NOp, COMPLEX>& prec, SolverSpace solverSpace)
    {
      enum SolverType {
	FGMRES,
	BICGSTAB,
	MR,
	GCR,
	MG,
	EO,
	HIE,
	DD,
	BJ,
	SHIFT,
	PROJ,
	SPINEO,
	IGD,
	DDAG,
	G5,
	BLOCKING,
	CASTING
      };
      static const std::map<std::string, SolverType> solverTypeMap{{"fgmres", FGMRES},
								   {"bicgstab", BICGSTAB},
								   {"mr", MR},
								   {"gcr", GCR},
								   {"mg", MG},
								   {"eo", EO},
								   {"hie", HIE},
								   {"dd", DD},
								   {"bj", BJ},
								   {"shift", SHIFT},
								   {"proj", PROJ},
								   {"spineo", SPINEO},
								   {"igd", IGD},
								   {"g5", G5},
								   {"blocking", BLOCKING},
								   {"casting", CASTING}};
      SolverType solverType = getOption<SolverType>(ops, "type", solverTypeMap);
      switch (solverType)
      {
      case FGMRES: // flexible GMRES
	return detail::getFGMRESSolver(op, ops, prec);
      case BICGSTAB: // bicgstab
	return detail::getBicgstabSolver(op, ops, prec);
      case MR: // minimal residual
	return detail::getMRSolver(op, ops, prec);
      case GCR: // generalized conjugate residual
	return detail::getGCRSolver(op, ops, prec);
      case MG: // Multigrid
	return detail::getMGPrec(op, ops, prec, solverSpace);
      case EO: // even-odd Schur preconditioner
	return detail::getEvenOddPrec(op, ops, prec, solverSpace);
      case HIE: // hierarchical Schur preconditioner
	return detail::getHierarchicalPrec(op, ops, prec);
      case DD: // domain decomposition with domains local to processes
	return detail::getDomainDecompositionPrec(op, ops, prec);
      case BJ: // block Jacobi
	return detail::getBlockJacobi(op, ops, prec);
      case SHIFT: // shift operator
	return detail::getShiftedOp(op, ops, prec);
      case PROJ: // projector preconditioner
	return detail::getProjPrec(op, ops, prec);
      case SPINEO: // even-odd Schur complement on the spin components
	return detail::getSpinEvenOddPrec(op, ops, prec);
      case IGD: // inexact Generalized Davidson
	return detail::getInexactGD(op, ops, prec);
      case DDAG: // return the operator conjugate transposed
	return detail::getDagger(op, ops, prec);
      case G5: // apply \gamma_5
	return detail::getG5(op, ops, prec);
      case BLOCKING: // reshape the operator
	return detail::getBlocking(op, ops, prec);
      case CASTING: // change the precision
	return detail::getCasting(op, ops, prec, solverSpace);
      }
      throw std::runtime_error("This shouldn't happen");
    }

    /// Return an Operator that wraps up a LinearOperator<LatticeFermion>
    inline Operator<Nd + 7, Complex> asOperatorView(const LinearOperator<LatticeFermion>& linOp,
						    bool use_kron_format = true)
    {
      LatticeFermion a;
      auto d = asTensorView(a).toComplex();
      auto blkd =
	d.template reshape_dimensions<Nd + 7>({{"x", "0x"}, {"y", "1y"}, {"z", "2z"}, {"t", "3t"}},
					      {{'0', 1}, {'1', 1}, {'2', 1}, {'3', 1}}, true)
	  .make_eg();
      const auto dim = d.kvdim();
      const auto blkdim = blkd.kvdim();
      return {
	[&, blkdim](Tensor<Nd + 8, Complex> x, Tensor<Nd + 8, Complex> y) {
	  Tracker _t("chroma's matvec ");
	  unsigned int n = x.kvdim().at('n');
	  _t.arity = n;
	  auto tx = x.template reshape_dimensions<Nd + 4>(
	    {{"0x", "x"}, {"1y", "y"}, {"2z", "z"}, {"3t", "t"}}, {}, true);
	  auto ty = tx.make_compatible();
	  LatticeFermion x0, y0;
	  for (unsigned int i = 0; i < n; ++i)
	  {
	    tx.kvslice_from_size({{'n', i}}, {{'n', 1}}).copyTo(asTensorView(x0));
	    y0 = zero;
	    linOp(y0, x0, PLUS /* I believe, it's ignored */);
	    asTensorView(y0).copyTo(ty.kvslice_from_size({{'n', i}}, {{'n', 1}}));
	  }
	  ty.template reshape_dimensions<Nd + 8>(
	      {{"x", "0x"}, {"y", "1y"}, {"z", "2z"}, {"t", "3t"}}, blkdim, true)
	    .copyTo(y);
	},
	blkd,	 // domain
	blkd,	 // image
	nullptr, // no conjugate
	"",	 // no order_t
	XEvenOddLayout,
	XEvenOddLayout,
	detail::getNeighbors(dim, 1 /* near-neighbors links only */, XEvenOddLayout),
	ColumnMajor,	// preferred ordering
	use_kron_format /* has a Kronecker form */
      };
    }

    //
    // High-level chroma operations
    //

    /// Constructor
    /// \param fermAction: XML for the fermion action
    /// \param invParam: XML for the quark propagator
    /// \param u: gauge fields

    ChimeraSolver::ChimeraSolver(const GroupXML_t& fermAction, const GroupXML_t& invParam,
				 const multi1d<LatticeColorMatrix>& u)
    {
      // Initialize fermion action and state
      std::istringstream xml_s(fermAction.xml);
      XMLReader fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << fermAction.id << std::endl;
      S = TheFermionActionFactory::Instance().createObject(fermAction.id, fermacttop,
							   fermAction.path);
      state = S->createState(u);

      // If the inverter is MGPROTON, use this infrastructure
      if (invParam.id == std::string("MGPROTON"))
      {
	QDPIO::cout << "Setting up MGPROTON invertor..." << std::endl;
	Tracker _t("setup mgproton");

	// Parse XML with the inverter options
	std::shared_ptr<Options> ops = getOptionsFromXML(broadcast(invParam.xml));

	// Clone the matvec
	LinearOperator<LatticeFermion>* fLinOp = S->genLinOp(state);
	ColOrdering co = getOption<ColOrdering>(*ops, "InvertParam/operator_ordering",
						getColOrderingMap(), ColumnMajor);
	ColOrdering co_blk = getOption<ColOrdering>(*ops, "InvertParam/operator_block_ordering",
						    getColOrderingMap(), RowMajor);
	Operator<Nd + 7, Complex> linOp = detail::cloneOperator(
	  asOperatorView(*fLinOp), co, co_blk, detail::ConsiderBlockingSparse, "chroma's operator");

	// Destroy chroma objects
	delete fLinOp;
	state = State();
	S = Action();

	// Construct the solver
	op = Operator<Nd + 7, Complex>(getSolver(linOp, getOptions(*ops, "InvertParam")));

	// Clean cache of operators
	detail::cleanOperatorsCache();

	QDPIO::cout << "MGPROTON invertor ready; setup time: "
		    << detail::tostr(_t.stopAndGetElapsedTime()) << " s" << std::endl;
      }
      else
      {
	PP = S->qprop(state, invParam);
      }
    }

    namespace detail
    {
      /// Apply the inverse to LatticeColorVec tensors for a list of spins
      /// \param PP: invertor
      /// \param chi: lattice color tensor on a t_slice, cxyzXn
      /// \param t_source: time-slice in chi
      /// \param Nt_forward: return the next Nt_forward time-slices after t_source
      /// \param Nt_backward: return the previous Nt_backward time-slices before t_source
      /// \param spin_sources: list of spins
      /// \param max_rhs: maximum number of vectors solved at once
      /// \param order_out: coordinate order of the output tensor, a permutation of cSxyztXns where
      ///        s is the spin source and S is the spin sink
      /// \return: tensor cSxyztXns where the first t_slice is the t_source-Nt_backward time-slice of
      ///        the vectors after the inversion, and goes increasingly until time-source t_source+Nt_forward

      template <typename COMPLEX_CHI, typename COMPLEX_OUT>
      Tensor<Nd + 5, COMPLEX_OUT> doInversion(const SystemSolver<LatticeFermion>& PP,
					      const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
					      int first_tslice_out, int n_tslice_out,
					      const std::vector<int>& spin_sources, int max_rhs,
					      const std::string& order_out = "cSxyztXns")
      {
	int num_vecs = chi.kvdim()['n'];
	Tensor<Nd + 5, COMPLEX_OUT> psi(
	  order_out,
	  latticeSize<Nd + 5>(
	    order_out,
	    {{'t', n_tslice_out}, {'S', Ns}, {'s', spin_sources.size()}, {'n', num_vecs}}),
	  chi.getDev());

	int max_step = std::max(num_vecs, max_rhs);
	std::vector<std::shared_ptr<LatticeFermion>> chis(max_step), quark_solns(max_step);
	for (int col = 0; col < max_step; col++)
	  chis[col].reset(new LatticeFermion);
	for (int col = 0; col < max_step; col++)
	  quark_solns[col].reset(new LatticeFermion);

	for (int spin_source : spin_sources)
	{
	  for (int n0 = 0, n_step = std::min(max_rhs, num_vecs); n0 < num_vecs;
	       n0 += n_step, n_step = std::min(n_step, num_vecs - n0))
	  {
	    for (int n = n0, col = 0; col < n_step; ++n, ++col)
	    {
	      // Put the colorvec sources for the t_source on chis for spin `spin_source`
	      // chis[col][s=spin_source] = chi[n=n0]
	      *chis[col] = zero;
	      chi.kvslice_from_size({{'n', n}}, {{'n', 1}})
		.copyTo(SB::asTensorView(*chis[col])
			  .kvslice_from_size({{'t', t_source}, {'s', spin_source}}));

	      *quark_solns[col] = zero;
	    }

	    // Solve
	    std::vector<SystemSolverResults_t> res =
	      PP(std::vector<std::shared_ptr<LatticeFermion>>(quark_solns.begin(),
							      quark_solns.begin() + n_step),
		 std::vector<std::shared_ptr<const LatticeFermion>>(chis.begin(),
								    chis.begin() + n_step));

	    for (int n = n0, col = 0; col < n_step; ++n, ++col)
	    {
	      // psi[n=n] = quark_solns[col][t=first_tslice+(0:n_tslice_out-1)]
	      asTensorView(*quark_solns[col])
		.kvslice_from_size({{'t', first_tslice_out}}, {{'t', n_tslice_out}})
		.rename_dims({{'s', 'S'}})
		.copyTo(psi.kvslice_from_size({{'n', n}, {'s', spin_source}}));
	    }
	  }
	}

	return psi;
      }

      /// Apply the inverse to LatticeFermion tensors
      /// \param sol: invertor, "linear" operator in cs0123xyztX
      /// \param chi: spin-color lattice tensor, csxyztXn
      /// \param max_rhs: maximum number of vectors solved at once
      /// \return: tensor with the same labels as the input

      template <std::size_t N, typename COMPLEX_CHI, typename COMPLEX_OUT>
      Tensor<N, COMPLEX_OUT> doInversion(const Operator<Nd + 7, COMPLEX_OUT>& op,
					 const Tensor<N, COMPLEX_CHI>& chi, int max_rhs)
      {
	Tracker _t(std::string("mgproton solver"));

	// Get the columns labels, which are the ones not contracted with the operator
	std::string order_cols = remove_dimensions(chi.order, op.i.order);

	// Create tensors with full support on the lattice
	auto x0 = chi.template reshape_dimensions<Nd + 8>(
	  {{"x", "0x"}, {"y", "1y"}, {"z", "2z"}, {"t", "3t"}, {order_cols, "n"}},
	  {{'0', 1}, {'1', 1}, {'2', 1}, {'3', 0}}, true);
	auto y0 = x0.make_compatible();
	foreachInChuncks(
	  x0, y0, max_rhs,
	  [=](Tensor<Nd + 8, COMPLEX_CHI> x, Tensor<Nd + 8, COMPLEX_CHI> y) { op(x, y); });
	return y0.template reshape_dimensions<N, COMPLEX_OUT>({{"n", order_cols}});
      }

      /// Apply the inverse to LatticeColorVec tensors for a list of spins
      /// \param sol: invertor
      /// \param chi: lattice color tensor on a t_slice, cxyzXn
      /// \param t_source: time-slice in chi
      /// \param Nt_forward: return the next Nt_forward time-slices after t_source
      /// \param Nt_backward: return the previous Nt_backward time-slices before t_source
      /// \param spin_sources: list of spins
      /// \param max_rhs: maximum number of vectors solved at once
      /// \param order_out: coordinate order of the output tensor, a permutation of cSxyztXns where
      ///        s is the spin source and S is the spin sink
      /// \return: tensor cSxyztXns where the first t_slice is the t_source-Nt_backward time-slice of
      ///        the vectors after the inversion, and goes increasingly until time-source t_source+Nt_forward

      template <typename COMPLEX_CHI, typename COMPLEX_OUT>
      Tensor<Nd + 5, COMPLEX_OUT> doInversion(const Operator<Nd + 7, COMPLEX_OUT>& op,
					      const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
					      int first_tslice_out, int n_tslice_out,
					      const std::vector<int>& spin_sources, int max_rhs,
					      const std::string& order_out = "cSxyztXns")
      {
	Tracker _t(std::string("mgproton solver"));

	int num_vecs = chi.kvdim()['n'];
	Tensor<Nd + 5, COMPLEX_OUT> psi(
	  order_out,
	  latticeSize<Nd + 5>(
	    order_out,
	    {{'t', n_tslice_out}, {'S', Ns}, {'s', spin_sources.size()}, {'n', num_vecs}}),
	  chi.getDev());

	// Create tensors with full support on the lattice
	int max_step = std::max(num_vecs, max_rhs);
	auto aux = chi.template make_compatible<Nd + 8>(
	  op.preferred_col_ordering == ColumnMajor ? "0123csxyztXn" : "0123ncsxyztX",
	  {{'n', max_step},
	   {'t', Layout::lattSize()[3]},
	   {'s', Ns},
	   {'0', 1},
	   {'1', 1},
	   {'2', 1},
	   {'3', 1}});

	for (int spin_source : spin_sources)
	{
	  for (int n0 = 0, n_step = std::min(max_rhs, num_vecs); n0 < num_vecs;
	       n0 += n_step, n_step = std::min(n_step, num_vecs - n0))
	  {
	    auto aux0 = aux.kvslice_from_size({}, {{'n', n_step}});
	    aux0.set_zero();
	    chi.kvslice_from_size({{'n', n0}}, {{'n', n_step}})
	      .copyTo(aux0.kvslice_from_size({{'t', t_source}, {'s', spin_source}}));

	    // Solve
	    op(aux0)
	      .kvslice_from_size({{'t', first_tslice_out}}, {{'t', n_tslice_out}})
	      .rename_dims({{'s', 'S'}})
	      .copyTo(psi.kvslice_from_size({{'n', n0}, {'s', spin_source}}));
	  }
	}

	return psi;
      }

      /// Apply the inverse to a list of LatticeFermions
      /// \param PP: invertor
      /// \param chi: lattice spin-color field tensor, csxyztX
      /// \param max_rhs: maximum number of vectors solved at once
      /// \return:
      template <typename COMPLEX, std::size_t N>
      Tensor<N, COMPLEX> doInversion(const SystemSolver<LatticeFermion>& PP,
				     const Tensor<N, COMPLEX>& chi, int max_rhs)
      {
	detail::check_order_contains(chi.order, "csxyztX");
	std::string n_order = detail::remove_dimensions(chi.order, "csxyztX");
	Coor<N - 7> n_dim = latticeSize<N - 7>(n_order, chi.kvdim());
	int n_vol = (N == 7 ? 1 : superbblas::detail::volume(n_dim));

	Tensor<N, COMPLEX> r = chi.make_compatible(); // output tensor
	int max_step = std::max(1, std::max(n_vol, max_rhs));

	// Quick exit
	if (n_vol == 0)
	  return r;

	if (N == 7)
	{
	  // For a single vector
	  LatticeFermion chi0, psi0;
	  chi.copyTo(asTensorView(chi0));
	  SystemSolverResults_t res = PP(psi0, chi0);
	  asTensorView(psi0).copyTo(r);
	}
	else
	{
	  // Auxiliary LatticeFermion
	  std::vector<std::shared_ptr<LatticeFermion>> chis(max_step), quark_solns(max_step);
	  for (int col = 0; col < max_step; col++)
	    chis[col].reset(new LatticeFermion);
	  for (int col = 0; col < max_step; col++)
	    quark_solns[col].reset(new LatticeFermion);

	  Coor<N - 7> n_strides = detail::get_strides<Index>(n_dim, superbblas::FastToSlow);
	  for (int n0 = 0, n_step = std::min(max_rhs, n_vol); n0 < n_vol;
	       n0 += n_step, n_step = std::min(n_step, n_vol - n0))
	  {
	    for (int n = n0, col = 0; col < n_step; ++n, ++col)
	    {
	      // Get the field to copy from the tensor chi
	      Coor<N - 7> ni = detail::index2coor<N - 7>(n, n_dim, n_strides);
	      std::map<char, int> from{}, size{};
	      for (int d = 0; d < N - 7; ++d)
		from[n_order[d]] = ni[d], size[n_order[d]] = 1;

	      // Copy the field into a LatticeFermion
	      chi.kvslice_from_size(from, size).copyTo(asTensorView(*chis[col]));

	      *quark_solns[col] = zero;
	    }

	    // Solve
	    std::vector<SystemSolverResults_t> res =
	      PP(std::vector<std::shared_ptr<LatticeFermion>>(quark_solns.begin(),
							      quark_solns.begin() + n_step),
		 std::vector<std::shared_ptr<const LatticeFermion>>(chis.begin(),
								    chis.begin() + n_step));

	    for (int n = n0, col = 0; col < n_step; ++n, ++col)
	    {
	      // Get the field to copy from the tensor chi
	      Coor<N - 7> ni = detail::index2coor<N - 7>(n, n_dim, n_strides);
	      std::map<char, int> from{}, size{};
	      for (int d = 0; d < N - 7; ++d)
		from[n_order[d]] = ni[d], size[n_order[d]] = 1;

	      // Copy from LatticeFermion to the output tensor
	      asTensorView(*quark_solns[col]).copyTo(r.kvslice_from_size(from, size));
	    }
	  }
	}

	return r;
      }

      /// Returns a projector onto the approximate right singular value space of the given operator
      /// for the smallest singular values.
      ///
      /// It can work as a projector on the right of inv(op):
      ///   P = V*inv(V'*g5*op*V)*V'*g5*op
      /// Then:
      ///   tr P*inv(op) = tr V*inv(V'*g5*op*V)*V'*g5 = \sum_i v_i'*g5*vi/lambda_i
      ///
      /// It computes a number of the smallest right singular vectors of op as the eigenvectors of
      /// inv(g5*op), which is Hermitian.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Projector<NOp, COMPLEX> getDeflationProj(Operator<NOp, COMPLEX> op, const Options& ops)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup defl proj ") + prefix);

	// Get options
	unsigned int rank = getOption<unsigned int>(ops, "rank");
	if (rank == 0)
	{
	  // Return trivial solution
	  auto almost_proj = Operator<NOp, COMPLEX>{
	    [=](const Tensor<NOp + 1, COMPLEX>&, Tensor<NOp + 1, COMPLEX> y) { y.set_zero(); },
	    op.d, op.i, nullptr, op};

	  auto VUfun = [=](unsigned int ifrom, unsigned int isize,
			   char col_label) -> Tensor<NOp + 1, COMPLEX> { return {}; };

	  std::vector<COMPLEX> lambdas(0);

	  return {almost_proj, VUfun, VUfun, lambdas, op};
	}
	double tol = getOption<double>(ops, "tol", 0.1);
	auto solver = getSolver(op, ops.getValue("solver"));
	auto default_eig_ops = DictionaryOption(ops);
	auto eigensolver = SB::getInexactEigensolverGD(
	  solver, ops.getValue("eigensolver", Maybe<const Option&>(default_eig_ops)));

	// Compute the eigenpairs of inv(op) * g5, the right singular vectors of op
	auto values_vectors = eigensolver(rank, tol);
	auto inv_values = std::get<0>(values_vectors);
	auto vectors = std::get<1>(values_vectors).rename_dims({{'n', 'i'}});

	// Return the projector
	int ns = op.d.kvdim().at('s');
	auto g5 = getGamma5<COMPLEX>(ns, op.d.getDev(), op.d.dist);
	auto mult_by_g5 = [=](Tensor<NOp + 1, COMPLEX> v) {
	  char i = detail::get_free_label(v.order);
	  return contract<NOp + 1>(g5.rename_dims({{'j', 's'}, {'i', i}}), v, "s")
	    .rename_dims({{i, 's'}});
	};
	auto Vt_g5_op_V = contract<2>(
	  vectors.conj(), mult_by_g5(op(vectors.rename_dims({{'i', 'j'}}))), op.d.order);
	auto inv_Vt_g5_op_V = inv(Vt_g5_op_V, "i", "j");
	auto almost_proj = Operator<NOp, COMPLEX>{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    // Do V'*g5*x, dims: i x n
	    auto Vt_g5_x = contract<2>(vectors.conj(), mult_by_g5(x), op.d.order);

	    // Do inv(V'*g5*op*V)*V'*g5*x, dims I x i * i x n = I x n -> i x n
	    auto lambda_inv_Vt_g5_x =
	      contract<2>(inv_Vt_g5_op_V.rename_dims({{'i', 'I'}, {'j', 'i'}}), Vt_g5_x, "i")
		.rename_dims({{'I', 'i'}});

	    // Do y = V*inv(V'*g5*op*V)*V'*g5*x, dims latt x i * i x n = latt x n
	    contract(vectors, lambda_inv_Vt_g5_x, "i", CopyTo, y);
	  },
	  op.d, op.i, nullptr, solver};

	// Return V*inv_Vt_g5_op_V[:,ifrom..ifrom+isize]
	auto Vfun = [=](unsigned int ifrom, unsigned int isize, char col_label) {
	  return contract<NOp + 1>(
		   vectors, inv_Vt_g5_op_V.kvslice_from_size({{'j', ifrom}}, {{'j', isize}}), "i")
	    .rename_dims({{'j', col_label}});
	};

	// Returns \gamma_5*V[ifrom..ifrom+isize-1]
	auto Ufun = [=](unsigned int ifrom, unsigned int isize, char col_label) {
	  return mult_by_g5(vectors.kvslice_from_size({{'i', ifrom}}, {{'i', isize}})
			      .rename_dims({{'i', col_label}}));
	};

	// Return all ones
	std::vector<COMPLEX> lambdas(inv_values.size(), COMPLEX{1});

	return {almost_proj, Vfun, Ufun, lambdas, op};
      }

      /// Returns a projector onto the approximate right singular value space of the given operator
      /// for the smallest singular values.
      ///
      /// It can work as a projector on the right of inv(op):
      ///   P = Q*V*inv(V'*g5*Q'*op*Q*V)*V'*g5*Q'*op
      /// Then:
      ///   tr P*inv(op) = tr Q*V*inv(V'*g5*Q'*op*Q*V)*V'*g5*Q' = \sum_i v_i'*g5*vi/lambda_i
      ///
      /// It computes a number of the smallest right singular vectors of op as the eigenvectors of
      /// inv(g5*op), which is Hermitian.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Projector<NOp, COMPLEX> getMGDeflationProj(Operator<NOp, COMPLEX> op, const Options& ops)
      {
	std::string prefix = getOption<std::string>(ops, "prefix", "");
	Tracker _t(std::string("setup mg proj ") + prefix);

	// Get prolongator and coarse operator
	auto prolongator_coarse_spin_splitting =
	  getProlongatorAndCoarse(op, ops.getValue("prolongator"), FullSpace);
	Operator<NOp, COMPLEX> Q = std::get<0>(prolongator_coarse_spin_splitting);
	Operator<NOp, COMPLEX> op_c = std::get<1>(prolongator_coarse_spin_splitting);

	// Return the projector for the coarse operator
	auto proj_c = getProjector(op_c, ops.getValue("proj"));

	/// Return Q*V*inv(U'*Q'*op*Q*V)*U'*Q' = Q proj_c * Q'
	auto almost_proj = Operator<NOp, COMPLEX>{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    Q(proj_c.V_inv_Ut(Q.tconj()(x)), y);
	  },
	  op.d,
	  op.i,
	  nullptr,
	  op.order_t,
	  op.imgLayout,
	  op.domLayout,
	  DenseOperator(),
	  op.preferred_col_ordering,
	  false /* no Kronecker blocking */};

	// Return Q*V[:,ifrom..ifrom+isize]
	auto Vfun = [=](unsigned int ifrom, unsigned int isize, char col_label) {
	  return Q(proj_c.V(ifrom, isize, col_label));
	};

	// Returns Q*U[ifrom..ifrom+isize-1]
	auto Ufun = [=](unsigned int ifrom, unsigned int isize, char col_label) {
	  return Q(proj_c.U(ifrom, isize, col_label));
	};

	return {almost_proj, Vfun, Ufun, proj_c.lambdas, op};
      }
    }

    /// Apply the inverse to LatticeColorVec tensors for a list of spins
    /// \param PP: invertor
    /// \param chi: lattice color tensor on a t_slice, cxyzXn
    /// \param t_source: time-slice in chi
    /// \param Nt_forward: return the next Nt_forward time-slices after t_source
    /// \param Nt_backward: return the previous Nt_backward time-slices before t_source
    /// \param spin_sources: list of spins
    /// \param max_rhs: maximum number of vectors solved at once
    /// \param order_out: coordinate order of the output tensor, a permutation of cSxyztXns where
    ///        s is the spin source and S is the spin sink
    /// \return: tensor cSxyztXns where the first t_slice is the t_source-Nt_backward time-slice of
    ///        the vectors after the inversion, and goes increasingly until time-source t_source+Nt_forward

    template <typename COMPLEX_CHI>
    Tensor<Nd + 5, SB::Complex>
    doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
		int first_tslice_out, int n_tslice_out, const std::vector<int>& spin_sources,
		int max_rhs, const std::string& order_out)
    {
      detail::check_order_contains(order_out, "cSxyztXns");
      if (chi.kvdim()['t'] != 1)
	throw std::runtime_error("Expected one time-slice");
      const int num_vecs = chi.kvdim()['n'];

      if (n_tslice_out > Layout::lattSize()[3])
	throw std::runtime_error("Too many tslices");

      StopWatch snarss1;
      snarss1.reset();
      snarss1.start();

      Tensor<Nd + 5, SB::Complex> r;
      if (sol.op)
	r = detail::doInversion<COMPLEX_CHI, SB::Complex>(
	  sol.op, chi, t_source, first_tslice_out, n_tslice_out, spin_sources, max_rhs, order_out);
      else
	r = detail::doInversion<COMPLEX_CHI, SB::Complex>(
	  *sol.PP, chi, t_source, first_tslice_out, n_tslice_out, spin_sources, max_rhs, order_out);

      snarss1.stop();
      QDPIO::cout << "Time to compute inversions for " << spin_sources.size()
		  << " spin sources and " << num_vecs
		  << " colorvecs : " << snarss1.getTimeInSeconds() << " secs" << std::endl;

      return r;
    }

    template Tensor<Nd + 5, SB::Complex>
    doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, SB::ComplexD> chi, int t_source,
		int first_tslice_out, int n_tslice_out, const std::vector<int>& spin_sources,
		int max_rhs, const std::string& order_out);
    template Tensor<Nd + 5, SB::Complex>
    doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, SB::ComplexF> chi, int t_source,
		int first_tslice_out, int n_tslice_out, const std::vector<int>& spin_sources,
		int max_rhs, const std::string& order_out);

    /// Apply the inverse to a list of LatticeFermions
    /// \param sol: Chimera invertor
    /// \param chi: lattice spin-color tensor, at least dimensions csxyztX
    /// \param max_rhs: maximum number of vectors solved at once
    /// \param conjugate: whether to apply the invertor transpose-conjugate
    /// \return: tensor with the same ordering as `chi`.

    // template <std::size_t N, typename COMPLEX_CHI>
    // Tensor<N, SB::Complex> doInversion(const ChimeraSolver& sol,
    //   				 const Tensor<N, COMPLEX_CHI>& chi, int max_rhs = 0,
    //   				 Conjugation conj = NotConjugate)
    // {
    //   detail::check_order_contains(chi.order, "csxyztX");
    //   const int num_vecs =
    //     (N == 7 ? 1
    //   	  : detail::volume(chi.kvdim(), detail::remove_dimensions(chi.order, "csxyztX")));

    //   StopWatch snarss1;
    //   snarss1.reset();
    //   snarss1.start();

    //   // Multiply the input by g5 if applied conjugate
    //   Tensor<2, COMPLEX_CHI> g5;
    //   if (conj == Conjugate)
    //   {
    //     g5 = Gamma(Ns * Ns - 1).cloneOn<COMPLEX_CHI>(chi.getDev());
    //     chi.contract(chi, {}, NotConjugate, g5, {{'j', 's'}}, NotConjugate, {{'s', 'i'}});
    //   }

    //   Tensor<N, SB::Complex> r;
    //   if (sol.op)
    //     r = detail::doInversion<N, COMPLEX_CHI, COMPLEX_OUT>(sol.op, chi, max_rhs);
    //   else
    //     r = detail::doInversion<N, COMPLEX_CHI, COMPLEX_OUT>(*sol.PP, chi, max_rhs);

    //   // Multiply the input by g5 if applied conjugate
    //   if (conj == Conjugate)
    //   {
    //     r.contract(r, {}, NotConjugate, g5, {{'j', 's'}}, NotConjugate, {{'s', 'i'}});
    //   }

    //   snarss1.stop();
    //   QDPIO::cout << "Time to compute " << num_vecs
    //   	    << " inversions: " << snarss1.getTimeInSeconds() << " secs" << std::endl;

    //   return r;
    // }

    /// Apply the inverse to a list of LatticeFermions
    /// \param sol: Chimera invertor
    /// \param psis: output lattice spin-color tensor, at least dimensions csxyztX
    /// \param chis: input lattice spin-color tensor, at least dimensions csxyztX
    /// \param max_rhs: maximum number of vectors solved at once

    void doInversion(const ChimeraSolver& sol, MultipleLatticeFermions& psis,
		     const ConstMultipleLatticeFermions& chis, int max_rhs)
    {
      StopWatch snarss1;
      snarss1.reset();
      snarss1.start();

      if (max_rhs <= 0)
	max_rhs = chis.size();

      // Do the inversion
      if (sol.op)
      {
	auto op = sol.op;
	auto tchi = op.make_compatible_dom<Nd + 8>("n", {{'n', max_rhs}});
	auto tpsi = op.make_compatible_img<Nd + 8>("n", {{'n', max_rhs}});
	for (int i = 0, n = std::min(max_rhs, (int)chis.size()); i < chis.size();
	     i += n, n = std::min((int)chis.size() - i, max_rhs))
	{
	  // Adjust the size of tchi and tpsi to n
	  auto this_tchi = tchi.kvslice_from_size({{'n', 0}}, {{'n', n}});
	  auto this_tpsi = tpsi.kvslice_from_size({{'n', 0}}, {{'n', n}});

	  // Copy chis into this_tchi
	  for (int j = 0; j < n; ++j)
	    asTensorView(*chis[i + j]).copyTo(this_tchi.kvslice_from_size({{'n', j}}, {{'n', 1}}));

	  // Do the inversion: this_tpsi = D^{-1} * this_tchi
	  op(this_tchi, this_tpsi);

	  // Copy the solution into psis
	  for (int j = 0; j < n; ++j)
	    this_tpsi.kvslice_from_size({{'n', j}}, {{'n', 1}}).copyTo(asTensorView(*psis[i + j]));
	}
      }
      else
      {
	for (int i = 0, n = std::min(max_rhs, (int)chis.size()); i < chis.size();
	     i += n, n = std::min((int)chis.size() - i, max_rhs))
	{
	  (*sol.PP)(MultipleLatticeFermions(psis.begin() + i, psis.begin() + i + n),
		    ConstMultipleLatticeFermions(chis.begin() + i, chis.begin() + i + n));
	}
      }

      snarss1.stop();
      QDPIO::cout << "Time to compute " << chis.size()
		  << " inversions: " << snarss1.getTimeInSeconds() << " secs" << std::endl;
    }

    /// Apply the inverse to a list of LatticeFermions
    /// \param sol: Chimera invertor
    /// \param psis: output lattice spin-color tensor, at least dimensions csxyztX
    /// \param chis: input lattice spin-color tensor, at least dimensions csxyztX
    /// \param max_rhs: maximum number of vectors solved at once

    Operator<Nd + 7, Complex> getOperator(const ChimeraSolver& sol, int max_rhs)
    {
      if (sol.op)
	return sol.op;

      LatticeFermion a;
      auto blkd =
	asTensorView(a)
	  .toComplex()
	  .template reshape_dimensions<Nd + 7>({{"x", "0x"}, {"y", "1y"}, {"z", "2z"}, {"t", "3t"}},
					       {{'0', 1}, {'1', 1}, {'2', 1}, {'3', 1}}, true)
	  .make_eg();
      auto PP_handle = sol.PP;
      return {
	[=](Tensor<Nd + 8, Complex> x, Tensor<Nd + 8, Complex> y) {
	  Tracker _t("chroma's inversion ");
	  detail::doInversion(*PP_handle, x, max_rhs).copyTo(y);
	},
	blkd,	 // domain
	blkd,	 // image
	nullptr, // no conjugate
	"",	 // no order_t
	XEvenOddLayout,
	XEvenOddLayout,
	DenseOperator(),
	ColumnMajor, // preferred ordering
	false	     /* has a Kronecker form */
      };
    }

    /// Returns a projector onto the smallest singular space on an operator
    /// \param op: operator to make the projector onto
    /// \param ops: options to select the projector and influence the projector construction

    template <std::size_t NOp, typename COMPLEX>
    Projector<NOp, COMPLEX> getProjector(const Operator<NOp, COMPLEX>& op, const Options& ops)
    {
      enum ProjectorType { DEFL, MGDEFL };
      static const std::map<std::string, ProjectorType> projectorTypeMap{{"defl", DEFL},
									 {"mg", MGDEFL}};
      ProjectorType projectorType = getOption<ProjectorType>(ops, "type", projectorTypeMap);
      switch (projectorType)
      {
      case DEFL: // plain deflation
	return detail::getDeflationProj(op, ops);
      case MGDEFL: // Multigrid deflation
	return detail::getMGDeflationProj(op, ops);
      }
      throw std::runtime_error("This shouldn't happen");
    }

    /// Constructor
    /// \param fermAction: XML for the fermion action
    /// \param projParam: XML for a projector onto the quark propagator
    /// \param u: gauge fields

    ChimeraProjector::ChimeraProjector(const GroupXML_t& fermAction, const GroupXML_t& projParam,
				       const multi1d<LatticeColorMatrix>& u)
    {
      // Initialize fermion action and state
      std::istringstream xml_s(fermAction.xml);
      XMLReader fermacttop(xml_s);
      QDPIO::cout << "FermAct = " << fermAction.id << std::endl;
      S = TheFermionActionFactory::Instance().createObject(fermAction.id, fermacttop,
							   fermAction.path);
      state = S->createState(u);

      // If the inverter is MGPROTON, use this infrastructure
      if (projParam.id == std::string("MGPROTON"))
      {
	QDPIO::cout << "Setting up MGPROTON projector..." << std::endl;
	Tracker _t("setup mgproton projector");

	// Parse XML with the inverter options
	std::shared_ptr<Options> ops = getOptionsFromXML(broadcast(projParam.xml));

	// Clone the matvec
	LinearOperator<LatticeFermion>* fLinOp = S->genLinOp(state);
	ColOrdering co = getOption<ColOrdering>(*ops, "Projector/operator_ordering",
						getColOrderingMap(), ColumnMajor);
	ColOrdering co_blk = getOption<ColOrdering>(*ops, "Projector/operator_block_ordering",
						    getColOrderingMap(), RowMajor);
	Operator<Nd + 7, Complex> linOp = detail::cloneOperator(
	  asOperatorView(*fLinOp), co, co_blk, detail::ConsiderBlockingSparse, "chroma's operator");

	// Destroy chroma objects
	delete fLinOp;
	state = State();
	S = Action();

	// Construct the projector
	op = getProjector(linOp, getOptions(*ops, "Projector"));

	// Clean cache of operators
	detail::cleanOperatorsCache();

	QDPIO::cout << "MGPROTON projector ready; setup time: "
		    << detail::tostr(_t.stopAndGetElapsedTime()) << " s" << std::endl;
      }
      else
      {
	chroma_proj = S->projector(state, projParam);
      }
    }

    /// Apply the oblique projector V*inv(U^H*op*V)*U^H*op
    /// \param proj: Chimera projector
    /// \param psis: output lattice spin-color tensor, at least dimensions csxyztX
    /// \param chis: input lattice spin-color tensor, at least dimensions csxyztX
    /// \param max_rhs: maximum number of vectors solved at once

    void doVUAObliqueProjector(const ChimeraProjector& proj, MultipleLatticeFermions& psis,
			       const ConstMultipleLatticeFermions& chis, int max_rhs)
    {
      if (chis.size() == 0)
	return;

      StopWatch snarss1;
      snarss1.reset();
      snarss1.start();

      if (max_rhs <= 0)
	max_rhs = chis.size();

      // Do the projection
      if (proj.op.op)
      {
	auto tchi = proj.op.op.make_compatible_dom<Nd + 8>("n", {{'n', max_rhs}});
	auto tpsi = proj.op.op.make_compatible_img<Nd + 8>("n", {{'n', max_rhs}});
	for (int i = 0, n = std::min(max_rhs, (int)chis.size()); i < chis.size();
	     i += n, n = std::min((int)chis.size() - i, max_rhs))
	{
	  // Adjust the size of tchi and tpsi to n
	  auto this_tchi = tchi.kvslice_from_size({{'n', 0}}, {{'n', n}});
	  auto this_tpsi = tpsi.kvslice_from_size({{'n', 0}}, {{'n', n}});

	  // Copy chis into this_tchi
	  for (int j = 0; j < n; ++j)
	    asTensorView(*chis[i + j]).copyTo(this_tchi.kvslice_from_size({{'n', j}}, {{'n', 1}}));

	  // Do the projection: this_tpsi = V * inv(U'*D*V) * U' * D * this_tchi
	  proj.op.V_inv_Ut(proj.op.op(this_tchi), this_tpsi);

	  // Copy the solution into psis
	  for (int j = 0; j < n; ++j)
	    this_tpsi.kvslice_from_size({{'n', j}}, {{'n', 1}}).copyTo(asTensorView(*psis[i + j]));
	}
      }
      else
      {
	for (int i = 0, n = std::min(max_rhs, (int)chis.size()); i < chis.size();
	     i += n, n = std::min((int)chis.size() - i, max_rhs))
	{
	  proj.chroma_proj->VUAObliqueProjector(
	    MultipleLatticeFermions(psis.begin() + i, psis.begin() + i + n),
	    ConstMultipleLatticeFermions(chis.begin() + i, chis.begin() + i + n));
	}
      }

      snarss1.stop();
      QDPIO::cout << "Time to compute " << chis.size()
		  << " projection: " << snarss1.getTimeInSeconds() << " secs" << std::endl;
    }

    /// Get the rank of the oblique projector V*inv(U^H*op*V)*U^H*op
    /// \param proj: Chimera projector

    unsigned int getProjectorRank(const ChimeraProjector& proj)
    {
      if (proj.op.op)
	return proj.op.lambdas.size();
      return proj.chroma_proj->rank();
    }

    /// Get the left basis of the oblique projector V*inv(U^H*op*V)*U^H*op
    /// \param proj: Chimera projector
    /// \param from: index of the first basis column
    /// \param psis: output lattice spin-color tensor, at least dimensions csxyztX

    void getV(const ChimeraProjector& proj, unsigned int from, MultipleLatticeFermions& psis)
    {
      if (!proj.op.op)
      {
	for (unsigned int i = 0; i < psis.size(); ++i)
	  proj.chroma_proj->V(from + i, *psis[i]);
      }
      else
      {
	auto V = proj.op.V(from, psis.size(), 'i');
	for (unsigned int i = 0; i < psis.size(); ++i)
	  V.kvslice_from_size({{'i', i}}, {{'i', 1}}).copyTo(asTensorView(*psis[i]));
      }
    }

    /// Get the right basis of the oblique projector V*inv(U^H*op*V)*U^H*op
    /// \param proj: Chimera projector
    /// \param from: index of the first basis column
    /// \param psis: output lattice spin-color tensor, at least dimensions csxyztX

    void getU(const ChimeraProjector& proj, unsigned int from, MultipleLatticeFermions& psis)
    {
      if (!proj.op.op)
      {
	for (unsigned int i = 0; i < psis.size(); ++i)
	  proj.chroma_proj->U(from + i, *psis[i]);
      }
      else
      {
	auto U = proj.op.U(from, psis.size(), 'i');
	for (unsigned int i = 0; i < psis.size(); ++i)
	  U.kvslice_from_size({{'i', i}}, {{'i', 1}}).copyTo(asTensorView(*psis[i]));
      }
    }

    /// Get U_i^H*op*V_i for the oblique projector V*inv(U^H*op*V)*U^H*op
    /// \param proj: Chimera projector
    /// \param index: index of the U and V basis

    DComplex getLambda(const ChimeraProjector& proj, unsigned int index)
    {
      DComplex r;
      if (!proj.op.op)
      {
	proj.chroma_proj->lambda(index, r);
      }
      else
      {
	auto l = proj.op.lambdas[index];
	r.elem().elem().elem() = QDP::RComplex<double>(std::real(l), std::imag(l));
      }
      return r;
    }
  }
}

#endif // BUILD_SB
