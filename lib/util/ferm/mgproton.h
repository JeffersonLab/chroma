// -*- C++ -*-
/*! \file                                                                    
 * \brief Multigrid prototype next
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_MGPROTON__
#define __INCLUDE_MGPROTON__

#include "chromabase.h"
#include "meas/hadron/greedy_coloring.h"
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

    /// Representation of an operator, function of type tensor -> tensor where the input and the
    /// output tensors have the same dimensions

    template <std::size_t NOp, typename COMPLEX>
    using OperatorFun =
      std::function<void(const Tensor<NOp + 1, COMPLEX>&, Tensor<NOp + 1, COMPLEX>)>;

    /// Representation of an operator together with a map to convert from domain labels (columns) to
    /// image labels (rows)

    template <std::size_t NOp, typename COMPLEX>
    struct Operator {
      /// Function that the operators applies
      OperatorFun<NOp, COMPLEX> fop;
      /// Example tensor for the input tensor (domain)
      Tensor<NOp, COMPLEX> d;
      /// Example tensor for the output tensor (image)
      Tensor<NOp, COMPLEX> i;
      /// Function to apply when conjugate transposed (optional)
      OperatorFun<NOp, COMPLEX> fop_tconj;
      /// Labels that distinguish different operator instantiations
      std::string order_t;
      /// Distance of the farthest direct neighbor in each non-blocking (sparse) direction
      unsigned int max_dist_neighbors;

      /// Constant used in `Operator<NOp,COMPLEX>::max_dist_neighbors` to indicate that the operator is dense
      static const unsigned int DenseOperator = std::numeric_limits<unsigned int>::max();

      /// Empty constructor
      Operator()
      {
      }

      /// Constructor
      Operator(const OperatorFun<NOp, COMPLEX>& fop, Tensor<NOp, COMPLEX> d, Tensor<NOp, COMPLEX> i,
	       const OperatorFun<NOp, COMPLEX>& fop_tconj = nullptr,
	       const std::string& order_t = "", unsigned int max_dist_neighbors = 1)
	: fop(fop),
	  d(d),
	  i(i),
	  fop_tconj(fop_tconj),
	  order_t(order_t),
	  max_dist_neighbors(max_dist_neighbors)
      {
      }

      /// Return the transpose conjugate of the operator
      Operator<NOp, COMPLEX> tconj() const
      {
	if (!fop_tconj)
	  throw std::runtime_error("Operator does not have conjugate transpose form");
	return {fop_tconj, i, d, fop, order_t, max_dist_neighbors};
      }

      /// Apply the operator
      template <std::size_t N, typename T>
      Tensor<N, T> operator()(const Tensor<N, T>& t) const
      {
	// The `t` labels that are not in `d` are the column labels
	std::string cols = detail::union_dimensions(t.order, "", d.order); // t.order - d.order

	auto x = t.template collapse_dimensions<NOp + 1>(cols, 'n').template make_sure<COMPLEX>();
	auto y = i.template like_this<NOp + 1>("%n", '%', "", {{'n', x.kvdim()['n']}});
	fop(x, y);
	return y.template split_dimension<N>('n', cols, t.kvdim()).template make_sure<T>();
      }

      /// Apply the operator
      template <std::size_t N, typename T>
      void operator()(const Tensor<N, T>& x, Tensor<N, T> y) const
      {
	// The `t` labels that are not in `d` are the column labels
	std::string cols = detail::union_dimensions(x.order, "", d.order); // t.order - d.order

	auto x0 = x.template collapse_dimensions<NOp + 1>(cols, 'n').template make_sure<COMPLEX>();
	auto y0 = y.template collapse_dimensions<NOp + 1>(cols, 'n').template make_sure<COMPLEX>();
	fop(x0, y0);
	if (y.data != y0.data)
	  y0.copyTo(y);
      }
    };

    namespace detail
    {
      /// Return the inverse map
      /// \param map: from domain labels to image labels

      inline remap reverse(const remap& map)
      {
	remap o;
	for (const auto& it : map)
	  o.insert({it.second, it.first});
	return o;
      }

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

      /// Return a string version of the number in scientific notation
      /// \param v: number to convert
      /// \param prec: number of digits to print

      inline std::string tostr(double v, unsigned int prec = 2)
      {
	std::stringstream ss;
	ss << std::scientific << std::setprecision(prec) << v;
	return ss.str();
      }
    }

    /// Orthonormalize W against V, W_out <- (I-V*V')*W_in*R, with R such that W_out'*W_out = I,
    /// for each combination of values of the order_t dimensions, or once if it is empty.
    /// \param V: orthogonal matrix to orthogonalize against
    /// \param W: matrix to orthogonalize
    /// \param order_t: dimension labels that do not participate in the orthogonalization
    /// \param order_rows: dimension labels that are the rows of the matrices V and W
    /// \param order_cols: dimension labels that are the columns of the matrices V and W

    template <std::size_t Nrows, std::size_t Ncols, std::size_t NV, typename COMPLEX,
	      std::size_t NW>
    void ortho(Tensor<NV, COMPLEX> V, Tensor<NW, COMPLEX> W, const std::string& order_t,
	       const std::string& order_rows, const std::string& order_cols,
	       unsigned int max_its = 3)
    {
      // Check Nrows
      if (order_rows.size() != Nrows)
	throw std::runtime_error("ortho: invalid template argument `Nrows`");
      if (order_cols.size() != Ncols)
	throw std::runtime_error("ortho: invalid template argument `Ncols`");

      // Check that V and W are compatible, excepting the column dimensions
      if (V)
	detail::check_compatible(V, W, none, order_cols);

      // Find the ordering
      std::string Wcorder = detail::remove_dimensions(W.order, order_t + order_rows);
      std::string Vcorder =
	V ? detail::remove_dimensions(V.order, order_t + order_rows) : std::string{};

      // Create an alternative view of W with different labels for the column
      remap Wac = detail::getNewLabels(Wcorder, (V ? V.order : std::string{}) + W.order);
      std::string Wacorder = detail::update_order(Wcorder, Wac);
      auto Wa = W.rename_dims(Wac);

      constexpr std::size_t Nt = NW - Nrows - Ncols;
      for (unsigned int i = 0; i < max_its; ++i)
      {
	// W = W - V*(V'*W)
	if (V)
	  contract(V.scale(-1), contract<NV + NW - Nrows * 2>(V.conj(), W, order_rows), Vcorder,
		   AddTo, W);

	// W = W/chol(Wa'*W), where Wa'*W has dimensions (rows,cols)=(Wacorder,Wcorder)
	cholInv<NW, Nt + 2 * Ncols, NW, COMPLEX>(contract<Nt + 2 * Ncols>(Wa.conj(), W, order_rows),
						 Wacorder, Wcorder, Wa, Wacorder, CopyTo, W);
      }
    }

    /// Orthonormalize W, W_out <- W_in*R, with R such that W_out'*W_out = I,
    /// for each combination of values of the order_t dimensions, or once if it is empty.
    /// \param W: matrix to orthogonalize
    /// \param order_t: dimension labels that do not participate in the orthogonalization
    /// \param order_rows: dimension labels that are the rows of the matrices V and W
    /// \param order_cols: dimension labels that are the columns of the matrices V and W

    template <std::size_t Nrows, std::size_t Ncols, std::size_t NW, typename COMPLEX>
    void ortho(Tensor<NW, COMPLEX> W, const std::string& order_t, const std::string& order_rows,
	       const std::string& order_cols, unsigned int max_its = 3)
    {
      ortho<Nrows, Ncols, NW, COMPLEX, NW>(Tensor<NW, COMPLEX>{}, W, order_t, order_rows,
					   order_cols, max_its);
    }

    /// Solve iteratively op * y = x using FGMRES
    /// \param op: problem matrix
    /// \param prec: preconditioner
    /// \param x: input right-hand-sides
    /// \param y: the solution vectors
    /// \param max_basis_size: maximum rank of the search subspace basis per t
    /// \param tol: maximum tolerance
    /// \param max_its: maximum number of iterations
    /// \param error_if_not_converged: throw an error if the tolerance was not satisfied
    /// \param passing_initial_guess: whether `y` contains a solution guess
    /// \param verb: verbosity level
    /// \param prefix: prefix printed before every line

    template <std::size_t NOp, typename COMPLEX>
    void fgmres(const Operator<NOp, COMPLEX>& op, Maybe<Operator<NOp, COMPLEX>> prec,
		const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX>& y,
		unsigned int max_basis_size, double tol, unsigned int max_its = 0,
		bool error_if_not_converged = true, bool passing_initial_guess = false,
		Verbosity verb = NoOutput, std::string prefix = "")
    {
      detail::log(1, prefix + " starting fgmres");

      const bool do_ortho = false;

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

      // Get an unused label for the search subspace columns
      char Vc = detail::get_free_label(x.order);
      char Vac = detail::get_free_label(x.order + std::string(1, Vc));
      std::string order_cols = detail::remove_dimensions(x.order, op.i.order);
      std::string order_rows = detail::remove_dimensions(op.d.order, op.order_t);
      std::size_t num_cols = x.volume(order_cols);

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
	r = y.like_this();
	r.set_zero();
	y.set_zero();
      }
      x.scale(-1).addTo(r);
      auto normr0 = norm<1>(r, op.order_t + order_cols); // type std::vector<real of T>
      if (max(normr0) == 0)
	return;

      // Allocate the search subspace U (onto the left singular space), and
      // Z (onto the right singular space)
      auto U = r.template like_this<NOp + 2>(std::string("%") + std::string(1, Vc), '%', "",
					     {{Vc, max_basis_size + 1}});

      // Allocate the search subspace Z (onto the right singular space); when no preconditioning
      // then Z = U

      auto Z = U.kvslice_from_size({}, {{Vc, max_basis_size}}); // set Z an alias of U
      if (prec.hasSome())
	Z = Z.like_this(); // create new space when using preconditioning

      auto normr = normr0.clone();
      unsigned int it = 0;
      double max_tol = HUGE_VAL;
      for (it = 0; it < max_its; ++it)
      {
	// U(:,0) = r;
	r.copyTo(U.kvslice_from_size({}, {{Vc, 1}}));

	// Expand the search subspace from residual
	for (unsigned int i = 0; i < max_basis_size; ++i)
	{
	  // Z(:,i) = prec * U(:,i)
	  if (prec.hasSome())
	  {
	    prec.getSome()(U.kvslice_from_size({{Vc, i}}, {{Vc, 1}}),
			   Z.kvslice_from_size({{Vc, i}}, {{Vc, 1}}));
	    nprecs += num_cols;
	  }

	  // U(:,i+1) = op * Z(:,i)
	  op(Z.kvslice_from_size({{Vc, i}}, {{Vc, 1}}),
	     U.kvslice_from_size({{Vc, i + 1}}, {{Vc, 1}}));
	  nops += num_cols;
	}

	// Orthogonalize U and put it into W: W = orth(U(:,2:end))
	// NOTE: for small max_basis_size and assuming r is far from a left singular vector,
	//       a light or none orthogonalization should be enough
	auto Up = U.kvslice_from_size({{Vc, 1}}, {{Vc, max_basis_size}});
	auto Uo = Up.rename_dims({{Vc, Vac}});
	if (do_ortho)
	{
	  //Uo = Uo.clone();
	  //ortho<NOp, NOp + 1, COMPLEX>(Uo, op.order_t + order_cols, order_rows, std::string(1, Vac),
	  //			       1 /* one iteration should be enough */);
	}

	// Restrict to Uo: [x_rt H_rt] = Uo'*U = Up'*[r Up]
	auto x_rt = contract<2>(Uo.conj(), r, order_rows);
	auto H_rt = contract<3>(Uo.conj(), Up, order_rows);

	// Solve the projected problem: y_rt = (Uo'*U(:2:end))\(Uo'*r);
	auto y_rt = solve<2>(H_rt, std::string(1, Vac), std::string(1, Vc),
			     x_rt.rename_dims({{Vac, Vc}}), std::string(1, Vc))
		      .rename_dims({{Vac, Vc}})
		      .make_sure(none, none, OnEveryoneReplicated);

	// Update solution: y += -Z*y_rt
	contract(Z.scale(-1), y_rt, std::string(1, Vc), AddTo, y);

	// Update residual: r += -U(2:end)*y_rt
	contract(Up.scale(-1), y_rt, std::string(1, Vc), AddTo, r);

	// Compute the norm
	auto normr = norm<1>(r, op.order_t + order_cols);

	// Get the worse tolerance
	max_tol = 0;
	for (int i = 0, vol = normr.volume(); i < vol; ++i)
	  max_tol = std::max(max_tol, normr.get({{i}}) / normr0.get({{i}}));

	// Report iteration
	if (verb >= Detailed)
	  QDPIO::cout << prefix << " MGPROTON FGMRES iteration it.: " << it
		      << " max rel. residual: " << detail::tostr(max_tol, 2) << std::endl;

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
	max_tol = 0;
	for (int i = 0, vol = normr.volume(); i < vol; ++i)
	  max_tol = std::max(max_tol, normr.get({{i}}) / normr0.get({{i}}));
	if (tol > 0 && max_tol > tol)
	  throw std::runtime_error("fmgres didn't converged and you ask for checking the error");
      }

      // Report iteration
      if (verb >= JustSummary)
	QDPIO::cout << prefix << " MGPROTON FGMRES summary its.: " << (it + 1)
		    << " max rel. residual: " << detail::tostr(max_tol, 2) << " matvecs: " << nops
		    << " precs: " << nprecs << std::endl;
    }

    template <std::size_t NOp, typename COMPLEX>
    Operator<NOp, COMPLEX> getSolver(const Operator<NOp, COMPLEX>& op, const Options& ops);

    namespace detail
    {
      /// Function that gets an operator and options and returns an operator that approximates the inverse
      template <std::size_t NOp, typename COMPLEX>
      using Solver =
	std::function<Operator<NOp, COMPLEX>(const Operator<NOp, COMPLEX>& op, const Options& ops)>;

      /// Returns an approximate inverse operator given an operator, options, and a list of solvers
      /// \param solvers: map of solvers
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getSolver(const std::map<std::string, Solver<NOp, COMPLEX>>& solvers,
				       const Operator<NOp, COMPLEX>& op, const Options& ops)
      {
	std::string type = getOption<std::string>(ops, "type");
	if (type.size() == 0)
	  ops.getValue("type").throw_error(
	    "Error constructing the solver: invalid `type', it's empty");
	if (solvers.count(type) == 0)
	{
	  std::string supported_tags;
	  for (const auto& it : solvers)
	    supported_tags += it.first + " ";
	  ops.getValue("type").throw_error(
	    "Error constructing the solver: unsupported `type' value `" + type +
	    "'; supported values: " + supported_tags);
	}

	return solvers.at(type)(op, ops);
      }

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
      Operator<NOp, COMPLEX> getFGMRESSolver(Operator<NOp, COMPLEX> op, const Options& ops)
      {
	// Get preconditioner
	Maybe<Operator<NOp, COMPLEX>> prec = none;
	Maybe<const Options&> precOps = getOptionsMaybe(ops, "prec");
	if (precOps)
	  prec = Maybe<Operator<NOp, COMPLEX>>(getSolver(op, precOps.getSome()));

	// Get the remainder options
	unsigned int max_basis_size = getOption<unsigned int>(ops, "max_basis_size", 0);
	double tol = getOption<double>(ops, "tol", 0.0);
	unsigned int max_its = getOption<unsigned int>(ops, "max_its", 0);
	if (max_its == 0 && tol <= 0)
	  ops.throw_error("set either `tol` or `max_its`");
	bool error_if_not_converged = getOption<bool>(ops, "error_if_not_converged", true);
	unsigned int max_simultaneous_rhs = getOption<unsigned int>(ops, "max_simultaneous_rhs", 0);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);
	std::string prefix = getOption<std::string>(ops, "prefix", "");

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  foreachInChuncks(
		    x, y, max_simultaneous_rhs,
		    [=](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		      fgmres(op, prec, x, y, max_basis_size, tol, max_its, error_if_not_converged,
			     false /* no init guess */, verb, prefix);
		    },
		    'n');
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.DenseOperator};
      }

      /// Returns the \gamma_5 for a given number of spins
      /// \param ns: number of spins

      template <typename COMPLEX = Complex>
      Tensor<2, COMPLEX> getGamma5(int ns, DeviceHost dev = OnDefaultDevice)
      {
	if (ns == 1)
	{
	  Tensor<2, COMPLEX> r("ij", {1, 1}, OnHost, OnEveryoneReplicated);
	  r.set({{0, 0}}, COMPLEX{1});
	  return r.make_sure(none, dev);
	}
	else if (ns == Ns)
	{
	  return SB::Gamma(Ns * Ns - 1);
	}
	else if (ns == 2)
	{
	  Tensor<2, COMPLEX> r("ij", {2, 2}, OnHost, OnEveryoneReplicated);
	  r.set_zero();
	  r.set({{0, 0}}, COMPLEX{1});
	  r.set({{1, 1}}, COMPLEX{-1});
	  return r.make_sure(none, dev);
	}
	else
	{
	  throw std::runtime_error("Error in getGamma5: Unsupported spin number");
	}
      }

      inline std::map<char, int> update_dims(const std::map<char, int>& m, const remap& rm)
      {
	std::map<char, int> r;
	for (const auto& it : m)
	  r[rm.at(it.first)] = it.second;
	return r;
      }

      /// Return a sparse tensor with the content of the given operator
      /// \param op: operator to extract the nonzeros from
      /// \param power: maximum distance to recover the nonzeros:
      ///               0, block diagonal; 1: near-neighbors...
      /// \return: a pair of a sparse tensor and a remap; the sparse tensor has the same image
      ///          labels as the given operator and domain labels are indicated by the returned remap.

      template <std::size_t NOp, typename COMPLEX>
      std::pair<SpTensor<NOp, NOp, COMPLEX>, remap>
      cloneOperatorToSpTensor(const Operator<NOp, COMPLEX>& op, unsigned int power = 1)
      {
	log(1, "starting cloneOperatorToSpTensor");

	// Unsupported explicitly colorized operators
	if (op.d.kvdim().count('X') == 0)
	  throw std::runtime_error(
	    "cloneOperatorToSpTensor: unsupported not explicitly colored operators");

	// Unsupported other powers than zero or one
	if (power != 0 && power != op.max_dist_neighbors)
	  throw std::runtime_error("cloneOperatorToSpTensor: unsupported value for `power`: "
				   "either zero or `max_dist_neighbors`");

	// Unsupported explicitly colorized operators
	if (op.d.kvdim().count('X') == 0)
	  throw std::runtime_error("cloneOperator: unsupported not explicitly colored operators");

	// If the operator is empty, just return itself
	if (op.d.volume() == 0 || op.i.volume() == 0)
	  return {{}, {}};

	// TODO: add optimizations for multiple operators
	if (op.order_t.size() > 0)
	  throw std::runtime_error("Not implemented");

	// Create the ordering for the domain and the image, moving the lattice dimensions as the slowest indices
	// NOTE: assuming that x,y,z,t are the only non-blocking dimensions
	std::set<char> lattice_dims{'x', 'y', 'z', 't', 'X'};
	remap rd = getNewLabels(op.d.order, op.i.order + "u~");
	auto d = op.d.reorder("%xyztX", '%').rename_dims(rd);
	auto i = op.i.reorder("%xyztX", '%');
	auto dim = i.kvdim();

	// Get the blocking for the domain and the image
	std::map<char, int> blkd, blki;
	remap rev_rd = reverse(rd);
	for (const auto& it : d.kvdim())
	  blkd[it.first] = (lattice_dims.count(rev_rd.at(it.first)) == 0 ? it.second : 1);
	for (const auto& it : i.kvdim())
	  blki[it.first] = (lattice_dims.count(it.first) == 0 ? it.second : 1);

	// Construct the probing vectors, which they have as the rows the domain labels and as
	// columns the domain blocking dimensions

	constexpr int Nblk = NOp - Nd - 1;
	auto t_bl = d.like_this(none, blkd, OnHost, OnEveryoneReplicated);
	t_bl.set_zero();
	auto t_bl_rows = i.template like_this<Nblk>("%", '%', "xyztX");
	t_bl_rows.set_zero();
	auto t_blbl =
	  contract<NOp + Nblk>(t_bl_rows, t_bl, "").make_sure(none, none, OnEveryoneReplicated);
	assert(!t_blbl.isSubtensor());
	COMPLEX* t_blbl_data = t_blbl.data.get();
	for (std::size_t i = 0, vol = t_bl.volume(); i < vol; ++i)
	  t_blbl_data[i * vol + i] = COMPLEX{1};

	std::map<char, int> nonblkd;
	for (const auto& it : dim)
	  nonblkd[it.first] = (lattice_dims.count(it.first) == 0 ? 1 : it.second);

	// Compute the coloring
	Coor<Nd> dims{{dim['x'] * dim['X'], dim['y'], dim['z'], dim['t']}};
	Coloring coloring{0, // zero displacement in coloring
			  power == 0 ? op.max_dist_neighbors
				     : op.max_dist_neighbors * 2 + 1, // k-distance coloring
			  dims};
	unsigned int num_colors = coloring.numColors();

	// Get the number of neighbors
	std::vector<Coor<Nd>> neighbors = Coloring::all_neighbors(power, dims);

	// Create the sparse tensor
	int maxX = dim['X'];
	auto d_sop =
	  (power == 0 ? d
		      : d.extend_support({{rd.at('x'), (op.max_dist_neighbors + maxX - 1) / maxX},
					  {rd.at('y'), op.max_dist_neighbors},
					  {rd.at('z'), op.max_dist_neighbors},
					  {rd.at('t'), op.max_dist_neighbors}}));
	SpTensor<NOp, NOp, COMPLEX> sop{d_sop, i, Nblk, Nblk, (unsigned int)neighbors.size()};

	int maxY = std::min(2, dims[1]);
	int maxZ = std::min(2, dims[2]);
	int maxT = std::min(2, dims[3]);
	for (unsigned int color = 0; color < num_colors; ++color)
	{
	  // Extracting the proving vector
	  auto t_l = fillLatticeField<5,
				      COMPLEX>("xyztX", nonblkd, OnDefaultDevice, [&](Coor<4> c) {
	    return coloring.getColor({{c[0], c[1], c[2], c[3]}}) == color ? COMPLEX{1} : COMPLEX{0};
	  });

	  // Contracting the proving vector with the blocking components
	  auto probs = contract<NOp * 2>(t_l, t_blbl, "");

	  // Compute the matvecs
	  auto mv = op(std::move(probs));

	  // Construct an indicator tensor where all blocking dimensions but only the nodes colored `color` are copied
	  auto ones = t_blbl.template like_this<NOp * 2 - 5, float>();
	  ones.set(1);
	  auto sel = contract<NOp * 2>(t_l.template transformWithCPUFun<float>(
					 [](const COMPLEX& t) { return (float)std::real(t); }),
				       ones, "");
	  auto mv_mask = mv.create_mask();
	  auto data_mask = sop.data.create_mask();

	  // Split mv and sel into subsets where the neighbors are at the same position
	  auto sel_eo = sel.split_dimension('y', "Yy", maxY)
			  .split_dimension('z', "Zz", maxZ)
			  .split_dimension('t', "Tt", maxT);
	  auto mv_mask_eo = mv_mask.split_dimension('y', "Yy", maxY)
			      .split_dimension('z', "Zz", maxZ)
			      .split_dimension('t', "Tt", maxT);
	  auto mv_eo = mv.split_dimension('y', "Yy", maxY)
			 .split_dimension('z', "Zz", maxZ)
			 .split_dimension('t', "Tt", maxT);
	  auto data_eo = sop.data.split_dimension('y', "Yy", maxY)
			   .split_dimension('z', "Zz", maxZ)
			   .split_dimension('t', "Tt", maxT);
	  auto data_mask_eo = data_mask.split_dimension('y', "Yy", maxY)
				.split_dimension('z', "Zz", maxZ)
				.split_dimension('t', "Tt", maxT);

	  // Populate the nonzeros by copying pieces from `mv` into sop.data. We want to copy only the
	  // nonzeros in `mv`, which are `neighbors` away from the nonzeros of `probs`. Assume `probs` having a
	  // single nonzero at the natural coordinates (X,Y,Z,T), and `mv` has some nonzero at direction `dir`, that is,
	  // at the natural coordinate (x+dirx,Y+diry,Z+dirz,T+dirt), which has even-odd coordinates (X,x,y,z,t) =
	  // ((X+Y+Z+T+dirx+diry+dirz+dirt)%maxX,(X+dirx)/maxX,Y+diry,Z+dirz,T+dirt)

	  std::map<char, int> XYZTsize{{'X', 1}, {'Y', 1}, {'Z', 1}, {'T', 1}};
	  for (int mu = 0; mu < neighbors.size(); ++mu)
	  {
	    const auto& dir = neighbors[mu];
	    int sumdir = std::accumulate(dir.begin(), dir.end(), int{0});
	    if ((dir[0] + dims[0]) % maxX == 0)
	    {
	      // Nat coor (0,diry,dirz,dirt) to even-odd coordinate
	      std::map<char, int> to{
		{'X', sumdir}, {'x', dir[0]}, {'y', dir[1]}, {'z', dir[2]}, {'t', dir[3]}};
	      auto mv_mask_mu = mv_mask.kvslice_from_size(to);
	      auto data_mask_mu =
		data_mask.kvslice_from_size(to).kvslice_from_size({{'u', mu}}, {{'u', 1}});
	      sel.copyTo(mv_mask_mu);
	      sel.copyTo(data_mask_mu);
	      mv.kvslice_from_size(to).copyToWithMask(
		sop.data.kvslice_from_size(to).kvslice_from_size({{'u', mu}}, {{'u', 1}}),
		mv_mask_mu, data_mask_mu);
	    }
	    else
	    {
	      for (int T = 0; T < maxT; ++T)
	      {
		for (int Z = 0; Z < maxZ; ++Z)
		{
		  for (int Y = 0; Y < maxY; ++Y)
		  {
		    for (int X = 0; X < maxX; ++X)
		    {
		      // Nat coor (X,Y,Z,T) to even-odd coordinate
		      std::map<char, int> XYZTfrom{
			{'X', X + Y + Z + T}, {'Y', Y}, {'Z', Z}, {'T', T}};
		      // Nat coor (x+dirx,Y+diry,Z+dirz,T+dirt) to even-odd coordinate
		      std::map<char, int> to{{'X', X + Y + Z + T + sumdir},
					     {'x', (dir[0] + dims[0] + X) / maxX},
					     {'Y', Y + dir[1]},
					     {'y', (dir[1] + dims[1] + Y) / maxY},
					     {'Z', Z + dir[2]},
					     {'z', (dir[2] + dims[2] + Z) / maxZ},
					     {'T', T + dir[3]},
					     {'t', (dir[3] + dims[3] + T) / maxT}};
		      auto sel_eo_mu = sel_eo.kvslice_from_size(XYZTfrom, XYZTsize);
		      auto mv_mask_eo_mu = mv_mask_eo.kvslice_from_size(to, XYZTsize);
		      auto data_mask_eo_mu = data_mask_eo.kvslice_from_size(to, XYZTsize)
					       .kvslice_from_size({{'u', mu}}, {{'u', 1}});
		      sel_eo_mu.copyTo(mv_mask_eo_mu);
		      sel_eo_mu.copyTo(data_mask_eo_mu);
		      mv_eo.kvslice_from_size(to, XYZTsize)
			.copyToWithMask(data_eo.kvslice_from_size(to, XYZTsize)
					  .kvslice_from_size({{'u', mu}}, {{'u', 1}}),
					mv_mask_eo_mu, data_mask_eo_mu);
		    } // X
		  }   // Y
		}     // Z
	      }	      // T
	    }
	  } // mu
	}   // color

	// Populate the coordinate of the columns, that is, to give the domain coordinates of first nonzero in each
	// BSR nonzero block. Assume that we are processing nonzeros block for the image coordinate `c` on the
	// direction `dir`, that is, the domain coordinates will be (cx-dirx,cy-diry,cz-dirz,dt-dirt) in natural
	// coordinates. But we get the image coordinate `c` in even-odd coordinate, (cX,cx,cy,cz,ct), which has the
	// following natural coordinates (cx*2+(cX+cy+cz+ct)%2,cy,cz,ct). After subtracting the direction we get the
	// natural coordinates (cx*2+(cX+cy+cz+ct)%2-dirx,cy-diry,cz-dirz,ct-dirt), which corresponds to the following
	// even-odd coordinates ((cX-dirx-diry-dirz-dirt)%2,(cx*2+(cX+cy+cz+ct)%2-dirx)/2,cy-diry,cz-dirz,ct-dirt).

	sop.jj.fillCpuFunCoor([&](const Coor<NOp + 2>& c) {
	  // c has order '~u%xyztX' where xyztX were remapped by ri
	  int domi = c[0];	  // the domain label to evaluate, label ~
	  int mu = c[1];	  // the direction, label u
	  int base = c[domi + 2]; // the image coordinate value for label `domi`

	  // Do nothing for a blocking direction
	  if (domi < Nblk)
	    return 0;

	  const auto& dir = neighbors[mu];

	  // For labels X and x
	  if (domi == Nblk || domi == Nblk + Nd)
	  {
	    int sumdir = std::accumulate(dir.begin(), dir.end(), int{0});
	    if (domi == Nblk + Nd)
	      return (base + sumdir + maxX * Nd) % maxX;
	    int sumyzt = std::accumulate(c.begin() + 2 + Nblk + 1, c.end() - 1, int{0});
	    const auto& cX = c[2 + Nblk + Nd];
	    return ((base * maxX + (cX + sumyzt) % maxX + dims[0] - dir[0]) / maxX) %
		   (dims[0] / maxX);
	  }

	  int latd = domi - Nblk;
	  return (base - dir[latd] + dims[latd]) % dims[latd];
	});

	// Construct the sparse operator
	sop.construct();

	// Return the sparse tensor and the remap from original operator to domain of the sparse tensor
	return {sop, rd};
      }

      /// Return a sparse tensor with the content of the given operator
      /// \param op: operator to extract the nonzeros from
      /// \param power: maximum distance to recover the nonzeros:
      ///               0, block diagonal; 1: near-neighbors...
      /// \return: a pair of a sparse tensor and a remap; the sparse tensor has the same image
      ///          labels as the given operator and domain labels are indicated by the returned remap.

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> cloneOperator(const Operator<NOp, COMPLEX>& op)
      {
	// If the operator is empty, just return itself
	if (op.d.volume() == 0 || op.i.volume() == 0)
	  return op;

	// Get a sparse tensor representation of the operator
	auto t = detail::cloneOperatorToSpTensor(op, op.max_dist_neighbors);
	auto sop = t.first;
	remap rd = t.second;

	// Construct the operator to return
	Operator<NOp, COMPLEX> rop{
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    auto x0 = x.reorder(sop.i.order + "%", '%');
	    auto y0 = (y.order == x0.order ? y : x0.like_this());
	    sop.contractWith(x0, rd, y0);
	    if (y.data != y0.data)
	      y0.copyTo(y);
	  },
	  sop.i,
	  sop.i,
	  nullptr,
	  op.order_t,
	  op.max_dist_neighbors};

	// Do a test
	auto x = op.d.template like_this<NOp + 1>("%n", '%', "", {{'n', 2}});
	urand(x, -1, 1);
	auto y = op(x);
	auto base_norm = norm<1>(y, "n");
	rop(x).scale(-1).addTo(y);
	auto error = norm<1>(y, "n");
	auto eps = std::sqrt(std::numeric_limits<typename real_type<COMPLEX>::type>::epsilon());
	for (int i = 0; i < base_norm.volume(); ++i)
	  if (error.get({{i}}) > eps * base_norm.get({{i}}))
	    throw std::runtime_error("cloneOperator: too much error on the cloned operator");

	return rop;
      }

      /// Return the block diagonal of an operator
      /// \param op: operator to extract the block diagonal
      /// \param block_labels: labels that compose the blocks
      /// \param m: map from the operator block labels to the new labels to represent the rows
      /// \return: tensor with the block labels as columns and the renamed labels as rows

      template <std::size_t N, std::size_t NOp, typename COMPLEX>
      Tensor<N, COMPLEX> getBlockDiag(const Operator<NOp, COMPLEX>& op,
				      const std::string& block_labels, const remap& m)
      {
	// Check
	// Get a sparse tensor representation of the operator
	auto t = cloneOperatorToSpTensor(op, 0 /* clone only the block diagonal */);
	auto sop = t.first;
	remap rd = t.second;

	// Create tensor to return
	std::string order = update_order(block_labels, m) + op.d.order;
	auto dims = op.d.kvdim();
	for (const auto& it : m)
	  dims[it.second] = dims[it.first];
	Tensor<N, COMPLEX> r = op.d.template like_this<N>(order, dims);

	// Copy the blocks
	remap m_blk;
	for (const auto& it : m)
	  m_blk[rd[it.first]] = it.second;
	sop.data.rename_dims(m_blk).copyTo(r);

	// Return tensor
	return r;
      }

      /// Returns the prolongator constructed from
      /// \param solvers: map of solvers
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getMGProlongator(const Operator<NOp, COMPLEX>& op,
					      unsigned int num_null_vecs,
					      const std::map<char, unsigned int>& blocking,
					      const Operator<NOp, COMPLEX>& null_solver)
      {
	detail::log(1, "starting getMGProlongator");

	// For now blocking on x should be divisible by X
	auto opdims = op.d.kvdim();
	int X = opdims.at('X');
	if (blocking.at('x') > 1 && blocking.at('x') % X != 0)
	  throw std::runtime_error("Unsupported blocking in x direction with isn't divisible by 2 "
				   "when using even-odd layout");

	// Create num_null_vecs random vectors and solve them
	auto b = op.d.template like_this<NOp + 1>("%n", '%', "", {{'n', num_null_vecs}});
	urand(b, -1, 1);
	auto nv = null_solver(std::move(b));

	// Do chirality splitting nv2 = [ nv * gpos, nv * gneg ]
	// TODO: experiment without chirality splitting
	int ns = opdims.at('s');
	if (ns != 2 && ns != Ns)
	  throw std::runtime_error("Error in getMGProlongator: Unsupported spin number");

	auto nv2 = nv.like_this(none, {{'n', num_null_vecs * 2}});
	auto g5 = getGamma5(ns, OnHost), g5pos = g5.clone(), g5neg = g5.clone();
	for (int i = 0; i < ns; ++i) // make diagonal entries of gpos all positive or zero
	  g5pos.set({{i, i}}, g5.get({{i, i}}) + COMPLEX{1});
	for (int i = 0; i < ns; ++i) // make diagonal entries of gneg all negative or zero
	  g5neg.set({{i, i}}, g5.get({{i, i}}) - COMPLEX{1});
	nv2.kvslice_from_size({}, {{'n', num_null_vecs}})
	  .contract(g5pos, {{'i', 's'}}, NotConjugate, nv, {{'s', 'j'}}, NotConjugate);
	nv2.kvslice_from_size({{'n', num_null_vecs}}, {{'n', num_null_vecs}})
	  .contract(g5neg, {{'i', 's'}}, NotConjugate, nv, {{'s', 'j'}}, NotConjugate);

	// Do the blocking
	// NOTE: Xx is transformed into WwXx, where W,X has size X, and w has size blocking(x)/X,
	//       and the remaining x has size x/blocking(x)
	int bx = blocking.at('x');
	int dimw = bx > 1 ? bx / X : 1;
	auto nv_blk = nv2.rename_dims({{'X', 'W'}})
			.template split_dimension<NOp + 1 + 3 - 1>(
			  'x', "wXx", {{'w', dimw}, {'X', 1}, {'x', opdims.at('x') / bx * X}})
			.split_dimension('y', "Yy", blocking.at('y'))
			.split_dimension('z', "Zz", blocking.at('z'))
			.split_dimension('t', "Tt", blocking.at('t'))
			.rename_dims({{'c', 'C'}, {'s', 'S'}})
			.split_dimension('n', "cs", num_null_vecs);

	// Do the orthogonalization on each block and chirality
	ortho<7, 1>(nv_blk, "Xxyzts", "WwYZTSC", "c");

	// Return the operator
	Tensor<NOp, COMPLEX> d = op.d.like_this(none, nv_blk.kvdim()), i = op.i;
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  auto y_blk =
		    y.rename_dims({{'X', 'W'}})
		      .template split_dimension<NOp + 1 + 3 - 1>(
			'x', "wXx", {{'w', dimw}, {'X', 1}, {'x', opdims.at('x') / bx * X}})
		      .split_dimension('y', "Yy", blocking.at('y'))
		      .split_dimension('z', "Zz", blocking.at('z'))
		      .split_dimension('t', "Tt", blocking.at('t'))
		      .rename_dims({{'c', 'C'}, {'s', 'S'}});
		  contract(nv_blk, x, "cs", CopyTo, y_blk);
		},
		d,
		i,
		[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  auto x_blk =
		    x.rename_dims({{'X', 'W'}})
		      .template split_dimension<NOp + 1 + 3 - 1>(
			'x', "wXx", {{'w', dimw}, {'X', 1}, {'x', opdims.at('x') / bx * X}})
		      .split_dimension('y', "Yy", blocking.at('y'))
		      .split_dimension('z', "Zz", blocking.at('z'))
		      .split_dimension('t', "Tt", blocking.at('t'))
		      .rename_dims({{'c', 'C'}, {'s', 'S'}});
		  contract<NOp + 1>(nv_blk.conj(), x_blk, "WwYZTSC", CopyTo, y);
		},
		op.order_t,
		op.DenseOperator};
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
      Operator<NOp, COMPLEX> getMGPrec(Operator<NOp, COMPLEX> op, const Options& ops)
      {
	// Get prolongator, V
	unsigned int num_null_vecs = getOption<unsigned int>(ops, "num_null_vecs");
	std::vector<unsigned int> blocking_v = getVectorOption<unsigned int>(ops, "blocking");
	if (blocking_v.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getMGPrec: the blocking should be a vector with four elements");
	std::map<char, unsigned int> blocking{
	  {'x', blocking_v[0]}, {'y', blocking_v[1]}, {'z', blocking_v[2]}, {'t', blocking_v[3]}};
	const Operator<NOp, COMPLEX> nullSolver =
	  getSolver(op, getOptions(ops, "solver_null_vecs"));
	const Operator<NOp, COMPLEX> V = getMGProlongator(op, num_null_vecs, blocking, nullSolver);

	// Compute the coarse operator, V' * op * V
	unsigned int create_coarse_max_rhs =
	  getOption<unsigned int>(ops, "create_coarse_max_rhs", 0);
	const Operator<NOp, COMPLEX> op_c = cloneOperator(Operator<NOp, COMPLEX>{
	  [&](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    foreachInChuncks(x, y, create_coarse_max_rhs,
			     [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
			       V.tconj()(op(V(x)), y);
			     });
	  },
	  V.d, V.d, nullptr, op.order_t, op.max_dist_neighbors});

	// Get the solver for the projector
	const Operator<NOp, COMPLEX> coarseSolver =
	  getSolver(op_c, getOptions(ops, "solver_coarse"));

	// Get the solver for the smoother
	const Operator<NOp, COMPLEX> opSolver = getSolver(op, getOptions(ops, "solver_smoother"));

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  // y0 = V*solver(V'*Op*V, V'x)
		  auto y0 = V(coarseSolver(V.tconj()(x)));

		  // x1 = x - op*y0
		  auto x1 = op(y0.scale(-1));
		  x.addTo(x1);

		  // y = y0 + solver(Op, x1)
		  opSolver(std::move(x1), y);
		  y0.addTo(y);
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.DenseOperator};
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
      /// (Op_ee-Op_eo*Op_oo^{-1}*Op_oe)^{-1} is approximate by an iterative solver. Note that
      /// the residual norm while solving A_ee=Op_ee-Op_eo*Op_oo^{-1}*Op_oe is the same as the original
      /// residual norm. To prove that, notice that the global residual will be zero for the odds,
      /// and that the global residual for a solution x'_e of A_ee is
      ///   r = L*A*R * R^{-1}*[x'_e x'_o] - b = [r_e; 0];
      /// therefore ||L^{-1}r|| = ||r|| because of the form of L, making
      ///   ||r|| = || A_e*x'_e - (b_e-Op_eo*Op_oo^{-1}*b_o) ||.
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getEvenOddPrec(Operator<NOp, COMPLEX> op, const Options& ops)
      {
	auto dims = op.d.kvdim();
	if (dims.count('X') == 0 || dims.at('X') != 2)
	  ops.throw_error(
	    "getEvenOddPrec: only supported on explicitly colored operators with two colors");

	// Get the block diagonal of the operator
	auto opDiag = getBlockDiag<NOp + 2>(op, "sc", {{'s', 'S'}, {'c', 'C'}});

	// Get an explicit form for A:=I-Op_eo*Op_oo^{-1}*Op_oe*Op_ee^{-1}
	unsigned int create_operator_max_rhs =
	  getOption<unsigned int>(ops, "create_operator_max_rhs", 0);
	const Operator<NOp, COMPLEX> opA = cloneOperator(Operator<NOp, COMPLEX>{
	  [&](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    foreachInChuncks(
	      x, y, create_operator_max_rhs,
	      [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
		// y0 = Op_ee^{-1} * x
		auto y0 = op.d.template like_this<NOp + 1>(none, {{'n', x.kvdim()['n']}});
		y0.set_zero();
		solve<NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		  opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}), "SC", "sc", x, "sc", CopyTo,
		  y0.kvslice_from_size({{'X', 0}}, {{'X', 1}}));

		// y1 = Op_oe * y0
		auto y1 = op(y0).kvslice_from_size({{'X', 1}}, {{'X', 1}});
		// y2 = Op_oo^{-1} * y1
		auto y2 = y0;
		y2.set_zero();
		solve<NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		  opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "SC", "sc", y1, "sc", CopyTo,
		  y2.kvslice_from_size({{'X', 1}}, {{'X', 1}}));

		// y3 = Op_eo * y2
		auto y3 = op(y2).kvslice_from_size({{'X', 0}}, {{'X', 1}});

		// y = x - y3
		x.copyTo(y);
		y3.scale(-1).addTo(y);
	      });
	  },
	  op.d.kvslice_from_size({{'X', 0}}, {{'X', 1}}),
	  op.d.kvslice_from_size({{'X', 0}}, {{'X', 1}}), nullptr, op.order_t,
	  op.max_dist_neighbors * 2});

	// Get solver on opA
	const Operator<NOp, COMPLEX> solver = getSolver(opA, getOptions(ops, "solver"));

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  // be = x_e - Op_eo*Op_oo^{-1}*x_o
		  auto be = solver.d.template like_this<NOp + 1>(none, {{'n', x.kvdim()['n']}});
		  x.kvslice_from_size({{'X', 0}}, {{'X', 1}}).copyTo(be);
		  auto ya = op.d.template like_this<NOp + 1>(none, {{'n', x.kvdim()['n']}});
		  ya.set_zero();
		  solve<NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		    opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "SC", "sc",
		    x.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "sc", CopyTo,
		    ya.kvslice_from_size({{'X', 1}}, {{'X', 1}}));
		  op(ya).kvslice_from_size({{'X', 0}}, {{'X', 1}}).scale(-1).addTo(be);

		  // Solve opA * ye0 = be
		  auto ye0 = solver(be);

		  // y_e = Op_ee^{-1} * ye0
		  y.set_zero();
		  solve<NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		    opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}), "SC", "sc", ye0, "sc", CopyTo,
		    y.kvslice_from_size({{'X', 0}}, {{'X', 1}}));

		  // y_o = Op_oo^{-1}*(-Op_oe*y_e + x_o)
		  op(y, ya);
		  auto yo0 = be;
		  x.kvslice_from_size({{'X', 1}}, {{'X', 1}}).copyTo(yo0);
		  ya.kvslice_from_size({{'X', 1}}, {{'X', 1}}).scale(-1).addTo(yo0);
		  solve<NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		    opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "SC", "sc", yo0, "sc", CopyTo,
		    y.kvslice_from_size({{'X', 1}}, {{'X', 1}}));
		},
		op.i,
		op.d,
		nullptr,
		op.order_t,
		op.DenseOperator};
      }
    }

    /// Returns an operator that approximate the inverse of a given operator
    /// \param op: operator to make the inverse of
    /// \param ops: options to select the solver from `solvers` and influence the solver construction

    template <std::size_t NOp, typename COMPLEX>
    Operator<NOp, COMPLEX> getSolver(const Operator<NOp, COMPLEX>& op, const Options& ops)
    {
      static const std::map<std::string, detail::Solver<NOp, COMPLEX>> solvers{
	{"fgmres", detail::getFGMRESSolver<NOp, COMPLEX>}, // flexible GMRES
	{"mg", detail::getMGPrec<NOp, COMPLEX>},	   // Multigrid
	{"eo", detail::getEvenOddPrec<NOp, COMPLEX>}	   // even-odd Schur preconditioner
      };

      return detail::getSolver(solvers, op, ops);
    }

    /// Return an Operator that wraps up a LinearOperator<LatticeFermion>
    inline Operator<Nd + 3, Complex> asOperatorView(const LinearOperator<LatticeFermion>& linOp)
    {
      LatticeFermion a;
      auto d = asTensorView(a).toComplex().make_eg();
      return {
	[&](const Tensor<Nd + 4, Complex>& x, Tensor<Nd + 4, Complex> y) {
	  LatticeFermion x0, y0;
	  unsigned int n = x.kvdim()['n'];
	  for (unsigned int i = 0; i < n; ++i)
	  {
	    x.kvslice_from_size({{'n', i}}, {{'n', 1}}).copyTo(asTensorView(x0));
	    y0 = zero;
	    linOp(y0, x0, PLUS /* I believe, it's ignored */);
	    asTensorView(y0).copyTo(y.kvslice_from_size({{'n', i}}, {{'n', 1}}));
	  }
	},
	d,	 // domain
	d,	 // image
	nullptr, // no conjugate
	"",	 // no order_t
	1	 // links with near-neighbors only
      };
    }

    //
    // High-level chroma operations
    //

    /// Either a Chroma solver or a superb solver
    struct ChimeraSolver {
      /// Action
      Handle<
	FermionAction<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>
	S;

      /// State
      Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>
	state;

      /// Chroma solver (optional)
      Handle<SystemSolver<LatticeFermion>> PP;

      /// Operator on scxyztX (optional)
      Maybe<Operator<Nd + 3, Complex>> op;

      /// Constructor
      /// \param fermAction: XML for the fermion action
      /// \param invParam: XML for the quark propagator
      /// \param u: gauge fields

      ChimeraSolver(const GroupXML_t& fermAction, const GroupXML_t& invParam,
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
	  // Parse XML with the inverter options
	  std::shared_ptr<Options> ops = getOptionsFromXML(broadcast(invParam.xml));

	  // Clone the matvec
	  LinearOperator<LatticeFermion>* fLinOp = S->genLinOp(state);
	  Operator<Nd + 3, Complex> linOp = detail::cloneOperator(asOperatorView(*fLinOp));

	  // Construct the solver
	  op = Maybe<Operator<Nd + 3, Complex>>{getSolver(linOp, getOptions(*ops, "InvertParam"))};
	}
	else
	{
	  PP = S->qprop(state, invParam);
	}
      }
    };

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
      Tensor<Nd + 5, COMPLEX_OUT> doInversion(const Operator<Nd + 3, COMPLEX_OUT>& op,
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

	// Create tensors with full support on the lattice
	int max_step = std::max(num_vecs, max_rhs);
	auto aux = chi.template like_this<Nd + 4>(
	  "csxyztXn", {{'n', max_step}, {'t', Layout::lattSize()[3]}, {'s', Ns}});

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

    template <typename COMPLEX_CHI, typename COMPLEX_OUT>
    Tensor<Nd + 5, COMPLEX_OUT>
    doInversion(const ChimeraSolver& sol, const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
		int first_tslice_out, int n_tslice_out, const std::vector<int>& spin_sources,
		int max_rhs, const std::string& order_out = "cSxyztXns")
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

      Tensor<Nd + 5, COMPLEX_OUT> r;
      if (sol.op.hasSome())
	r = detail::doInversion<COMPLEX_CHI, COMPLEX_OUT>(sol.op.getSome(), chi, t_source,
							  first_tslice_out, n_tslice_out,
							  spin_sources, max_rhs, order_out);
      else
	r = detail::doInversion<COMPLEX_CHI, COMPLEX_OUT>(
	  *sol.PP, chi, t_source, first_tslice_out, n_tslice_out, spin_sources, max_rhs, order_out);

      snarss1.stop();
      QDPIO::cout << "Time to compute inversions for " << spin_sources.size()
		  << " spin sources and " << num_vecs
		  << " colorvecs : " << snarss1.getTimeInSeconds() << " secs" << std::endl;

      return r;
    }
  }
}

#endif // BUILD_SB

#endif // __INCLUDE_MGPROTON__
