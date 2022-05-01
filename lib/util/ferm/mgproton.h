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
#include "meas/hadron/greedy_coloring.h"

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
						      {"JustSummary", JustSummary},
						      {"Detailed", Detailed},
						      {"VeryDetailed", VeryDetailed},
						      {"false", NoOutput}};
      return m;
    }

    /// Representation of an operator, function of type tensor -> tensor where the input and the
    /// output tensors have the same dimensions

    template <std::size_t NOp, typename COMPLEX>
    using OperatorFun = std::function<Tensor<NOp + 1, COMPLEX>(const Tensor<NOp + 1, COMPLEX>&)>;

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
      /// Map to rename image labels into domain labels (optional)
      remap map;

      /// Constructor
      Operator(const OperatorFun<NOp, COMPLEX>& fop, Tensor<NOp, COMPLEX> d, Tensor<NOp, COMPLEX> i,
	       const OperatorFun<NOp, COMPLEX>& fop_tconj = nullptr, const remap& map = {})
	: fop(fop), d(d), i(i), fop_tconj(fop_tconj), map(map)
      {
      }

      /// Return the transpose conjugate of the operator
      Operator<NOp, COMPLEX> tconj() const
      {
	if (!fop_tconj)
	  throw std::runtime_error("Operator does not have conjugate transpose form");
	return {fop_tconj, i, d, fop, reverse(map)};
      }

      /// Apply the operator
      template <std::size_t N, typename T>
      Tensor<N, T> operator()(const Tensor<N, T>& t)
      {
	// The `t` labels that are not in `d` are the column labels
	std::string cols = detail::union_dimensions(t.order, "", d.order); // t.order - d.order

	return fop(t.colapse_dimensions<NOp + 1>(cols, "n"))
	  .rename_dims(map)
	  .split_dimensions<N>("n", cols, t.kvdim());
      }
    };

    namespace detail
    {
      /// Return an order with all domain labels
      /// \param map: from domain labels to image labels
      /// \param exclude: do not return the given labels

      std::string getDomain(const remap& map, const std::string& exclude = "")
      {
	std::string o;
	for (const auto& it : map)
	  if (std::find(exclude.begin(), exclude.end(), it.first) == exclude.end())
	    o += it.first;
	return o;
      }

      /// Return an order with all domain labels
      /// \param map: from domain labels to image labels
      /// \param exclude: do not return the given labels

      std::string getImage(const remap& map, const std::string& exclude = "")
      {
	std::string o;
	for (const auto& it : map)
	  if (std::find(exclude.begin(), exclude.end(), it.second) == exclude.end())
	    o += it.second;
	return o;
      }

      /// Return the inverse map
      /// \param map: from domain labels to image labels

      remap reverse(const remap& map)
      {
	remap o;
	for (const auto& it : map)
	  o.insert({it.second, it.first});
	return o;
      }
    }

    namespace detail
    {

    }

    /// Orthonormalize W against V, W_out <- (I-V*V')*W_in*R, with R such that W_out'*W_out = I,
    /// for each combination of values of the order_t dimensions, or once if it is empty.
    /// \param V: orthogonal matrix to orthogonalize against
    /// \param W: matrix to orthogonalize
    /// \param order_t: dimension labels that do not participate in the orthogonalization
    /// \param order_rows: dimension labels that are the rows of the matrices V and W
    /// \param order_cols: dimension labels that are the columns of the matrices V and W

    template <std::size_t NV, std::size_t NW, typename COMPLEX>
    void ortho(Maybe<Tensor<NV, COMPLEX>> V, Tensor<NW, COMPLEX>& W, Maybe<std::string> order_t,
	       Maybe<std::string> order_rows, Maybe<std::string> order_cols,
	       unsigned int max_its = 3)
    {
      // Check that V and W are compatible, excepting the column dimensions
      if (V)
	detail::check_compatible(V, W, none, order_cols);

      // Find the ordering
      std::string torder, rorder, Wcorder, Vcorder;
      detail::get_mat_ordering(W, order_t, order_rows, order_cols, torder, rorder, Wcorder);
      if (V)
	detail::get_mat_ordering(V, order_t, order_rows, order_cols, torder, rorder, Vcorder);

      // Create an alternative view of W with different labels for the column
      remap Wac = getNewLabels(Wcorder, toString(V.kvdim()) + toString(W.kvdim()));
      std::string Wacorder = detail::update_order(Wcorder, Wac);
      auto Wa = W.rename_dims(Wac);

      for (unsigned int i = 0; i < max_its; ++i)
      {
	// W = W - V*(V'*W)
	if (V)
	  contract(V.scale(-1), contract(V.conj(), W, rorder), Vcorder, AddTo, W);

	// W = W/chol(Wa'*W), where Wa'*W has dimensions (rows,cols)=(Wacorder,Wcorder)
	cholInv(W, contract(Wa.conj(), W, rorder), torder, Wacorder, Wcorder, CopyTo, W);
      }
    }

    /// Solve iteratively op * y = x using FGMRES
    /// \param op: problem matrix
    /// \param prec: preconditioner
    /// \param x: input right-hand-sides
    /// \param y: the solution vectors
    /// \param order_t: op and prec dimension labels that do not participate in the solver
    /// \param max_basis_size: maximum rank of the search subspace basis per t
    /// \param tol: maximum tolerance
    /// \param max_its: maximum number of iterations
    /// \param error_if_not_converged: throw an error if the tolerance was not satisfied
    /// \param max_simultaneous_rhs: maximum number of right-hand-sides solved simultaneously; all by default
    /// \param verb: verbosity level
    /// \param prefix: prefix printed before every line

    template <std::size_t NOp, typename COMPLEX>
    void fgmres(const Operator<NOp, COMPLEX>& op, Maybe<const Operator<NOp, COMPLEX>&> prec,
		const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX>& y,
		const std::string& order_t, unsigned int max_basis_size, double tol,
		unsigned int max_its = 0, bool error_if_not_converged = true,
		unsigned int max_simultaneous_rhs = 0, Verbosity verb = NoOutput,
		std::string prefix = "")
    {
      const bool do_ortho = false;

      // Get an unused label for the search subspace columns
      char Vc = get_free_label(x.order);
      char Vac = get_free_label(x.order + std::string(1, Vc));
      std::string order_cols = detail::union_dimensions(x.order, "", op.i.order);
      std::string order_rows = detail::union_dimensions(op.d.order, "", order_t);

      // Compute residual, r = op * y - x
      auto r = op(y) x.scale(-1).copyTo(r);
      auto normr0 = norm(r, order_t + order_cols); // type std::vector<real of T>

      // Allocate the search subspace U (onto the left singular space), and
      // Z (onto the right singular space)
      auto U = r.like_this<Nd + 1>(std::string("%") + std::string(1, Vc), '%', "",
				   {{Vc, max_basis_size + 1}});

      // Allocate the search subspace Z (onto the right singular space); when no preconditioning
      // then Z = U

      auto Z = U.kvslice_from_size({}, {{Vc, max_basis_size}}); // set Z an alias of U
      if (prec.hasSome())
	Z = Z.like_this(); // create new space when using preconditioning

      auto normr = normr0;
      unsigned int it = 0;
      double max_tol = HUGE_VAL;
      for (it = 0; it < max_its; ++it)
      {
	// U(:,0) = r;
	r.copyTo(U.kvslice_from_size({}, {{Vc, 1}}));

	// Expand the search subspace from residual
	for (unsigned int i = 0; i < max_basis_size; ++i) {
	  // Z(:,i) = prec * U(:,i)
	  if (prec.hasSome())
	    prec(U.kvslice_from_size({{Vc, i}}, {{Vc, 1}}))
	      .copyTo(Z.kvslice_from_size({{Vc, i}}, {{Vc, 1}}));

	  // U(:,i+1) = op * Z(:,i)
	  op(Z.kvslice_from_size({{Vc, i}}, {{Vc, 1}}))
	    .copyTo(U.kvslice_from_size({{Vc, i + 1}}, {{Vc, 1}}));
	}

	// Orthogonalize U and put it into W: W = orth(U(:,2:end))
	// NOTE: for small max_basis_size and assuming r is far from a left singular vector,
	//       a light or none orthogonalization should be enough
	auto Up = U.kvslice_from_size({{Vc, 1}}, {{Vc, max_basis_size}}).rename_dims({{Vc, Vac}});
	auto Uo = Up;
	if (do_ortho)
	{
	  Uo = Uo.clone();
	  ortho(none, Uo, order_t + order_cols, order_rows, std::string(1, Vac),
		1 /* one iteration should be enough */);
	}

	// Restrict to Uo: [x_rt H_rt] = Uo'*U = Up'*[r Up]
	auto x_rt = contract<NRhs - Nr + 1>(Uo.conj(), r, order_rows);
	auto H_rt = contract<NRhs - Nr + 2>(Uo.conj(), Up, order_rows);

	// Solve the projected problem: y_rt = (Uo'*U(:2:end))\(Uo'*r);
	auto y_rt = solve<NRhs - Nr + 1>(H_rt, x_rt, none, std::string(1, Vac), std::string(1, Vc));

	// Update solution: x += -Z*y_rt
	contract(Z.scale(-1), y_rt, order_columns, AddTo, x);

	// Update residual: r += -U(2:end)*y_rt
	contract(Up.scale(-1), y_rt, order_columns, AddTo, r);

	// Compute the norm
	auto normr = norm(r, order_t + order_cols);

	// Get the worse tolerance
	max_tol = 0;
	for (unsigned int i = 0; i < normr.size(); ++i)
	  max_tol = std::max(max_tol, normr[i] / normr0[i]);

	// Report iteration
	if (Verbosity >= Detailed)
	  QDPIO::cout << prefix << " MGPROTON FGMRES iteration it.: " << it
		      << " max rel. residual: " << max_tol << std::endl;

	if (max_tol <= tol) break;
      }

      // Report iteration
      if (Verbosity >= Summary)
	QDPIO::cout << prefix << " MGPROTON FGMRES summary its.: " << it
		    << " max rel. residual: " << max_tol << std::endl;
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
	std::string type = ops.getValue("type", StringOption{""}).getString();
	if (type.size() == 0)
	  throw std::runtime_error("Error getting a solver from option `" + ops.prefix +
				   "': missing tag `type'");
	if (solvers.count(type) == 0)
	{
	  std::string supported_tags;
	  for (const auto& it : solvers)
	    supported_tags += it->first + " ";
	  throw std::runtime_error("Error getting a solver from option `" + ops.prefix +
				   "': unsupported tag value `" + type +
				   "'; supported values: " + supported_tags);
	}

	return solvers[type](op, ops);
      }

      /// Returns a FGMRES solver
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getFGMRESSolver(const Operator<NOp, COMPLEX>& op, const Options& ops)
      {
	// Get preconditioner
	Maybe<const Operator<NOp, COMPLEX>> prec = none;
	const Options& precOps = ops.getValue("prec", NoneOption{}, Dictionary);
	if (precOps) prec = Maybe<const Operator<NOp, COMPLEX>>(getSolver(op, precOps));

	// Get the remainder options
	std::string order_t = getOption<std::string>(ops, "order_t", "");
	std::string order_rows =
	  getOption<std::string>(ops, "order_rows", getImage(op.map, order_t));
	std::string order_cols =
	  getOption<std::string>(ops, "order_cols", getDomain(op.map, order_t));
	unsigned int max_basis_size = getOption<unsigned int>(ops, "max_basis_size", 0);
	double tol = getOption<double>(ops, "tol", 0.0);
	unsigned int max_its = getOption<unsigned int>(ops, "max_its", 0);
	bool error_if_not_converged = getOption<bool>(ops, "max_basis_size", true);
	unsigned int max_simultaneous_rhs = getOption<unsigned int>(ops, "max_simultaneous_rhs", 0);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);
	std::string prefix = getOption<std::string>(ops, "prefix", "");

	// Return the solver
	const Operator<NOp, COMPLEX> op_ = op; // copy by value of the operator
	return {[=](const Tensor<NOp + 1, COMPLEX>& x) {
		  auto y = x.like_this();
		  fgmres(op_, prec, x, y, order_t, order_rows, order_cols, max_basis_size, tol,
			 max_its, error_if_not_converged, max_simultaneous_rhs, verb, prefix);
		  return y;
		},
		op.i, op.d};
      }

      /// Returns the \gamma_5 or a projection for a given number of spins
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
	  Tensor<1, COMPLEX> ones("j", {Ns}, OnHost, OnEveryoneReplicated);
	  for (int i = 0; i < Ns; ++i)
	    ones.set({{i}}, COMPLEX{1});
	  Tensor<1, COMPLEX> p("i", {Ns}, OnHost, OnEveryoneReplicated);
	  p.contract(getGamma5(Ns), {}, NotConjugate, ones, {}, NotConjugate);
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

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> cloneOperator(const Operator<NOp, COMPLEX>& op,
					   unsigned int max_dist_neigbors = 1)
      {
	// Unsupported explicitly colorized operators
	if (op.domain_eg.kvdim().count('X') == 0 || op.domain_eg.kvdim()['X'] != 2)
	  throw std::runtime_error("cloneEOOperator: unsupported not explicitly colored operators");

	// If the operator is empty, just return itself
	if (op.d.volume() == 0 || op.i.volume() == 0)
	  return op;

	// Create the ordering for the domain and the image, moving the lattice dimensions as the slowest indices
	// NOTE: assuming that x,y,z,t are the only non-blocking dimensions
	std::set<char> lattice_dims{'x', 'y', 'z', 't', 'X'};
	remap ri = getNewLabels(op.i.order, op.d.order + std::string{"u~"});
	auto d = op.d.reorder("%Xxyzt", '%');
	auto i = op.i.reorder("%Xxyzt", '%').rename_dims(ri);
	const auto dim = d.kvdim();

	// Get the blocking for the domain and the image
	std::map<char, int> blkd, blki;
	for (const auto& it : d.kvdim())
	  blkd[it.first] = (lattice_dims.count(it.first) == 0 ? it.second : 1);
	for (const auto& it : i.kvdim())
	  blki[it.first] = (lattice_dims.count(it.first) == 0 ? it.second : 1);

	// Generate a remap for dimensions that are going to be just one element
	remap rd1 = getNewLabels(op.d.order, op.d.order);
	remap rd1_latt, rd1_nonlatt;
	for (const auto& it : dim)
	{
	  if (lattice_dims.count(it.first) == 0)
	    rd1_nonlatt[it.first] = it.second;
	  else
	    rd1_latt[it.first] = it.second;
	}
	std::string rd1_labels;
	for (const auto& it : dim)
	  rd1_labels.push_back(it.first);

	// Generate a remap for the columns of the probing vectors, which has
	// one column for each element of the non-lattice dimensions
	remap rdc = getNewLabels(op.d.order, op.d.order + rd1_labels);

	// Construct the probing vectors, which they have as the rows the domain and as columns the domain blocking
	// dimensions
	
	auto t_bl = op.d.like_this(none, blkd, OnHost, OnEveryoneReplicated);
	t_bl.set_zero();
	auto t_blbl = contract<NOp * 2>(t_bl.rename_dims(rd1_latt), t_bl.rename_dims(rdc), "");
	COMPLEX* t_blbl_data = t_blbl.data.get();
	for (std::size_t i = 0, vol = t_bl.volume(); i < vol; ++i)
	  t_blbl_data[i * vol + i] = COMPLEX{1};

	std::map<char, int> nonblkd;
	for (const auto& it : kvdim)
	  nonblkd[it.first] = (lattice_dims.count(it.first) == 0 ? 1 : it.second);

	// Compute the coloring
	Coloring coloring{
	  0, max_dist_neigbors * 2 + 1, {{dim['x'] * dim['X'], dim['y'], dim['z'], dim['t']}}};
	unsigned int num_colors = coloring.numColors();

	// Get the number of neighbors
	unsigned int num_neighbors = 1;
	for (const auto& it : kvdim)
	  if (lattice_dims.count(it.first) == 1)
	    num_neighbors += (dim[it.first] <= 1 ? 0 : (dim[it.first] <= 2 ? 1 : 2));

	// Create the sparse tensor
	SpTensor<NOp, NOp, COMPLEX> sop{d, i, blkd, blki, num_neighbors};

	for (unsigned int color = 0; color < num_colors; ++color)
	{
	  // Extracting the proving vector
	  auto t_l = fillLatticeField<NOp, T>(
	    d.reorder("xyztX%", '%').order, nonblkd, OnDefaultDevice, [&](Coor<NOp - 1> c) {
	      return coloring.getColor({{c[0], c[1], c[2], c[3]}}) == color ? T{1} : T{0};
	    });

	  // Contracting the proving vector with the blocking components
	  auto probs = contract<NOp * 3>(t_l, t_blbl, "");

	  // Compute the matvecs
	  auto mv = op(probs).rename_dims(ri).rename_dims(reverse(rdc));

	  // Construct an indicator tensor where all blocking dimensions but only the nodes colored `color` are copied
	  auto sel = contract<NOp * 3, float>(
	    fillLatticeField<NOp, float>(
	      d.reorder("xyztX%", '%').order, nonblkd, OnDefaultDevice,
	      [&](Coor<NOp - 1> c) {
		return coloring.getColor({{c[0], c[1], c[2], c[3]}}) == color ? 1 : 0;
	      }),
	    op.d.like_this<float>(none, blkd, OnHost, OnEveryoneReplicated).set(1), "");

	  // Split mv and sel into subsets where the neighbors are at the same position
	  int maxX = dim['X'];
	  int maxY = std::min(2, dim['y']);
	  int maxZ = std::min(2, dim['z']);
	  int maxT = std::min(2, dim['t']);
	  auto sel_eo = sel.split_dimension('y', "Yy", maxY)
			  .split_dimension('z', "Zz", maxZ)
			  .split_dimension('t', "Tt", maxT);
	  auto mv_eo = mv.split_dimension('y', "Yy", maxY)
			 .split_dimension('z', "Zz", maxZ)
			 .split_dimension('t', "Tt", maxT);

	  // Populate the nonzeros
	  int dimx2 = dim['x'] * 2;
	  int mu = 0; // 0: self; 1: neg. x; 2: pos. x; 3: neg. y; 4: pos. y...
	  for (char ldir : std::vector<char>{0, 'x', 'y', 'z', 't'})
	  {
	    // Process lattice dimensions in this operator only
	    if (ldir != 0 && d.kvdim().count(ldir) == 0)
	      continue;

	    for (unsigned int dir = -1; dir < 1; ++dir)
	    {
	      // Skip invalid combinations:
	      //  a) ldir=0 is the self direction and does not have any other neighbors
	      //  b) the lattice directions, x,y,z,t may have a neighbor on the negative or on the positive direction
	      //  c) a lattice direction with size <= 1 does not have any neighbors
	      //  d) a lattice direction with size == 2 only have one neighbor
	      if ((ldir == 0 && dir != 0) ||		     // a)
		  (ldir != 0 && dir == 0) ||		     // b)
		  (ldir != 0 && dim[ldir] <= 1) ||	     // c)
		  (ldir != 0 && dim[ldir] == 2 && dir >= 0)) // d)
		continue;

	      // Copy the nonzeros from the proving vector results to the sparse tensor data
	      if (ldir != 'x' || maxX == 1)
	      {
		// Case where the direction isn't 'x' or not using even-odd layout (dim['X'] == 1)
		std::map<char, int> to{{'X', ldir == 0 ? 0 : 1}, {'u', mu}};
		if (ldir != 0)
		  to[ldir] = dir;

		mv.copyToWithMask(sop.data.kvslice_from_size(to), sel);
	      }
	      else
	      {
		// Case where the direction is 'x' and we are using even-odd layout
		std::map<char, int> XYZTsize{{'X', 1}, {'Y', 1}, {'Z', 1}, {'T', 1}};
		std::map<char, int> to{{'X', 1}, {'u', mu}};
		int dir0 = (dir < 0 ? dir + dimx2 : dir); // avoid modulus with negative numbers
		for (int T = 0; T < maxT; ++T)
		{
		  for (int Z = 0; Z < maxZ; ++Z)
		  {
		    for (int Y = 0; Y < maxY; ++Y)
		    {
		      for (int X = 0; X < maxX; ++X)
		      {
			std::map<char, int> XYZTfrom{{'X', X}, {'Y', Y}, {'Z', Z}, {'T', T}};
			to['x'] = (dir0 + ((X + Y + Z + T) % 2)) / 2;
			mv.kvslice_from_size(XYZTfrom, XYZTsize)
			  .copyToWithMask(
			    sop.data.kvslice_from_size(to).kvslice_from_size(XYZTfrom, XYZTsize),
			    sel.kvslice_from_size(XYZTfrom, XYZTsize));
		      }
		    }
		  }
		}
	      }

	      // Increase mu
	      mu++;
	    } // dir
	  }   // ldir
	}     // color

	// Populate the coordinate of the columns
	Coor<Nd> dims{{dim['x'] * dim['X'], dim['y'], dim['z'], dim['t']}}}
	sop.jj.fillCpuFunCoor([](const Coor<NOp + 2>& c) {
	  // c has order '~uXxyzt...' where Xxyzt... were remapped by ri
	  int base = c[c[0] + 2];
	  if (c[1] == 0)
	  {
	    // self direction, do nothing
	  }
	  else if ((c[1] == 1 || c[1] == 2) && c[0] == 0)
	  {
	    // Move one step behind or forward on the x direction and update X
	    base = (base + 1) % maxX;
	  }
	  else if ((c[1] == 1 || c[1] == 2) && c[0] == 1)
	  {
	    // Move one step behind or forward on the x direction and update x
	    int dir0 = (c[1] == 1 ? dimx2 - 1 : 1); // avoid modulus with negative numbers
	    base = (dir0 + ((c[2] / 2 + c[3] / 2 + c[4] / 2 + c[5] / 2) % 2)) / 2;
	  }
	  else if ((c[1] - 1) / 2 == c[0] - 1)
	  {
	    // Move on the same direction as the coordinate
	    base += (c[1] % 2 == 1 ? dims[(c[1] - 1) / 2] - 1 : 1);
	  }
	  return base;
	});

	// Construct the sparse operator
	sop.construct();

	// Return the operator
	return getOperator(sop);
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
	// Create num_null_vecs random vectors and solve them
	auto b = op.d.like_this<NOp + 1>("%n", '%', {{'n', num_null_vecs}});
	urand(b, -1, 1);
	auto nv = null_solver(std::move(b));

	// Do chirality splitting nv2 = [ nv * gpos, nv * gneg ]
	// TODO: experiment without chirality splitting
	int ns = nv.kvdims()['s'];
	if (ns != 2 && ns != Ns)
	  throw std::runtime_error("Error in getMGProlongator: Unsupported spin number");

	auto nv2 = nv.like_this(none, {'n', num_null_vecs * 2});
	auto g5 = getGamma5(ns, OnHost), gpos = g5.clone(), gneg = g5.clone();
	for (int i = 0; i < Ns; ++i) // make diagonal entries of gpos all positive or zero
	  gpos.set({{i, i}}, g5.get({{i, i}} + 1));
	for (int i = 0; i < Ns; ++i) // make diagonal entries of gneg all negative or zero
	  gneg.set({{i, i}}, g5.get({{i, i}} - 1));
	nv2.kvslice_from_size({}, {{'n', num_null_vecs}})
	  .contract(g5pos, {}, NotConjugate, nv, {}, NotConjugate);
	nv2.kvslice_from_size({{'n', num_null_vecs}}, {{'n', num_null_vecs}})
	  .contract(g5neg, {}, NotConjugate, nv, {}, NotConjugate);

	// Do the blocking
	auto nv_blk = nv.split_dimension('x', "Wx", blocking['x'])
			.split_dimension('y', "Yy", blocking['y'])
			.split_dimension('z', "Zz", blocking['z'])
			.split_dimension('t', "Tt", blocking['t'])
			.split_dimension('n', "CS", num_null_vecs);

	// Do the orthogonalization on each block and chirality
	ortho(Tensor<1, COMPLEX>(), nv_blk, "xyztn", "WYZTXsc", "N");

	// Return the operator
	auto dimd = nv_blk.kvdim();
	dimd['c'] = num_null_vecs;
	dimd['s'] = 2;
	Tensor<NOp, COMPLEX> d = d.like_this(none, dimd), i = op.i;
	return {[=](Tensor<NOp + 1, COMPLEX> t) {
		  auto out = i.like_this<NOp + 1>("%n", '%', {{'n', t.kvdim()['n']}});
		  auto out_blk = out.split_dimension('x', "Wx", blocking['x'])
				  .split_dimension('y', "Yy", blocking['y'])
				  .split_dimension('z', "Zz", blocking['z'])
				  .split_dimension('t', "Tt", blocking['t']);
		  contract(nv_blk, t, "", CopyTo, out_blk);
		  return out;
		},
		d, i,
		[=](Tensor<NOp + 1, COMPLEX> t) {
		  auto t_blk = t.split_dimension('x', "Wx", blocking['x'])
				 .split_dimension('y', "Yy", blocking['y'])
				 .split_dimension('z', "Zz", blocking['z'])
				 .split_dimension('t', "Tt", blocking['t']);
		  return contract<NOp + 1>(nv_blk.conj(), t_blk, "WYZT");
		}};
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
	std::vector<unsigned int> blocking_v =
	  getOption<std::vector<unsigned int>>(ops, "blocking");
	if (blocking_v.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getMGPrec: the blocking should be a vector with four elements");
	std::map<char, unsigned int> blocking{
	  {'x', blocking_v[0]}, {'y', blocking_v[1]}, {'z', blocking_v[2]}, {'t', blocking_v[3]}};
	const Options& nvInvOps = ops.getValue("solver_null_vecs", none, Dictionary);
	const Operator<NOp, COMPLEX> nullSolver = getSolver(op, nvInvOps);
	const Operator<NOp, COMPLEX> V = getMGProlongator(op, num_null_vecs, blocking, nullSolver);

	// Compute the coarse operator, V' * op * V
	const Operator<NOp, COMPLEX> op_c = cloneOperator(Operator<NOp, COMPLEX>{
	  [](Tensor<NOp + 1, COMPLEX> t) { return V.tconj()(op(V(x))); }, V.d, V.d});

	// Get the solver for the projector
	const Options& coarseSolverOps = ops.getValue("solver_coarse", none, Dictionary);
	const Operator<NOp, COMPLEX> coarseSolver = getSolver(op_c, coarseSolverOps);

	// Get the solver for the smoother
	const Options& opSolverOps = ops.getValue("solver_smoother", none, Dictionary);
	const Operator<NOp, COMPLEX> opSolver = getSolver(op, opSolverOps);

	// Return the solver
	return {[=](Tensor<NOp + 1, COMPLEX> x) {
		  // y0 = V*solver(V'*Op*V, V'x)
		  auto y0 = V(coarseSolver(V.tconj()(x)));

		  // x1 = x - op*y0
		  auto x1 = op(y0.scale(-1));
		  x1 += std::move(x);

		  // y = y0 + solver(Op, x1)
		  auto y = opSolver(std::move(x1));
		  y += y0;

		  return y;
		},
		op.i, op.d, nullptr, reverse(op.map)};
      }
    }

    /// Returns an operator that approximate the inverse of a given operator
    /// \param op: operator to make the inverse of
    /// \param ops: options to select the solver from `solvers` and influence the solver construction

    template <std::size_t NOp, typename COMPLEX>
    Operator<NOp, COMPLEX> getSolver(const Operator<NOp, COMPLEX>& op, const Options& ops)
    {
      static const std::map<std::string, Solver<NOp, COMPLEX>> solvers{
	{"fgmres", detail::getFGMRESSolver<NOp, COMPLEX>}, // flexible GMRES
	{"mg", detail::getMGPrec<NOp, COMPLEX>}		   // Multigrid
      };

      return detail::getSolver(solvers, op, ops);
    }

    /// Return an Operator that wraps up a LinearOperator<LatticeFermion>
    inline Operator<Nd + 3, Complex> asOperatorView(const LinearOperator<LatticeFermion>& linOp)
    {
      LatticeFermion a;
      auto d = asTensorView(a).toComplex().make_eg();
      return {[](Tensor<Nd + 4, Complex> t) {
		auto r = t.like_this();
		LatticeFermion x, y;
		unsigned int n = t.kvdim()['n'];
		for (unsigned int i = 0; i < n; ++i)
		{
		  t.kvslice_from_size({{'n', i}}, {{'n', 1}}).copyTo(x);
		  y = zero();
		  linOp(x, y, PLUS /* I believe, it's ignored */);
		  asTensorView(y).copyTo(r.kvslice_from_size({{'n', i}}, {{'n', 1}}));
		}
		return r;
	      },
	      d, d};
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
      Handle<
	FermionState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>
	state;

      /// Chroma solver (optional)
      Handle<SystemSolver<LatticeFermion>> PP;

      /// Operator on scxyztX (optional)
      Maybe<Operator<Nd + 3>> op;

      /// Constructor
      /// \param fermAction: XML for the fermion action
      /// \param invParam: XML for the quark propagator
      /// \param u: gauge fields

      ChimeraSolver(const GroupXML_t& fermAction, const GroupXML_t& invParam,
		    const multi1d<LatticeColorMatrix>& u)
      {
	// Initialize fermion action and state
	std::istringstream  xml_s(fermAction.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << fermAction.id << std::endl;
	S = TheFermionActionFactory::Instance().createObject(fermAction.id, fermacttop,
							     fermAction.path);
	state = S->createState(u);

	// If the inverter is MGPROTON, use this infrastructure
	if (invParam.id == "MGPROTON")
	{
	  // Parse XML with the inverter options
	  Options ops = getOptionsFromXML(invParam.xml);

	  // Clone the matvec
	  LinearOperator<LatticeFermion>* fLinOp = S->genLinOp(state);
	  Operator<Nd + 3, Complex> linOp = cloneOperator(asOperatorView(*fLinOp));

	  // Construct the solver
	  op = Maybe{getSolver(linOp, ops)};
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
      Tensor<Nd + 5, COMPLEX_OUT> doInversion(const Operator<Nd + 5, COMPLEX_OUT>& op,
					      const Tensor<Nd + 3, COMPLEX_CHI> chi, int t_source,
					      int first_tslice_out, int n_tslice_out,
					      const std::vector<int>& spin_sources, int max_rhs,
					      const std::string& order_out = "cSxyztXns")
      {
	Tensor<Nd + 5, COMPLEX_OUT> psi(
	  order_out,
	  latticeSize<Nd + 5>(
	    order_out,
	    {{'t', n_tslice_out}, {'S', Ns}, {'s', spin_sources.size()}, {'n', num_vecs}}),
	  chi.getDev());

	// Create tensors with full support on the lattice
	int max_step = std::max(num_vecs, max_rhs);
	auto aux = chi.like_this(none, {{'n', max_step}, {Layout::lattSize()[3]}});

	for (int spin_source : spin_sources)
	{
	  for (int n0 = 0, n_step = std::min(max_rhs, num_vecs); n0 < num_vecs;
	       n0 += n_step, n_step = std::min(n_step, num_vecs - n0))
	  {
	    auto aux0 = aux.kvslice_from_size({}, {{'n', n_step}});
	    aux0.set_zero();
	    chi.kvslice_from_size({{'n', n}}, {{'n', n_step}})
	      .copyTo(aux0.kvslice_from_size({{'t', t_source}, {'s', spin_source}}));

	    // Solve
	    op(aux0)
	      .kvslice_from_size({{'t', first_tslice_out}}, {{'t', n_tslice_out}})
	      .rename_dims({{'s', 'S'}})
	      .copyTo(psi.kvslice_from_size({{'n', n}, {'s', spin_source}}));
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

      if (sol.op.hasSome())
	detail::doInversion<Nd + 5, COMPLEX_OUT>(sol.op.getSome(), chi, t_source, first_tslice_out,
						 n_tslice_out, spin_sources, max_rhs, order_out);
      else
	detail::doInversion<Nd + 5, COMPLEX_OUT>(sol.PP, chi, t_source, first_tslice_out,
						 n_tslice_out, spin_sources, max_rhs, order_out);

      snarss1.stop();
      QDPIO::cout << "Time to compute inversions for " << spin_sources.size()
		  << " spin sources and " << num_vecs
		  << " colorvecs : " << snarss1.getTimeInSeconds() << " secs" << std::endl;
    }
  }
}

#endif // BUILD_SB

#endif // __INCLUDE_MGPROTON__
