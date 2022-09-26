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

    /// Ordering of matrices

    enum ColOrdering {
      RowMajor,	   ///< row-major ordering, the fastest index is the column
      ColumnMajor, ///< row-major ordering, the fastest index is the row
    };

    /// Operator's layout
    enum OperatorLayout {
      NaturalLayout,  ///< natural ordering
      XEvenOddLayout, ///< X:(x+y+z+t)%2, x:x/2, y:y, z:z, t:t
      EvensOnlyLayout ///< x:x/2, y:y, z:z, t:t for all (x+y+z+t)%2==0
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

    /// Displacements of each site nonzero edge for every operator's site

    using NaturalNeighbors = std::vector<std::map<char, int>>;

    /// Special value to indicate that the operator is dense

    inline const NaturalNeighbors& DenseOperator()
    {
      static const NaturalNeighbors dense{{{(char)0, 0}}};
      return dense;
    }

    /// Representation of an operator, function of type tensor -> tensor where the input and the
    /// output tensors have the same dimensions

    template <std::size_t NOp, typename COMPLEX>
    using OperatorFun =
      std::function<void(const Tensor<NOp + 1, COMPLEX>&, Tensor<NOp + 1, COMPLEX>)>;

    /// Representation of an operator, function of type tensor -> tensor where the output tensor
    /// has powers of the operator applied to the input tensor

    template <std::size_t NOp, typename COMPLEX>
    using OperatorPowerFun =
      std::function<void(const Tensor<NOp + 2, COMPLEX>&, Tensor<NOp + 2, COMPLEX>, char)>;

    /// Representation of a function that takes and returns tensors with the same labels, although the
    /// dimensions may be different.

    template <std::size_t NOp, typename COMPLEX>
    struct Operator {
      /// Function that the operators applies (optional)
      OperatorFun<NOp, COMPLEX> fop;
      /// Example tensor for the input tensor (domain)
      Tensor<NOp, COMPLEX> d;
      /// Example tensor for the output tensor (image)
      Tensor<NOp, COMPLEX> i;
      /// Function to apply when conjugate transposed (optional)
      OperatorFun<NOp, COMPLEX> fop_tconj;
      /// Labels that distinguish different operator instantiations
      std::string order_t;
      /// Operator's domain space layout
      OperatorLayout domLayout;
      /// Operator's image space layout
      OperatorLayout imgLayout;
      /// Neighbors for each site in this operator
      NaturalNeighbors neighbors;
      /// Preferred ordering
      ColOrdering preferred_col_ordering;
      /// Operator based on sparse tensor (optional)
      SpTensor<NOp, NOp, COMPLEX> sp;
      /// Sparse tensor map from image labels to domain labels (optional)
      remap rd;
      /// Operator maximum power support (optional)
      unsigned int max_power;
      /// Identity
      std::shared_ptr<std::string> id;

      /// Empty constructor
      Operator()
      {
      }

      /// Constructor
      Operator(const OperatorFun<NOp, COMPLEX>& fop, Tensor<NOp, COMPLEX> d, Tensor<NOp, COMPLEX> i,
	       const OperatorFun<NOp, COMPLEX>& fop_tconj, const std::string& order_t,
	       OperatorLayout domLayout, OperatorLayout imgLayout, NaturalNeighbors neighbors,
	       ColOrdering preferred_col_ordering, const std::string& id = "")
	: fop(fop),
	  d(d),
	  i(i),
	  fop_tconj(fop_tconj),
	  order_t(order_t),
	  domLayout(domLayout),
	  imgLayout(imgLayout),
	  neighbors(neighbors),
	  preferred_col_ordering(preferred_col_ordering),
	  sp{},
	  rd{},
	  max_power{0},
	  id(std::make_shared<std::string>(id))
      {
      }

      /// Constructor for a power-supported function
      Operator(const SpTensor<NOp, NOp, COMPLEX>& sp, const remap& rd, unsigned int max_power,
	       Tensor<NOp, COMPLEX> d, Tensor<NOp, COMPLEX> i, const std::string& order_t,
	       OperatorLayout domLayout, OperatorLayout imgLayout, NaturalNeighbors neighbors,
	       ColOrdering preferred_col_ordering, const std::string& id = "")
	: fop{},
	  d(d),
	  i(i),
	  fop_tconj{},
	  order_t(order_t),
	  domLayout(domLayout),
	  imgLayout(imgLayout),
	  neighbors(neighbors),
	  preferred_col_ordering(preferred_col_ordering),
	  sp{sp},
	  rd{rd},
	  max_power{max_power},
	  id(std::make_shared<std::string>(id))
      {
      }

      /// Constructor from other operator
      Operator(const OperatorFun<NOp, COMPLEX>& fop, Tensor<NOp, COMPLEX> d, Tensor<NOp, COMPLEX> i,
	       const OperatorFun<NOp, COMPLEX>& fop_tconj, const Operator<NOp, COMPLEX>& op,
	       const std::string& id = "")
	: fop(fop),
	  d(d),
	  i(i),
	  fop_tconj(fop_tconj),
	  order_t(op.order_t),
	  domLayout(op.domLayout),
	  imgLayout(op.imgLayout),
	  neighbors(op.neighbors),
	  preferred_col_ordering(op.preferred_col_ordering),
	  sp{},
	  rd{},
	  max_power{0},
	  id(std::make_shared<std::string>(id))
      {
      }

      /// Return whether the operator is not empty
      explicit operator bool() const noexcept
      {
	return (bool)d;
      }

      /// Return the transpose conjugate of the operator
      Operator<NOp, COMPLEX> tconj() const
      {
	if (sp || !fop_tconj)
	  throw std::runtime_error("Operator does not have conjugate transpose form");
	return {
	  fop_tconj, i, d, fop, order_t, imgLayout, domLayout, neighbors, preferred_col_ordering};
      }

      /// Return whether the operator has transpose conjugate
      bool has_tconj() const
      {
	return !sp && fop_tconj;
      }

      /// Apply the operator
      template <std::size_t N, typename T>
      Tensor<N, T> operator()(const Tensor<N, T>& t) const
      {
	// The `t` labels that are not in `d` are the column labels
	std::string cols = detail::remove_dimensions(t.order, d.order); // t.order - d.order

	if (sp)
	{
	  remap mcols = detail::getNewLabels(cols, sp.d.order + sp.i.order);
	  auto y = i.template like_this<N, T>(
	    preferred_col_ordering == ColumnMajor ? std::string("%") + cols : cols + "%", '%', "",
	    t.kvdim());
	  sp.contractWith(t.rename_dims(mcols), rd, y.rename_dims(mcols), {});
	  return y;
	}
	else
	{
	  auto x =
	    t.template collapse_dimensions<NOp + 1>(cols, 'n', true).template make_sure<COMPLEX>();
	  auto y = i.template like_this<NOp + 1>(
	    preferred_col_ordering == ColumnMajor ? "%n" : "n%", '%', "", {{'n', x.kvdim()['n']}});
	  fop(x, y);
	  return y.template split_dimension<N>('n', cols, t.kvdim()).template make_sure<T>();
	}
      }

      /// Apply the operator
      template <std::size_t N, typename T>
      void operator()(const Tensor<N, T>& x, Tensor<N, T> y, char power_label = 0) const
      {
	// The `x` labels that are not in `d` are the column labels
	std::string cols_and_power =
	  detail::remove_dimensions(x.order, d.order); // x.order - d.order
	std::string cols =
	  power_label == 0 ? cols_and_power
			   : detail::remove_dimensions(cols_and_power, std::string(1, power_label));
	int power = power_label == 0 ? 1 : y.kvdim().at(power_label);
	if (power <= 0)
	  return;

	if (sp)
	{
	  remap mcols = detail::getNewLabels(cols_and_power, sp.d.order + sp.i.order);
	  if (power == 1)
	  {
	    sp.contractWith(x.rename_dims(mcols), rd, y.rename_dims(mcols), {});
	  }
	  else
	  {
	    char power_label0 = mcols.count(power_label) == 1 ? mcols.at(power_label) : power_label;
	    auto x0 = x.rename_dims(mcols);
	    auto y0 = y.rename_dims(mcols);
	    int power0 = (max_power == 0 ? power : std::min(power, (int)max_power));
	    sp.contractWith(x0, rd, y0.kvslice_from_size({}, {{power_label0, power0}}), {},
			    power_label0);
	    for (int i = power0, ni = std::min((int)max_power, power - i); i < power;
		 i += ni, ni = std::min((int)max_power, power - i))
	    {
	      sp.contractWith(y0.kvslice_from_size({{power_label0, i - 1}}, {{power_label0, 1}}),
			      rd, y0.kvslice_from_size({{power_label0, i}}, {{power_label0, ni}}),
			      {}, power_label0);
	    }
	  }
	}
	else
	{
	  if (power_label == 0)
	  {
	    auto x0 = x.template collapse_dimensions<NOp + 1>(cols_and_power, 'n', true)
			.template cast<COMPLEX>();
	    auto y0 = y.template collapse_dimensions<NOp + 1>(cols_and_power, 'n', true)
			.template cast_like<COMPLEX>();
	    fop(x0, y0);
	    y0.copyTo(y);
	  }
	  else if (power > 0)
	  {
	    operator()(x.kvslice_from_size({}, {{power_label, 1}}),
		       y.kvslice_from_size({}, {{power_label, 1}}));
	    for (int i = 1, p = power; i < p; ++i)
	      operator()(y.kvslice_from_size({{power_label, i - 1}}, {{power_label, 1}}),
			 y.kvslice_from_size({{power_label, i}}, {{power_label, 1}}));
	  }
	}
      }

      /// Return this operator with an implicit different type
      /// \tparam T: new implicit precision

      template <typename T = COMPLEX,
		typename std::enable_if<std::is_same<T, COMPLEX>::value, bool>::type = true>
      Operator<NOp, T> cast() const
      {
	return *this;
      }

      template <typename T,
		typename std::enable_if<!std::is_same<T, COMPLEX>::value, bool>::type = true>
      Operator<NOp, T> cast() const
      {
	if (!*this)
	  return {};

	if (sp)
	{
	  return Operator<NOp, T>(sp.template cast<T>(), rd, max_power, d.template cast<T>(),
				  i.template cast<T>(), order_t, domLayout, imgLayout, neighbors,
				  preferred_col_ordering);
	}
	else
	{
	  const Operator<NOp, COMPLEX> op = *this,
				       op_tconj = has_tconj() ? tconj() : Operator<NOp, COMPLEX>{};
	  return Operator<NOp, T>(
	    [=](const Tensor<NOp + 1, T>& x, Tensor<NOp + 1, T> y) { op(x, y); },
	    d.template cast<T>(), i.template cast<T>(),
	    op_tconj ? [=](const Tensor<NOp + 1, T>& x, Tensor<NOp + 1, T> y) { op_tconj(x, y); }
		     : OperatorFun<NOp, T>{},
	    order_t, domLayout, imgLayout, neighbors, preferred_col_ordering);
	}
      }

      /// Return a slice of the tensor starting at coordinate `dom_kvfrom`, `img_kvfrom` and taking
      /// `dom_kvsize`, `img_kvsize` elements in each direction. The missing dimensions in `*_kvfrom`
      /// are set to zero and the missing directions in `*_kvsize` are set to the size of the tensor.
      ///
      /// \param dom_kvfrom: dictionary with the index of the first element in each domain direction
      /// \param dom_kvsize: dictionary with the number of elements in each domain direction
      /// \param img_kvfrom: dictionary with the index of the first element in each domain direction
      /// \param img_kvsize: dictionary with the number of elements in each domain direction
      /// \return: a copy of the tensor or an implicit operator

      Operator<NOp, COMPLEX> kvslice_from_size(const std::map<char, int>& dom_kvfrom = {},
					       const std::map<char, int>& dom_kvsize = {},
					       const std::map<char, int>& img_kvfrom = {},
					       const std::map<char, int>& img_kvsize = {}) const
      {
	if (!*this)
	  return {};

	// Update the eg layouts and the data layouts
	auto new_d = d.kvslice_from_size(dom_kvfrom, dom_kvsize);
	auto new_i = i.kvslice_from_size(img_kvfrom, img_kvsize);
	OperatorLayout new_domLayout =
	  (new_d.kvdim().at('X') != d.kvdim().at('X') && domLayout == XEvenOddLayout
	     ? EvensOnlyLayout
	     : domLayout);
	OperatorLayout new_imgLayout =
	  (new_i.kvdim().at('X') != i.kvdim().at('X') && imgLayout == XEvenOddLayout
	     ? EvensOnlyLayout
	     : imgLayout);
	if (sp)
	{
	  return Operator<NOp, COMPLEX>(sp.kvslice_from_size(detail::update_kvcoor(dom_kvfrom, rd),
							     detail::update_kvcoor(dom_kvsize, rd),
							     img_kvfrom, img_kvsize),
					rd, max_power, new_d, new_i, order_t, new_domLayout,
					new_imgLayout, neighbors, preferred_col_ordering);
	}
	else
	{
	  const Operator<NOp, COMPLEX> op = *this,
				       op_tconj = has_tconj() ? tconj() : Operator<NOp, COMPLEX>{};
	  return Operator<NOp, COMPLEX>(
	    [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = op.d.template like_this<NOp + 1>(
		op.preferred_col_ordering == ColumnMajor ? "%n" : "n%", '%', "",
		{{'n', x.kvdim().at('n')}});
	      x0.set_zero();
	      x.copyTo(x0.kvslice_from_size(dom_kvfrom, dom_kvsize));
	      auto y0 = op.i.template like_this<NOp + 1>(
		op.preferred_col_ordering == ColumnMajor ? "%n" : "n%", '%', "",
		{{'n', x.kvdim().at('n')}});
	      op(x0, y0);
	      y0.kvslice_from_size(img_kvfrom, img_kvsize).copyTo(y);
	    },
	    new_d, new_i,
	    op_tconj ? [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      auto x0 = op_tconj.d.template like_this<NOp + 1>(
		op_tconj.preferred_col_ordering == ColumnMajor ? "%n" : "n%", '%', "",
		{{'n', x.kvdim().at('n')}});
	      x0.set_zero();
	      x.copyTo(x0.kvslice_from_size(img_kvfrom, img_kvsize));
	      auto y0 = op_tconj.i.template like_this<NOp + 1>(
		op_tconj.preferred_col_ordering == ColumnMajor ? "%n" : "n%", '%', "",
		{{'n', x.kvdim().at('n')}});
	      op_tconj(x0, y0);
	      y0.kvslice_from_size(dom_kvfrom, dom_kvsize).copyTo(y);
            } : OperatorFun<NOp, COMPLEX>{},
	    order_t, new_domLayout, new_imgLayout, neighbors, preferred_col_ordering);
	}
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

      /// Returns max ||I - C||_F, which is a heuristic of the orthogonality level of C=V^\dagger*V
      /// \param C: matrix to test
      /// \param order_t: dimension labels that do not participate in the orthogonalization
      /// \param order_rows: dimension labels that are the rows of the matrices V and W
      /// \param order_cols: dimension labels that are the columns of the matrices V and W

      template <std::size_t Nrows, std::size_t Ncols, std::size_t N, typename COMPLEX>
      double ortho_level(Tensor<N, COMPLEX> C, const std::string& order_t,
			 const std::string& order_rows, const std::string& order_cols)
      {
	// Check Nrows
	if (order_rows.size() != Nrows)
	  throw std::runtime_error("ortho_level: invalid template argument `Nrows`");
	if (order_cols.size() != Ncols)
	  throw std::runtime_error("ortho_level: invalid template argument `Ncols`");

	// Check that the matrix is square
	auto dim = C.kvdim();
	std::size_t m = volume(dim, order_rows);
	if (m != volume(dim, order_cols))
	  throw std::runtime_error("ortho_level: the input tensor is not square");

	// Construct the identity
	auto t = C.template like_this<Nrows + Ncols>(order_rows + order_cols, {}, OnHost, OnMaster);
	t.set_zero();
	if (t.getLocal())
	{
	  COMPLEX* t_data = t.data.get();
	  for (std::size_t i = 0; i < m; ++i)
	    t_data[i * m + i] = COMPLEX{1};
	}
	auto tones = C.template like_this<N - Nrows - Ncols>(order_t, {}, OnHost, OnMaster);
	tones.set(COMPLEX{1});
	auto I = contract<N>(t.make_sure(none, OnDefaultDevice),
			     tones.make_sure(none, OnDefaultDevice), "");

	// Compute ||I - C||_F^2
	C.scale(-1).addTo(I);
	auto fnorm2 = contract<N - Nrows - Ncols>(I.conj(), I, order_rows + order_cols);

	// Return the square root
	double tol = 0;
	fnorm2.make_sure(none, OnHost, OnEveryoneReplicated)
	  .foreachWithCPUFun(
	    [&](const COMPLEX& t) { tol = std::max(tol, (double)std::sqrt(std::real(t))); });
	return tol;
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
	       unsigned int max_its = 4, Verbosity verb = NoOutput, const std::string& prefix = "")
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
      double l = 0;
      unsigned int i = 0;
      for (; i <= max_its;)
      {
	// W = W - V*(V'*W)
	if (V)
	  contract(V.scale(-1), contract<NV + NW - Nrows * 2>(V.conj(), W, order_rows), Vcorder,
		   AddTo, W);

	// Compute Wa'*W and the orthogonality level of the basis
	auto C = contract<Nt + 2 * Ncols>(Wa.conj(), W, order_rows);
	l = detail::ortho_level<Ncols, Ncols>(C, order_t, Wacorder, Wcorder);
	if (verb >= Detailed)
	  QDPIO::cout << prefix << " ortho #its: " << i << " |I-V'*V|_F: " << detail::tostr(l)
		      << std::endl;

	// If ||I-C|| < 1, then the basis has no linear dependencies; the conditioning of the basis is
	// also related to that norm, but we choose an arbitrary close but smaller value than one.
	if (l < 0.7)
	  break;

	if (i >= max_its)
	  throw std::runtime_error("ortho: failing in orthonormalizing the basis");

	// W = W/chol(Wa'*W), where Wa'*W has dimensions (rows,cols)=(Wacorder,Wcorder)
	C = C.make_sure(none, none, OnEveryoneReplicated);
	cholInv<NW, Nt + 2 * Ncols, NW, COMPLEX>(std::move(C), Wacorder, Wcorder, Wa, Wacorder,
						 CopyTo, W);
	++i;
      }

      if (verb >= JustSummary)
	QDPIO::cout << prefix << " ortho summary rank: " << detail::volume(W.kvdim(), Wcorder)
		    << " #its: " << i << " |I-V'*V|_F: " << detail::tostr(l) << std::endl;
    }

    /// Orthonormalize W, W_out <- W_in*R, with R such that W_out'*W_out = I,
    /// for each combination of values of the order_t dimensions, or once if it is empty.
    /// \param W: matrix to orthogonalize
    /// \param order_t: dimension labels that do not participate in the orthogonalization
    /// \param order_rows: dimension labels that are the rows of the matrices V and W
    /// \param order_cols: dimension labels that are the columns of the matrices V and W

    template <std::size_t Nrows, std::size_t Ncols, std::size_t NW, typename COMPLEX>
    void ortho(Tensor<NW, COMPLEX> W, const std::string& order_t, const std::string& order_rows,
	       const std::string& order_cols, unsigned int max_its = 4, Verbosity verb = NoOutput,
	       const std::string& prefix = "")
    {
      ortho<Nrows, Ncols, NW, COMPLEX, NW>(Tensor<NW, COMPLEX>{}, W, order_t, order_rows,
					   order_cols, max_its, verb, prefix);
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
    void fgmres(const Operator<NOp, COMPLEX>& op, Maybe<Operator<NOp, COMPLEX>> prec,
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
      if (!op.d.is_compatible(x) || !op.d.is_compatible(y) ||
	  (prec && !prec.getSome().d.is_compatible(op.d)))
	throw std::runtime_error("Either the input or the output vector isn't compatible with the "
				 "operator or de the preconditioner");

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
	r = op.i.template like_this<NOp + 1>(
	  op.preferred_col_ordering == ColumnMajor ? std::string("%") + order_cols
						   : order_cols + std::string("%"),
	  '%', "", {{order_cols[0], x.kvdim().at(order_cols[0])}});

	r.set_zero();
	y.set_zero();
      }
      x.scale(-1).addTo(r);
      auto normr0 = norm<1>(r, op.order_t + order_cols); // type std::vector<real of T>
      if (max(normr0) == 0)
	return;

      // Allocate the search subspace U (onto the left singular space), and
      // Z (onto the right singular space)
      auto U = r.template like_this<NOp + 2>(op.preferred_col_ordering == ColumnMajor
					       ? std::string("%") + std::string(1, Vc)
					       : std::string(1, Vc) + std::string("%"),
					     '%', "", {{Vc, max_basis_size}});

      // Allocate the search subspace Z (onto the right singular space)
      auto Z = U.like_this();

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
	if (prec.hasSome())
	{
	  for (unsigned int i = 0; i < expansion_size; ++i)
	  {
	    // Z(:,ires+i) = prec * U(:,ires+i-1)
	    prec.getSome()(i == 0 ? r_Vc : U.kvslice_from_size({{Vc, ires + i - 1}}, {{Vc, 1}}),
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
	Uo = Uo.clone();
	ortho<NOp, 1>(Uo, op.order_t + order_cols, order_rows, std::string(1, Vac), 4, verb,
		      prefix);

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
	  auto rd = r.like_this();
	  op(y, rd); // rd = op(y)
	  nops += num_cols;
	  x.scale(-1).addTo(rd);
	  r.scale(-1).addTo(rd);
	  auto normrd = norm<1>(rd, op.order_t + order_cols);
	  double max_tol_d = 0;
	  for (int i = 0, vol = normr.volume(); i < vol; ++i)
	    max_tol_d = std::max(max_tol_d, (double)normrd.get({{i}}) / normr.get({{i}}));
	  QDPIO::cout << prefix
		      << " MGPROTON FGMRES error in residual vector: " << detail::tostr(max_tol_d)
		      << std::endl;
	}

	// Get the worse tolerance
	max_tol = 0;
	for (int i = 0, vol = normr.volume(); i < vol; ++i)
	  max_tol = std::max(max_tol, (double)normr.get({{i}}) / normr0.get({{i}}));

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
	max_tol = 0;
	for (int i = 0, vol = normr.volume(); i < vol; ++i)
	  max_tol = std::max(max_tol, (double)normr.get({{i}}) / normr0.get({{i}}));
	if (tol > 0 && max_tol > tol)
	  throw std::runtime_error("fgmres didn't converged and you ask for checking the error");
      }

      // Report iteration
      if (verb >= JustSummary)
	QDPIO::cout << prefix << " MGPROTON FGMRES summary #its.: " << it
		    << " max rel. residual: " << detail::tostr(max_tol, 2) << " matvecs: " << nops
		    << " precs: " << nprecs << std::endl;
    }

    template <std::size_t NOp, typename COMPLEX>
    Operator<NOp, COMPLEX> getSolver(const Operator<NOp, COMPLEX>& op, const Options& ops,
				     const Operator<NOp, COMPLEX>& prec = Operator<NOp, COMPLEX>(),
				     SolverSpace solverSpace = FullSpace);

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
	// Get preconditioner
	Maybe<Operator<NOp, COMPLEX>> prec = none;
	Maybe<const Options&> precOps = getOptionsMaybe(ops, "prec");
	if (precOps && prec_)
	  throw std::runtime_error(
	    "getFGMRESSolver: invalid `prec` tag: the solver got a preconditioner already");
	if (precOps)
	  prec = Maybe<Operator<NOp, COMPLEX>>(getSolver(op, precOps.getSome()));
	else if (prec_)
	  prec = Maybe<Operator<NOp, COMPLEX>>(prec_);

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
	std::string prefix = getOption<std::string>(ops, "prefix", "");

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  Tracker _t(std::string("fgmres ") + prefix);
		  _t.cost = x.kvdim().at('n');
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
		op.preferred_col_ordering};
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
	  return SB::Gamma(Ns * Ns - 1).template make_sure<COMPLEX>(none, dev);
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

      /// Return the natural lattice dimensions
      /// \param dim: dimension for each label
      /// \param layout: operator's layout

      inline std::map<char, int> getNatLatticeDims(const std::map<char, int>& dim,
						   OperatorLayout layout)
      {
	int nX = (layout == EvensOnlyLayout ? 2 : (dim.count('X') == 1 ? dim.at('X') : 1));
	return std::map<char, int>{
	  {'x', dim.at('x') * (dim.count('0') == 1 ? dim.at('0') : 1) * nX},
	  {'y', dim.at('y') * (dim.count('1') == 1 ? dim.at('1') : 1)},
	  {'z', dim.at('z') * (dim.count('2') == 1 ? dim.at('2') : 1)},
	  {'t', dim.at('t') * (dim.count('3') == 1 ? dim.at('3') : 1)}};
      }

      /// Return the neighbors as displacements from origin in natural coordinates
      /// \param blocking: blocking for each natural direction
      /// \param dim: operator dimensions
      /// \param neighbors: operator's neighbors in natural coordinates
      /// \param layout: operator's layout

      inline NaturalNeighbors
      getNeighborsAfterBlocking(const std::map<char, unsigned int>& blocking,
				const std::map<char, int>& dim, const NaturalNeighbors& neighbors,
				OperatorLayout layout)
      {
	using superbblas::detail::operator/;
	using superbblas::detail::operator+;

	// Get the natural dimensions of the lattice
	Coor<Nd> blk = kvcoors<Nd>("xyzt", blocking, 1);
	Coor<Nd> nat_dims = kvcoors<Nd>("xyzt", getNatLatticeDims(dim, layout));

	// Filter out odd neighbors if `Xsubrange` and block them
	std::set<Index> idx_neighbors;
	Coor<Nd> blk_strides = superbblas::detail::get_strides<Index>(blk, superbblas::FastToSlow);
	Coor<Nd> blk_nat_dims = nat_dims / blk;
	Coor<Nd> blk_nat_dims_strides =
	  superbblas::detail::get_strides<Index>(blk_nat_dims, superbblas::FastToSlow);
	std::size_t blk_vol = superbblas::detail::volume(blk);
	for (const auto& kvcoor : neighbors)
	{
	  Coor<Nd> c = kvcoors<Nd>("xyzt", kvcoor);
	  for (Index i = 0; i < blk_vol; ++i)
	  {
	    Coor<Nd> blk_c =
	      normalize_coor(c + superbblas::detail::index2coor(i, blk, blk_strides), nat_dims) /
	      blk;
	    Index idx_blk_c =
	      superbblas::detail::coor2index(blk_c, blk_nat_dims, blk_nat_dims_strides);
	    idx_neighbors.insert(idx_blk_c);
	  }
	}

	// Convert the indices into maps
	NaturalNeighbors r;
	for (Index idx : idx_neighbors)
	{
	  Coor<Nd> c = superbblas::detail::index2coor(idx, blk_nat_dims, blk_nat_dims_strides);
	  r.push_back(std::map<char, int>{{{'x', c[0]}, {'y', c[1]}, {'z', c[2]}, {'t', c[3]}}});
	}

	return r;
      }

      /// Return the Manhattan distance of the furthest neighbor
      /// \param neighbors: operator's neighbors in natural coordinates
      /// \param dim: operator dimensions
      /// \param layout: operator's layout

      inline unsigned int getFurthestNeighborDistance(const NaturalNeighbors& neighbors,
						      const std::map<char, int>& dim,
						      OperatorLayout layout)
      {
	const auto natdim = getNatLatticeDims(dim, layout);
	unsigned int max_distance = 0;
	for (const auto& kvcoor : neighbors)
	{
	  unsigned int dist = 0;
	  for (const auto& it : kvcoor)
	  {
	    int label_dim = natdim.at(it.first);
	    int label_coor = normalize_coor(it.second, label_dim);
	    dist += std::min(label_coor, label_dim - label_coor);
	  }
	  max_distance = std::max(max_distance, dist);
	}

	return max_distance;
      }

      /// Return the Manhattan distance of the furthest neighbor in an operator
      /// \param op: given operator

      template <std::size_t NOp, typename COMPLEX>
      unsigned int getFurthestNeighborDistance(Operator<NOp, COMPLEX> op)
      {
	return getFurthestNeighborDistance(op.neighbors, op.i.kvdim(), op.imgLayout);
      }

      /// Return the neighbors as displacements from origin in natural coordinates
      /// \param dim: operator dimensions
      /// \param max_dist_neighbors: the distance of the farthest neighbor for each site
      /// \param layout: operator's layout

      inline NaturalNeighbors getNeighbors(const std::map<char, int>& dim,
					   unsigned int max_dist_neighbors, OperatorLayout layout)
      {
	// Get the natural dimensions of the lattice and all the neighbors up to distance `max_dist_neighbors`
	Coor<Nd> nat_dims = kvcoors<Nd>("xyzt", getNatLatticeDims(dim, layout));
	std::vector<Coor<Nd>> neighbors = Coloring::all_neighbors(max_dist_neighbors, nat_dims);

	// Filter out odd neighbors if the layout is `EvensOnlyLayout`
	std::set<std::size_t> idx_neighbors;
	Stride<Nd> strides =
	  superbblas::detail::get_strides<std::size_t>(nat_dims, superbblas::FastToSlow);
	for (const auto& c : neighbors)
	{
	  if (layout == EvensOnlyLayout && std::accumulate(c.begin(), c.end(), Index{0}) % 2 != 0)
	    continue;
	  idx_neighbors.insert(superbblas::detail::coor2index(c, nat_dims, strides));
	}

	// Convert the indices into maps
	NaturalNeighbors r;
	for (std::size_t idx : idx_neighbors)
	{
	  Coor<Nd> c = superbblas::detail::index2coor(idx, nat_dims, strides);
	  r.push_back(std::map<char, int>{{'x', c[0]}, {'y', c[1]}, {'z', c[2]}, {'t', c[3]}});
	}

	return r;
      }

      /// Return the color for each site
      /// \param dim: operator dimensions
      /// \param layout: operator's layout
      /// \param neighbors: operator's neighbors in natural coordinates
      /// \param power: maximum distance to recover the nonzeros:
      ///               0, block diagonal; 1: near-neighbors...

      template <std::size_t N>
      std::pair<Tensor<N, float>, std::size_t>
      getColors(const std::map<char, int>& dim, OperatorLayout layout,
		const NaturalNeighbors& neighbors, unsigned int power)
      {
	// Unsupported other powers than zero or one
	unsigned int max_dist_neighbors = getFurthestNeighborDistance(neighbors, dim, layout);
	if (power != 0 && power != max_dist_neighbors)
	  throw std::runtime_error("getColors: unsupported value for `power`: either zero or the "
				   "distance to the furthest neighbor");

	// Compute the coloring
	Coor<Nd> nat_dims = kvcoors<Nd>("xyzt", getNatLatticeDims(dim, layout));
	Coloring coloring{0, // zero displacement in coloring
			  power == 0 ? max_dist_neighbors + 1
				     : max_dist_neighbors * 2 + 1, // k-distance coloring
			  nat_dims};

	// Create a field with the color of each site
	std::string order("xyztX");
	for (const auto& it : dim)
	  if (std::find(order.begin(), order.end(), it.first) == order.end())
	    order.push_back(it.first);
	auto real_dim = dim;
	real_dim['X'] = (layout == EvensOnlyLayout ? 2 : (dim.count('X') == 1 ? dim.at('X') : 1));
	auto t =
	  fillLatticeField<N, float>(order, {}, real_dim, real_dim, OnDefaultDevice,
				     [&](Coor<N - 1> c) {
				       return (float)coloring.getColor({{c[0], c[1], c[2], c[3]}});
				     })
	    .kvslice_from_size({}, {{'X', dim.at('X')}});

	return {t, coloring.numColors()};
      }

      /// Return a mask for the sites with even or odd x coordinate
      /// \param xoddity: 0 for even, 1 for odd x coordinates
      /// \param t: return a tensor with this distribution
      /// \param layout: tensor's layout
      /// NOTE: this implementation may be too slow for large tensors

      template <std::size_t N, typename T>
      Tensor<N, float> getXOddityMask_aux(int xoddity, const Tensor<N, T>& t, OperatorLayout layout)
      {
	if (xoddity != 0 && xoddity != 1)
	  throw std::runtime_error("getXOddityMask: invalid input argument `xoddity`");
	auto dim = t.kvdim();
	auto r = t.template make_compatible<N, float>();
	if (layout == NaturalLayout)
	{
	  if (dim.at('X') != 1 && dim.at('X') != 2)
	    throw std::runtime_error("getXOddityMask: invalid dimension size `X`");
	  Index xLabelPos = -1;
	  char xLabel = (dim.at('X') == 1 ? 'x' : 'X');
	  for (std::size_t i = 0; i < N; ++i)
	    if (r.order[i] == xLabel)
	      xLabelPos = i;
	  if (xLabelPos < 0)
	    throw std::runtime_error("getXOddityMask: given tensor doesn't have `x` dimension");
	  r.fillCpuFunCoor([&](const Coor<N>& coor) {
	    return coor[xLabelPos] % 2 == xoddity ? float{1} : float{0};
	  });
	}
	else if (layout == XEvenOddLayout)
	{
	  if (dim.at('X') != 2)
	    throw std::runtime_error("getXOddityMask: invalid dimension size `X`");
	  Coor<4> c;
	  for (std::size_t i = 0, j = 0; i < N; ++i)
	    if (is_in("Xyzt", r.order[i]))
	      c[j++] = i;
	  r.fillCpuFunCoor([&](const Coor<N>& coor) {
	    Index s = 0;
	    for (Index i : c)
	      s += coor[i];
	    return s % 2 == xoddity ? float{1} : float{0};
	  });
	}
	else if (layout == EvensOnlyLayout)
	{
	  if (dim.at('X') != 1)
	    throw std::runtime_error("getXOddityMask: invalid dimension size `X`");
	  r.set(xoddity == 0 ? float{1} : float{0});
	}
	else
	  throw std::runtime_error("getXOddityMask: unsupported layout");

	return r;
      }

      /// Return a mask for the sites with even or odd x coordinate
      /// \param xoddity: 0 for even, 1 for odd x coordinates
      /// \param t: return a tensor with this distribution
      /// \param layout: tensor's layout

      template <std::size_t N, typename T>
      Tensor<N, float> getXOddityMask(int xoddity, const Tensor<N, T>& t, OperatorLayout layout)
      {
	// Create the mask on the lattice components
	auto r_lat = t.template make_compatible<5, float>("Xxyzt");
	r_lat = getXOddityMask_aux(xoddity, r_lat, layout);

	// Create a matrix of ones to extend the mask onto the other components
	auto dim_dense = t.kvdim();
	for (auto& it : dim_dense)
	  if (is_in("xyztX", it.first))
	    it.second = 1;
	auto r_dense = t.template like_this<N - 5, float>("%", '%', "xyztX", dim_dense, none,
							  OnEveryoneReplicated);
	r_dense.set(1);

	// Contract both to create the output tensor
	auto r = t.template make_compatible<N, float>();
	contract(r_lat, r_dense, "", CopyTo, r);

	return r;
      }

      /// Return a copy of the given tensor in natural ordering into an even-odd ordering.
      ///
      /// \param v: origin tensor

      template <std::size_t N, typename T>
      Tensor<N, T> toEvenOddOrdering(const Tensor<N, T>& v)
      {
	// If the tensor is already in even-odd ordering, return it
	auto vdim = v.kvdim();
	if (vdim.at('X') == 2)
	  return v;

	// Check that the tensor can be reordered in even-ordering compressed on the x-direction
	if (!(vdim.at('x') % 2 == 0 && //
	      (vdim.at('y') == 1 || vdim.at('y') % 2 == 0) &&
	      (vdim.at('z') == 1 || vdim.at('z') % 2 == 0) &&
	      (vdim.at('t') == 1 || vdim.at('t') % 2 == 0)))
	  throw std::runtime_error("toEvenOddOrdering: invalid tensor dimensions");

	// All even/odd coordinate x elements cannot be selected with slicing at once for even-odd
	// ordering, so we use arbitrary selection of elements: masks. The approach to convert between
	// orderings is to mask all elements with even/odd x coordinate and copy them to a new tensor
	// with the target layout. The copying with mask superbblas operation wasn't design to support
	// different mask on the origin and destination tensor. But it's going to produce the desired
	// effect if the following properties match:
	//   a) the origin and destination masks are active in all dimensions excepting some dimensions,
	//      only X in this case;
	//   b) for all coordinates only one element is active on the excepting dimensions, the even or the odd
	//      x coordinates in the X dimension in this case;
	//   c) the excepting dimensions are fully supported on all processes; and
	//   d) the excepting dimensions are the fastest index and have the same ordering in the origin and
	//      the destination tensors.

	auto v0 = v.reshape_dimensions({{"Xx", "Xx"}}, {{'X', 2}}).reorder("Xxyzt%", '%');
	auto r = v0.make_compatible();
	for (int oddity = 0; oddity < 2; ++oddity)
	{
	  auto nat_mask = getXOddityMask<N>(oddity, v0, NaturalLayout);
	  auto eo_mask = getXOddityMask<N>(oddity, r, XEvenOddLayout);
	  v0.copyToWithMask(r, nat_mask, eo_mask);
	}
	return r;
      }

      /// Return a copy of the given tensor in even-odd ordering into a natural ordering.
      ///
      /// \param v: origin tensor

      template <std::size_t N, typename T>
      Tensor<N, T> toNaturalOrdering(const Tensor<N, T>& v)
      {
	// If the tensor is already in natural ordering, return it
	auto vdim = v.kvdim();
	if (vdim.at('X') == 1)
	  return v;

	// All even/odd coordinate x elements cannot be selected with slicing at once for even-odd
	// ordering, so we use arbitrary selection of elements: masks. The approach to convert between
	// orderings is to mask all elements with even/odd x coordinate and copy them to a new tensor
	// with the target layout. The copying with mask superbblas operation wasn't design to support
	// different mask on the origin and destination tensor. But it's going to produce the desired
	// effect if the following properties match:
	//   a) the origin and destination masks are active in all dimensions excepting some dimensions,
	//      only X in this case;
	//   b) for all coordinates only one element is active on the excepting dimensions, the even or the odd
	//      x coordinates in the X dimension in this case;
	//   c) the excepting dimensions are fully supported on all processes; and
	//   d) the excepting dimensions are the fastest index and have the same ordering in the origin and
	//      the destination tensors.

	auto v0 = v.reorder("Xxyzt%", '%');
	auto r = v0.make_compatible();
	for (int oddity = 0; oddity < 2; ++oddity)
	{
	  auto eo_mask = getXOddityMask<N>(oddity, v0, XEvenOddLayout);
	  auto nat_mask = getXOddityMask<N>(oddity, r, NaturalLayout);
	  v0.copyToWithMask(r, eo_mask, nat_mask);
	}
	return r.reshape_dimensions({{"Xx", "Xx"}}, {{'X', 1}, {'x', vdim.at('x') * 2}});
      }

      /// Return a sparse tensor with the content of the given operator
      /// \param op: operator to extract the nonzeros from
      /// \param power: maximum distance to recover the nonzeros:
      ///               0, block diagonal; 1: near-neighbors...
      /// \param coBlk: ordering of the nonzero blocks of the sparse operator
      /// \return: a pair of a sparse tensor and a remap; the sparse tensor has the same image
      ///          labels as the given operator and domain labels are indicated by the returned remap.

      template <std::size_t NOp, typename COMPLEX>
      std::pair<SpTensor<NOp, NOp, COMPLEX>, remap>
      cloneUnblockedOperatorToSpTensor(const Operator<NOp, COMPLEX>& op, unsigned int power = 1,
				       ColOrdering coBlk = RowMajor, const std::string& prefix = "")
      {
	log(1, "starting cloneUnblockedOperatorToSpTensor");

	Tracker _t(std::string("clone unblocked operator ") + prefix);

	// Unsupported explicitly colorized operators
	if (op.d.kvdim().count('X') == 0)
	  throw std::runtime_error(
	    "cloneUnblockedOperatorToSpTensor: unsupported not explicitly colored operators");

	// If the operator is empty, just return itself
	if (op.d.volume() == 0 || op.i.volume() == 0)
	  return {{}, {}};

	// TODO: add optimizations for multiple operators
	if (op.order_t.size() > 0)
	  throw std::runtime_error("Not implemented");

	// Create the ordering for the domain and the image where the dense dimensions indices run faster than the sparse dimensions
	// NOTE: assuming that x,y,z,t are the only sparse dimensions; X remains sparse
	std::string sparse_labels("xyztX");
	std::string dense_labels = remove_dimensions(op.i.order, sparse_labels);
	remap rd = getNewLabels(op.d.order, op.i.order + "u~0123");
	auto i = op.i.reorder("%xyztX", '%').like_this(none, {}, OnDefaultDevice).make_eg();

	// Get the blocking for the domain and the image
	std::map<char, int> blkd, blki;
	for (const auto& it : i.kvdim())
	  blki[it.first] = (is_in(dense_labels, it.first) ? it.second : 1);
	for (const auto& it : blki)
	  blkd[rd.at(it.first)] = it.second;

	// Construct the probing vectors, which they have as the rows the domain labels and as
	// columns the domain blocking dimensions

	constexpr int Nblk = NOp - Nd - 1;
	auto t_blk = identity<Nblk, COMPLEX>(blki, rd, dense_labels)
		       .make_sure(none, none, OnEveryoneReplicated);

	// Compute the coloring
	auto t_ = getColors<NOp>(i.kvdim(), op.imgLayout, op.neighbors, power);
	Tensor<NOp, float> colors = t_.first;
	unsigned int num_colors = t_.second;

	// The first half of the colors are for even nodes
	int maxX = op.i.kvdim().at('X');
	int real_maxX = (op.imgLayout == EvensOnlyLayout ? 2 : maxX);

	// Get the neighbors
	unsigned int max_dist_neighbors = getFurthestNeighborDistance(op);
	std::vector<Coor<Nd>> neighbors;
	if (power == 0)
	{
	  neighbors.push_back(Coor<Nd>{{}});
	}
	else if (power == max_dist_neighbors)
	{
	  for (const auto& it : op.neighbors)
	    neighbors.push_back(kvcoors<Nd>("xyzt", it, 0));
	}

	// Create masks for the elements with even natural x coordinate and with odd natural x coordinate
	auto even_x_mask = getXOddityMask<NOp>(0, i, op.imgLayout);
	auto odd_x_mask = getXOddityMask<NOp>(1, i, op.imgLayout);
	auto ones_blk = t_blk.template like_this<Nblk * 2, float>();
	ones_blk.set(1);

	// Create the sparse tensor
	auto d_sop =
	  (power == 0 ? i
		      : i.extend_support({{'x', (max_dist_neighbors + real_maxX - 1) / real_maxX},
					  {'y', max_dist_neighbors},
					  {'z', max_dist_neighbors},
					  {'t', max_dist_neighbors}}))
	    .rename_dims(rd);
	SpTensor<NOp, NOp, COMPLEX> sop{
	  d_sop, i, Nblk, Nblk, (unsigned int)neighbors.size(), coBlk == ColumnMajor};

	// Extract the nonzeros with probing
	sop.data.set_zero(); // all values may not be populated when using blocking
	for (unsigned int color = 0; color < num_colors; ++color)
	{
	  // Extracting the proving vector
	  auto t_l = colors.template transformWithCPUFun<COMPLEX>(
	    [&](float site_color) { return site_color == color ? COMPLEX{1} : COMPLEX{0}; });

	  // Skip empty masks
	  // NOTE: the call to split_dimension add a fake dimension that acts as columns
	  if (std::norm(norm<1>(t_l.split_dimension('X', "Xn", maxX), "n").get({{0}})) == 0)
	    continue;

	  // Contracting the proving vector with the blocking components
	  auto probs = contract<NOp + Nblk>(t_l, t_blk, "");

	  // Compute the matvecs
	  auto mv = op(std::move(probs));

	  // Construct an indicator tensor where all blocking dimensions but only the nodes colored `color` are copied
	  auto color_mask = t_l.template transformWithCPUFun<float>(
	    [](const COMPLEX& t) { return (float)std::real(t); });
	  auto sel_x_even =
	    contract<NOp + Nblk>(contract<NOp>(color_mask, even_x_mask, ""), ones_blk, "");
	  auto sel_x_odd =
	    contract<NOp + Nblk>(contract<NOp>(color_mask, odd_x_mask, ""), ones_blk, "");

	  // Populate the nonzeros by copying pieces from `mv` into sop.data. We want to copy only the
	  // nonzeros in `mv`, which are `neighbors` away from the nonzeros of `probs`.
	  latticeCopyToWithMask(mv, sop.data, 'u', neighbors, {{'X', real_maxX}}, sel_x_even,
				sel_x_odd);
	}

	// Populate the coordinate of the columns, that is, to give the domain coordinates of first nonzero in each
	// BSR nonzero block. Assume that we are processing nonzeros block for the image coordinate `c` on the
	// direction `dir`, that is, the domain coordinates will be (cx-dirx,cy-diry,cz-dirz,dt-dirt) in natural
	// coordinates. But we get the image coordinate `c` in even-odd coordinate, (cX,cx,cy,cz,ct), which has the
	// following natural coordinates (cx*2+(cX+cy+cz+ct)%2,cy,cz,ct). After subtracting the direction we get the
	// natural coordinates (cx*2+(cX+cy+cz+ct)%2-dirx,cy-diry,cz-dirz,ct-dirt), which corresponds to the following
	// even-odd coordinates ((cX-dirx-diry-dirz-dirt)%2,(cx*2+(cX+cy+cz+ct)%2-dirx)/2,cy-diry,cz-dirz,ct-dirt).

	Coor<Nd> real_dims = kvcoors<Nd>("xyzt", getNatLatticeDims(i.kvdim(), op.imgLayout));
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
	      return (base + sumdir + real_maxX * Nd) % real_maxX;
	    int sumyzt = std::accumulate(c.begin() + 2 + Nblk + 1, c.end() - 1, int{0});
	    const auto& cX = c[2 + Nblk + Nd];
	    return ((base * real_maxX + (cX + sumyzt) % real_maxX + real_dims[0] - dir[0]) /
		    real_maxX) %
		   (real_dims[0] / real_maxX);
	  }

	  int latd = domi - Nblk;
	  return (base - dir[latd] + real_dims[latd]) % real_dims[latd];
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
      /// \param coBlk: ordering of the nonzero blocks of the sparse operator
      /// \return: a pair of a sparse tensor and a remap; the sparse tensor has the same image
      ///          labels as the given operator and domain labels are indicated by the returned remap.

      template <std::size_t NOp, typename COMPLEX>
      std::pair<SpTensor<NOp, NOp, COMPLEX>, remap>
      cloneOperatorToSpTensor(const Operator<NOp, COMPLEX>& op, unsigned int power,
			      ColOrdering coBlk = RowMajor, const std::string& prefix = "")
      {
	Tracker _t(std::string("clone blocked operator ") + prefix);

	// Unblock the given operator, the code of `cloneUnblockedOperatorToSpTensor` is too complex as it is
	auto unblki = op.i
			.template reshape_dimensions<NOp - 4>(
			  {{"0x", "x"}, {"1y", "y"}, {"2z", "z"}, {"3t", "t"}}, {}, true)
			.make_eg();
	auto opdim = op.i.kvdim();
	Operator<NOp - 4, COMPLEX> unblocked_op{
	  [=](const Tensor<NOp - 4 + 1, COMPLEX>& x, Tensor<NOp - 4 + 1, COMPLEX> y) {
	    op(x.template reshape_dimensions<NOp + 1>(
		 {{"x", "0x"}, {"y", "1y"}, {"z", "2z"}, {"t", "3t"}}, opdim, true))
	      .template reshape_dimensions<NOp - 4 + 1>(
		{{"0x", "x"}, {"1y", "y"}, {"2z", "z"}, {"3t", "t"}}, {}, true)
	      .copyTo(y);
	  },
	  unblki,
	  unblki,
	  nullptr,
	  op.order_t,
	  op.domLayout,
	  op.imgLayout,
	  op.neighbors,
	  coBlk};

	// Get a sparse tensor representation of the operator
	unsigned int op_dist = getFurthestNeighborDistance(op);
	if (op_dist > 1 && power % op_dist != 0)
	  throw std::runtime_error("cloneOperatorToSpTensor: invalid power value, it isn't "
				   "divisible by the furthest neighbor distance");
	auto opdims = op.i.kvdim();
	auto t =
	  cloneUnblockedOperatorToSpTensor(unblocked_op, std::min(power, op_dist), coBlk, prefix);
	remap rd = t.second;
	for (const auto& it :
	     detail::getNewLabels("0123", op.d.order + op.i.order + "0123" + t.first.d.order))
	  rd[it.first] = it.second;
	int max_op_power = (op_dist == 0 ? 0 : std::max(power / op_dist, 1u) - 1u);
	std::map<char, int> m_power{};
	for (char c : std::string("xyzt") + update_order("xyzt", rd))
	  m_power[c] = max_op_power;
	auto sop = t.first.extend_support(m_power)
		     .split_dimension(rd['x'], update_order("0x", rd), opdim.at('0'), 'x', "0x",
				      opdim.at('0'))
		     .split_dimension(rd['y'], update_order("1y", rd), opdim.at('1'), 'y', "1y",
				      opdim.at('1'))
		     .split_dimension(rd['z'], update_order("2z", rd), opdim.at('2'), 'z', "2z",
				      opdim.at('2'))
		     .split_dimension(rd['t'], update_order("3t", rd), opdim.at('3'), 't', "3t",
				      opdim.at('3'))
		     .reorder(std::string("%") + update_order("0123xyztX", rd), "%0123xyztX", '%');

	return {sop, rd};
      }

      /// Return an efficient operator application
      /// \param op: operator to extract the nonzeros from
      /// \param power: maximum distance to recover the nonzeros:
      ///               0, block diagonal; 1: near-neighbors...
      /// \param co: preferred ordering of dense input and output tensors
      /// \param coBlk: ordering of the nonzero blocks of the sparse operator

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> cloneOperator(const Operator<NOp, COMPLEX>& op, unsigned int power,
					   ColOrdering co, ColOrdering coBlk,
					   const std::string& prefix = "")
      {
	// If the operator is empty, just return itself
	if (op.d.volume() == 0 || op.i.volume() == 0)
	  return op;

	// Get a sparse tensor representation of the operator
	unsigned int op_dist = getFurthestNeighborDistance(op);
	if (op_dist > 1 && power % op_dist != 0)
	  throw std::runtime_error("cloneOperator: invalid power value");
	auto t = cloneOperatorToSpTensor(op, power, coBlk, prefix);
	auto sop = t.first;
	remap rd = t.second;

	// Construct the operator to return
	Operator<NOp, COMPLEX> rop{sop,	       rd,	     power,	   sop.i,	 sop.i,
				   op.order_t, op.domLayout, op.imgLayout, op.neighbors, co};

	// Do a test
	Tracker _t(std::string("clone blocked operator (testing) ") + prefix);
	for (const auto& test_order : std::vector<std::string>{"%n", "n%"})
	{
	  auto x = op.d.template like_this<NOp + 1>(test_order, '%', "", {{'n', 2}});
	  auto y_rop = op.d.template like_this<NOp + 1>(test_order, '%', "", {{'n', 2}});
	  urand(x, -1, 1);
	  auto y_op = op(x);
	  for (int nfrom = 0; nfrom < 2; ++nfrom)
	  {
	    for (int nsize = 1; nsize <= 2; ++nsize)
	    {
	      y_rop.set(detail::NaN<COMPLEX>::get());
	      auto x0 = x.kvslice_from_size({{'n', nfrom}, {'n', nsize}});
	      auto y_op0 = y_op.kvslice_from_size({{'n', nfrom}, {'n', nsize}});
	      auto y_rop0 = y_rop.kvslice_from_size({{'n', nfrom}, {'n', nsize}});
	      auto base_norm0 = norm<1>(y_op0, "n");
	      rop(x0, y_rop0); // y_rop0 = rop(x0)
	      y_op0.scale(-1).addTo(y_rop0);
	      auto error = norm<1>(y_rop0, "n");
	      auto eps =
		std::sqrt(std::numeric_limits<typename real_type<COMPLEX>::type>::epsilon());
	      for (int i = 0; i < base_norm0.volume(); ++i)
		if (error.get({{i}}) > eps * base_norm0.get({{i}}))
		  throw std::runtime_error("cloneOperator: too much error on the cloned operator");
	    }
	  }
	}

	// Test for powers
	for (const auto& test_order : std::vector<std::string>{"%n^", "n%^"})
	{
	  const int max_power = 3;
	  auto x = op.d.template like_this<NOp + 2>(test_order, '%', "", {{'n', 2}, {'^', 1}});
	  auto y_op =
	    op.d.template like_this<NOp + 2>(test_order, '%', "", {{'n', 2}, {'^', max_power}});
	  urand(x, -1, 1);
	  op(x).copyTo(y_op.kvslice_from_size({}, {{'^', 1}}));
	  for (unsigned int i = 1; i < max_power; ++i)
	    op(y_op.kvslice_from_size({{'^', i - 1}}, {{'^', 1}}))
	      .copyTo(y_op.kvslice_from_size({{'^', i}}, {{'^', 1}}));
	  auto y_rop = y_op.like_this();
	  for (int nfrom = 0; nfrom < 2; ++nfrom)
	  {
	    for (int nsize = 1; nsize <= 2; ++nsize)
	    {
	      y_rop.set(detail::NaN<COMPLEX>::get());
	      auto x0 = x.kvslice_from_size({{'n', nfrom}, {'n', nsize}});
	      auto y_op0 = y_op.kvslice_from_size({{'n', nfrom}, {'n', nsize}});
	      auto y_rop0 = y_rop.kvslice_from_size({{'n', nfrom}, {'n', nsize}});
	      auto base_norm0 = norm<2>(y_op0, "n^").template collapse_dimensions<1>("n^", 'n');
	      rop(x0, y_rop0, '^'); // y_rop0 = {rop(x0), rop(rop(x0)), ...}
	      y_op0.scale(-1).addTo(y_rop0);
	      auto error = norm<2>(y_rop0, "n^").template collapse_dimensions<1>("n^", 'n');
	      auto eps =
		std::sqrt(std::numeric_limits<typename real_type<COMPLEX>::type>::epsilon());
	      for (int i = 0; i < base_norm0.volume(); ++i)
		if (error.get({{i}}) > eps * base_norm0.get({{i}}))
		  throw std::runtime_error(
		    "cloneOperator: too much error on the cloned operator for the power");
	    }
	  }
	}

	return rop;
      }

      /// Return an efficient operator application
      /// \param op: operator to extract the nonzeros from
      /// \param co: preferred ordering of dense input and output tensors
      /// \param coBlk: ordering of the nonzero blocks of the sparse operator

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> cloneOperator(const Operator<NOp, COMPLEX>& op, ColOrdering co,
					   ColOrdering coBlk, const std::string& prefix = "")
      {
	return cloneOperator(op, getFurthestNeighborDistance(op), co, coBlk, prefix);
      }

      /// Return the block diagonal of an operator
      /// \param op: operator to extract the block diagonal
      /// \param block_labels: labels that compose the blocks (and will be the rows)
      /// \param m: map from the operator block labels to the new labels for the columns
      /// \return: tensor with the block labels as rows and the renamed labels as columns

      template <std::size_t N, std::size_t NOp, typename COMPLEX>
      Tensor<N, COMPLEX> getBlockDiag(const Operator<NOp, COMPLEX>& op,
				      const std::string& block_labels, const remap& m)
      {
	// Check
	// Get a sparse tensor representation of the operator
	auto t = cloneOperatorToSpTensor(op, 0 /* clone only the block diagonal */, RowMajor,
					 "block diag");
	const auto sop = t.first;
	const remap rd = t.second;

	// Create tensor to return
	const std::string cols = update_order(block_labels, m);
	std::string order =
	  block_labels + cols + detail::remove_dimensions(op.d.order, block_labels);
	auto dims = op.d.kvdim();
	for (const auto& it : m)
	  dims[it.second] = dims[it.first];
	Tensor<N, COMPLEX> r = op.d.template like_this<N>(order, dims);

	// Copy the blocks
	remap m_cols = getNewLabels(cols, sop.data.order);
	remap m_blk;
	for (const auto& it : m)
	  m_blk[rd.at(it.first)] = m_cols.at(it.second);
	sop.data.rename_dims(m_blk).copyTo(r.rename_dims(m_cols));

	// Return tensor
	return r;
      }

      /// Returns the prolongator constructed from
      /// \param solvers: map of solvers
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver from `solvers` and influence the solver construction

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX>
      getMGProlongator(const Operator<NOp, COMPLEX>& op, unsigned int num_null_vecs,
		       const std::map<char, unsigned int>& mg_blocking,
		       const std::map<char, unsigned int>& layout_blocking,
		       bool do_chirality_splitting, const Operator<NOp, COMPLEX>& null_solver,
		       SolverSpace solverSpace)
      {
	detail::log(1, "starting getMGProlongator");

	// For now blocking on x should be divisible by X
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

	// Solve Ax=0 with random initial guesses
	auto b = op.d.template like_this<NOp + 1>("%n", '%', "", {{'n', num_null_vecs}});
	if (solverSpace == FullSpace)
	{
	  nrand(b);
	}
	else
	{
	  b.set_zero();
	  nrand(b.kvslice_from_size({}, {{'X', 1}}));
	}
	auto nv = null_solver(op(b));
	b.scale(-1).addTo(nv);

	// Do chirality splitting nv2 = [ nv * gpos, nv * gneg ]
	int ns = opdims.at('s');
	if (ns != 1 && ns != 2 && ns != Ns)
	  throw std::runtime_error("Error in getMGProlongator: Unsupported spin number");
	auto nv2 = nv;
	if (ns > 1 && do_chirality_splitting)
	{
	  nv2 = nv.like_this(none, {{'n', num_null_vecs * 2}});
	  auto g5 = getGamma5<COMPLEX>(ns, OnHost), g5pos = g5.cloneOn(OnHost),
	       g5neg = g5.cloneOn(OnHost);
	  for (int i = 0; i < ns; ++i) // make diagonal entries of gpos all positive or zero
	    g5pos.set({{i, i}}, g5.get({{i, i}}) + COMPLEX{1});
	  for (int i = 0; i < ns; ++i) // make diagonal entries of gneg all negative or zero
	    g5neg.set({{i, i}}, g5.get({{i, i}}) - COMPLEX{1});
	  nv2.kvslice_from_size({}, {{'n', num_null_vecs}})
	    .contract(g5pos, {{'i', 's'}}, NotConjugate, nv, {{'s', 'j'}}, NotConjugate);
	  nv2.kvslice_from_size({{'n', num_null_vecs}}, {{'n', num_null_vecs}})
	    .contract(g5neg, {{'i', 's'}}, NotConjugate, nv, {{'s', 'j'}}, NotConjugate);
	}

	// If the blocking in x isn't divisible by X, then work on natural ordering
	if (not x_blocking_divide_X)
	{
	  nv2 = toNaturalOrdering(nv2);
	  opdims['x'] *= X;
	  opdims['X'] = 1;
	}

	// Do the blocking, which encompasses the following transformations:
	//  X0x -> WX0x, 1y -> Y1y, 2z -> Z2z, 3t -> T3t,
	// where output X is a singlet dimension, and W,Y,Z, and T have size mg_blocking,
	// and output's 0,1,2, and 3 have size layout_blocking, and the output's x,y,z, and t
	// have the remaining

	auto nv_blk = nv2.template reshape_dimensions<NOp + 1 + 5>(
	  {{"X0x", "WX0x"},
	   {"1y", "Y1y"},
	   {"2z", "Z2z"},
	   {"3t", "T3t"},
	   {"n", "cs"},
	   {"c", "C"},
	   {"s", "S"}},
	  {{'X', 1}, // we don't do even-odd layout on the coarse operator space
	   {'W', mg_blocking.at('x') / (op.imgLayout == EvensOnlyLayout ? 2 : 1)},
	   {'Y', mg_blocking.at('y')},
	   {'Z', mg_blocking.at('z')},
	   {'T', mg_blocking.at('t')},
	   {'0', layout_blocking.at('x')},
	   {'1', layout_blocking.at('y')},
	   {'2', layout_blocking.at('z')},
	   {'3', layout_blocking.at('t')},
	   {'c', num_null_vecs},
	   {'s', nv2.kvdim()['n'] / num_null_vecs}},
	  true);

	// User even-odd ordering for nv_blk
	if (nv_blk.kvdim().at('0') != 1)
	  throw std::runtime_error("getMGProlongator: unsupported blocking on the x direction");

	// Do the orthogonalization on each block and chirality
	ortho<6, 1>(nv_blk, "X0123xyzts", "WYZTSC", "c", 4, JustSummary, "prolongator");

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
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  auto y0 = contract<NOp + 1 + 4>(nv_blk, toNaturalOrdering(x), "cs")
			      .template reshape_dimensions<NOp + 1>({{"WX0x", "X0x"},
								     {"Y1y", "1y"},
								     {"Z2z", "2z"},
								     {"T3t", "3t"},
								     {"C", "c"},
								     {"S", "s"}},
								    opdims, true);
		  if (not x_blocking_divide_X)
		    toEvenOddOrdering(y0).copyTo(y);
		  else
		    y0.copyTo(y);
		},
		d,
		i,
		[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  auto x0 = x;
		  if (not x_blocking_divide_X)
		    x0 = toNaturalOrdering(x);
		  auto x_blk = x0.template reshape_dimensions<NOp + 1 + 4>({{"X0x", "WX0x"},
									    {"1y", "Y1y"},
									    {"2z", "Z2z"},
									    {"3t", "T3t"},
									    {"c", "C"},
									    {"s", "S"}},
									   nv_blk.kvdim(), true);
		  if (nv_blk_eo_dim.at('X') == 1)
		    contract<NOp + 1>(nv_blk.conj(), x_blk, "WYZTSC", CopyTo, y);
		  else
		    toEvenOddOrdering(contract<NOp + 1>(nv_blk.conj(), x_blk, "WYZTSC")).copyTo(y);
		},
		op.order_t,
		nv_blk_eo_dim.at('X') == 2 ? XEvenOddLayout : NaturalLayout,
		op.imgLayout,
		getNeighborsAfterBlocking(mg_blocking, op.d.kvdim(), op.neighbors, op.imgLayout),
		op.preferred_col_ordering};
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

	// Get prolongator, V
	unsigned int num_null_vecs = getOption<unsigned int>(ops, "num_null_vecs");
	std::vector<unsigned int> mg_blocking_v = getVectorOption<unsigned int>(ops, "blocking");
	if (mg_blocking_v.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getMGPrec: the blocking should be a vector with four elements");
	std::map<char, unsigned int> mg_blocking{{'x', mg_blocking_v[0]},
						 {'y', mg_blocking_v[1]},
						 {'z', mg_blocking_v[2]},
						 {'t', mg_blocking_v[3]}};
	std::vector<unsigned int> layout_blocking_v = getVectorOption<unsigned int>(
	  ops, "coarse_layout_blocking", std::vector<unsigned int>{{1, 1, 1, 1}});
	if (layout_blocking_v.size() != Nd)
	  ops.getValue("coarse_layout_blocking")
	    .throw_error("getMGPrec: the blocking should be a vector with four elements");
	std::map<char, unsigned int> layout_blocking{{'x', layout_blocking_v[0]},
						     {'y', layout_blocking_v[1]},
						     {'z', layout_blocking_v[2]},
						     {'t', layout_blocking_v[3]}};
	const Operator<NOp, COMPLEX> nullSolver =
	  getSolver(op, getOptions(ops, "solver_null_vecs"));
	bool do_chirality_splitting = getOption<bool>(ops, "chirality_splitting", true);
	const Operator<NOp, COMPLEX> V =
	  getMGProlongator(op, num_null_vecs, mg_blocking, layout_blocking, do_chirality_splitting,
			   nullSolver, solverSpace);

	// Compute the coarse operator, either V' * op * V or V' * op * g5 * V
	unsigned int create_coarse_max_rhs =
	  getOption<unsigned int>(ops, "create_coarse_max_rhs", 0);
	ColOrdering co = getOption<ColOrdering>(ops, "operator_ordering", getColOrderingMap(),
						op.preferred_col_ordering);
	ColOrdering co_blk =
	  getOption<ColOrdering>(ops, "operator_block_ordering", getColOrderingMap(), RowMajor);
	int ns = op.d.kvdim().at('s');
	auto g5 = getGamma5<COMPLEX>(ns);
	const Operator<NOp, COMPLEX> op_c = cloneOperator(
	  Operator<NOp, COMPLEX>{
	    [&](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	      foreachInChuncks(x, y, create_coarse_max_rhs,
			       [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
				 if (do_chirality_splitting || ns == 1)
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
	    op.preferred_col_ordering},
	  co, co_blk, "coarse");

	// Get the solver for the projector
	const Operator<NOp, COMPLEX> coarseSolver =
	  getSolver(op_c, getOptions(ops, "solver_coarse"));

	// Get the solver for the smoother
	const Operator<NOp, COMPLEX> opSolver = getSolver(op, getOptions(ops, "solver_smoother"));

	// Return the solver
	return {
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
	    // x0 = g5 * x if !do_chirality_splitting && ns > 1
	    auto x0 = x;
	    if (!do_chirality_splitting && ns > 1)
	    {
	      x0 =
		contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), x, "s").rename_dims({{'i', 's'}});
	    }

	    // y0 = V*solver(V'*Op*V, V'x0)
	    auto y0 = V(coarseSolver(V.tconj()(x0)));

	    // x1 = x - op*y0
	    auto x1 = op(y0.scale(-1));
	    x.addTo(x1);

	    // y = y1 + solver(Op, x1)
	    opSolver(std::move(x1), y);
	    y0.addTo(y);
	  },
	  op.i,
	  op.d,
	  nullptr,
	  op.order_t,
	  op.imgLayout,
	  op.domLayout,
	  DenseOperator(),
	  op.preferred_col_ordering};
      }

      /// Destroy function
      using DestroyFun = std::function<void()>;

      /// Return a list of destroy callbacks after setting a solver

      inline std::vector<DestroyFun>& getEvenOddOperatorsCacheDestroyList()
      {
	static std::vector<DestroyFun> list;
	return list;
      }

      /// Call the destroy callbacks set up in `getEvenOddOperatorsCacheDestroyList`
      inline void cleanEvenOddOperatorsCache()
      {
	// Destroy in reverse order to the allocated one
	while (getEvenOddOperatorsCacheDestroyList().size() > 0)
	{
	  getEvenOddOperatorsCacheDestroyList().back()();
	  getEvenOddOperatorsCacheDestroyList().pop_back();
	}
      }

      /// Tuple storing the operator even-odd and odd-even parts and the block diagonal

      template <std::size_t NOp, typename COMPLEX>
      using EvenOddOperatorParts =
	std::tuple<Operator<NOp, COMPLEX>, Operator<NOp, COMPLEX>, Tensor<NOp + 2, COMPLEX>>;

      /// Return the cache of block diagonals for even-odd operators generated by getEvenOddPrec

      template <std::size_t NOp, typename COMPLEX>
      std::map<void*, EvenOddOperatorParts<NOp, COMPLEX>>& getEvenOddOperatorsPartsCache()
      {
	static std::map<void*, EvenOddOperatorParts<NOp, COMPLEX>> m = []() {
	  getEvenOddOperatorsCacheDestroyList().push_back(
	    []() { getEvenOddOperatorsPartsCache<NOp, COMPLEX>().clear(); });
	  return std::map<void*, EvenOddOperatorParts<NOp, COMPLEX>>{};
	}();
	return m;
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
      /// Notice that A_ee=Op_ee-Op_eo*Op_oo^{-1}*Op_oe is Op^{-1}_ee. So a way for preconditioning
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

	// Get the block diagonal of the operator with rows cs and columns CS
	remap m_sc{{'s', 'S'}, {'c', 'C'}};
	Tensor<NOp + 2, COMPLEX> opDiag;
	Operator<NOp, COMPLEX> op_eo, op_oe;
	if (getEvenOddOperatorsPartsCache<NOp, COMPLEX>().count(op.id.get()) == 0)
	{
	  opDiag = getBlockDiag<NOp + 2>(op, "cs", m_sc); // rows cs, cols CS
	  op_eo = op.kvslice_from_size({{'X', 1}}, {{'X', 1}}, {{'X', 0}}, {{'X', 1}});
	  op_oe = op.kvslice_from_size({{'X', 0}}, {{'X', 1}}, {{'X', 1}}, {{'X', 1}});
	  getEvenOddOperatorsPartsCache<NOp, COMPLEX>()[op.id.get()] = {op_eo, op_oe, opDiag};
	}
	else
	{
	  auto t = getEvenOddOperatorsPartsCache<NOp, COMPLEX>().at(op.id.get());
	  op_eo = std::get<0>(t);
	  op_oe = std::get<1>(t);
	  opDiag = std::get<2>(t);
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
		foreachInChuncks(x, y, create_operator_max_rhs,
				 [&](Tensor<NOp + 1, COMPLEX> x, Tensor<NOp + 1, COMPLEX> y) {
				   // y = x
				   x.copyTo(y);

				   // y0 = Op_ee^{-1} * x
				   auto y0 = solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
				     opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}), "cs", "CS",
				     x.rename_dims(m_sc), "CS");

				   // y1 = Op_oe * y0
				   auto y1 = op_oe(std::move(y0));

				   // y2 = Op_oo^{-1} * y1
				   auto y2 = solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
				     opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "cs", "CS",
				     std::move(y1).rename_dims(m_sc), "CS");

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
				   // y = Op_ee * x
				   contract(opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}),
					    x.rename_dims(m_sc), "CS", CopyTo, y);

				   // y1 = Op_oe * x
				   auto y1 = op_oe(x);

				   // y2 = Op_oo^{-1} * y1
				   auto y2 = solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
				     opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "cs", "CS",
				     std::move(y1).rename_dims(m_sc), "CS");

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
	  op.preferred_col_ordering};

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
		// y = Op_ee * prec_ee0(x)
		contract(opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}),
			 prec_ee0(x).rename_dims(m_sc), "CS", CopyTo, y);
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
	    // be = x_e - Op_eo*Op_oo^{-1}*x_o
	    Tensor<NOp + 1, COMPLEX> be;
	    if (solverSpace == FullSpace)
	    {
	      be = solver.d.template like_this<NOp + 1>(
		solver.preferred_col_ordering == ColumnMajor ? "%n" : "n%", '%', "",
		{{'n', x.kvdim().at('n')}});
	      x.kvslice_from_size({{'X', 0}}, {{'X', 1}}).copyTo(be);
	      op_eo(solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		      opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "cs", "CS",
		      x.kvslice_from_size({{'X', 1}}, {{'X', 1}}).rename_dims(m_sc), "CS"))
		.scale(-1)
		.addTo(be);
	    }
	    else
	    {
	      be = x.kvslice_from_size({{'X', 0}}, {{'X', 1}});
	    }

	    // Solve opA * y_e = be
	    auto ye = y.kvslice_from_size({{'X', 0}}, {{'X', 1}});
	    solver(be, ye);

	    // Do y_e = Op_ee^{-1} y_e if use_Aee_prec
	    if (use_Aee_prec)
	    {
	      solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		opDiag.kvslice_from_size({{'X', 0}}, {{'X', 1}}), "cs", "CS", ye.rename_dims(m_sc),
		"CS")
		.copyTo(ye);
	    }

	    // y_o = Op_oo^{-1}*(-Op_oe*y_e + x_o)
	    if (solverSpace == FullSpace)
	    {
	      auto yo0 = be;
	      x.kvslice_from_size({{'X', 1}}, {{'X', 1}}).copyTo(yo0);
	      op_oe(ye).scale(-1).addTo(yo0);
	      solve<2, NOp + 1, NOp + 2, NOp + 1, COMPLEX>(
		opDiag.kvslice_from_size({{'X', 1}}, {{'X', 1}}), "cs", "CS", yo0.rename_dims(m_sc),
		"CS", CopyTo, y.kvslice_from_size({{'X', 1}}, {{'X', 1}}));
	    }
	  },
	  op.i,
	  op.d,
	  nullptr,
	  op.order_t,
	  op.imgLayout,
	  op.domLayout,
	  DenseOperator(),
	  op.preferred_col_ordering};

	// Do a test
	if (superbblas::getDebugLevel() > 0)
	{
	  auto x = op.d.template like_this<NOp + 1>("%n", '%', "", {{'n', 2}});
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

	// Get the block diagonal of the operator with rows cs and columns CS
	const std::string blk_rows = "cs0123"; // order of the block of rows to invert
	remap m_blk = getNewLabels(blk_rows, op.d.order + op.i.order); // column labels
	const std::string blk_cols =
	  update_order(blk_rows, m_blk); // order of the block of columns to invert
	Tensor<NOp + 6, COMPLEX> opDiag = getBlockDiag<NOp + 6>(op, blk_rows, m_blk);

	// Return the solver
	return {[=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
		  // y = Op_diag^{-1} * x
		  solve<6, NOp + 1, NOp + 6, NOp + 1, COMPLEX>(
		    opDiag, blk_rows, blk_cols, x.rename_dims(m_blk), blk_cols, CopyTo, y);
		},
		op.d, op.i, nullptr, op};
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
	  getVectorOption<unsigned int>(ops, "blocking", default_blocking);
	if (blocking.size() != Nd)
	  ops.getValue("blocking")
	    .throw_error("getBlocking: the blocking should be a vector with four elements");
	std::map<char, unsigned int> mblk{
	  {'x', blocking[0]}, {'y', blocking[1]}, {'z', blocking[2]}, {'t', blocking[3]}};

	// Check that blocking is congruent with the lattice dimensions
	const auto nat_dim = getNatLatticeDims(op.d.kvdim(), op.domLayout);
	for (const auto& it : mblk)
	  if (nat_dim.at(it.first) % it.second != 0)
	    ops.getValue("blocking")
	      .throw_error("getBlocking: one of the given blocking does not divide the operator "
			   "dimension size");

	// Construct the blocked operator
	ColOrdering co = getOption<ColOrdering>(ops, "operator_ordering", getColOrderingMap(),
						op.preferred_col_ordering);
	ColOrdering co_blk =
	  getOption<ColOrdering>(ops, "operator_block_ordering", getColOrderingMap(), RowMajor);
	unsigned int power = getOption<unsigned int>(ops, "power", 1);

	int nX = (dim.count('X') == 0 ? 1 : dim.at('X'));
	auto blkd = op.d
		      .template reshape_dimensions<NOp>(
			{{"0x", "0x"}, {"1y", "1y"}, {"2z", "2z"}, {"3t", "3t"}},
			{{'0', std::max(nX, (int)mblk.at('x')) / nX},
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
	  getFurthestNeighborDistance(op) * power, co, co_blk, "blocking");

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

      // Auxiliary structure passed to PRIMME's matvec

      template <std::size_t NOp, typename COMPLEX>
      struct GDOperatorAux {
	const DeviceHost primme_dev; // where primme allocations are
	const Operator<NOp, COMPLEX>& op;

	// Wrapper for PRIMME of `LaplacianOperator`
	/// \param x: pointer to input vector
	/// \param ldx: leading dimension for `x`
	/// \param y: pointer to output vector
	/// \param ldy: leading dimension for `y`
	/// \param blockSize: number of input/output vectors
	/// \param ierr: output error state (zero means ok)

	static void GDPrimmeMatvec(void* x, PRIMME_INT* ldx, void* y, PRIMME_INT* ldy,
				   int* blockSize, primme_params* primme, int* ierr)
	{
	  *ierr = -1;
	  try
	  {
	    // The implementation assumes that ldx and ldy is nLocal
	    if (*blockSize > 1 && (*ldx != primme->nLocal || *ldy != primme->nLocal))
	      throw std::runtime_error("We cannot play with the leading dimensions");

	    GDOperatorAux<NOp, COMPLEX>& aux = *(GDOperatorAux<NOp, COMPLEX>*)primme->matrix;
	    auto dim = aux.op.d.kvdim();
	    dim['n'] = *blockSize;
	    std::string order = aux.op.d.order + "n";
	    Coor<NOp + 1> size = latticeSize<NOp + 1>(order, dim);
	    Tensor<NOp + 1, ComplexD> tx(order, size, aux.primme_dev, OnEveryone,
					 std::shared_ptr<ComplexD>((ComplexD*)x, [](ComplexD*) {}));
	    Tensor<NOp + 1, ComplexD> ty(order, size, aux.primme_dev, OnEveryone,
					 std::shared_ptr<ComplexD>((ComplexD*)y, [](ComplexD*) {}));

	    // ty = op * g5 * x
	    auto g5 = getGamma5<ComplexD>(dim['s'], aux.primme_dev);
	    aux.op(
	      contract<NOp + 1>(g5.rename_dims({{'j', 's'}}), tx, "s").rename_dims({{'i', 's'}}),
	      ty);

	    *ierr = 0;
	  } catch (...)
	  {
	  }
	}
      };

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

	// Get eigensolver properties
	unsigned int max_basis_size = getOption<unsigned int>(ops, "max_basis_size", 0);
	unsigned int max_block_size = getOption<unsigned int>(ops, "max_block_size", 1);
	Verbosity verb = getOption<Verbosity>(ops, "verbosity", getVerbosityMap(), NoOutput);
	double tol = getOption<double>(getOptions(ops, "solver"), "tol", 0.0) * 3.0;
	const Operator<NOp, COMPLEX> solver = getSolver(op, getOptions(ops, "solver"));

	// Return the solver
	return {
	  [=](const Tensor<NOp + 1, COMPLEX>& x, Tensor<NOp + 1, COMPLEX> y) {
#  if defined(SUPERBBLAS_USE_CUDA) && defined(BUILD_MAGMA)
	    DeviceHost primme_dev = OnDefaultDevice;
#  else
	    DeviceHost primme_dev = OnHost;
#  endif

	    // Create an auxiliary struct for the PRIMME's matvec
	    // NOTE: Please keep 'n' as the slowest index; the rows of vectors taken by PRIMME's matvec has dimensions 'cxyztX',
	    // and 'n' is the dimension for the columns.
	    GDOperatorAux<NOp, COMPLEX> opaux{primme_dev, solver};

	    // Make a bigger structure holding
	    primme_params primme;
	    primme_initialize(&primme);

	    // Primme solver setup
	    unsigned int numEvals = x.kvdim()['n'];
	    primme.numEvals = numEvals;
	    primme.printLevel =
	      (verb == NoOutput ? 0 : verb == JustSummary ? 1 : verb == Detailed ? 3 : 5);
	    primme.n = solver.d.volume();
	    primme.eps = tol;
	    primme.target = primme_largest_abs;
	    double zeros = 0;
	    primme.targetShifts = &zeros;
	    primme.numTargetShifts = 1;

	    // Set parallel settings
	    primme.nLocal = solver.d.getLocal().volume();
	    primme.numProcs = QDP::Layout::numNodes();
	    primme.procID = QDP::Layout::nodeNumber();
	    primme.globalSumReal = ns_getColorvecs::primmeGlobalSum;

	    // No preconditioner for my matrix
	    primme.matrixMatvec = opaux.GDPrimmeMatvec;
	    primme.matrix = &opaux;

	    // Set block size
	    primme.maxBasisSize = max_basis_size;
	    primme.maxBlockSize = max_block_size;
	    primme.ldOPs = primme.nLocal;

	    // Should set lots of defaults
	    if (primme_set_method(PRIMME_DEFAULT_MIN_TIME, &primme) < 0)
	    {
	      QDPIO::cerr << __func__ << ": invalid preset method\n";
	      QDP_abort(1);
	    }

	    // Allocate space for converged Ritz values and residual norms
	    std::vector<double> evals(primme.numEvals);
	    std::vector<double> rnorms(primme.numEvals);
	    Tensor<Nd + 4, ComplexD> evecs = solver.d.template like_this<Nd + 4, ComplexD>(
	      "%n", '%', "", {{'n', (int)numEvals}}, primme_dev, OnEveryone);
#  if defined(SUPERBBLAS_USE_CUDA) && defined(BUILD_MAGMA)
	    primme.queue = &*detail::getMagmaContext();
#  endif

	  // Call primme
#  if defined(SUPERBBLAS_USE_CUDA) && defined(BUILD_MAGMA)
	    int ret = magma_zprimme(evals.data(), evecs.data.get(), rnorms.data(), &primme);
#  else
	    int ret = zprimme(evals.data(), evecs.data.get(), rnorms.data(), &primme);
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

	    // Cleanup
	    primme_free(&primme);

	    // Copy evecs into y
	    evecs.copyTo(y);
	  },
	  solver.i,
	  solver.d,
	  nullptr,
	  solver.order_t,
	  solver.imgLayout,
	  solver.domLayout,
	  DenseOperator(),
	  solver.preferred_col_ordering};
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
	auto g5 = getGamma5<COMPLEX>(ns);

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
	auto g5 = getGamma5<COMPLEX>(ns);

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
		op.preferred_col_ordering};
      }

      /// Returns a solver with possible different precision than the operator's
      ///
      /// \param op: operator to make the inverse of
      /// \param ops: options to select the solver and the null-vectors creation

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX> getCasting(Operator<NOp, COMPLEX> op, const Options& ops,
					Operator<NOp, COMPLEX> prec_, SolverSpace solverSpace)
      {
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
      enum SolverType { FGMRES, MG, EO, BJ, IGD, DDAG, G5, BLOCKING, CASTING };
      static const std::map<std::string, SolverType> solverTypeMap{{"fgmres", FGMRES},
								   {"mg", MG},
								   {"eo", EO},
								   {"bj", BJ},
								   {"igd", IGD},
								   {"g5", G5},
								   {"blocking", BLOCKING},
								   {"casting", CASTING}};
      SolverType solverType = getOption<SolverType>(ops, "type", solverTypeMap);
      switch (solverType)
      {
      case FGMRES: // flexible GMRES
	return detail::getFGMRESSolver(op, ops, prec);
      case MG: // Multigrid
	return detail::getMGPrec(op, ops, prec, solverSpace);
      case EO: // even-odd Schur preconditioner
	return detail::getEvenOddPrec(op, ops, prec, solverSpace);
      case BJ: // block Jacobi
	return detail::getBlockJacobi(op, ops, prec);
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
    inline Operator<Nd + 7, Complex> asOperatorView(const LinearOperator<LatticeFermion>& linOp)
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
	  _t.cost = n;
	  auto tx = x.template reshape_dimensions<Nd + 4>(
	    {{"0x", "x"}, {"1y", "y"}, {"2z", "z"}, {"3t", "t"}}, {}, true);
	  auto ty = tx.like_this();
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
	ColumnMajor // preferred ordering
      };
    }

    //
    // High-level chroma operations
    //

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
      Maybe<Operator<Nd + 7, Complex>> op;

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
	  Operator<Nd + 7, Complex> linOp =
	    detail::cloneOperator(asOperatorView(*fLinOp), co, co_blk, "chroma's operator");

	  // Destroy chroma objects
	  delete fLinOp;
	  state = State();
	  S = Action();

	  // Construct the solver
	  op = Maybe<Operator<Nd + 7, Complex>>(getSolver(linOp, getOptions(*ops, "InvertParam")));

	  // Clean cache of operators
	  detail::cleanEvenOddOperatorsCache();

	  QDPIO::cout << "MGPROTON invertor ready; setup time: "
		      << detail::tostr(_t.stopAndGetElapsedTime()) << " s" << std::endl;
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
      Tensor<Nd + 5, COMPLEX_OUT> doInversion(const Operator<Nd + 7, COMPLEX_OUT>& op,
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
	auto aux = chi.template like_this<Nd + 8>(
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
