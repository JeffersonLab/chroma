// -*- C++ -*-
/*! \file                                                                    
 * \brief Copy, permuting, contracting tensors with superbblas
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_SUPERB_CONTRACTIONS__
#define __INCLUDE_SUPERB_CONTRACTIONS__

#include "chromabase.h"

#ifdef BUILD_SB

// Activate the MPI support in Superbblas
#  define SUPERBBLAS_USE_MPI
//#  define SUPERBBLAS_USE_MPIIO

#  include "actions/ferm/fermacts/fermact_factory_w.h"
#  include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#  include "meas/hadron/greedy_coloring.h"
#  include "meas/smear/link_smearing_factory.h"
#  include "qdp.h"
#  include "qdp_map_obj_disk_multiple.h"
#  include "superbblas.h"
#  include "util/ferm/key_timeslice_colorvec.h"
#  include <algorithm>
#  include <array>
#  include <chrono>
#  include <cmath>
#  include <cstring>
#  include <iomanip>
#  include <map>
#  include <memory>
#  include <random>
#  include <set>
#  include <sstream>
#  include <stdexcept>
#  include <string>
#  include <type_traits>

#  ifndef M_PI
#    define M_PI                                                                                   \
      3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068L
#  endif

#  ifdef BUILD_PRIMME
#    include <primme.h>
#  endif

namespace Chroma
{

  namespace SB
  {

    using Index = superbblas::IndexType;   ///< Default index type, `int` for now
    using Complex = std::complex<REAL>;	   ///< Default chroma complex precision
    using ComplexD = std::complex<REAL64>; ///< Complex double
    using ComplexF = std::complex<REAL32>; ///< Complex single

    /// Implicit complex type, there's a 2-size dimension for representing
    /// the real and the imaginary part, which usually has the label `.`
    template <typename T>
    struct DIYComplex {
      using value_type = T;
    };

    /// Type to represent coordinates
    template <std::size_t N>
    using Coor = superbblas::Coor<N>;

    /// Type to represent tensor layouts, its coordinate has the number of
    /// elements to jump to the element with the next coordinate
    template <std::size_t N>
    using Stride = superbblas::Coor<N, std::size_t>;

    /// Type of checksum, used by StorageTensor
    using checksum_type = superbblas::checksum_type;

    /// Where to store the tensor (see class Tensor)
    enum DeviceHost {
      OnHost,	      ///< on cpu memory
      OnDefaultDevice ///< on GPU memory if possible
    };

    /// How to distribute the tensor (see class Tensor)
    using Distribution = std::string;
    /// Fully supported on node with index zero
    static const Distribution OnMaster("__OnMaster__");
    /// Distribute the lattice dimensions (x, y, z, t)
    static const Distribution OnEveryone("tzyx");
    /// Distribute the lattice dimensions (x, y, z, t) as chroma does
    static const Distribution OnEveryoneAsChroma("__OnEveryonAsChroma__");
    /// All nodes have a copy of the tensor
    static const Distribution OnEveryoneReplicated("__OnEveryoneReplicated__");
    /// Local (single process) and non-collective
    static const Distribution Local("");
    /// Only the local process has support and non-collective 
    static const Distribution Glocal("__glocal__");

    /// Whether complex conjugate the elements before contraction (see Tensor::contract)
    enum Conjugation { NotConjugate, Conjugate };

    /// Whether the tensor is dense or sparse (see StorageTensor)
    enum Sparsity { Dense, Sparse };

    /// Whether to copy or add the values into the destination tensor (see Tensor::doAction)
    enum Action { CopyTo, AddTo };

    /// Auxiliary class for initialize Maybe<T> with no value
    struct None {
    };

    template <typename T, bool B = std::is_reference<T>::value>
    struct Maybe;

    template <typename T>
    struct is_maybe {
      static constexpr bool value = false;
    };

    template <>
    struct is_maybe<None> {
      static constexpr bool value = true;
    };

    template <typename T>
    struct is_maybe<Maybe<T>> {
      static constexpr bool value = true;
    };

    /// Class for optional values
    template <typename T>
    struct Maybe<T, false> {
      /// Whether the value is set
      bool has_value;

      /// The value
      T value;

      /// Constructor without a value
      Maybe() : has_value{false}, value{}
      {
      }

      /// Constructor without a value
      Maybe(None) : Maybe()
      {
      }

      /// Constructor with a value
      template <typename Q,
		typename std::enable_if<(!is_maybe<Q>::value && std::is_convertible<Q, T>::value),
					bool>::type = true>
      Maybe(const Q& t) : has_value{true}, value{T(t)}
      {
      }

      /// Return whether it has been initialized with a value
      bool hasSome() const
      {
	return has_value;
      }

      /// Return whether it has been initialized with a value
      explicit operator bool() const noexcept
      {
	return has_value;
      }

      /// Return the value if it has been initialized with some
      T getSome() const
      {
	if (has_value)
	  return value;
	throw std::runtime_error("Maybe::getSome: value isn't set");
      }

      /// Return the value if it has been initialized with some; otherwise return `def`
      T getSome(T def) const
      {
	if (has_value)
	  return getSome();
	else
	  return def;
      }
    };

    /// Class for optional values
    template <typename T>
    struct Maybe<T, true> {
      /// Whether the value is set
      bool has_value;

      /// The value type: change references by pointers
      using Tvalue = typename std::remove_reference<T>::type*;

      /// The value
      Tvalue value;

      /// Constructor without a value
      Maybe() : has_value{false}, value{}
      {
      }

      /// Constructor without a value
      Maybe(None) : Maybe()
      {
      }

      /// Constructor with a value
      Maybe(const T& t) : has_value{true}, value{&t}
      {
      }

      /// Return whether it has been initialized with a value
      bool hasSome() const
      {
	return has_value;
      }

      /// Return whether it has been initialized with a value
      explicit operator bool() const noexcept
      {
	return has_value;
      }

      /// Return the value if it has been initialized with some
      T getSome() const
      {
	if (has_value)
	  return *value;
	throw std::runtime_error("Maybe::getSome: value isn't set");
      }

      /// Return the value if it has been initialized with some; otherwise return `def`
      T getSome(T def) const
      {
	if (has_value)
	  return getSome();
	else
	  return def;
      }
    };

    /// Initialize Maybe<T> without value
    constexpr None none = None{};

    /// Return the number of seconds from some start
    inline double w_time()
    {
      return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch())
	.count();
    }

    /// Wrapper around superbblas time tracking
    struct Tracker : public superbblas::detail::tracker<superbblas::detail::Cpu> {
      Tracker(const std::string& funcName)
	: superbblas::detail::tracker<superbblas::detail::Cpu>(funcName, superbblas::detail::Cpu{},
							       true /* time it anyway */)
      {
      }
    };

    namespace detail
    {
      /// Throw an error if it is not a valid order, that is, if some label is repeated
      template <std::size_t N>
      void check_order(const std::string& order)
      {
	if (order.size() != N)
	{
	  std::stringstream ss;
	  ss << "The length of the dimension labels `" << order
	     << "` should match the template argument N `" << N << "`";
	  throw std::runtime_error(ss.str());
	}

	char s[256];
	for (unsigned int i = 0; i < 256; ++i)
	  s[i] = 0;
	for (unsigned int i = 0; i < N; ++i)
	{
	  if (s[order[i]] != 0)
	  {
	    std::stringstream ss;
	    ss << "Invalid order: some label names are repeated `" << order << "`";
	    throw std::runtime_error(ss.str());
	  }
	  s[order[i]] = 1;
	}
      }
    }

    enum Throw_kvcoors { NoThrow, ThrowOnUnmatchLabel, ThrowOnMissing };

    template <std::size_t N, typename INT = int>
    Coor<N> kvcoors(const std::string& order, const std::map<char, INT>& m, Index missing = 0,
		    Throw_kvcoors t = ThrowOnUnmatchLabel)
    {
      detail::check_order<N>(order);

      Coor<N> r;
      unsigned int found = 0;
      for (std::size_t i = 0; i < N; ++i)
      {
	auto it = m.find(order[i]);
	if (it != m.end())
	{
	  r[i] = it->second;
	  ++found;
	}
	else if (t == ThrowOnMissing)
	{
	  std::stringstream ss;
	  ss << "kvcoors: Missing value for label `" << order[i] << "` on dimension labels `"
	     << order << "`";
	  throw std::runtime_error(ss.str());
	}
	else
	{
	  r[i] = missing;
	}
      }

      if (found != m.size() && t == ThrowOnUnmatchLabel)
      {
	std::stringstream ss;
	ss << "kvcoors: Some dimension label on the given map m does not correspond to a dimension "
	      "label `"
	   << order << "`.";
	throw std::runtime_error(ss.str());
      }

      return r;
    }

    template <std::size_t N>
    Coor<N> latticeSize(const std::string& order, const std::map<char, int>& m = {})
    {
#  if QDP_USE_LEXICO_LAYOUT
      // No red-black ordering
      std::map<char, int> m0 = {{'x', Layout::lattSize()[0]},
				{'y', Layout::lattSize()[1]},
				{'z', Layout::lattSize()[2]},
				{'t', Layout::lattSize()[3]},
				{'X', 1},
				{'s', Ns},
				{'c', Nc},
				{'.', 2}};
#  elif QDP_USE_CB2_LAYOUT
      // Red-black ordering
      assert(Layout::lattSize()[0] % 2 == 0);
      std::map<char, int> m0 = {{'x', Layout::lattSize()[0] / 2},
				{'y', Layout::lattSize()[1]},
				{'z', Layout::lattSize()[2]},
				{'t', Layout::lattSize()[3]},
				{'X', 2},
				{'s', Ns},
				{'c', Nc},
				{'.', 2}};
#  else
      throw std::runtime_error("Unsupported layout");
#  endif
      for (const auto& it : m)
	m0[it.first] = it.second;
      return kvcoors<N>(order, m0, 0, NoThrow);
    }

    // Replace a label by another label
    using remap = std::map<char, char>;

    // Return the equivalent value of the coordinate `v` in the interval [0, dim[ for a periodic
    // dimension with length `dim`.

    inline int normalize_coor(int coor, int dim)
    {
      return (dim == 0 ? 0 : (coor + dim * (coor < 0 ? -coor / dim + 1 : 0)) % dim);
    }

    // Return the equivalent value of the coordinate `v` in the interval [0, dim[ for a periodic
    // dimension with length `dim`.

    template <std::size_t N>
    Coor<N> normalize_coor(Coor<N> v, Coor<N> dim)
    {
      Coor<N> r;
      for (std::size_t i = 0; i < N; ++i)
	r[i] = normalize_coor(v[i], dim[i]);
      return r;
    }

    namespace detail
    {
      using namespace superbblas::detail;

      /// Return whether a character is in a string
      /// \param s: the string
      /// \param c: the character

      inline bool is_in(const std::string& s, char c)
      {
	return std::find(s.begin(), s.end(), c) != s.end();
      }

      inline std::map<char, int> update_kvcoor(const std::map<char, int>& kvcoor, const remap& m)
      {
	std::map<char, int> r;
	for (auto const& it : kvcoor)
	  r[m.at(it.first)] = it.second;
	return r;
      }

      // Throw an error if `order` does not contain a label in `should_contain`
      inline void check_order_contains(const std::string& order, const std::string& should_contain)
      {
	for (char c : should_contain)
	{
	  if (order.find(c) == std::string::npos)
	  {
	    std::stringstream ss;
	    ss << "The input order `" << order << "` is missing the label `" << c << "`";
	    throw std::runtime_error(ss.str());
	  }
	}
      }

      /// Concatenate the two given orders without repeating characters and removing some other characters
      /// \param order0: input string
      /// \param order1: input string
      /// \param remove_dims: remove these characters

      inline std::string union_dimensions(const std::string& order0, const std::string& order1,
					  const std::string& remove_dims = "")
      {
	std::string out;
	out.reserve(order0.size() + order1.size());
	for (char c : order0)
	  if (out.find(c) == std::string::npos && remove_dims.find(c) == std::string::npos)
	    out.push_back(c);
	for (char c : order1)
	  if (out.find(c) == std::string::npos && remove_dims.find(c) == std::string::npos)
	    out.push_back(c);
	return out;
      }

      /// Remove some characters from a given string
      /// \param order: input string
      /// \param remove_dims: remove these characters from `order`

      inline std::string remove_dimensions(const std::string& order, const std::string& remove_dims)
      {
	return union_dimensions(order, "", remove_dims);
      }

      inline std::string update_order(std::string order, const remap& m)
      {
	for (std::size_t i = 0; i < order.size(); ++i)
	{
	  auto it = m.find(order[i]);
	  if (it != m.end())
	    order[i] = it->second;
	}
	return order;
      }

      template <std::size_t N>
      std::string update_order_and_check(std::string order, const remap& m)
      {
	std::string new_order = update_order(order, m);
	check_order<N>(new_order);
	return new_order;
      }

      template <std::size_t N>
      Coor<N - 1> remove_coor(Coor<N> v, std::size_t pos)
      {
	assert(pos < N);
	Coor<N - 1> r;
	for (std::size_t i = 0, j = 0; i < N; ++i)
	  if (i != pos)
	    r[j++] = v[i];
	return r;
      }

      inline std::string remove_coor(const std::string& v, std::size_t pos)
      {
	std::string r = v;
	r.erase(pos, 1);
	return r;
      }

      template <std::size_t N>
      Coor<N + 1> insert_coor(Coor<N> v, std::size_t pos, Index value)
      {
	assert(pos <= N);
	Coor<N + 1> r;
	for (std::size_t i = 0, j = 0; j < N + 1; ++j)
	{
	  if (j != pos)
	    r[j] = v[i++];
	  else
	    r[j] = value;
	}
	return r;
      }

      template <std::size_t N>
      Coor<N> replace_coor(Coor<N> v, std::size_t pos, Index value)
      {
	assert(pos <= N);
	v[pos] = value;
	return v;
      }

      inline std::string insert_coor(std::string v, std::size_t pos, char value)
      {
	assert(pos <= v.size());
	v.insert(pos, 1, value);
	return v;
      }

      inline std::string replace_coor(std::string v, std::size_t pos, char value)
      {
	assert(pos <= v.size());
	v[pos] = value;
	return v;
      }

      /// Return a character not given
      /// \param used_labels: labels not to use

      inline char get_free_label(const std::string& used_labels)
      {
	for (char c = '0'; true; ++c)
	{
	  if (used_labels.find(c) == std::string::npos)
	    return c;
	  if (c == std::numeric_limits<char>::max())
	    break;
	}
	throw std::runtime_error("get_free_labels: out of labels");
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

      /// Return a map to transform given labels into another ones
      /// \param labels: labels to remap
      /// \param used_labels: labels not to use, besides `labels`

      inline remap getNewLabels(const std::string& labels, std::string used_labels)
      {
	remap r;
	used_labels += labels;
	for (unsigned int i = 0; i < labels.size(); ++i)
	{
	  char c = get_free_label(used_labels);
	  r[labels[i]] = c;
	  used_labels.push_back(c);
	}
	return r;
      }

      // Return an array with an order, used by superbblas

      template <std::size_t N>
      std::array<char, N> to_sb_order(const std::string& order)
      {
	if (order.size() != N)
	  throw std::runtime_error(
	    "to_sb_order: the given string doesn't match the template parameter");
	std::array<char, N> c;
	std::copy_n(order.begin(), N, c.begin());
	return c;
      }

      /// Return the volume associated to an order
      /// \param m: size of each dimension
      /// \param labels: labels to consider

      inline std::size_t volume(const std::map<char, int>& m, const std::string& labels)
      {
	if (labels.size() == 0)
	  return 0;
	std::size_t vol = 1;
	for (char c : labels)
	  vol *= (std::size_t)m.at(c);
	return vol;
      }

      enum CoorType { From, Size };

      /// Return whether two from-size ranges are compatible on common dimensions
      /// \param o0: order for the first coordinates
      /// \param from0: first coordinate on the first range
      /// \param size0: range size on the first range
      /// \param o0: order for the second coordinates
      /// \param from1: first coordinate on the second range
      /// \param size1: range size on the second range
      /// \param labelsToCompare: labels to compare

      template <std::size_t N0, std::size_t N1>
      bool compatibleRanges(const std::string& o0, const Coor<N0>& from0, const Coor<N0>& size0,
			    const std::string& o1, const Coor<N1>& from1, const Coor<N1>& size1,
			    const std::string& labelsToCompare)
      {
	if (o0.size() != N0 || o1.size() != N1)
	  throw std::runtime_error("compatibleRanges: invalid size of input ordering");
	std::map<char, std::array<int, 2>> mfs0;
	for (unsigned int i = 0; i < N0; ++i)
	  if (std::find(labelsToCompare.begin(), labelsToCompare.end(), o0[i]) !=
	      labelsToCompare.end())
	    mfs0[o0[i]] = {{from0[i], size0[i]}};
	for (unsigned int i = 0; i < N1; ++i)
	{
	  if (mfs0.count(o1[i]) == 0)
	    continue;
	  auto fs0 = mfs0.at(o1[i]);
	  if (fs0[0] != from1[i] || fs0[1] != size1[i])
	    return false;
	}
	return true;
      }

      /// Coarse a range given a blocking
      /// \param fs: range to block
      /// \param blocking: blocking on each coordinate

      template <std::size_t N>
      std::array<Coor<N>, 2> coarse_range(const std::array<Coor<N>, 2>& fs, const Coor<N>& blocking)
      {
	std::array<Coor<N>, 2> r;
	for (unsigned int i = 0; i < N; ++i)
	{
	  r[0][i] = fs[0][i] / blocking[i] * blocking[i];
	  r[1][i] = (fs[0][i] + fs[1][i] + blocking[i] - 1) / blocking[i] * blocking[i] - r[0][i];
	}
	return r;
      }

      /// Split a dimension into another dimensions
      /// \param pos: dimension to split
      /// \param c: coordinate to transform
      /// \param new_dim: dimensions of the new tensor
      /// \param t: either `From` (first element) or `Size` (number of elements in each dimension)

      template <std::size_t Nout, std::size_t N>
      Coor<Nout> split_dimension(std::size_t pos, const Coor<N>& c, const Coor<Nout>& new_dim,
				 CoorType t)
      {
	constexpr std::size_t Nnew = Nout + 1 - N;
	Coor<Nout> r;
	std::copy_n(c.begin(), pos, r.begin());
	Index stride = 1;
	for (unsigned int k = 0; k < Nnew; ++k)
	{
	  if (!(t == From || c[pos] < stride || c[pos] % stride == 0))
	    throw std::runtime_error("split_dimension: Not supporting for this partition");
	  if (t == From)
	    r[pos + k] = (c[pos] / stride) % new_dim[pos + k];
	  else
	    r[pos + k] = std::min((c[pos] + stride - 1) / stride, new_dim[pos + k]);
	  stride *= new_dim[pos + k];
	}
	if (t == Size && new_dim[pos + Nnew - 1] != std::numeric_limits<int>::max() &&
	    c[pos] > stride)
	  throw std::runtime_error("split_dimension: dimension shorter than it should be");
	std::copy_n(c.begin() + pos + 1, N - pos - 1, r.begin() + pos + Nnew);
	return r;
      }

      /// Collapse several dimensions into another dimension
      /// \param pos: first dimension to collapse
      /// \param c: coordinate to transform
      /// \param old_dim: dimensions of the old tensor
      /// \param t: either `From` (first element) or `Size` (number of elements in each dimension)

      template <std::size_t Nout, std::size_t N>
      Coor<Nout> collapse_dimensions(std::size_t pos, const Coor<N>& c, const Coor<N>& old_dim,
				     CoorType t)
      {
	constexpr std::size_t Ncol = N + 1 - Nout; // number of dimensions to collapse
	Coor<Nout> r;
	std::copy_n(c.begin(), pos, r.begin());
	Index stride = 1, i = (t == From ? 0 : 1);
	bool odd_dim_watched = false;
	for (unsigned int k = 0; k < Ncol; ++k)
	{
	  if (t == Size && c[pos + k] > 0 && c[pos + k] != old_dim[pos + k])
	  {
	    if (odd_dim_watched)
	      throw std::runtime_error(
		"collapse_dimensions: unsupported to collapse a range with holes");
	    odd_dim_watched = true;
	  }
	  if (t == From)
	    i += c[pos + k] * stride;
	  else
	    i *= c[pos + k];
	  stride *= old_dim[pos + k];
	}
	r[pos] = i;
	std::copy_n(c.begin() + pos + Ncol, N - pos - Ncol, r.begin() + pos + 1);
	return r;
      }

      enum ReshapeDimensionsError {
	Success,		     ///< success
	CollapseRangeWithHoles,	     ///< unsupported to collapse a range with holes
	SizeNotDivisibleByPartition, ///< unsupported size for the partition
	NewDimensionIsTooShort	     ///< new dimension shorter than it should be
      };

      /// Reshape several dimensions into another dimension
      /// \param ncollapse: number of dimensions to collapse starting from each old dimension
      /// \param nsplit: number of dimensions to split starting from each old dimension
      /// \param old_dim: dimensions of the old tensor
      /// \param new_dim: maximum dimension size for the new tensor
      /// \param t: either `From` (first element) or `Size` (number of elements in each dimension)
      /// \param c: coordinate to transform

      template <std::size_t Nout, std::size_t N>
      std::pair<ReshapeDimensionsError, Coor<Nout>>
      reshape_dimensions(const Coor<N>& ncollapse, const Coor<N>& nsplit, const Coor<N>& old_dim,
			 const Coor<Nout>& new_dim, CoorType t, const Coor<N>& c)
      {
	Coor<Nout> r;
	unsigned int ri = 0;
	for (unsigned int ci = 0; ci < N; ++ci)
	{
	  if (ncollapse[ci] == 1 && nsplit[ci] == 1)
	  {
	    r[ri++] = c[ci];
	  }
	  else
	  {
	    // Collapse the dimensions from it[0] up to it[0]+it[1]-1
	    Index idx = 0;
	    {
	      Index stride = 1;
	      idx = (t == From ? 0 : 1);
	      bool odd_dim_watched = false;
	      for (unsigned int k = 0; k < ncollapse[ci]; ++k)
	      {
		if (t == Size && c[ci + k] > 0 && c[ci + k] != old_dim[ci + k])
		{
		  if (odd_dim_watched)
		    return {CollapseRangeWithHoles, {}};
		  odd_dim_watched = true;
		}
		if (t == From)
		  idx += c[ci + k] * stride;
		else
		  idx *= c[ci + k];
		stride *= old_dim[ci + k];
	      }
	    }

	    // Split the new dimension into it[2] new dimensions
	    {
	      Index stride = 1;
	      for (unsigned int k = 0; k < nsplit[ci]; ++k)
	      {
		if (!(t == From || idx < stride || idx % stride == 0))
		  return {SizeNotDivisibleByPartition, {}};
		if (t == From)
		  r[ri + k] = (idx / stride) % new_dim[ri + k];
		else
		  r[ri + k] = std::min((idx + stride - 1) / stride, new_dim[ri + k]);
		stride *= new_dim[ri + k];
	      }
	      if (t == Size && new_dim[ri + nsplit[ci] - 1] != std::numeric_limits<int>::max() &&
		  idx > stride)
		return {NewDimensionIsTooShort, {}};
	    }

	    ri += nsplit[ci];
	    ci += ncollapse[ci] - 1;
	  }
	}

	// Return the new coordinates
	return {Success, r};
      }

      /// Destroy function
      using DestroyFun = std::function<void()>;

      /// Return a list of destroy callbacks to execute just before finishing chroma

      inline std::vector<DestroyFun>& getDestroyList()
      {
	static std::vector<DestroyFun> list;
	return list;
      }

      // Get the cpu context
      inline std::shared_ptr<superbblas::Context>& getCpuContext()
      {
	static std::shared_ptr<superbblas::Context> cpuctx;
	if (!cpuctx)
	{
	  cpuctx = std::make_shared<superbblas::Context>(superbblas::createCpuContext());
	  getDestroyList().push_back([] { getCpuContext().reset(); });
	}
	return cpuctx;
      }

      // Return a context on either the host or the device
      inline std::shared_ptr<superbblas::Context>& getGpuContext()
      {
	// Creating GPU context can be expensive; so do it once
	static std::shared_ptr<superbblas::Context> cudactx;

#  ifdef SUPERBBLAS_USE_GPU
	if (!cudactx)
	{

	  int dev = -1;
#    if defined(QDP_IS_QDPJIT)
	  // When using QDP-JIT, the GPU device to use is already selected
#      ifdef SUPERBBLAS_USE_CUDA
	  superbblas::detail::gpuCheck(cudaGetDevice(&dev));
#      elif defined(SUPERBBLAS_USE_HIP)
	  superbblas::detail::gpuCheck(hipGetDevice(&dev));
#      else
#	error unsupported GPU platform
#      endif
#    else
	  // When not using QDP-JIT, select the GPU device based on either the local
	  // MPI rank or the global MPI rank and assuming that consecutive MPI ranks
	  // tends to be on the same node.
	  const char* l = std::getenv("SB_NUM_GPUS_ON_NODE");
	  if (l)
	  {
	    dev = Layout::nodeNumber() % std::atoi(l);
	  }
	  else
	  {
	    const char* l = std::getenv("SLURM_LOCALID");
	    if (l)
	    {
	      dev = std::atoi(l);
	    }
	    else
	    {
	      QDPIO::cerr << "Please set SB_NUM_GPUS_ON_NODE or SLURM_LOCALID" << std::endl;
	      QDP_abort(1);
	    }
	  }
#    endif

	  // Workaround on a potential issue in qdp-jit: avoid passing through the pool allocator
#    if defined(QDP_IS_QDPJIT)
	  if (jit_config_get_max_allocation() != 0)
	  {
	    // Make superbblas use the same memory allocator for gpu as any other qdp-jit lattice object
	    superbblas::getCustomAllocator() = [](std::size_t size,
						  superbblas::platform plat) -> void* {
	      if (size == 0)
		return nullptr;
	      if (plat == superbblas::CPU)
		return malloc(size);
	      void* ptr = nullptr;
	      QDP_get_global_cache().addDeviceStatic(&ptr, size, true);
	      assert(superbblas::detail::getPtrDevice(ptr) >= 0);
	      return ptr;
	    };

	    // The corresponding deallocator
	    superbblas::getCustomDeallocator() = [](void* ptr, superbblas::platform plat) {
	      if (ptr == nullptr)
		return;
	      if (plat == superbblas::CPU)
		free(ptr);
	      else
		QDP_get_global_cache().signoffViaPtr(ptr);
	    };
	  }
#    endif // defined(QDP_IS_QDPJIT)
	  cudactx = std::make_shared<superbblas::Context>(superbblas::createGpuContext(dev));
	  getDestroyList().push_back([] { getGpuContext().reset(); });
	}
	return cudactx;
#  else	 // SUPERBBLAS_USE_GPU
	return getCpuContext();
#  endif // SUPERBBLAS_USE_GPU
      }

      // Return a context on either the host or the device
      inline const std::shared_ptr<superbblas::Context>& getContext(DeviceHost dev)
      {
	return dev == OnHost ? getCpuContext() : getGpuContext();
      }

      /// Return if two devices are the same

      inline bool is_same(DeviceHost a, DeviceHost b)
      {
#  ifdef SUPERBBLAS_USE_GPU
	return a == b;
#  else
	// Without gpus, OnHost and OnDefaultDevice means on cpu.
	return true;
#  endif
      }

      /// is_complex<T>::value is true if `T` is complex

      template <typename T>
      struct is_complex : std::false_type {
      };

      template <typename T>
      struct is_complex<std::complex<T>> : std::true_type {
      };

      /// real_type<T>::type is T::value_type if T is complex or DIYComplex; otherwise it is T

      template <typename T>
      struct real_type {
	using type = T;
      };

      template <typename T>
      struct real_type<std::complex<T>> {
	using type = T;
      };

      template <typename T>
      struct real_type<DIYComplex<T>> {
	using type = T;
      };

      /// is_diycomplex<T>::value is true if `T` is DIYComplex

      template <typename T>
      struct is_diycomplex : std::false_type {
      };

      template <typename T>
      struct is_diycomplex<DIYComplex<T>> : std::true_type {
      };

      /// base_type<T>::type is T excepting for base_type<DIYComplex<T>>::type that is T

      template <typename T>
      struct base_type {
	using type = T;
      };

      template <typename T>
      struct base_type<DIYComplex<T>> {
	using type = T;
      };

      /// Return x if conjugate is false and conj(x) otherwise
      /// \param x: value to conjugate

      template <typename T, typename std::enable_if<is_complex<T>::value, bool>::type = true>
      T cond_conj(bool conjugate, const T& x)
      {
	return (conjugate ? std::conj(x) : x);
      }

      template <typename T, typename std::enable_if<!is_complex<T>::value, bool>::type = true>
      T cond_conj(bool, const T& x)
      {
	return x;
      }

      /// Return an ordering with labels 0, 1, ...
      inline std::string getTrivialOrder(std::size_t N)
      {
	std::string r(N, 0);
	for (std::size_t i = 0; i < N; ++i)
	  r[i] = (i + 1) % 128;
	return r;
      }

      /// Return a map with pairs from elements withdraw from two lists
      /// ita_begin, first element for the first pair
      /// ita_end, first element not to include
      /// itb_begin, second element for the first pair

      template <typename ITA_BEGIN, typename ITA_END, typename ITB_BEGIN>
      std::map<char, int> zip(ITA_BEGIN ita_begin, ITA_END ita_end, ITB_BEGIN itb_begin)
      {
	std::map<char, int> m;
	while (ita_begin != ita_end)
	{
	  m[*ita_begin] = *itb_begin;
	  ita_begin++;
	  itb_begin++;
	}
	return m;
      }

      /// Return whether the tensor isn't local or on master or replicated
      inline bool isDistributedOnEveryone(const Distribution& dist)
      {
	return dist != OnMaster && dist != OnEveryoneReplicated && dist != Local && dist != Glocal;
      }

      /// Stores the subtensor supported on each node (used by class Tensor)
      template <std::size_t N>
      struct TensorPartition {
      public:
	using PartitionStored = std::vector<superbblas::PartitionItem<N>>;
	Coor<N> dim;	   ///< Dimensions of the tensor
	PartitionStored p; ///< p[i] = {first coordinate, size} of tensor on i-th node
	bool isLocal;	   ///< Whether the partition is non-collective

	/// Constructor
	/// \param order: dimension labels (use x, y, z, t for lattice dimensions)
	/// \param dim: dimension size for the tensor
	/// \param dist: how to distribute the tensor among the nodes

	TensorPartition(const std::string& order, Coor<N> dim, Distribution dist) : dim(dim)
	{
	  detail::check_order<N>(order);
	  isLocal = false;
	  if (dist == OnMaster)
	  {
	    p = all_tensor_on_master(dim);
	  }
	  else if (dist == OnEveryoneReplicated)
	  {
	    p = all_tensor_replicated(dim);
	  }
	  else if (dist == OnEveryoneAsChroma)
	  {
	    p = partitioning_chroma_compatible(order, dim);
	  }
	  else if (dist == Local)
	  {
	    p = local(dim);
	    isLocal = true;
	  }
	  else if (dist == Glocal)
	  {
	    throw std::runtime_error("TensorPartition: unsupported distribution");
	  }
	  else
	  {
	    p = partitioning_distributed(order, dim, dist);
	  }
	}

	/// Constructor for `insert_dimension`
	/// \param dim: dimension size for the tensor
	/// \param p: partition
	/// \praam isLocal: whether the tensor is local

	TensorPartition(Coor<N> dim, const PartitionStored& p, bool isLocal)
	  : dim(dim), p(p), isLocal(isLocal)
	{
	}

	/// Empty constructor
	TensorPartition() : dim{}, p{}, isLocal{false}
	{
	}

	/// Return the volume of the tensor supported on this node
	std::size_t localVolume() const
	{
	  return superbblas::detail::volume(p[MpiProcRank()][1]);
	}

	/// Return the first coordinate store locally
	Coor<N> localFrom() const
	{
	  return p[MpiProcRank()][0];
	}

	/// Return the number of elements store locally in each dimension
	Coor<N> localSize() const
	{
	  return p[MpiProcRank()][1];
	}

	/// Return how many processes have support on this tensor
	/// Note that it may differ from MPI's numProcs if the tensor does not have support on all processes.

	unsigned int numProcs() const
	{
	  unsigned int numprocs = 0;
	  for (const auto& i : p)
	    if (superbblas::detail::volume(i[1]) > 0)
	      numprocs++;
	  return numprocs;
	}

	/// Return the process rank on this tensor or -1 if this process does not have support on the tensor
	/// Note that it may differ from MPI's rank if the tensor does not have support on all processes.

	int procRank() const
	{
	  // Return -1 if this process does not have support on the tensor
	  int mpi_rank = MpiProcRank();
	  if (superbblas::detail::volume(p[mpi_rank][1]) == 0)
	    return -1;

	  // Return as rank how many processes with MPI rank smaller than this process have support
	  int this_rank = 0;
	  for (int i = 0; i < mpi_rank; ++i)
	    if (superbblas::detail::volume(p[i][1]) > 0)
	      this_rank++;
	  return this_rank;
	}

	/// Return the MPI process rank
	int MpiProcRank() const
	{
	  return (isLocal ? 0 : Layout::nodeNumber());
	}

	/// Return the maximum local volume supported supported by a process

	std::size_t maxLocalVolume() const
	{
	  std::size_t maxLocalVol = 0;
	  for (const auto& it : p)
	    maxLocalVol = std::max(maxLocalVol, superbblas::detail::volume(it[1]));
	  return maxLocalVol;
	}

	/// Return whether other partition is compatible with this one

	template <std::size_t N0>
	bool is_compatible(const std::string& o0, const TensorPartition<N0>& t,
			   const std::string& o1, const std::string& labelToCompare) const
	{
	  if (t.p.size() != p.size())
	    return false;
	  for (unsigned int i = 0; i < p.size(); ++i)
	    if (!compatibleRanges(o0, p[i][0], p[i][1], o1, t.p[i][0], t.p[i][1], labelToCompare))
	      return false;

	  return true;
	}

	/// Make a partition compatible with a given one

	template <std::size_t N1>
	TensorPartition<N1> make_compatible(const std::string& o, const std::string& o1,
					    const Coor<N1>& new_dim,
					    const std::string& labelsToCompare) const
	{
	  if (o.size() != N || o1.size() != N1)
	    throw std::runtime_error(
	      "make_compatible: one the given orders does not match the expected length");

	  typename TensorPartition<N1>::PartitionStored r(
	    p.size(), std::array<Coor<N1>, 2>{Coor<N1>{{}}, new_dim});
	  for (unsigned int i = 0; i < p.size(); ++i)
	  {
	    for (unsigned int j = 0; j < N1; ++j)
	    {
	      if (std::find(labelsToCompare.begin(), labelsToCompare.end(), o1[j]) !=
		  labelsToCompare.end())
	      {
		int j0 = std::find(o.begin(), o.end(), o1[j]) - o.begin();
		r[i][0][j] = p[i][0][j0];
		r[i][1][j] = p[i][1][j0];
	      }
	      if (superbblas::detail::volume(r[i][1]) == 0)
		r[i] = std::array<Coor<N1>, 2>{Coor<N1>{{}}, Coor<N1>{{}}};
	    }
	  }
	  TensorPartition<N1> new_t{new_dim, r, isLocal};

	  // TODO: Fusing different tensor partitions isn't trivial. If both partitions differ in the number of
	  // active processes (processes with nonzero support), then the current implementation will fail
	  if (!is_compatible(o, new_t, o1, labelsToCompare))
	    throw std::runtime_error("make_compatible is broken and you hit a corner case");

	  return new_t;
	}

	/// Insert a new non-distributed dimension

	TensorPartition<N + 1> insert_dimension(std::size_t pos, std::size_t dim_size) const
	{
	  typename TensorPartition<N + 1>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back({insert_coor(i[0], pos, 0), insert_coor(i[1], pos, dim_size)});
	  return TensorPartition<N + 1>{insert_coor(dim, pos, dim_size), r, isLocal};
	}

	/// Remove a non-distributed dimension

	TensorPartition<N - 1> remove_dimension(std::size_t pos) const
	{
	  typename TensorPartition<N - 1>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back({remove_coor(i[0], pos), remove_coor(i[1], pos)});
	  return TensorPartition<N - 1>{remove_coor(dim, pos), r, isLocal};
	}

	/// Split a dimension into a non-distributed dimension and another dimension

	template <std::size_t Nout, typename std::enable_if<(N > 0), bool>::type = true>
	TensorPartition<Nout> split_dimension(std::size_t pos, const Coor<Nout>& new_dim) const
	{
	  typename TensorPartition<Nout>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back({detail::split_dimension(pos, i[0], new_dim, From),
			 detail::split_dimension(pos, i[1], new_dim, Size)});
	  return TensorPartition<Nout>{detail::split_dimension(pos, dim, new_dim, Size), r,
				       isLocal};
	}

	/// Coarse the ranges on each process

	template <typename std::enable_if<(N > 0), bool>::type = true>
	TensorPartition<N> coarse_support(const Coor<N>& blocking) const
	{
	  typename TensorPartition<N>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back(detail::coarse_range(i, blocking));
	  return TensorPartition<N>{dim, r, isLocal};
	}

	/// Collapse several dimensions into another dimension

	template <std::size_t Nout, typename std::enable_if<(Nout > 0), bool>::type = true>
	TensorPartition<Nout> collapse_dimensions(std::size_t pos) const
	{
	  typename TensorPartition<Nout>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back({detail::collapse_dimensions<Nout>(pos, i[0], dim, From),
			 detail::collapse_dimensions<Nout>(pos, i[1], dim, Size)});
	  return TensorPartition<Nout>{detail::collapse_dimensions<Nout>(pos, dim, dim, Size), r,
				       isLocal};
	}

	/// Reshape several dimensions into another dimension
	/// \param ncollapse: number of dimensions to collapse starting from each old dimension
	/// \param nsplit: number of dimensions to split starting from each old dimension
	/// \param new_dim: maximum dimension size for the new tensor

	template <std::size_t Nout, typename std::enable_if<(N > 0 && Nout > 0), bool>::type = true>
	std::pair<ReshapeDimensionsError, TensorPartition<Nout>>
	reshape_dimensions(const Coor<N>& ncollapse, const Coor<N>& nsplit,
			   const Coor<Nout>& new_dim) const
	{
	  auto new_dim_aux =
	    detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, new_dim, Size, dim);
	  if (new_dim_aux.first != Success)
	    return {new_dim_aux.first, {}};
	  typename TensorPartition<Nout>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	  {
	    auto new_from =
	      detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, new_dim, From, i[0]);
	    auto new_size =
	      detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, new_dim, Size, i[1]);
	    if (new_from.first != Success)
	      return {new_from.first, {}};
	    if (new_size.first != Success)
	      return {new_size.first, {}};
	    r.push_back({new_from.second, new_size.second});
	  }
	  return {Success, TensorPartition<Nout>{new_dim_aux.second, r, isLocal}};
	}

	/// Extend the support of distributed dimensions by one step in each direction

	TensorPartition<N> extend_support(Coor<N> m) const
	{
	  typename TensorPartition<N>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	  {
	    superbblas::PartitionItem<N> fs;
	    for (unsigned int j = 0; j < N; ++j)
	    {
	      fs[1][j] = std::min(i[1][j] + 2 * m[j], dim[j]);
	      fs[0][j] = (fs[1][j] < dim[j] ? (i[0][j] - m[j] + dim[j]) % dim[j] : 0);
	    }
	    r.push_back(fs);
	  }
	  return TensorPartition<N>{dim, r, isLocal};
	}

	/// Return a subpartition given a range

	TensorPartition<N> get_subpartition(const Coor<N>& from, const Coor<N>& size) const
	{
	  typename TensorPartition<N>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	  {
	    Coor<N> lfrom, lsize;
	    superbblas::detail::intersection(i[0], i[1], from, size, dim, lfrom, lsize);
	    r.push_back(superbblas::detail::volume(lsize) == 0
			  ? std::array<Coor<N>, 2>{Coor<N>{{}}, Coor<N>{{}}}
			  : std::array<Coor<N>, 2>{normalize_coor(lfrom - from, size), lsize});
	  }
	  return TensorPartition<N>{size, r, isLocal};
	}

	/// Return a partition with the local portion of the tensor

	TensorPartition<N> get_local_partition() const
	{
	  return TensorPartition<N>{
	    localSize(), PartitionStored(1, superbblas::PartitionItem<N>{{{}, localSize()}}), true};
	}

	/// Return a partition with the local portion of the tensor

	TensorPartition<N> get_glocal_partition() const
	{
	  PartitionStored r(p.size());
	  r[MpiProcRank()] = p[MpiProcRank()];
	  return TensorPartition<N>{dim, r, isLocal};
	}

	/// Return a copy of this tensor with a compatible distribution to be contracted with the given tensor
	/// \param order: labels for this distribution
	/// \param t: given tensor distribution
	/// \param ordert: labels for the given distribution

	template <std::size_t Nt>
	TensorPartition<N> make_suitable_for_contraction(const std::string& order,
							 const TensorPartition<Nt>& t,
							 const std::string& ordert) const
	{
	  PartitionStored r(p.size());
	  std::map<char, Index> mf, ms;
	  for (std::size_t i = 0; i < N; ++i)
	    mf[order[i]] = 0;
	  for (std::size_t i = 0; i < N; ++i)
	    ms[order[i]] = dim[i];
	  for (std::size_t pi = 0; pi < p.size(); ++pi)
	  {
	    std::map<char, Index> mfrom = mf;
	    for (std::size_t i = 0; i < Nt; ++i)
	      mfrom[ordert[i]] = t.p[pi][0][i];
	    for (std::size_t i = 0; i < N; ++i)
	      r[pi][0][i] = mfrom[order[i]];

	    std::map<char, Index> msize = ms;
	    for (std::size_t i = 0; i < Nt; ++i)
	      msize[ordert[i]] = t.p[pi][1][i];
	    for (std::size_t i = 0; i < N; ++i)
	      r[pi][1][i] = msize[order[i]];

	    if (superbblas::detail::volume(r[pi][1]) == 0)
	      r[pi] = std::array<Coor<N>, 2>{Coor<N>{{}}, Coor<N>{{}}};
	  }
	  return TensorPartition<N>{dim, r, isLocal};
	}

      private:
	/// Return a partitioning for a non-collective tensor
	/// \param dim: dimension size for the tensor

	static PartitionStored local(Coor<N> dim)
	{
	  return PartitionStored(1, superbblas::PartitionItem<N>{{{}, dim}});
	}

	/// Return a partitioning where the root node has support for the whole tensor
	/// \param dim: dimension size for the tensor

	static PartitionStored all_tensor_on_master(Coor<N> dim)
	{
	  int nprocs = Layout::numNodes();
	  // Set the first coordinate and size of tensor supported on each proc to zero excepting
	  // on proc 0, where the size is set to dim
	  PartitionStored fs(nprocs);
	  if (1 <= nprocs)
	    fs[0][1] = dim;
	  return fs;
	}

	/// Return a partitioning where all nodes have support for the whole tensor
	/// \param dim: dimension size for the tensor

	static PartitionStored all_tensor_replicated(Coor<N> dim)
	{
	  int nprocs = Layout::numNodes();
	  // Set the first coordinate of the tensor supported on each prop to zero and the size
	  // to dim
	  PartitionStored fs(nprocs);
	  for (auto& it : fs)
	    it[1] = dim;
	  return fs;
	}

	/// Return a partitioning for a tensor of `dim` dimension onto a grid of processes
	/// \param order: dimension labels (use x, y, z, t for lattice dimensions)
	/// \param dim: dimension size for the tensor

	static PartitionStored partitioning_chroma_compatible(const std::string& order, Coor<N> dim)
	{
	  // Find a dimension label in `order` that is going to be distributed
	  const char dist_labels[] = "xyzt"; // distributed dimensions
	  int first_dist_label = -1;
	  for (unsigned int i = 0; i < std::strlen(dist_labels); ++i)
	  {
	    const auto& it = std::find(order.begin(), order.end(), dist_labels[i]);
	    if (it != order.end())
	    {
	      first_dist_label = it - order.begin();
	      break;
	    }
	  }

	  // If no dimension is going to be distributed, the whole tensor will have support only on node zero
	  if (first_dist_label < 0)
	    return all_tensor_on_master(dim);

	  // Get the number of procs use in each dimension; for know we put as many as chroma
	  // put onto the lattice dimensions
	  multi1d<int> procs_ = Layout::logicalSize();
	  Coor<N> procs = kvcoors<N>(
	    order, {{'x', procs_[0]}, {'y', procs_[1]}, {'z', procs_[2]}, {'t', procs_[3]}}, 1,
	    NoThrow);

	  // For each proc, get its coordinate in procs (logical coordinate) and compute the
	  // fair range of the tensor supported on the proc
	  int num_procs = Layout::numNodes();
	  PartitionStored fs(num_procs);
	  for (int rank = 0; rank < num_procs; ++rank)
	  {
	    multi1d<int> cproc_ = Layout::getLogicalCoordFrom(rank);
	    Coor<N> cproc = kvcoors<N>(
	      order, {{'x', cproc_[0]}, {'y', cproc_[1]}, {'z', cproc_[2]}, {'t', cproc_[3]}}, 0,
	      NoThrow);
	    for (unsigned int i = 0; i < N; ++i)
	    {
	      // First coordinate in process with rank 'rank' on dimension 'i'
	      fs[rank][0][i] = dim[i] / procs[i] * cproc[i] + std::min(cproc[i], dim[i] % procs[i]);
	      // Number of elements in process with rank 'cproc[i]' on dimension 'i'
	      fs[rank][1][i] =
		dim[i] / procs[i] + (dim[i] % procs[i] > cproc[i] ? 1 : 0) % procs[i];
	    }

	    // Avoid replicating parts of tensor if some of the lattice dimensions does not participate on this tensor
	    for (unsigned int i = 0; i < std::strlen(dist_labels); ++i)
	    {
	      if (std::find(order.begin(), order.end(), dist_labels[i]) != order.end())
		continue;
	      if (cproc_[i] > 0)
	      {
		fs[rank][1][first_dist_label] = 0;
		break;
	      }
	    }

	    // Normalize
	    if (superbblas::detail::volume(fs[rank][1]) == 0)
	      fs[rank] = std::array<Coor<N>, 2>{Coor<N>{{}}, Coor<N>{{}}};
	  }
	  return fs;
	}

	/// Return a partitioning for a tensor of `dim` dimension onto a grid of processes
	/// \param order: dimension labels
	/// \param dim: dimension size for the tensor
	/// \param order: labels to distribute

	static PartitionStored partitioning_distributed(const std::string& order,
							const Coor<N>& dim,
							const std::string& dist_labels)
	{
	  Coor<N> dist_dim = dim;

	  // Update the dimension x with the even-odd label X
	  {
	    const auto& itX = std::find(order.begin(), order.end(), 'X');
	    const auto& itx = std::find(order.begin(), order.end(), 'x');
	    if (itX != order.end() && itx != order.end())
	    {
	      dist_dim[itx - order.begin()] *= dim[itX - order.begin()];
	    }
	  }

	  // Avoid splitting even and odds components for xyzt
	  const std::string even_odd_labels = "xyzt";
	  for (unsigned int i = 0; i < even_odd_labels.size(); ++i)
	  {
	    const auto& it = std::find(order.begin(), order.end(), even_odd_labels[i]);
	    if (it != order.end() && dist_dim[it - order.begin()] % 2 == 0)
	      dist_dim[it - order.begin()] /= 2;
	  }

	  int num_procs = Layout::numNodes();
	  auto procs = superbblas::partitioning_distributed_procs(order.c_str(), dist_dim,
								  dist_labels.c_str(), num_procs);
	  return superbblas::basic_partitioning(order.c_str(), dim, procs, dist_labels.c_str(),
						num_procs);
	}
      };

      template <typename T>
      struct WordType {
	using type = T;
      };

      template <typename T>
      struct WordType<std::complex<T>> {
	using type = T;
      };

      /// Return a Nan for float, double, and complex variants
      template <typename T>
      struct NaN;

      /// Specialization for int
      template <>
      struct NaN<int> {
	static int get()
	{
	  return std::numeric_limits<int>::min();
	}
      };

      /// Specialization for float
      template <>
      struct NaN<float> {
	static float get()
	{
	  return std::nanf("");
	}
      };

      /// Specialization for double
      template <>
      struct NaN<double> {
	static double get()
	{
	  return std::nan("");
	}
      };

      /// Specialization for std::complex
      template <typename T>
      struct NaN<std::complex<T>> {
	static std::complex<T> get()
	{
	  return std::complex<T>{NaN<T>::get(), NaN<T>::get()};
	}
      };

      /// Specialization for DIYComplex
      template <typename T>
      struct NaN<DIYComplex<T>> {
	static T get()
	{
	  return NaN<T>::get();
	}
      };

      /// Return if a float, double, and std::complex is finite
      template <typename T>
      struct IsFinite {
	static bool get(T v)
	{
	  return std::isfinite(v);
	}
      };

      /// Specialization for std::complex
      template <typename T>
      struct IsFinite<std::complex<T>> {
	static bool get(std::complex<T> v)
	{
	  return std::isfinite(v.real()) && std::isfinite(v.imag());
	}
      };

      namespace repr
      {
	template <typename Ostream, std::size_t N>
	Ostream& operator<<(Ostream& s, Coor<N> o)
	{
	  s << "[";
	  if (N > 0)
	    s << o[0];
	  for (unsigned int i = 1; i < N; ++i)
	    s << "," << o[i];
	  s << "]";
	  return s;
	}

	template <typename Ostream, typename T>
	Ostream& operator<<(Ostream& s, std::complex<T> o)
	{
	  s << std::real(o) << "+" << std::imag(o) << "i";
	  return s;
	}

	template <typename Ostream, typename T,
		  typename std::enable_if<std::is_floating_point<T>::value, bool>::type = true>
	Ostream& operator<<(Ostream& s, T o)
	{
	  s.operator<<(o);
	  return s;
	}

	template <typename Ostream, typename T, std::size_t N>
	Ostream& operator<<(Ostream& s, const std::array<T, N>& o);

	template <typename Ostream, typename T>
	Ostream& operator<<(Ostream& s, const std::vector<T>& o)
	{
	  s << "{";
	  for (const auto& i : o)
	    s << i;
	  s << "}";
	  return s;
	}

	template <typename Ostream, typename T, std::size_t N>
	Ostream& operator<<(Ostream& s, const std::array<T, N>& o)
	{
	  s << "{";
	  for (const auto& i : o)
	    s << i;
	  s << "}";
	  return s;
	}
      }

      inline void log(int level, const std::string& s)
      {
	static int log_level = []() {
	  const char* l = std::getenv("SB_LOG");
	  if (l)
	    return std::atoi(l);
	  return 0;
	}();
	if (log_level < level)
	  return;
	QDPIO::cout << s << std::endl;
	QDPIO::cout.flush();
      }

      inline void log_mem()
      {
	if (!superbblas::getTrackingMemory())
	  return;
	std::stringstream ss;
	ss << "mem usage, CPU: " << std::fixed << std::setprecision(0)
	   << superbblas::getCpuMemUsed(0) / 1024 / 1024
	   << " MiB   GPU: " << superbblas::getGpuMemUsed(0) / 1024 / 1024 << " MiB";
	log(1, ss.str());
      }

      template <typename T, typename A, typename B,
		typename std::enable_if<!is_complex<T>::value, bool>::type = true>
      T safe_div(A a, B b)
      {
	if (std::fabs(std::imag(a)) != 0 || std::fabs(std::imag(b)) != 0)
	  throw std::runtime_error("Invalid division");
	return std::real(a) / std::real(b);
      }

      template <typename T, typename A, typename B,
		typename std::enable_if<is_complex<T>::value, bool>::type = true>
      T safe_div(A a, B b)
      {
	return (T)a / (T)b;
      }

      inline bool is_default_device_gpu()
      {
	static bool v = []() {
	  const char* l = std::getenv("SB_DEFAULT_DEVICE_GPU");
	  if (l)
	    return std::atoi(l) != 0;
	  return true;
	}();
	return v;
      }

      /// Data allocation
      template <typename T>
      struct Allocation {
	/// Allocation
	T* ptr;

	/// Context of the allocation
	std::shared_ptr<superbblas::Context> ctx;

	/// Unfinished operations on the allocation
	std::vector<superbblas::Request> pending_operations;

	/// Deallocate the pointer on destruction
	bool destroy_ptr;

	/// Allocate n elements of type T with the given context
	/// \param n: number of elements to allocate
	/// \param ctx: context of the allocation

	Allocation(std::size_t n, const std::shared_ptr<superbblas::Context>& ctx)
	  : ptr(superbblas::allocate<T>(n, *ctx)), ctx(ctx), destroy_ptr(true)
	{
	}

	/// User given pointer that will not be deallocated automatically
	/// \param ptr: pointer to an allocation
	/// \param ctx: context of the allocation

	Allocation(T* ptr, const std::shared_ptr<superbblas::Context>& ctx)
	  : ptr(ptr), ctx(ctx), destroy_ptr(false)
	{
	}

	/// Destructor

	~Allocation()
	{
	  finish_pending_operations();
	  if (destroy_ptr)
	    superbblas::deallocate(ptr, *ctx);
	}

	/// Return the pointer

	T* data()
	{
	  finish_pending_operations();
	  return ptr;
	}

	/// Append a pending operation
	/// \param req: pending operation to finish later
	/// NOTE: finish the operation for pointers not managed; those come from Chroma objects
	///       and Chroma users may get unexpected results if they access the Chroma object
	///       while the proxy Tensor object created with `asTensorView` is still alive.

	void append_pending_operation(const superbblas::Request& req)
	{
	  if (destroy_ptr)
	  {
	    if (req)
	      pending_operations.push_back(req);
	  }
	  else
	  {
	    superbblas::wait(req);
	  }
	}

	/// Finish all pending operations
	void finish_pending_operations()
	{
	  for (const auto& i : pending_operations)
	    superbblas::wait(i);
	  pending_operations.clear();
	}
      };
    }

    enum CopyingTrash { dontCopyingTrash, doCopyingTrash };

    /// Class for operating dense tensors

    template <std::size_t N, typename T>
    struct Tensor {
      using value_type = typename detail::base_type<T>::type;
      static_assert(superbblas::supported_type<value_type>::value, "Not supported type");

      /// Allocation type
      /// NOTE: the complex types decay to the base type, which is needed by `toFakeReal`
      using Allocation = detail::Allocation<typename detail::real_type<T>::type>;

    public:
      std::string order;		      ///< Labels of the tensor dimensions
      Coor<N> dim;			      ///< Length of the tensor dimensions
      std::shared_ptr<Allocation> allocation; ///< Tensor storage
      std::shared_ptr<detail::TensorPartition<N>>
	p;		      ///< Distribution of the tensor among the processes
      Distribution dist;      ///< Whether the tensor is stored on the cpu or a device
      Coor<N> from;	      ///< First active coordinate in the tensor
      Coor<N> size;	      ///< Number of active coordinates on each dimension
      Stride<N> strides;      ///< Displacement for the next element along every direction
      value_type scalar;      ///< Scalar factor of the tensor for unconjugated
      bool conjugate;	      ///< Whether the values are implicitly conjugated
      bool eg;		      ///< Whether this tensor is an example
      bool unordered_writing; ///< Whether to allow execution of writing operations
			      /// in a different order than the one issue
      char complexLabel;      /// label for the dimension having the real and the imaginary part
			      /// (only used when `T` is DIYComplex

      /// Return a string describing the tensor
      /// \param ptr: pointer to the memory allocation
      /// \return: the string representing the tensor

      std::string repr(value_type* ptr = nullptr) const
      {
	using namespace detail::repr;
	std::stringstream ss;
	ss << "Tensor{";
	if (ptr)
	  ss << "data:" << ptr << ", ";
	std::size_t sizemb = p->localVolume() * sizeof(T) / 1024 / 1024;
	ss << "order:" << order << ", from:" << from << ", size:" << size << ", dim:" << dim
	   << ", dist:" << dist << ", local_storage:" << sizemb << " MiB}";
	return ss.str();
      }

      /// Constructor
      /// \param order: dimension labels of the tensor
      /// \param dim: size for each dimension
      /// \param dev: where to allocate the content on the GPU if available (`OnDefaultDevice`)
      ///        or on CPU always (`OnHost`)
      /// \param dist: how to distribute the tensor, see `Distribution`

      Tensor(const std::string& order, Coor<N> dim, DeviceHost dev = OnDefaultDevice,
	     Distribution dist = OnEveryone, char complexLabel = 0)
	: Tensor(order, dim, dev, dist,
		 std::make_shared<detail::TensorPartition<N>>(
		   detail::TensorPartition<N>(order, dim, dist)),
		 false /*= unordered_writing */, complexLabel)
      {
      }

      /// Empty constructor

      Tensor()
	: order(detail::getTrivialOrder(N)),
	  dim{{}},
	  p(std::make_shared<detail::TensorPartition<N>>(
	    detail::TensorPartition<N>(detail::getTrivialOrder(N), {{}}, OnEveryoneReplicated))),
	  dist(OnEveryoneReplicated),
	  from{{}},
	  size{{}},
	  strides{{}},
	  scalar{0},
	  conjugate{false},
	  eg{false},
	  unordered_writing{false},
	  complexLabel{0}
      {
      }

      /// Constructor for bringing the memory allocation (see `asTensorView`)
      /// \param order: dimension labels of the tensor
      /// \param dim: size for each dimension
      /// \param dev: where to allocate the content on the GPU if available (`OnDefaultDevice`)
      ///        or on CPU always (`OnHost`)
      /// \param dist: how to distribute the tensor, see `Distribution`
      /// \param ptr: pointer to the first element

      Tensor(const std::string& order, Coor<N> dim, DeviceHost dev, Distribution dist,
	     value_type* ptr, char complexLabel = 0)
	: order(order),
	  dim(dim),
	  allocation(std::make_shared<Allocation>((typename detail::real_type<T>::type*)ptr,
						  detail::getContext(dev))),
	  dist(dist),
	  from{{}},
	  size(dim),
	  strides(detail::get_strides<std::size_t, N>(dim, superbblas::FastToSlow)),
	  scalar{1},
	  conjugate{false},
	  eg{false},
	  unordered_writing{false},
	  complexLabel{complexLabel}
      {
	checkOrder();

	// For now, TensorPartition creates the same distribution as chroma for tensor with
	// dimensions divisible by chroma logical dimensions
	p = std::make_shared<detail::TensorPartition<N>>(
	  detail::TensorPartition<N>(order, dim, dist));
      }

      /// Internal constructor, used by `toFakeReal`
      /// \param order: dimension labels of the tensor
      /// \param dim: size for each dimension
      /// \param allocation: allocation
      /// \param p: partition of the tensor among the processes
      /// \param dist: how to distribute the tensor, see `Distribution`
      /// \param from: coordinate of the first element in this view
      /// \param size: elements in each direction in this view
      /// \param scalar: scalar factor of this view
      /// \param conjugate: whether the elements are implicitly conjugated
      /// \param eg: whether the tensor is an example gratia
      /// \param unordered_writing: whether its allow to apply the writings in different order
      /// \param complexLabel: complexity label

      Tensor(const std::string& order, Coor<N> dim, std::shared_ptr<Allocation> allocation,
	     std::shared_ptr<detail::TensorPartition<N>> p, Distribution dist, Coor<N> from,
	     Coor<N> size, value_type scalar, bool conjugate, bool eg, bool unordered_writing,
	     char complexLabel)
	: order(order),
	  dim(dim),
	  allocation(allocation),
	  p(p),
	  dist(dist),
	  from(normalize_coor(from, dim)),
	  size(size),
	  strides(detail::get_strides<std::size_t, N>(dim, superbblas::FastToSlow)),
	  scalar(scalar),
	  conjugate(conjugate),
	  eg(eg),
	  unordered_writing(unordered_writing),
	  complexLabel(complexLabel)
      {
	checkOrder();
      }

      /// Internal constructor, used by `make_suitable_for_contraction`
      /// \param order: dimension labels of the tensor
      /// \param dim: size for each dimension
      /// \param dist: how to distribute the tensor, see `Distribution`
      /// \param p: partition of the tensor among the processes

      Tensor(const std::string& order, Coor<N> dim, DeviceHost dev, Distribution dist,
	     std::shared_ptr<detail::TensorPartition<N>> p, bool unordered_writing,
	     char complexLabel)
	: order(order),
	  dim(dim),
	  /// NOTE: the extra two factor shouldn't apply for DIYComplex
	  allocation(std::make_shared<Allocation>(
	    p->localVolume() * (detail::is_complex<T>::value ? 2u : 1u), detail::getContext(dev))),
	  p(p),
	  dist(dist),
	  from{{}},
	  size(dim),
	  strides(detail::get_strides<std::size_t, N>(dim, superbblas::FastToSlow)),
	  scalar{1},
	  conjugate{false},
	  eg{false},
	  unordered_writing{unordered_writing},
	  complexLabel{complexLabel}
      {
	checkOrder();
	detail::log_mem();
      }

    protected:
      /// Internal constructor, used by functions making slices, eg. `kvslice_from_size`
      /// \param order: dimension labels of the tensor
      /// \param from: coordinate of the first element in this view
      /// \param size: elements in each direction in this view

      Tensor(const Tensor& t, const std::string& order, Coor<N> from, Coor<N> size)
	: order(order),
	  dim(t.dim),
	  allocation(t.allocation),
	  p(t.p),
	  dist(t.dist),
	  from(normalize_coor(from, t.dim)),
	  size(size),
	  strides(t.strides),
	  scalar{t.scalar},
	  conjugate{t.conjugate},
	  eg{t.eg},
	  unordered_writing{t.unordered_writing},
	  complexLabel{t.complexLabel}
      {
	checkOrder();
      }

      /// Internal constructor, used by `scale` and `conj`
      /// \param scalar: scalar factor of this view
      /// \param conjugate: whether the elements are implicitly conjugated

      Tensor(const Tensor& t, value_type scalar, bool conjugate)
	: order(t.order),
	  dim(t.dim),
	  allocation(t.allocation),
	  p(t.p),
	  dist(t.dist),
	  from(t.from),
	  size(t.size),
	  strides(t.strides),
	  scalar{scalar},
	  conjugate{conjugate},
	  eg{false},
	  unordered_writing{t.unordered_writing},
	  complexLabel{t.complexLabel}
      {
	checkOrder();
      }

    public:
      /// Return whether the tensor is not empty
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   Tensor<2,Complex> q("cs", {{Nc,0}});
      ///   if (t) std::cout << "t is not empty";  // print this
      ///   if (q) std::cout << "q is not empty";  // doesn't print this

      explicit operator bool() const noexcept
      {
	return volume() > 0;
      }

      /// Return whether the view doesn't start at the origin or doesn't encompass the whole original tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.isSubtensor(); // is false
      ///   t.slice_from_size({{0,1}}).isSubtensor(); // is true
      ///   t.slice_from_size({}, {{0,1}}).isSubtensor(); // is true

      bool isSubtensor() const
      {
	return (from != Coor<N>{{}} || size != dim);
      }

      /// Return the first coordinate supported by the tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.kvfrom(); // is {{'c',0},{'s',0}}
      ///   t.kvslice_from_size({{'s',1}}).kvfrom(); // is {{'c',0},{'s',1}}

      std::map<char, int> kvfrom() const
      {
	std::map<char, int> d;
	for (unsigned int i = 0; i < N; ++i)
	  d[order[i]] = from[i];
	return d;
      }

      /// Return the dimensions of the tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.kvdim(); // is {{'c',Nc},{'s',Ns}}
      ///   t.kvslice_from_size({}, {{'s',1}}).kvdim(); // is {{'c',Nc},{'s',1}}

      std::map<char, int> kvdim() const
      {
	std::map<char, int> d;
	for (unsigned int i = 0; i < N; ++i)
	  d[order[i]] = size[i];
	return d;
      }

      /// Return the allocated dimensions of the tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.kvslice_from_size({}, {{'s',1}}).kvdim(); // is {{'c',Nc},{'s',1}}
      ///   t.kvslice_from_size({}, {{'s',1}}).alloc_kvdim(); // is {{'c',Nc},{'s',2}}

      std::map<char, int> alloc_kvdim() const
      {
	std::map<char, int> d;
	for (unsigned int i = 0; i < N; ++i)
	  d[order[i]] = dim[i];
	return d;
      }

      /// Return the number of the elements in the tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.volume(); // is Nc*Ns
      ///   t.kvslice_from_size({}, {{'s',1}}).volume(); // is Nc*1

      std::size_t volume() const
      {
	return superbblas::detail::volume(size);
      }

      /// Return the number of local elements in the tensor to this process
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.localVolume(); // is Nc*Ns if t replicated among all processes

      std::size_t localVolume() const
      {
	return p->localVolume();
      }

      /// Return the product of the size for each given label
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.volume("c"); // is Nc

      std::size_t volume(const std::string& labels) const
      {
	const auto d = kvdim();
	std::size_t vol = 1;
	for (char l : labels)
	  vol *= d.at(l);
	return vol;
      }

      /// Return whether the tensor is an example
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.volume(); // is Nc*Ns
      ///   t.kvslice_from_size({}, {{'s',1}}).volume(); // is Nc*1

      bool is_eg() const
      {
	return eg;
      }

      /// Return if the given tensor has the same distribution as this
      /// \param w: tensor to compare with

      template <std::size_t Nw, typename Tw>
      bool isDistributedAs(Tensor<Nw, Tw> w, Maybe<std::string> labels = none) const
      {
	return p->is_compatible(order, *w.p, w.order, labels.getSome(order));
      }

      /// Return if the given tensor the same length for the shared dimensions
      /// \param w: tensor to compare with
      /// \param labels: dimension labels to compare if given; all labels otherwise

      template <std::size_t Nw, typename Tw>
      bool is_compatible(Tensor<Nw, Tw> w, Maybe<std::string> labels = none) const
      {
	std::string labels_to_compare = labels.getSome(order);
	auto dims = kvdim();
	auto wdims = w.kvdim();
	for (char c : labels_to_compare)
	  if (wdims.count(c) == 1 && wdims.at(c) != dims.at(c))
	    return false;
	return true;
      }

      /// Return the allocation is managed by superbblas

      bool is_managed() const
      {
	if (!allocation)
	  return true;
	return allocation->destroy_ptr;
      }

      /// Return the pointer to the first local element
      /// NOTE: there will be no pending writing operations

      value_type* data() const
      {
	if (!allocation)
	  return nullptr;

	// If the pointer isn't managed by supperbblas, it may be managed by Chroma
	// and we make sure that all operations from Chroma side are finished
	if (!is_managed())
	  superbblas::syncLegacyStream(ctx());

	return (value_type*)allocation->data();
      }

      /// Return the pointer to the first local element
      /// NOTE: there may be pending writing operations if `unordered_writing` is true

      value_type* data_for_writing() const
      {
	if (!allocation)
	  return nullptr;

	// If the pointer isn't managed by supperbblas, it may be managed by Chroma
	// and we make sure that all operations from Chroma side are finished
	if (!allocation->destroy_ptr)
	  superbblas::syncLegacyStream(ctx());

	if (unordered_writing)
	  return (value_type*)allocation->ptr;

	return (value_type*)allocation->data();
      }

      /// Return the allocation context

      superbblas::Context& ctx() const
      {
	if (!allocation)
	  return *detail::getContext(OnHost);
	return *allocation->ctx;
      }

      /// Get an element of the tensor
      /// \param coor: coordinates of the element to get
      /// \return: the value of the element at the coordinate
      ///
      /// NOTE:
      /// - operation allowed only for tensors supported on the CPU and replicated on every process (or local)
      /// - the operation is slow, avoid in critical performance operations
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}, OnHost, OnEveryoneReplicated);
      ///   t.set({0,1}, 1.0); // set the element with c=0 and s=1 to 1.0
      ///   t.get({0,1}); // get the element with c=0 and s=1
      ///
      ///   Tensor<5,double> q("xyztX", latticeSize<5>("xyztX"), OnHost);
      ///   q.getLocal().set({0,0,0,0,0}, 1.0); // set the first local element in this process to 1.0
      ///   q.getLocal().get({0,0,0,0,0}); // get the first local element in this process

      value_type get(Coor<N> coor) const
      {
	if (ctx().plat != superbblas::CPU)
	  throw std::runtime_error(
	    "Unsupported to `get` elements from tensors not stored on the host");
	if (dist != OnEveryoneReplicated && dist != Local)
	  throw std::runtime_error(
	    "Unsupported to `get` elements on a distributed tensor; change the distribution to "
	    "`OnEveryoneReplicated` or local");
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// coor[i] = coor[i] + from[i]
	for (unsigned int i = 0; i < N; ++i)
	  coor[i] = normalize_coor(normalize_coor(coor[i], size[i]) + from[i], dim[i]);

	return detail::cond_conj(conjugate,
				 data()[detail::coor2index<N>(coor, dim, strides)] * scalar);
      }

      /// Set an element of the tensor
      /// \param coor: coordinates of the element to set
      /// \param v: the new value of the element
      ///
      /// NOTE:
      /// - operation allowed only for tensors supported on the CPU and replicated on every process (or local)
      /// - the operation is slow, avoid in critical performance operations
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}, OnHost, OnEveryoneReplicated);
      ///   t.set({0,1}, 1.0); // set the element with c=0 and s=1 to 1.0
      ///
      ///   Tensor<5,double> q("xyztX", latticeSize<5>("xyztX"), OnHost);
      ///   q.getLocal().set({0,0,0,0,0}, 1.0); // set the first local element in this process to 1.0

      void set(Coor<N> coor, value_type v)
      {
	if (dist != OnEveryoneReplicated && dist != Local)
	  throw std::runtime_error(
	    "Unsupported to `set` elements on a distributed tensor; change the distribution to "
	    "`OnEveryoneReplicated` or local");
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	if (ctx().plat == superbblas::CPU)
	{
	  // coor[i] = coor[i] + from[i]
	  for (unsigned int i = 0; i < N; ++i)
	    coor[i] = normalize_coor(normalize_coor(coor[i], size[i]) + from[i], dim[i]);

	  data_for_writing()[detail::coor2index<N>(coor, dim, strides)] =
	    detail::cond_conj(conjugate, v) / scalar;
	}
	else
	{
	  Tensor<1, T>(std::string(1, order[0]), Coor<1>{1}, OnHost, OnEveryoneReplicated, &v)
	    .copyTo(this->slice_from_size(coor, superbblas::detail::ones<N>()));
	}
      }

      /// Modify the content this tensor with the result of a function on each element
      /// \param func: function () -> value_type
      /// \param threaded: whether to run threaded

      template <typename Func>
      void fillWithCPUFuncNoArgs(Func func, bool threaded = true)
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	std::size_t vol = t.getLocal().volume();
	value_type* ptr = t.data_for_writing();

	if (threaded)
	{
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < vol; ++i)
	    ptr[i] = func();
	}
	else
	{
	  for (std::size_t i = 0; i < vol; ++i)
	    ptr[i] = func();
	}

	t.copyTo(*this);
      }

      /// Fill the tensor with the value of the function applied to each element
      /// \param func: function (Coor<N>) -> value_type
      /// \param threaded: whether to run threaded

      template <typename Func>
      void fillCpuFunCoor(Func func, bool threaded = true)
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	using superbblas::detail::operator+;

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	std::size_t vol = t.getLocal().volume();
	value_type* ptr = t.data_for_writing();
	/// Number of elements in each direction for the local part
	Coor<N> local_size = t.getLocal().size;
	/// Stride for the local volume
	Stride<N> local_stride =
	  superbblas::detail::get_strides<std::size_t>(local_size, superbblas::FastToSlow);
	/// Coordinates of first elements stored locally
	Coor<N> local_from = t.p->localFrom();

	if (threaded)
	{
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < vol; ++i)
	  {
	    // Get the global coordinates
	    Coor<N> c = normalize_coor(
	      superbblas::detail::index2coor(i, local_size, local_stride) + local_from, t.dim);
	    ptr[i] = func(c);
	  }
	}
	else
	{
	  for (std::size_t i = 0; i < vol; ++i)
	  {
	    Coor<N> c = normalize_coor(
	      superbblas::detail::index2coor(i, local_size, local_stride) + local_from, t.dim);
	    ptr[i] = func(c);
	  }
	}

	t.copyTo(*this);
      }

      /// Return a new tensor with the value of the function applied to each element
      /// \param func: function value_type -> value_type for Tr
      /// \param threaded: whether to run threaded

      template <typename Tr, typename Func>
      Tensor<N, Tr> transformWithCPUFun(Func func, bool threaded = true) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	auto r = t.template make_compatible<N, Tr>();
	assert(!r.isSubtensor() && !t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	value_type* tptr = t.data();
	typename Tensor<N, Tr>::value_type* rptr = r.data();

	if (threaded)
	{
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < vol; ++i)
	    rptr[i] = func(tptr[i]);
	}
	else
	{
	  for (std::size_t i = 0; i < vol; ++i)
	    rptr[i] = func(tptr[i]);
	}

	return r.make_sure(none, getDev());
      }

      /// Return a new tensor with the value of the function applied to each element
      /// \param func: function (Coor<N>, value_type) -> value_type for Tr
      /// \param threaded: whether to run threaded

      template <typename Tr, typename FuncWithCoor>
      Tensor<N, Tr> transformWithCPUFunWithCoor(FuncWithCoor func, bool threaded = true) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	using superbblas::detail::operator+;

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	auto r = t.template make_compatible<N, Tr>();
	assert(!r.isSubtensor() && !t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	value_type* tptr = t.data();
	typename Tensor<N, Tr>::value_type* rptr = r.data();
	/// Number of elements in each direction for the local part
	Coor<N> local_size = t.getLocal().size;
	/// Stride for the local volume
	Stride<N> local_stride =
	  superbblas::detail::get_strides<std::size_t>(local_size, superbblas::FastToSlow);
	/// Coordinates of first elements stored locally
	Coor<N> local_from = t.p->localFrom();

	if (threaded)
	{
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < vol; ++i)
	  {
	    // Get the global coordinates
	    Coor<N> c = normalize_coor(
	      superbblas::detail::index2coor(i, local_size, local_stride) + local_from, t.dim);
	    rptr[i] = func(c, tptr[i]);
	  }
	}
	else
	{
	  for (std::size_t i = 0; i < vol; ++i)
	  {
	    // Get the global coordinates
	    Coor<N> c = normalize_coor(
	      superbblas::detail::index2coor(i, local_size, local_stride) + local_from, t.dim);
	    rptr[i] = func(c, tptr[i]);
	  }
	}

	return r.make_sure(none, getDev());
      }

      /// Return the coordinates of the first element returning true by the given function
      /// \param func: function (value_type) -> bool

      template <typename Func>
      Maybe<Coor<N>> find(Func func) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	using superbblas::detail::operator+;

	auto t_ = make_sure(none, OnHost, OnMaster);
	auto t = isSubtensor() ? t_.cloneOn(OnHost) : t_;
	assert(!t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	value_type* tptr = t.data();
	/// Number of elements in each direction for the local part
	Coor<N> local_size = t.getLocal().size;
	/// Stride for the local volume
	Stride<N> local_stride =
	  superbblas::detail::get_strides<std::size_t>(local_size, superbblas::FastToSlow);
	/// Coordinates of first elements stored locally
	Coor<N> local_from = t.p->localFrom();

	Maybe<Coor<N>> r = none;
	for (std::size_t i = 0; i < vol; ++i)
	{
	  if (func(detail::cond_conj(t.conjugate, tptr[i])))
	  {
	    // Get the global coordinates
	    Coor<N> c = normalize_coor(
	      superbblas::detail::index2coor(i, local_size, local_stride) + local_from, t.dim);
	    r = Maybe<Coor<N>>(c);
	    break;
	  }
	}

	return broadcast(r);
      }

      /// Apply the function to each tensor element
      /// \param func: function value_type -> void
      /// \param threaded: whether to run threaded

      template <typename Func>
      void foreachWithCPUFun(Func func, bool threaded = true) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	assert(!t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	value_type* tptr = t.data();

	if (threaded)
	{
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < vol; ++i)
	    func(detail::cond_conj(t.conjugate, tptr[i]));
	}
	else
	{
	  for (std::size_t i = 0; i < vol; ++i)
	    func(detail::cond_conj(t.conjugate, tptr[i]));
	}
      }

      /// Set all elements with the given value
      /// \param v: the new value for all the elements
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}, OnHost, OnEveryoneReplicated);
      ///   t.set({0,1}, 1.0); // set the element with c=0 and s=1 to 1.0
      ///
      ///   Tensor<5,double> q("xyztX", latticeSize<5>("xyztX"), OnHost);
      ///   q.getLocal().set({0,0,0,0,0}, 1.0); // set the first local element in this process to 1.0

      void set(value_type v)
      {
	if (std::norm(v) == 0)
	  set_zero();
	else
	  fillWithCPUFuncNoArgs([=]() { return v; });
      }

      /// Return a new tensors with the dimension labels renamed
      /// \param m: dictionary with the dimensions to rename
      /// \return: new view of the tensor with the dimension labels renamed
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   Tensor<2,Complex> q = t.rename_dims({{'s','S'}});
      ///   q.order; // is "cS"

      Tensor<N, T> rename_dims(const SB::remap& m) const
      {
	return Tensor<N, T>(*this, detail::update_order_and_check<N>(order, m), this->from,
			    this->size);
      }

      /// Return a slice of the tensor starting at coordinate `kvfrom` and taking `kvsize` elements
      /// in each direction. The missing dimensions in `kvfrom` are set to zero and the missing
      /// directions in `kvsize` are set to the size of the tensor.
      ///
      /// \param kvfrom: dictionary with the index of the first element in each direction
      /// \param kvsize: dictionary with the number of elements in each direction
      /// \return: new view of the tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.kvslice_from_size({{'s',1}}); // is a view where the origin element is t(c=0,s=1)
      ///   t.slice_from_size({{0,1}},{{Nc,Ns}}); // equivalent view

      Tensor<N, T> kvslice_from_size(const std::map<char, int>& kvfrom = {},
				     const std::map<char, int>& kvsize = {}) const
      {
	std::map<char, int> updated_kvsize = this->kvdim();
	for (const auto& it : kvsize)
	  updated_kvsize[it.first] = it.second;
	return slice_from_size(kvcoors<N>(order, kvfrom), kvcoors<N>(order, updated_kvsize));
      }

      /// Return a slice of the tensor starting at coordinate `from` and taking `size` elements
      /// in each direction.
      ///
      /// \param from: first coordinate in the view
      /// \param size: number of elements in each direction
      /// \return: new view of the tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   t.slice_from_size({{0,1}},{{Nc,Ns}}); // view of the tensor starting at element (0,1)
      ///   t.slice_from_size({{0,1}},{{1,1}}); // view of a single element at (0,1)
      ///   t.slice_from_size({{0,1}},{{Nc,1}}); // view of all elements with s=1

      Tensor<N, T> slice_from_size(Coor<N> from, Coor<N> size) const
      {
	for (unsigned int i = 0; i < N; ++i)
	{
	  if (size[i] > this->size[i])
	    throw std::runtime_error(
	      "The size of the slice cannot be larger than the original tensor");
	  if (normalize_coor(from[i], this->size[i]) + size[i] > this->size[i] &&
	      this->size[i] != this->dim[i])
	    throw std::runtime_error(
	      "Unsupported to make a view on a non-contiguous range on the tensor");
	}

	using superbblas::detail::operator+;
	return Tensor<N, T>(*this, order, this->from + from, size);
      }

      /// Return a similar tensor keeping the same distribution
      /// \param new_order: dimension labels of the new tensor
      /// \param kvsize: override the length of the given dimensions
      /// \param new_dev: device
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   // Create a new tensor as a collection of three `t` tensors
      ///   Tensor<3,Complex> q = t.make_compatible<3>("csn", {{'n',3}});
      ///   // Create a tensor like q but with allocation on host
      ///   Tensor<3,Complex> v = q.make_compatible(none, {}, OnHost);

      template <std::size_t Nn = N, typename Tn = T>
      Tensor<Nn, Tn> make_compatible(const Maybe<std::string>& new_order = none,
				     const std::map<char, int>& kvsize = {},
				     Maybe<DeviceHost> new_dev = none) const
      {
	std::map<char, int> new_kvdim = kvdim();
	for (const auto& it : kvsize)
	  new_kvdim[it.first] = it.second;
	std::string new_order_ = new_order.getSome(order);
	auto new_dim = kvcoors<Nn>(new_order_, new_kvdim, 0, ThrowOnMissing);
	std::string same_dim_labels;
	auto dim_ = kvdim();
	for (char c : new_order_)
	  if (dim_.count(c) == 1 && dim_.at(c) == new_kvdim.at(c))
	    same_dim_labels.push_back(c);
	return Tensor<Nn, Tn>(new_order_, new_dim, new_dev.getSome(getDev()), dist,
			      std::make_shared<detail::TensorPartition<Nn>>(
				p->get_subpartition(from, size)
				  .make_compatible(order, new_order_, new_dim, same_dim_labels)),
			      false /* unordered_writing */, complexLabel);
      }

      /// Return a tensor on the same device and following the same distribution
      /// \param new_order: dimension labels of the new tensor
      /// \param remaining_char: placeholder for the remaining dimensions
      /// \param kvsize: override the length of the given dimensions
      /// \param new_dev: device
      /// \param new_dist: distribution
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   // Create a new tensor as a collection of three `t` tensors
      ///   Tensor<3,Complex> q = t.make_compatible<3>("%n", '%', "", {{'n',3}});
      ///   // Create a tensor like q but without the dimension c
      ///   Tensor<2,Complex> v = q.make_compatible<2>("%", '%', "c");

      template <std::size_t Nn = N, typename Tn = T>
      Tensor<Nn, Tn>
      make_compatible(const std::string& new_order, char remaining_char,
		const std::string& remove_dims = "", const std::map<char, int>& kvsize = {},
		Maybe<DeviceHost> new_dev = none) const
      {
	return make_compatible<Nn, Tn>(
	  detail::remove_dimensions(get_order_for_reorder(new_order, remaining_char), remove_dims),
	  kvsize, new_dev);
      }

      /// Return a tensor on the same device and following the same distribution
      /// \param new_order: dimension labels of the new tensor
      /// \param kvsize: override the length of the given dimensions
      /// \param new_dev: device
      /// \param new_dist: distribution
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   // Create a new tensor as a collection of three `t` tensors
      ///   Tensor<3,Complex> q = t.like_this<3>("csn", {{'n',3}});
      ///   // Create a tensor like q but with allocation on host
      ///   Tensor<3,Complex> v = q.like_this(none, {}, OnHost);

      template <std::size_t Nn = N, typename Tn = T>
      Tensor<Nn, Tn>
      like_this(const Maybe<std::string>& new_order = none, const std::map<char, int>& kvsize = {},
		Maybe<DeviceHost> new_dev = none, Maybe<Distribution> new_dist = none) const
      {
	std::map<char, int> new_kvdim = kvdim();
	for (const auto& it : kvsize)
	  new_kvdim[it.first] = it.second;
	std::string new_order_ = new_order.getSome(order);
	auto new_dim = kvcoors<Nn>(new_order_, new_kvdim, 0, ThrowOnMissing);
	return Tensor<Nn, Tn>(new_order_, new_dim, new_dev.getSome(getDev()),
			      new_dist.getSome(dist), complexLabel);
      }

      /// Return a tensor on the same device and following the same distribution
      /// \param new_order: dimension labels of the new tensor
      /// \param remaining_char: placeholder for the remaining dimensions
      /// \param kvsize: override the length of the given dimensions
      /// \param new_dev: device
      /// \param new_dist: distribution
      ///
      /// Example:
      ///
      ///   Tensor<2,Complex> t("cs", {{Nc,Ns}});
      ///   // Create a new tensor as a collection of three `t` tensors
      ///   Tensor<3,Complex> q = t.like_this<3>("%n", '%', "", {{'n',3}});
      ///   // Create a tensor like q but without the dimension c
      ///   Tensor<2,Complex> v = q.like_this<2>("%", '%', "c");

      template <std::size_t Nn = N, typename Tn = T>
      Tensor<Nn, Tn>
      like_this(const std::string& new_order, char remaining_char,
		const std::string& remove_dims = "", const std::map<char, int>& kvsize = {},
		Maybe<DeviceHost> new_dev = none, Maybe<Distribution> new_dist = none) const
      {
	return like_this<Nn, Tn>(
	  detail::remove_dimensions(get_order_for_reorder(new_order, remaining_char), remove_dims),
	  kvsize, new_dev, new_dist);
      }

      /// Return a copy of this tensor, possibly with a new precision `nT`
      ///
      /// Example:
      ///
      ///   Tensor<2,std::complex<double>> t("cs", {{Nc,Ns}});
      ///   Tensor<2,std::complex<double>> v = t.clone(); // copy of t
      ///   Tensor<2,std::complex<float>> q = t.clone<std::complex<float>>(); // copy of t in single prec.

      template <typename Tn = T>
      Tensor<N, Tn> clone() const
      {
	return cloneOn<Tn>(getDev());
      }

      /// Return a copy of this tensor on device `new_dev`, possibly with a new precision `nT`
      /// \param new_dev: device that will hold the new tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,std::complex<double>> t("cs", {{Nc,Ns}});
      ///   Tensor<2,std::complex<double>> v = t.cloneOn(OnHost); // copy of t on host
      ///   Tensor<2,std::complex<float>> q = t.cloneOn<std::complex<float>>(OnHost); // copy of t in single prec. on host

      template <typename Tn = T>
      Tensor<N, Tn> cloneOn(DeviceHost new_dev) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	Tensor<N, Tn> r = dist != Glocal ? like_this<N, Tn>(none, {}, new_dev)
					 : make_compatible<N, Tn>(none, {}, new_dev);
	r.conjugate = conjugate;
	copyTo(r);
	return r;
      }

      /// Return a template of this tensor
      ///
      /// Example:
      ///
      ///   Tensor<2,std::complex<double>> t("cs", {{Nc,Ns}});
      ///   Tensor<2,std::complex<double>> t_eg = t.make_eg();
      ///   Tensor<2,std::complex<double>> q = t_eg.like_this(); // create a new tensor like t

      Tensor<N, T> make_eg() const
      {
	return Tensor<N, T>(
	  order, size,
	  std::make_shared<Allocation>(Allocation(std::size_t(0), detail::getContext(getDev()))),
	  std::make_shared<detail::TensorPartition<N>>(p->get_subpartition(from, size)), dist, {{}},
	  size, value_type{1}, false /* not conjugate */, true /* is eg */,
	  false /* ordered writing */, complexLabel);
      }

      /// Return this tensor but allowing consecutive writing operations (from `copyTo`, `contract`...)
      /// to apply non-atomically and in different order than issued. In return, this may reduce the
      /// latency impact by overlapping communications with other operations.

      Tensor<N, T> make_writing_nonatomic() const
      {
	Tensor<N, T> t = *this;
	t.unordered_writing = true;
	return t;
      }

      /// Return the new ordering based on a partial reordering
      /// \param new_order: new dimension labels order
      /// \param remaining_char: if it isn't the null char, placeholder for the dimensions not given
      ///
      /// If the dimension labels order does not match the current order, return a copy of this
      /// tensor with that ordering. If the given order does not contain all dimensions, only the
      /// involved dimensions are permuted.

      std::string get_order_for_reorder(const std::string& new_order, char remaining_char = 0) const
      {
	std::string new_order1;
	if (remaining_char != 0)
	{
	  std::string::size_type rem_pos = new_order.find(remaining_char);
	  if (rem_pos == std::string::npos)
	  {
	    new_order1 = new_order;
	  }
	  else
	  {
	    new_order1 = new_order.substr(0, rem_pos) +
			 detail::remove_dimensions(order, new_order) +
			 new_order.substr(rem_pos + 1, new_order.size() - rem_pos - 1);
	  }
	}
	else
	{
	  new_order1 = order;
	  unsigned int j = 0;
	  for (unsigned int i = 0; i < N; ++i)
	    if (new_order.find(order[i]) != std::string::npos)
	      new_order1[i] = new_order[j++];
	  if (j < new_order.size())
	    throw std::runtime_error("Unknown labels in the given order");
	}

	return new_order1;
      }

      /// Return a copy of this tensor with the given ordering
      /// \param new_order: new dimension labels order
      /// \param remaining_char: if it isn't the null char, placeholder for the dimensions not given
      ///
      /// If the dimension labels order does not match the current order, return a copy of this
      /// tensor with that ordering. If the given order does not contain all dimensions, only the
      /// involved dimensions are permuted.

      Tensor<N, T> reorder(const std::string& new_order, char remaining_char = 0) const
      {
	std::string new_order1 = get_order_for_reorder(new_order, remaining_char);
	if (order == new_order1)
	  return *this;
	Tensor<N, T> r = make_compatible(new_order1);
	if (is_eg())
	  r = r.make_eg();
	else
	  copyTo(r);
	r.unordered_writing = unordered_writing;
	return r;
      }

      /// Return whether the tensor has complex components although being stored with a non-complex type `T`

      static constexpr bool isFakeReal()
      {
	return detail::is_diycomplex<T>::value;
      }

      /// Check that the dimension labels are valid

      void checkOrder() const
      {
	// Check that all labels are different there are N
	detail::check_order<N>(order);

	/// Throw exception if this a diycomplex tensor but there's no complexity label
	if (isFakeReal() && order.find(complexLabel) == std::string::npos)
	  throw std::runtime_error("checkOrder: DIYComplex tensor missing the complexity label");

	for (auto s : size)
	  if (s < 0)
	    std::runtime_error("Invalid tensor size: it should be positive");
      }

      /// Return a view of this tensor with an extra label for the real and the imaginary parts
      /// \param complexLabel: label to represent the real and the imaginary part

      template <typename U = T>
      typename std::enable_if<!detail::is_diycomplex<U>::value && detail::is_complex<U>::value,
			      Tensor<N + 1, DIYComplex<typename U::value_type>>>::type
      toFakeReal(char complexLabel = '.') const
      {
	std::string new_order = std::string(1, complexLabel) + order;
	Coor<N + 1> new_from = {0};
	std::copy_n(from.begin(), N, new_from.begin() + 1);
	Coor<N + 1> new_size = {2};
	std::copy_n(size.begin(), N, new_size.begin() + 1);
	Coor<N + 1> new_dim = {2};
	std::copy_n(dim.begin(), N, new_dim.begin() + 1);
	if (std::fabs(std::imag(scalar)) != 0)
	  throw std::runtime_error(
	    "Unsupported conversion to fake real tensors with an implicit complex scale");
	typename T::value_type new_scalar = std::real(scalar);
	auto new_p = std::make_shared<detail::TensorPartition<N + 1>>(p->insert_dimension(0, 2));

	return Tensor<N + 1, DIYComplex<typename T::value_type>>(
	  new_order, new_dim, allocation, new_p, dist, new_from, new_size, new_scalar, conjugate,
	  eg, unordered_writing, complexLabel);
      }

      template <typename U = T>
      typename std::enable_if<!detail::is_diycomplex<U>::value && !detail::is_complex<U>::value,
			      Tensor<N, DIYComplex<T>>>::type
      toFakeReal(char complexLabel = '.') const
      {
	return Tensor<N, DIYComplex<T>>(order, dim, allocation, p, dist, from, size, scalar,
					conjugate, eg, unordered_writing, complexLabel);
      }

      template <typename U = T>
      typename std::enable_if<detail::is_diycomplex<U>::value, Tensor<N, T>>::type
      toFakeReal(char newComplexLabel = 0) const
      {
	if (newComplexLabel == 0 || newComplexLabel == complexLabel)
	  return *this;
	return rename_dims({{complexLabel, newComplexLabel}});
      }

      /// Return a view or a copy of this tensor where the real and the imaginary parts are together in an std::complex
      /// \param allow_cloning: whether to allow reordering to put the complexity label first

      template <typename U = T>
      typename std::enable_if<detail::is_diycomplex<U>::value,
			      Tensor<N - 1, std::complex<typename U::value_type>>>::type
      toComplex(bool allow_cloning = true) const
      {
	std::size_t dot_pos = order.find(complexLabel);
	std::string new_order = detail::remove_coor(order, dot_pos);

	if (dot_pos != 0)
	{
	  if (allow_cloning || is_eg())
	    return reorder(std::string(1, complexLabel) + new_order).toComplex(false);
	  else
	    throw std::runtime_error("Not allow to create a new tensor in `toComplex`");
	}

	Coor<N - 1> new_from = detail::remove_coor(from, dot_pos);
	Coor<N - 1> new_size = detail::remove_coor(size, dot_pos);
	Coor<N - 1> new_dim = detail::remove_coor(dim, dot_pos);
	std::complex<typename T::value_type> new_scalar{scalar};
	auto new_p = std::make_shared<detail::TensorPartition<N - 1>>(p->remove_dimension(dot_pos));

	return Tensor<N - 1, std::complex<typename T::value_type>>(
	  new_order, new_dim, allocation, new_p, dist, new_from, new_size, new_scalar, conjugate,
	  eg, unordered_writing, 0 /* no complex label */);
      }

      template <typename U = T>
      typename std::enable_if<!detail::is_diycomplex<U>::value, Tensor<N, T>>::type
      toComplex(bool allow_cloning = true) const
      {
	(void)allow_cloning;
	return *this;
      }

      /// Return a copy of the tensor with the values conjugated if the tensor is implicitly conjugated

      template <typename U = T>
      typename std::enable_if<!detail::is_diycomplex<U>::value && detail::is_complex<U>::value,
			      Tensor<N, T>>::type
      make_conjugate_explicit() const
      {
	if (!conjugate) return *this;
	return toFakeReal().make_conjugate_explicit().toComplex();
      }

      template <typename U = T>
      typename std::enable_if<detail::is_diycomplex<U>::value && !detail::is_complex<U>::value,
			      Tensor<N, T>>::type
      make_conjugate_explicit() const
      {
	if (!conjugate) return *this;
	auto t = make_compatible();
	auto this_unconj = conj();
	auto old_scalar = this_unconj.scalar;
	this_unconj.scalar = value_type{1};
	this_unconj.kvslice_from_size({{complexLabel, 0}}, {{complexLabel, 1}})
	  .copyTo(t.kvslice_from_size({{complexLabel, 0}}, {{complexLabel, 1}}));
	this_unconj.kvslice_from_size({{complexLabel, 1}}, {{complexLabel, 1}}).scale(-1)
	  .copyTo(t.kvslice_from_size({{complexLabel, 1}}, {{complexLabel, 1}}));
	return t.scale(old_scalar);
      }

      template <typename U = T>
      typename std::enable_if<!detail::is_diycomplex<U>::value && !detail::is_complex<U>::value,
			      Tensor<N, T>>::type
      make_conjugate_explicit() const
      {
	if (!conjugate) return *this;
	return conj();
      }

      /// Split a dimension into another dimensions
      /// \param dim_label: dimension to split
      /// \param new_labels: the labels of the new dimensions
      /// \param new_dim: number of elements in each new labels

      template <std::size_t Nout, typename std::enable_if<(N > 0), bool>::type = true>
      Tensor<Nout, T> split_dimension(char dim_label, std::string new_labels,
				      const std::map<char, int>& new_dim) const
      {
	using namespace detail;

	// Find the position of dim_label in order
	std::string::size_type pos = order.find(dim_label);
	if (pos == std::string::npos)
	{
	  std::stringstream ss;
	  ss << "Not found label `" << dim_label << "` in this tensor with dimension labels `"
	     << order;
	  throw std::runtime_error(ss.str());
	}

	// Check the length of the output tensor
	if (N + new_labels.size() - 1 != Nout)
	  throw std::runtime_error(
	    "split_dimension: `new_labels` doesn't match the output tensor dimensions!");

	// Check that the size is divisible by the new partition
	if (new_labels.size() == 0)
	{
	  if (size[pos] != 1)
	    throw std::runtime_error("Invalid operation: removing a dimension that isn't singlet");
	  if (dim[pos] != 1)
	    throw std::runtime_error("Unsupported remove a dimension that isn't singlet; clone "
				     "this object before doing the operator");
	}

	// Set the new characteristics of the tensor
	std::string new_order = std::string(order.begin(), order.begin() + pos) + new_labels +
				std::string(order.begin() + pos + 1, order.end());
	Coor<Nout> d;
	std::copy_n(dim.begin(), pos, d.begin());
	for (unsigned int i = 0; i < new_labels.size(); ++i)
	  d[pos + i] = new_dim.at(new_labels[i]);
	std::copy_n(dim.begin() + pos + 1, N - pos - 1, d.begin() + pos + new_labels.size());

	// Transform the partition
	auto new_p = std::make_shared<detail::TensorPartition<Nout>>(p->split_dimension(pos, d));

	return Tensor<Nout, T>(new_order, new_p->dim, allocation, new_p, dist,
			       detail::split_dimension(pos, from, d, From),
			       detail::split_dimension(pos, size, d, Size), scalar, conjugate, eg,
			       unordered_writing, complexLabel);
      }

      /// Split a dimension into another dimensions
      /// \param dim_label: dimension to split
      /// \param new_labels: the labels of the new dimensions
      /// \param step: length of the first label in `new_labels`

      Tensor<N + 1, T> split_dimension(char dim_label, const std::string& new_labels,
				       Index step) const
      {
	if (new_labels.size() != 2)
	  throw std::runtime_error(
	    "split_dimension: invalid `new_labels`, it should have size two");
	if (kvdim().at(dim_label) == 1)
	  step = 1;
	if (kvdim().at(dim_label) % step != 0)
	  throw std::runtime_error(
	    "split_dimension: invalid `step`, it should divide the dimension size");
	return split_dimension<N + 1>(
	  dim_label, new_labels,
	  {{new_labels[0], step}, {new_labels[1], kvdim().at(dim_label) / step}});
      }

      /// Collapse several dimensions into a new one
      /// \param dim_label: dimension to split
      /// \param new_labels: the labels of the new dimensions
      /// \param new_dim: number of elements in each new labels

      template <std::size_t Nout, typename std::enable_if<(N > 0), bool>::type = true>
      Tensor<Nout, T> collapse_dimensions(std::string labels, char new_label,
					  bool allow_copy = false) const
      {
	using namespace detail;

	// Check that all `labels` are together in `order`
	auto s_labels = std::search(order.begin(), order.end(), labels.begin(), labels.end());
	if (s_labels == order.end())
	{
	  if (!allow_copy)
	    throw std::runtime_error(
	      "collapse_dimensions: invalid labels to collapse or they are "
	      " not appear together in the same ordering and copying is not allow");

	  // Find the position of the first label
	  std::string new_order;
	  for (char c : order)
	  {
	    if (std::find(labels.begin(), labels.end(), c) != labels.end())
	      break;
	    new_order.push_back(c);
	  }
	  new_order += labels + "%";
	  return reorder(new_order, '%').template collapse_dimensions<Nout>(labels, new_label);
	}

	// Check the length of the output tensor
	if (N - labels.size() + 1 != Nout)
	  throw std::runtime_error(
	    "collapse_dimensions: `labels` doesn't match the output tensor dimensions!");

	// Lets put the new dimension on the first dimension to collapse
	std::size_t pos = s_labels - order.begin();

	// Set the new characteristics of the tensor
	std::string new_order = std::string(order.begin(), order.begin() + pos) +
				std::string(1, new_label) +
				std::string(order.begin() + pos + labels.size(), order.end());

	// Transform the partition
	auto new_p = std::make_shared<detail::TensorPartition<Nout>>(
	  p->template collapse_dimensions<Nout>(pos));

	return Tensor<Nout, T>(new_order, new_p->dim, allocation, new_p, dist,
			       detail::collapse_dimensions<Nout>(pos, from, dim, From),
			       detail::collapse_dimensions<Nout>(pos, size, dim, Size), scalar,
			       conjugate, eg, unordered_writing, complexLabel);
      }

      /// Rearrange several dimensions into new ones
      /// \param m: maps from suborder of the current tensor to new orders
      /// \param allow_copy: whether to allow to return a reordered copy of the current tensor

      template <std::size_t Nout = N, typename std::enable_if<(N > 0), bool>::type = true>
      Tensor<Nout, T> reshape_dimensions(const std::map<std::string, std::string>& m,
					 const std::map<char, int>& new_dim,
					 bool allow_copy = true) const
      {
	using namespace detail;

	// Check that all suborders in `m` are together in `order`
	std::string old_order = order;
	for (const auto& it : m)
	{
	  if (it.first.size() == 0 || it.second.size() == 0)
	    throw std::runtime_error("reshape_dimensions: invalid map element with empty string");

	  auto s_labels = std::search(order.begin(), order.end(), it.first.begin(), it.first.end());
	  if (s_labels == order.end())
	  {
	    if (!allow_copy)
	      throw std::runtime_error(
		"reshape_dimensions: invalid labels to reshape or they do "
		" not appear together in the same way as in the tensor and copying is not allow");

	    // Find the position of the first label to reshape and enforce the given subordering
	    std::string old_order0;
	    for (char c : old_order)
	    {
	      if (std::find(it.first.begin(), it.first.end(), c) != it.first.end())
		break;
	      old_order0.push_back(c);
	    }
	    old_order = old_order0 + it.first + remove_dimensions(old_order, old_order0 + it.first);
	  }
	}
	if (old_order != order)
	  return reorder(old_order).template reshape_dimensions<Nout>(m, new_dim, true);

	// Check the length of the output tensor
	int nout = N;
	for (const auto& it : m)
	  nout += (int)it.second.size() - (int)it.first.size();
	if (nout != Nout)
	  throw std::runtime_error("reshape_dimensions: the resulting tensor after the changes "
				   "given in `m` doesn't match the output tensor's dimensions!");

	// Compute the new order
	std::string new_order = order;
	for (const auto& it : m)
	{
	  auto s_first = std::find(new_order.begin(), new_order.end(), it.first.front());
	  new_order = std::string(new_order.begin(), s_first) + it.second +
		      std::string(s_first + it.first.size(), new_order.end());
	}

	// Compute the dimensions of the new tensor
	auto new_dim0 = kvdim();
	for (const auto& it : new_dim)
	  new_dim0[it.first] = it.second;

	// The last label on the new subordering is optional
	for (const auto& it : m)
	  if (new_dim.count(it.second.back()) == 0)
	    new_dim0[it.second.back()] = std::numeric_limits<int>::max();

	// Compute the number of dimensions to collapse and to split
	std::map<char, int> m_ncollapse, m_nsplit;
	for (const auto& it : m)
	{
	  m_ncollapse[it.first.front()] = it.first.size();
	  m_nsplit[it.first.front()] = it.second.size();
	}
	Coor<N> ncollapse = kvcoors<N>(order, m_ncollapse, 1);
	Coor<N> nsplit = kvcoors<N>(order, m_nsplit, 1);
	auto d_aux = detail::reshape_dimensions<Nout>(
	  ncollapse, nsplit, dim, kvcoors<Nout>(new_order, new_dim0, 0, ThrowOnMissing), Size, dim);
	if (d_aux.first != Success)
	  throw std::runtime_error(
	    "reshape_dimensions: invalid reshape, most likely some new dimension is too short");
	auto d = d_aux.second; // new dimensions

	// Transform the partition
	auto new_p_aux = p->template reshape_dimensions<Nout>(ncollapse, nsplit, d);
	auto new_from = detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, d, From, from);
	auto new_size = detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, d, Size, size);

	// Whether a compatible partition can be made that doesn't require a copy of the tensors' data
	bool success =
	  (new_p_aux.first == Success && new_from.first == Success && new_size.first == Success);

	// Return the new tensor
	if (!allow_copy && !success)
	{
	  throw std::runtime_error("reshape_dimensions: unsupported reshape without copying");
	}
	else if (success)
	{
	  // Return a tensor with this data but a different shape
	  auto new_p = std::make_shared<detail::TensorPartition<Nout>>(new_p_aux.second);
	  return Tensor<Nout, T>(new_order, new_p->dim, allocation, new_p, dist, new_from.second,
				 new_size.second, scalar, conjugate, eg, unordered_writing,
				 complexLabel);
	}
	else if (new_size.first != Success)
	{
	  // This shouldn't happen
	  throw std::runtime_error("reshape_dimensions: something is wrong...");
	}
	else
	{
	  // Try the other way around
	  Tensor<Nout, T> r(new_order, new_size.second, getDev(), dist, complexLabel);
	  r.scalar = scalar;
	  r.conjugate = conjugate;
	  r.unordered_writing = unordered_writing;
	  if (eg)
	    return r.make_eg();
	  std::map<std::string, std::string> reverse_m;
	  for (const auto& it : m)
	    reverse_m[it.second] = it.first;
	  copyTo(r.template reshape_dimensions<N>(reverse_m, kvdim(), false));
	  return r;
	}
      }

      /// Append a dimension with size one
      /// \param new_label: label for the new dimension

      template <typename std::enable_if<(N > 0), bool>::type = true>
      Tensor<N + 1, T> append_dimension(char new_label) const
      {
	std::string last_label{order.back()};
	return reshape_dimensions<N + 1>({{last_label, last_label + std::string(1, new_label)}},
					 {{new_label, 1}}, false);
      }

      /// Coarse the support range of the tensor on each process
      /// \param blocking: blocking for each dimension

      template <typename std::enable_if<(N > 0), bool>::type = true>
      Tensor<N, T> coarse_support(const std::map<char, int>& blocking) const
      {
	// Get the blocking and check that it divides each diension
	auto c_blk = kvcoors<N>(order, blocking, 1);
	for (std::size_t i = 0; i < N; ++i)
	  if (dim[i] % c_blk[i] != 0)
	    throw std::runtime_error(
	      "coarse_support: the given blocking isn't dividing the tensor dimensions");

	// Transform the partition
	auto new_p = std::make_shared<detail::TensorPartition<N>>(p->coarse_support(c_blk));

	// Create output tensor
	Tensor<N, T> r(order, dim, getDev(), dist, new_p, unordered_writing, complexLabel);
	r.from = from;
	r.size = size;
	r.conjugate = conjugate;

	// Return it
	if (is_eg())
	  return r.make_eg();
	copyTo(r);
	return r;
      }

      /// Copy/add this tensor into the given one
      /// NOTE: if this tensor or the given tensor is fake real, force both to be fake real

      template <std::size_t Nw, typename Tw, std::size_t Nm = N, typename Tm = float,
		std::size_t Nwm = Nw, typename Twm = float,
		typename std::enable_if<
		  detail::is_complex<T>::value != detail::is_complex<Tw>::value, bool>::type = true>
      void doAction(Action action, Tensor<Nw, Tw> w, Tensor<Nm, Tm> m = {},
		    Tensor<Nwm, Twm> wm = {}, const std::string& uneven_mask_labels = "",
		    CopyingTrash copying_trash = dontCopyingTrash) const
      {
	if (m || wm)
	  throw std::runtime_error(
	    "doAction: unsupported mixing real and complex types with masks");
	toFakeReal().doAction(action, w.toFakeReal(), {}, {}, uneven_mask_labels, dontCopyingTrash);
      }

      /// Return the local support of this tensor
      Tensor<N, T> getLocal() const
      {
	// Shortcut for empty and local tensors
	if (!*this || dist == Local)
	  return *this;

	// Finish writing operations: local tensor will not be able to finish pending writing operations
	data();

	// Compute the size of the intersection of the current view and the local support
	Coor<N> lfrom, lsize;
	superbblas::detail::intersection(p->localFrom(), p->localSize(), from, size, dim, lfrom,
					 lsize);

	// If the current process has no support, return the empty tensor
	if (superbblas::detail::volume(lsize) == 0)
	  return Tensor<N, T>{};

	using superbblas::detail::operator-;
	return Tensor<N, T>(order, p->localSize(), allocation,
			    std::make_shared<detail::TensorPartition<N>>(p->get_local_partition()),
			    Local, normalize_coor(from - p->localFrom(), dim), lsize, scalar,
			    conjugate, eg, false /* ordered writing */, complexLabel);
      }

      /// Return the local support of this tensor as a subset of the global tensor
      Tensor<N, T> getGlocal() const
      {
	// Shortcut for empty and local tensors
	if (!*this || dist == Glocal || dist == Local)
	  return *this;

	// Finish writing operations: local tensor will not be able to finish pending writing operations
	data();

	return Tensor<N, T>(order, dim, allocation,
			    std::make_shared<detail::TensorPartition<N>>(p->get_glocal_partition()),
			    Glocal, from, size, scalar, conjugate, eg,
			    false /* ordered writing */, complexLabel);
      }

      /// Set zero
      void set_zero()
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	value_type* ptr = data_for_writing();
	MPI_Comm comm =
	  (dist == OnMaster || dist == Local || dist == Glocal ? MPI_COMM_SELF : MPI_COMM_WORLD);
	auto p_disp = (dist == Glocal ? p->MpiProcRank() : 0);
	if (dist != OnMaster || Layout::nodeNumber() == 0)
	{
	  superbblas::copy<N, N>(value_type{0}, p->p.data() + p_disp, 1, order.c_str(), from, size,
				 dim, (const value_type**)&ptr, nullptr, &ctx(),
				 p->p.data() + p_disp, 1, order.c_str(), from, dim, &ptr, nullptr,
				 &ctx(), comm, superbblas::FastToSlow, superbblas::Copy);
	  // Force synchronization in superbblas stream if the allocation isn't managed by superbblas
	  if (!is_managed())
	    superbblas::sync(ctx());
	}
      }

      /// Return whether the given tensor has the same distribution as this one
      template <typename Tv>
      bool is_distributed_like(Tensor<N, Tv> v) const
      {
	return order == v.order && from == v.from && size == v.size && dim == v.dim &&
	       p->p == v.p->p;
      }

      /// Return whether the given tensor has the same distribution as this one
      template <std::size_t Nv, typename Tv, typename std::enable_if<N != Nv, bool>::type = true>
      bool is_distributed_like(Tensor<Nv, Tv>) const
      {
	return false;
      }

      /// Return whether the given tensor has the same shape, distribution, type, and implicit scalar
      /// and conjugacy.
      /// \param v: tensor to compare

      bool is_like(Tensor<N, T> v) const
      {
	return order == v.order && from == v.from && size == v.size && dim == v.dim &&
	       scalar == v.scalar && conjugate == v.conjugate && dist == v.dist && p->p == v.p->p;
      }

      template <
	std::size_t Nv, typename Tv,
	typename std::enable_if<(N != Nv || !std::is_same<T, Tv>::value), bool>::type = true>
      bool is_like(Tensor<Nv, Tv>) const
      {
	return false;
      }

      /// Return whether the given tensor has the same memory allocation as this one
      /// \param v: tensor to compare

      template <std::size_t Nv>
      bool has_same_allocation(Tensor<Nv, T> v) const
      {
	// Compare the allocation pointer, not the actual allocation.ptr; we are making sure that
	// the two allocations are on the same device in this way
	return allocation == v.allocation;
      }

      template <std::size_t Nv, typename Tv,
		typename std::enable_if<!std::is_same<typename detail::real_type<T>::type,
						      typename detail::real_type<Tv>::type>::value,
					bool>::type = true>
      bool has_same_allocation(Tensor<Nv, Tv>) const
      {
	return false;
      }

      /// Return whether the given tensor has the same distribution as this one
      Tensor<N, float> create_mask() const
      {
	Tensor<N, float> m{order, dim, getDev(), dist, p, false /* ordered writing */,
			   0 /* no complex label */};
	m.set_zero();
	m.from = from;
	m.size = size;
	m.conjugate = conjugate;
	return m;
      }

      /// Copy/Add this tensor into the given one; and copy only where the values of the mask are nonzero if given
      template <std::size_t Nw, typename Tw, std::size_t Nm = N, typename Tm = float,
		std::size_t Nwm = Nw, typename Twm = float,
		typename std::enable_if<
		  detail::is_complex<T>::value == detail::is_complex<Tw>::value, bool>::type = true>
      void doAction(Action action, Tensor<Nw, Tw> w, Tensor<Nm, Tm> m = {},
		    Tensor<Nwm, Twm> wm = {}, const std::string& uneven_mask_labels = "",
		    CopyingTrash copying_trash = dontCopyingTrash) const
      {
	if (is_eg() || w.is_eg() || (m && m.is_eg()) || (wm && wm.is_eg()))
	  throw std::runtime_error("Invalid operation from an example tensor");

	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > wsize[i] && !detail::is_in(uneven_mask_labels, order[i]))
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");
	if (m || wm)
	  for (unsigned int i = 0; i < N; ++i)
	    if (size[i] != wsize[i] && !detail::is_in(uneven_mask_labels, order[i]))
	      throw std::runtime_error("copying with masks tensor with different dimensions");

	if (action == AddTo && w.scalar != decltype(w.scalar){1})
	  throw std::runtime_error(
	    "Not allowed to add to a tensor whose implicit scalar factor is not one");

	if (conjugate != w.conjugate &&
	    (detail::is_complex<T>::value || detail::is_diycomplex<T>::value))
	{
	  auto this_conj = (conjugate ? *this : conj()).make_conjugate_explicit();
	  auto new_this = (conjugate ? this_conj : this_conj.conj());
	  new_this.doAction(action, w, m, wm, uneven_mask_labels);
	  return;
	}

	bool some_is_local =
	  dist == Local || w.dist == Local || (m && m.dist == Local) || (wm && wm.dist == Local);
	bool some_isnt_local =
	  dist != Local || w.dist != Local || (m && m.dist != Local) || (wm && wm.dist != Local);
	if (some_is_local && some_isnt_local)
	  throw std::runtime_error(
	    "Not allowed to copy or add a non-local tensor into a local tensor or vice versa");

	// Transform to a local copy when both tensors have the same shape and distribution
	if (action == CopyTo && is_like(w))
	{
	  if (some_isnt_local)
	  {
	    auto this_local = getLocal();
	    auto w_local = w.getLocal();
	    if (w_local && this_local)
	      this_local.doAction(action, w_local, m.getLocal(), wm.getLocal(), uneven_mask_labels);
	    return;
	  }
	  if (some_is_local && has_same_allocation(w) && !m && !wm)
	    return;
	}

	// Check if some dimension size doesn't match
	if (m || wm)
	{
	  std::map<char, int> new_size;
	  bool v_has_new_size = false, w_has_new_size = false;
	  for (unsigned int i = 0; i < N; ++i)
	  {
	    new_size[order[i]] = std::max(size[i], wsize[i]);
	    if (size[i] != wsize[i])
	    {
	      v_has_new_size |= new_size[order[i]] != size[i];
	      w_has_new_size |= new_size[order[i]] != wsize[i];
	    }
	  }
	  for (const auto& it : m.kvdim())
	    if (new_size.count(it.first) == 0)
	      new_size[it.first] = it.second;
	  for (const auto& it : w.kvdim())
	    if (new_size.count(it.first) == 0)
	      new_size[it.first] = it.second;
	  for (const auto& it : wm.kvdim())
	    if (new_size.count(it.first) == 0)
	      new_size[it.first] = it.second;
	  if (v_has_new_size || w_has_new_size)
	  {
	    auto v0 = *this;
	    auto w0 = w;
	    auto m0 = m;
	    auto wm0 = wm;
	    if (v_has_new_size)
	    {
	      v0 = w_has_new_size ? make_compatible(none, new_size)
				  : w0.template make_compatible<N, T>(order, new_size);
	      copyTo(v0, doCopyingTrash);
	      if (m)
	      {
		m0 = v0.template make_compatible<Nm, Tm>(m.order, new_size);
		m0.set_zero();
		m.copyTo(m0);
	      }
	    }
	    if (w_has_new_size)
	    {
	      w0 = v0.template make_compatible<Nw, Tw>(w.order, new_size);
	      w.copyTo(w0, doCopyingTrash);
	      if (wm)
	      {
		wm0 = w0.template make_compatible<Nwm, Twm>(wm.order, new_size);
		wm0.set_zero();
		wm.copyTo(wm0);
	      }
	    }
	    v0.doAction(action, w0, m0, wm0);
	    if (w_has_new_size)
	      w0.kvslice_from_size({}, w.kvdim()).doAction(CopyTo, w);
	    return;
	  }
	}

	// Compute masks
	float *m0ptr = nullptr, *m1ptr = nullptr;
	Tensor<N, float> m0;
	Tensor<Nw, float> m1;
	if (m || wm)
	{
	  if (m)
	  {
	    if (is_distributed_like(m))
	    {
	      m0 = m;
	    }
	    else
	    {
	      m0 = create_mask();
	      m.copyTo(m0);
	    }
	  }

	  if (wm)
	  {
	    if (w.is_distributed_like(wm))
	    {
	      m1 = wm;
	    }
	    else
	    {
	      m1 = w.create_mask();
	      wm.copyTo(m1);
	    }
	  }

	  if (m && !wm)
	    m0.copyTo(m1);
	  if (!m && wm)
	    m1.copyTo(m0);

	  m0ptr = m0.data();
	  m1ptr = m1.data();
	}

	// Get the pointers to data
	value_type* ptr = data();
	typename decltype(w)::value_type* w_ptr = w.data_for_writing();

	// Shortcuts for who is involved in the operation
	bool do_operation = true;
	int p_disp = 0;
	MPI_Comm comm = MPI_COMM_WORLD;
	// a) if the origin and destination tensors have full support on the master node
	// and the destination tensor is only supported on the master node, the operation
	// only happens on the master node
	if ((dist == OnMaster || dist == OnEveryoneReplicated) && w.dist == OnMaster)
	{
	  do_operation = Layout::nodeNumber() == 0;
	  comm = MPI_COMM_SELF;
	}
	// b) if the origin and destination tensors are replicated on every node or they are
	// local, the operation happens locally on each node
	else if ((dist == OnEveryoneReplicated && w.dist == OnEveryoneReplicated) || dist == Local)
	{
	  comm = MPI_COMM_SELF;
	}
	// c) any is glocal
	if (dist == Glocal || w.dist == Glocal) {
	  comm = MPI_COMM_SELF;
	  p_disp = p->MpiProcRank();
	}

	if (do_operation)
	{
	  superbblas::Request req;
	  superbblas::copy<N, Nw>(detail::safe_div<value_type>(scalar, w.scalar),
				  p->p.data() + p_disp, 1, order.c_str(), from, size, dim,
				  (const value_type**)&ptr, (const float**)&m0ptr, &ctx(),
				  w.p->p.data() + p_disp, 1, w.order.c_str(), w.from, w.dim, &w_ptr,
				  (const float**)&m1ptr, &w.ctx(), comm, superbblas::FastToSlow,
				  action == CopyTo ? superbblas::Copy : superbblas::Add,
				  &req /*, copying_trash == doCopyingTrash*/);
	  w.allocation->append_pending_operation(req);
	  // Force synchronization in superbblas stream if the destination allocation isn't managed by superbblas
	  if (!w.is_managed())
	    superbblas::sync(w.ctx());
	}
      }

      /// Copy this tensor into the given one
      template <std::size_t Nw, typename Tw>
      void copyTo(Tensor<Nw, Tw> w, CopyingTrash copying_trash = dontCopyingTrash) const
      {
	doAction(CopyTo, w, {}, {}, "", copying_trash);
      }

      /// Copy this tensor into the given one but only the elements where the mask is nonzero
      template <std::size_t Nw, typename Tw, std::size_t Nm, typename Tm, std::size_t Nwm,
		typename Twm>
      void copyToWithMask(Tensor<Nw, Tw> w, Tensor<Nm, Tm> m, Tensor<Nwm, Twm> wm,
			  const std::string uneven_mask_labels = "") const
      {
	doAction<Nw, Tw, Nm, Tm, Nwm, Twm>(CopyTo, w, m, wm, uneven_mask_labels);
      }

      // Add `this` tensor into the given one
      template <std::size_t Nw, typename Tw>
      void addTo(Tensor<Nw, Tw> w) const
      {
	doAction(AddTo, w);
      }

      /// Return a copy of this tensor with a compatible distribution to be contracted with the given tensor
      /// \param v: given tensor

      template <std::size_t Nv, typename Tv,
		typename std::enable_if<std::is_same<T, Tv>::value, bool>::type = true>
      Tensor<N, T> make_suitable_for_contraction(Tensor<Nv, Tv> v) const
      {
	if (dist != OnEveryoneReplicated)
	  throw std::runtime_error("Invalid tensor distribution for this function");

	Coor<N> vsize = kvcoors<N>(order, v.kvdim(), 0, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (vsize[i] != 0 && vsize[i] != size[i])
	    throw std::runtime_error(
	      "Invalid tensor contractions: one of the dimensions does not match");

	auto new_p = std::make_shared<detail::TensorPartition<N>>(
	  p->make_suitable_for_contraction(order, *v.p, v.order));

	Tensor<N, T> r(order, dim, getDev(), OnEveryone, new_p, unordered_writing, complexLabel);
	copyTo(r);
	return r;
      }

      // Contract the dimensions with the same label in `v` and `w` than do not appear on `this` tensor.
      template <std::size_t Nv, std::size_t Nw>
      void contract(Tensor<Nv, T> v, const remap& mv, Conjugation conjv, Tensor<Nw, T> w,
		    const remap& mw, Conjugation conjw, const remap& mr = {},
		    value_type beta = value_type{0})
      {
	if (is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// NOTE: Superbblas tensor contraction is shit and does not deal with contracting a host and
	// device tensor (for now)
	// a) If either v or w is on OnDevice, force both to be on device
	if (v.ctx().plat != w.ctx().plat)
	{
	  if (v.getDev() != OnDefaultDevice)
	    v = v.cloneOn(OnDefaultDevice);
	  if (w.getDev() != OnDefaultDevice)
	    w = w.cloneOn(OnDefaultDevice);
	}

	// b) Do arrangements if the input tensors are on a different device than the result tensor
	if (getDev() != v.getDev())
	{
	  Tensor<N, T> aux =
	    std::norm(beta) == 0 ? like_this(none, {}, v.getDev()) : cloneOn(v.getDev());
	  aux.contract(v, mv, conjv, w, mw, conjw, mr, beta);
	  aux.copyTo(*this);
	  return;
	}

	if ((v.dist == Local) != (w.dist == Local) || (w.dist == Local) != (dist == Local) ||
	    (v.dist == Glocal) != (w.dist == Glocal) || (w.dist == Glocal) != (dist == Glocal))
	  throw std::runtime_error("contract: one of the contracted tensors or the output tensor "
				   "is local/glocal and others are not!");

	MPI_Comm comm = (dist == Local || dist == Glocal ? MPI_COMM_SELF : MPI_COMM_WORLD);
	auto p_disp = (dist == Glocal ? p->MpiProcRank() : 0);
	
	value_type* v_ptr = v.data();
	value_type* w_ptr = w.data();
	value_type* ptr = std::norm(beta) == 0 ? data_for_writing() : data();
	std::string orderv_ = detail::update_order_and_check<Nv>(v.order, mv);
	std::string orderw_ = detail::update_order_and_check<Nw>(w.order, mw);
	std::string order_ = detail::update_order_and_check<N>(order, mr);
	bool conjv_ = (((conjv == Conjugate) xor v.conjugate) xor conjugate);
	bool conjw_ = (((conjw == Conjugate) xor w.conjugate) xor conjugate);
	superbblas::contraction<Nv, Nw, N>(
	  detail::cond_conj(conjv_, v.scalar) * detail::cond_conj(conjw_, w.scalar) / scalar, //
	  v.p->p.data() + p_disp, v.from, v.size, v.dim, 1, orderv_.c_str(), conjv_,
	  (const value_type**)&v_ptr, &v.ctx(), //
	  w.p->p.data() + p_disp, w.from, w.size, w.dim, 1, orderw_.c_str(), conjw_,
	  (const value_type**)&w_ptr, &w.ctx(), //
	  detail::cond_conj(conjugate, beta), p->p.data() + p_disp, from, size, dim, 1,
	  order_.c_str(), &ptr, &ctx(), comm, superbblas::FastToSlow);

	// Force synchronization in superbblas stream if the destination allocation isn't managed by superbblas
	if (!is_managed())
	  superbblas::sync(ctx());
      }

      /// Compute the Cholesky factor of `v' and contract its inverse with `w`
      /// \param v: tensor to compute the Cholesky factor
      /// \param order_rows: labels that are rows of the matrices to factor
      /// \param order_cols: labels that are columns of the matrices to factor
      /// \param w: the other tensor to contract

      template <std::size_t Nv, std::size_t Nw>
      void cholInv(Tensor<Nv, T> v, const std::string& order_rows, const std::string& order_cols,
		   Tensor<Nw, T> w)
      {
	if (is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// Conjugacy isn't supported
	if (v.conjugate || w.conjugate || conjugate)
	  throw std::runtime_error("cholInv: Unsupported implicit conjugate tensors");

	// If either v or w is on OnDevice, force both to be on device
	if (v.ctx().plat != w.ctx().plat)
	{
	  if (v.getDev() != OnDefaultDevice)
	    v = v.cloneOn(OnDefaultDevice);
	  if (w.getDev() != OnDefaultDevice)
	    w = w.cloneOn(OnDefaultDevice);
	}

	// Superbblas tensor contraction is shit and those not deal with subtensors or contracting a host and
	// device tensor (for now)
	if (v.isSubtensor())
	  v = v.clone();
	if (w.isSubtensor())
	  w = w.clone();
	if (isSubtensor() || getDev() != v.getDev())
	{
	  Tensor<N, T> aux = make_compatible(none, {}, v.getDev());
	  aux.cholInv(v, order_rows, order_cols, w);
	  aux.copyTo(*this);
	  return;
	}

	// v is going to be modified and is reference, make a clone
	if (v.allocation.use_count() > 1)
	  v = v.clone();

	if ((v.dist == Local) != (w.dist == Local) || (w.dist == Local) != (dist == Local))
	  throw std::runtime_error(
	    "One of the contracted tensors or the output tensor is local and others are not!");

	if (detail::isDistributedOnEveryone(v.dist) && w.dist == OnEveryoneReplicated)
	  w = w.make_suitable_for_contraction(v);

	if (v.dist == OnEveryoneReplicated && detail::isDistributedOnEveryone(w.dist))
	  v = v.make_suitable_for_contraction(w);

	if (std::fabs(std::imag(v.scalar)) != 0 || std::real(v.scalar) < 0)
	  throw std::runtime_error("cholInv: unsupported a negative or imaginary scale");

	value_type* v_ptr = v.data();
	value_type* w_ptr = w.data();
	value_type* ptr = data_for_writing();
	superbblas::cholesky<Nv>(v.p->p.data(), v.dim, 1, v.order.c_str(), &v_ptr,
				 order_rows.c_str(), order_cols.c_str(), &v.ctx(), MPI_COMM_WORLD,
				 superbblas::FastToSlow);
	superbblas::trsm<Nv, Nw, N>(
	  w.scalar / std::sqrt(v.scalar) / scalar, //
	  v.p->p.data(), v.dim, 1, v.order.c_str(), (const value_type**)&v_ptr, order_rows.c_str(),
	  order_cols.c_str(),
	  &v.ctx(),									  //
	  w.p->p.data(), w.dim, 1, w.order.c_str(), (const value_type**)&w_ptr, &w.ctx(), //
	  p->p.data(), dim, 1, order.c_str(), &ptr, &ctx(), MPI_COMM_WORLD, superbblas::FastToSlow);

	// Force synchronization in superbblas stream if the destination allocation isn't managed by superbblas
	if (!is_managed())
	  superbblas::sync(ctx());
      }

      /// Compute the inverse of `v'
      /// \param v: tensor to compute the Cholesky factor
      /// \param order_rows: labels that are rows of the matrices to factor
      /// \param order_cols: labels that are columns of the matrices to factor
      /// \param w: the other tensor to contract

      template <std::size_t Nv>
      void inv(Tensor<Nv, T> v, const std::string& order_rows, const std::string& order_cols)
      {
	if (is_eg() || v.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	if (isSubtensor() || scalar != T{1} || conjugate)
	{
	  Tensor<N, T> aux = make_compatible();
	  aux.inv(v, order_rows, order_cols);
	  aux.copyTo(*this);
	  return;
	}

	v.copyTo(*this);
	value_type* ptr = data_for_writing();
	superbblas::inversion<Nv>(p->p.data(), dim, 1, order.c_str(), &ptr, order_rows.c_str(),
				  order_cols.c_str(), &ctx(), MPI_COMM_WORLD,
				  superbblas::FastToSlow);

	// Force synchronization in superbblas stream if the destination allocation isn't managed by superbblas
	if (!is_managed())
	  superbblas::sync(ctx());
      }

      /// Solve the linear systems within tensor `v' and right-hand-sides `w`
      /// \param v: tensor to compute the Cholesky factor
      /// \param order_rows: labels that are rows of the matrices to factor
      /// \param order_cols: labels that are columns of the matrices to factor
      /// \param w: the other tensor to contract

      template <std::size_t Nv, std::size_t Nw>
      void solve(Tensor<Nv, T> v, const std::string& order_rows, const std::string& order_cols,
		 Tensor<Nw, T> w)
      {
	if (is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// Conjugacy isn't supported
	if (v.conjugate || w.conjugate || conjugate)
	  throw std::runtime_error("solve: Unsupported implicit conjugate tensors");

	// If either v or w is on OnDevice, force both to be on device
	if (v.ctx().plat != w.ctx().plat)
	{
	  if (v.getDev() != OnDefaultDevice)
	    v = v.cloneOn(OnDefaultDevice);
	  if (w.getDev() != OnDefaultDevice)
	    w = w.cloneOn(OnDefaultDevice);
	}

	// Superbblas tensor contraction is shit and those not deal with subtensors or contracting a host and
	// device tensor (for now)
	if (v.isSubtensor())
	  v = v.clone();
	if (w.isSubtensor())
	  w = w.clone();
	if (isSubtensor() || getDev() != v.getDev())
	{
	  Tensor<N, T> aux = make_compatible(none, {}, v.getDev());
	  aux.solve(v, order_rows, order_cols, w);
	  aux.copyTo(*this);
	  return;
	}

	if ((v.dist == Local) != (w.dist == Local) || (w.dist == Local) != (dist == Local))
	  throw std::runtime_error("solve: One of the contracted tensors or the output tensor is "
				   "local and others are not!");

	// Help superbblas to get the same verbatim value in all processes for the same tensor element in all
	// replicated copies
	// TODO: check whether superbblas does this already
	if ((v.dist == OnMaster || v.dist == OnEveryoneReplicated) ||
	    (w.dist == OnMaster && w.dist == OnEveryoneReplicated))
	{
	  v = v.make_sure(none, none, OnMaster);
	  w = w.make_sure(none, none, OnMaster);
	}

	if (detail::isDistributedOnEveryone(v.dist) && w.dist == OnEveryoneReplicated)
	  w = w.make_suitable_for_contraction(v);

	if (v.dist == OnEveryoneReplicated && detail::isDistributedOnEveryone(w.dist))
	  v = v.make_suitable_for_contraction(w);

	value_type* v_ptr = v.data();
	value_type* w_ptr = w.data();
	value_type* ptr = data_for_writing();
	superbblas::gesm<Nv, Nw, N>(
	  w.scalar / v.scalar / scalar, //
	  v.p->p.data(), v.dim, 1, v.order.c_str(), (const value_type**)&v_ptr, order_rows.c_str(),
	  order_cols.c_str(),
	  &v.ctx(),									  //
	  w.p->p.data(), w.dim, 1, w.order.c_str(), (const value_type**)&w_ptr, &w.ctx(), //
	  p->p.data(), dim, 1, order.c_str(), &ptr, &ctx(), MPI_COMM_WORLD, superbblas::FastToSlow);

	// Force synchronization in superbblas stream if the destination allocation isn't managed by superbblas
	if (!is_managed())
	  superbblas::sync(ctx());
      }

      /// Return a view of this tensor where the elements are scaled by the given argument
      /// \param s: scaling factor
      /// \return: a new view (it doesn't create a copy of the tensor)

      Tensor<N, T> scale(value_type s) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");
	return Tensor<N, T>(*this, scalar * detail::cond_conj(conjugate, s), conjugate);
      }

      /// Return a view of this tensor where the elements are conjugated
      /// \return: a new view (it doesn't create a copy of the tensor)

      Tensor<N, T> conj() const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// NOTE: don't conjugate `scalar`: it's a value associated to the allocation, NOT the view
	return Tensor<N, T>(*this, scalar, !conjugate);
      }

      void release()
      {
	dim = {{}};
	allocation.reset();
	p.reset();
	from = {{}};
	size = {{}};
	strides = {{}};
	scalar = value_type{0};
	conjugate = false;
	eg = false;
	unordered_writing = false;
	complexLabel = 0;
      }

      // Return whether the current view is contiguous in memory
      bool isContiguous() const
      {
	// Meaningless for tensors not been fully supported on a single node
	if (dist != OnMaster && dist != Local)
	  return false;

	if (volume() > 0 && N > 1)
	{
	  bool non_full_dim = false; // some dimension is not full
	  for (unsigned int i = 0; i < N - 1; ++i)
	  {
	    if (from[i] != 0 || size[i] != dim[i])
	    {
	      if (non_full_dim && size[i] != 1)
		return false;
	      non_full_dim = true;
	    }
	  }
	}
	return true;
      }

      /// Return a copy of this tensor if it does not have the same type, the same order, or is not on the given device or distribution
      /// \param new_order: dimension labels of the new tensor
      /// \param new_dev: device
      /// \param new_dist: distribution
      /// \tparam Tn: new precision

      template <typename Tn = T,
		typename std::enable_if<std::is_same<T, Tn>::value, bool>::type = true>
      Tensor<N, Tn> make_sure(const Maybe<std::string>& new_order = none,
			      Maybe<DeviceHost> new_dev = none,
			      Maybe<Distribution> new_dist = none) const
      {
	if (new_order.getSome(order) != order ||
	    !detail::is_same(new_dev.getSome(getDev()), getDev()) || new_dist.getSome(dist) != dist)
	{
	  Tensor<N, Tn> r = new_dist.getSome(dist) != dist
			      ? like_this(new_order, {}, new_dev, new_dist)
			      : make_compatible(new_order, {}, new_dev);
	  if (is_eg())
	  {
	    r = r.make_eg();
	  }
	  else
	  {
	    r.conjugate = conjugate;
	    r.unordered_writing = unordered_writing;
	    copyTo(r);
	  }
	  return r;
	}
	else
	{
	  return *this;
	}
      }

      template <typename Tn = T,
		typename std::enable_if<!std::is_same<T, Tn>::value, bool>::type = true>
      Tensor<N, Tn> make_sure(const Maybe<std::string>& new_order = none,
			      Maybe<DeviceHost> new_dev = none,
			      Maybe<Distribution> new_dist = none) const
      {
	Tensor<N, Tn> r = like_this<N, Tn>(new_order, {}, new_dev, new_dist);
	if (is_eg())
	{
	  r = r.make_eg();
	}
	else
	{
	  r.conjugate = conjugate;
	  r.unordered_writing = unordered_writing;
	  copyTo(r);
	}
	return r;
      }

      /// Return a copy of this tensor in a different type or this tensor if the type coincides
      /// \tparam Tn: new precision

      template <typename Tn = T,
		typename std::enable_if<std::is_same<T, Tn>::value, bool>::type = true>
      Tensor<N, Tn> cast() const
      {
	return *this;
      }

      template <typename Tn = T, typename std::enable_if<!std::is_same<T, Tn>::value &&
							   detail::is_diycomplex<T>::value ==
							     detail::is_diycomplex<Tn>::value,
							 bool>::type = true>
      Tensor<N, Tn> cast() const
      {
	auto r = make_compatible<N, Tn>();
	if (is_eg())
	{
	  r = r.make_eg();
	}
	else
	{
	  r.conjugate = conjugate;
	  r.unordered_writing = unordered_writing;
	  copyTo(r);
	}
	return r;
      }

      template <typename Tn = T, typename std::enable_if<
				   !std::is_same<T, Tn>::value && detail::is_diycomplex<T>::value &&
				     !detail::is_diycomplex<Tn>::value &&
				     std::is_same<typename detail::base_type<T>::type, Tn>::value,
				   bool>::type = true>
      Tensor<N, Tn> cast() const
      {
	return Tensor<N, Tn>(order, dim, allocation, p, dist, from, size, scalar, conjugate, eg,
			     unordered_writing, 0 /* no complexity label */);
      }

      /// Return a compatible tensor in a different type or this tensor if the type coincides
      /// \tparam Tn: new precision

      template <typename Tn = T,
		typename std::enable_if<std::is_same<T, Tn>::value, bool>::type = true>
      Tensor<N, Tn> cast_like() const
      {
	return *this;
      }

      template <typename Tn = T,
		typename std::enable_if<!std::is_same<T, Tn>::value, bool>::type = true>
      Tensor<N, Tn> cast_like() const
      {
	auto r = make_compatible<N, Tn>();
	if (is_eg())
	{
	  r = r.make_eg();
	}
	else
	{
	  r.conjugate = conjugate;
	  r.unordered_writing = unordered_writing;
	}
	return r;
      }

      /// Extend the support of each dimension by the given amount in each direction
      /// \param m: amount to extend the support for each process
      /// \return a new tensor with the extension

      Tensor<N, T> extend_support(const std::map<char, int>& m) const
      {
	Tensor<N, T> r{
	  order,
	  dim,
	  getDev(),
	  dist,
	  std::make_shared<detail::TensorPartition<N>>(p->extend_support(kvcoors<N>(order, m, 0))),
	  unordered_writing,
	  complexLabel};
	r.from = from;
	r.size = size;
	r.strides = strides;
	r.scalar = scalar;
	r.conjugate = conjugate;
	if (is_eg())
	  r = r.make_eg();
	else
	  copyTo(r);
	return r;
      }

      /// Get where the tensor is stored

      DeviceHost getDev() const
      {
#  ifdef SUPERBBLAS_USE_GPU
	return (ctx().plat != superbblas::CPU ? OnDefaultDevice : OnHost);
#  else
	return OnDefaultDevice;
#  endif
      }

      void binaryRead(BinaryReader& bin)
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");
	if (ctx().plat != superbblas::CPU)
	  throw std::runtime_error("Only supported to read on `OnHost` tensors");
	if (dist != OnMaster)
	  throw std::runtime_error("Only supported to read on `OnMaster` tensors");
	if (!isContiguous())
	  throw std::runtime_error("Only supported contiguous views in memory");
	if (scalar != value_type{1} || conjugate)
	  throw std::runtime_error(
	    "Not allowed for tensor with a scale not being one or implicitly conjugated");

	// Only on primary node read the data
	std::size_t vol = volume();
	std::size_t disp = detail::coor2index<N>(from, dim, strides);
	std::size_t word_size = sizeof(typename detail::WordType<T>::type);
	bin.readArrayPrimaryNode((char*)&data_for_writing()[disp], word_size,
				 sizeof(T) / word_size * vol);
      }

      void binaryWrite(BinaryWriter& bin) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// If the writing is collective, the root process needs to hold the whole tensor
	if (!bin.isLocal() && dist != OnMaster)
	  throw std::runtime_error("For collective writing, the tensor should be `OnMaster`");

	// If the writing is non-collective, the tensor should be local
	if (bin.isLocal() && dist != Local)
	  throw std::runtime_error("For non-collective writing, the tensor should be `Local`");

	if (conjugate)
	  throw std::runtime_error("Not allowed for tensors implicitly conjugated");

	// If the tensor has an implicit scale, view, or is not on host, make a copy
	if (scalar != value_type{1} || isSubtensor() || ctx().plat != superbblas::CPU)
	{
	  cloneOn(OnHost).binaryWrite(bin);
	  return;
	}

	// Write the local data
	std::size_t vol = p->localVolume();
	std::size_t word_size = sizeof(typename detail::WordType<T>::type);
	bin.writeArrayPrimaryNode((char*)data(), word_size, sizeof(T) / word_size * vol);
      }

      void print(const std::string& name) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	std::stringstream ss;
	auto t = toComplex();
	auto t_host = t.like_this(none, {}, OnHost, OnMaster);
	t.copyTo(t_host);
	if (Layout::nodeNumber() == 0)
	{
	  using namespace detail::repr;
	  ss << "% " << repr(data()) << std::endl;
	  ss << "% dist=" << p->p << std::endl;
	  ss << name << "=reshape([";
	  assert(!t_host.isSubtensor());
	  std::size_t vol = volume();
	  for (std::size_t i = 0; i < vol; ++i)
	  {
	    //using detail::repr::operator<<;
	    ss << " ";
	    detail::repr::operator<<(ss, t_host.data()[i]);
	  }
	  ss << "], [" << size << "]);" << std::endl;
	}
	detail::log(1, ss.str());
      }
#  if 0
      /// Get where the tensor is stored
      void setNan() 
      {
#    ifndef QDP_IS_QDPJIT
	T nan = detail::NaN<T>::get();
	std::size_t vol = superbblas::detail::volume<N>(dim);
	T* p = data.get();
	for (std::size_t i = 0; i < vol; i++)
	  p[i] = nan;
#    endif
    }

    void
    checkNan() const
    {
#    ifndef QDP_IS_QDPJIT
	std::size_t vol = superbblas::detail::volume<N>(dim);
	T* p = data.get();
	for (std::size_t i = 0; i < vol; i++)
	  assert(detail::IsFinite<T>::get(p[i]));
#    endif
      }
#  endif
    };

    /// Copy v.kvslice_from_size(dir) into w.kvslice_from_size(dir) for every dir in disps
    /// assuming that v and w have even-odd layout (if 'X' == 2) and the displacements `disps`
    /// are given in natural coordinates.
    ///
    /// \param v: origin tensor
    /// \param w: destination tensor
    /// \param label_mu: if given, copy each displacement into a separate coordinate.
    /// \param disps: displacements in natural coordinates
    /// \param real_dims: even-odd dimension of the original lattice
    /// \param even_mask: mask into the elements with even x natural coordinate
    /// \param odd_mask: mask into the elements with odd x natural coordinate

    template <std::size_t Nv, typename Tv, std::size_t Nw, typename Tw, std::size_t Nm = Nv,
	      typename Tm = float>
    void latticeCopyToWithMask(const Tensor<Nv, Tv>& v, const Tensor<Nw, Tw>& w, char label_mu,
			       const std::vector<Coor<Nd>>& disps,
			       const std::map<char, int>& real_dims,
			       const Tensor<Nm, Tm>& mask_even, const Tensor<Nm, Tm>& mask_odd)
    {
      // Shortcuts
      if (disps.size() == 0)
	return;

      // Make sure that v is distributed as w
      if (!w.isDistributedAs(v, "xyztX"))
      {
	auto v_ = w.template make_compatible<Nv, Tv>(v.order, v.kvdim());
	v.copyTo(v_);
	latticeCopyToWithMask(v_, w, label_mu, disps, real_dims, mask_even, mask_odd);
	return;
      }

      // Get the number of colors on the original lattice
      const auto dim = v.kvdim();
      int real_maxX = real_dims.count('X') == 1 ? real_dims.at('X') : dim.at('X');

      // Preallocate the masks for v and w
      Tensor<Nv, float> v_mask = v.create_mask();
      Tensor<Nw, float> w_mask = w.create_mask();

      for (unsigned int mu = 0; mu < disps.size(); ++mu)
      {
	const auto& dir = disps[mu];
	int sumdir = std::accumulate(dir.begin(), dir.end(), int{0});

	for (int x = 0; x < 2; ++x)
	{
	  // Nat coor (x+dirx,Y+diry,Z+dirz,T+dirt) to even-odd coordinate
	  std::map<char, int> from{{'x', x / real_maxX}};
	  std::map<char, int> to{{'X', sumdir},
				 {'x', (dir[0] + dim.at('x') * real_maxX + x) / real_maxX},
				 {'y', dir[1]},
				 {'z', dir[2]},
				 {'t', dir[3]}};

	  // Restrict the destination tensor to label_mu if given
	  auto w_mu = (label_mu == 0 ? w : w.kvslice_from_size({{label_mu, mu}}, {{label_mu, 1}}));

	  auto mask_mu = (x == 0 ? mask_even : mask_odd).kvslice_from_size(from, {});
	  auto v_mask_mu = v_mask.kvslice_from_size(to, {});
	  auto w_mask_mu =
	    (label_mu == 0 ? w_mask : w_mask.kvslice_from_size({{label_mu, mu}}, {{label_mu, 1}}));
	  w_mask_mu = w_mask_mu.kvslice_from_size(to, {});
	  mask_mu.copyTo(v_mask_mu);
	  mask_mu.copyTo(w_mask_mu);
	  v.kvslice_from_size(to, {}).copyToWithMask(w_mu.kvslice_from_size(to, {}), v_mask_mu,
						     w_mask_mu);
	} // x
      }	  // mu
    }

    /// Return an identity matrix
    /// \param dim: length for each of the row dimensions
    /// \param m: labels map from the row to the column dimensions and other dimensions

    template <std::size_t N, typename T>
    Tensor<N, T> identity(const std::map<char, int>& dim, const remap& m,
			  const Distribution& dist = OnEveryone)
    {
      using value_type = typename detail::base_type<T>::type;

      // Get the order for the rows
      std::string orows;
      for (const auto& it : m)
	orows.push_back(it.first);

      // Get the order for the columns
      std::string ocols = detail::update_order(orows, m);

      // Get the extra dimensions
      std::string ot;
      for (const auto& it : dim)
	if (m.count(it.first) == 0)
	  ot.push_back(it.first);
      ot = detail::remove_dimensions(ot, ocols);

      // Get the dimensions of the identity tensor
      std::map<char, int> iden_dim;
      for (const auto& it : dim) {
	iden_dim[it.first] = (detail::is_in(ot, it.first) ? 1 : it.second);
	if (detail::is_in(orows, it.first)) iden_dim[m.at(it.first)] = it.second;
      }

      // Create the identity tensor
      const std::string order = orows + ocols + ot;
      Tensor<N, T> iden{order, kvcoors<N>(order, iden_dim, 0, ThrowOnMissing), OnHost,
		     OnEveryoneReplicated};
      iden.set_zero();
      if (iden.getLocal())
      {
	value_type* p = iden.getLocal().data();
	for (unsigned int i = 0, vol = detail::volume(dim, orows); i < vol; ++i)
	  p[vol * i + i] = value_type{1};
      }

      // Get the dimensions of the returned tensor
      std::map<char, int> t_dim;
      for (const auto& it : dim) {
	t_dim[it.first] = (!detail::is_in(ot, it.first) ? 1 : it.second);
	if (detail::is_in(orows, it.first)) t_dim[m.at(it.first)] = 1;
      }
      Tensor<N, T> t{order, kvcoors<N>(order, t_dim, 0, ThrowOnMissing), OnDefaultDevice, dist};
      t.set(1);

      std::map<char, int> r_dim = dim;
      for (const auto& it : dim)
	if (detail::is_in(orows, it.first)) r_dim[m.at(it.first)] = it.second;
      Tensor<N, T> r{order, kvcoors<N>(order, r_dim, 0, ThrowOnMissing), OnDefaultDevice, dist};

      kronecker(t, iden, r);
      return r;
    }

    /// Contract some dimension of the given tensors
    /// \param v: one tensor to contract
    /// \param w: the other tensor to contract
    /// \param labels_to_contract: labels dimensions to contract from `v` and `w`
    /// \param action: either to copy or add to the given output tensor if given
    /// \param r: optional given tensor where to put the resulting contraction
    /// \param mr: map from the given `r` to the labels of the contraction
    /// \param beta: scale on `r` if the `action` in `AddTo`
    /// \param dev: device for the resulting tensor if `action` isn't given
    /// \param dist: distribution of the resulting tensor if `action` isn't given
    ///
    /// Example:
    ///
    ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}), q("Ss", {{Ns,Ns}});
    ///   Tensor<2,Complex> r0 = contract<2>(t, q, "s"); // r0 dims are "cS"
    ///   Tensor<3,Complex> r1 = contract<3>(t, q, ""); // r1 dims are "csS"
    ///   Tensor<2,Complex> r2("cS", {{Nc,Ns}});
    ///   contract<2>(t, q, "s", CopyTo, r2); // r2 = q * s
    ///   Tensor<2,Complex> r3("cs", {{Nc,Ns}});
    ///   contract<2>(t, q, "s", CopyTo, r3, {{'s','S'}}); // r2 = q * s
    ///   contract<2>(t, q.rename_dims({{'s','S'},{'S','s'}}).conj(), "s", CopyTo, r3, {{'s','S'}}); // r2 = q * s^*

    template <std::size_t Nr, std::size_t Nv, std::size_t Nw, typename T>
    Tensor<Nr, T>
    contract(const Tensor<Nv, T>& v, Tensor<Nw, T> w, const std::string& labels_to_contract,
	     Maybe<Action> action = none, Tensor<Nr, T> r = Tensor<Nr, T>{}, const remap& mr = {},
	     typename detail::base_type<T>::type beta = typename detail::base_type<T>::type{1},
	     Maybe<DeviceHost> dev = none, Maybe<Distribution> dist = none)
    {
      // Check arguments
      if (action.hasSome() != (bool)r)
	throw std::runtime_error(
	  "contract: invalid argument, if `action` is given, `r` should be given also");
      if ((dev.hasSome() || dist.hasSome()) && action.hasSome())
	throw std::runtime_error(
	  "contract: invalid argument, if `action` is given, `dev` and `dist` shouldn't be given");

      // Compute the labels of the output tensor: v.order + w.order - labels_to_contract
      std::string rorder = detail::union_dimensions(v.order, w.order, labels_to_contract);
      if (Nr != rorder.size())
	throw std::runtime_error(
	  "contract: The dimension of the output tensor does not match the template argument");
      if ((bool)r && detail::union_dimensions(rorder, r.order) != rorder)
	throw std::runtime_error("contract: The given output tensor has an unexpected ordering");

      // If any of the input tensors is glocal, make sure both are
      if ((v.dist == Glocal) != (w.dist == Glocal))
      {
	Tensor<Nv, T> v0 = v;
	Tensor<Nw, T> w0 = w;
	if (v.dist != Glocal)
	  v0 = v.getGlocal();
	if (w.dist != Glocal)
	  w0 = w.getGlocal();
	return contract(v0, w0, labels_to_contract, action, r, mr, beta, dev, dist);
      }

      // If the output tensor is not given create a new one
      Tensor<Nr, T> r0;
      if (!r)
      {
	r0 = (v.dist != Glocal && (dev.hasSome() || dist.hasSome()))
	       ? v.template like_this<Nr>(rorder, w.kvdim(), dev, dist)
	       : v.template make_compatible<Nr>(rorder, w.kvdim());
	beta = 0;
      }
      else
      {
	r0 = r;
      }

      // Correct beta for the action
      if (action.hasSome() && action.getSome() == CopyTo)
	beta = 0.0;

      // Do the contraction
      r0.contract(v, {}, NotConjugate, w, {}, NotConjugate, mr, beta);

      return r0;
    }

    /// Contract some dimension of the given tensors
    /// \param v: one tensor to contract
    /// \param w: the other tensor to contract
    /// \param labels_to_contract: labels dimensions to contract from `v` and `w`
    /// \param dev: device for the resulting tensor
    /// \param dist: distribution of the resulting tensor
    ///
    /// Example:
    ///
    ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}), q("Ss", {{Ns,Ns}});
    ///   Tensor<2,Complex> r0 = contract<2>(t, q, "s"); // r0 dims are "cS"
    ///   Tensor<3,Complex> r1 = contract<3>(t, q, ""); // r1 dims are "csS"
    ///   Tensor<3,Complex> r2 = contract<3>(t, q, "", OnMaster); // r2 is supported on master

    template <std::size_t Nr, std::size_t Nv, std::size_t Nw, typename T>
    Tensor<Nr, T> contract(const Tensor<Nv, T>& v, Tensor<Nw, T> w,
			   const std::string& labels_to_contract, Maybe<DeviceHost> dev,
			   Maybe<Distribution> dist)
    {
      return contract<Nr>(v, w, labels_to_contract, none, {}, {}, 0, dev, dist);
    }

    /// Contract some dimension of the given tensors
    /// \param v: one tensor to contract
    /// \param w: the other tensor to contract
    /// \param labels_to_contract: map of labels dimensions to contract from `v` to `w`
    /// \param action: either to copy or add to the given output tensor if given
    /// \param r: optional given tensor where to put the resulting contraction
    /// \param mr: map from the given `r` to the labels of the contraction
    /// \param beta: scale on `r` if the `action` in `AddTo`
    /// \param dev: device for the resulting tensor if `action` isn't given
    /// \param dist: distribution of the resulting tensor if `action` isn't given
    ///
    /// Example:
    ///
    ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}), q("St", {{Ns,Ns}});
    ///   Tensor<2,Complex> r0 = contract<2>(t, q, {{'s','t'}}); // r0 dims are "cS"
    ///   Tensor<2,Complex> r2("cS", {{Nc,Ns}});
    ///   contract<2>(t, q, {{'s','t'}}, CopyTo, r2); // r2 = q * s
    ///   Tensor<2,Complex> r3("cs", {{Nc,Ns}});
    ///   contract<2>(t, q, {{'s','t'}}, CopyTo, r3, {{'s','S'}}); // r2 = q * s
    ///   contract<2>(t, q.rename_dims({{'s','S'},{'S','s'}}).conj(), {{'s','t'}}, CopyTo, r3, {{'s','S'}}); // r2 = q * s^*

    template <std::size_t Nr, std::size_t Nv, std::size_t Nw, typename T>
    Tensor<Nr, T>
    contract(const Tensor<Nv, T>& v, Tensor<Nw, T> w, const remap& labels_to_contract,
	     Maybe<Action> action = none, Tensor<Nr, T> r = Tensor<Nr, T>{}, const remap& mr = {},
	     typename detail::base_type<T>::type beta = typename detail::base_type<T>::type{1},
	     Maybe<DeviceHost> dev = none, Maybe<Distribution> dist = none)
    {
      // Remap the labels to contract from v and w
      std::string labels_to_contract_v;
      for (const auto& it : labels_to_contract)
	labels_to_contract_v.push_back(it.first);
      remap mv = detail::getNewLabels(labels_to_contract_v, v.order + w.order);
      remap mw;
      for (const auto it : mv)
	mw[labels_to_contract.at(it.first)] = it.second;
      std::string labels_to_contract_str;
      for (const auto it : mv)
	labels_to_contract_str.push_back(it.second);
      return contract(v.rename_dims(mv), w.rename_dims(mw), labels_to_contract_str, action, r, mr,
		      beta, dev, dist);
    }

    /// Do the Kronecker product of two tensors
    /// \param v: one tensor to contract
    /// \param w: the other tensor to contract
    /// \param r: optional given tensor where to put the resulting contraction

    template <std::size_t Nr, std::size_t Nv, std::size_t Nw, typename T>
    Tensor<Nr, T> kronecker(const Tensor<Nv, T>& v, const Tensor<Nw, T>& w,
			    Tensor<Nr, T> r = Tensor<Nr, T>{})
    {
      // Make sure that no dimension in common has size larger than one in both tensors
      auto v_kvdim = v.kvdim();
      auto w_kvdim = w.kvdim();
      for (const auto& it : v_kvdim)
	if (it.second > 1 && w_kvdim.count(it.first) > 0 && w_kvdim.at(it.first) > 1)
	  throw std::runtime_error(
	    "kronecker: input tensors have a common dimension with size larger than one");

      // Renamed all dimensions in w to avoid a common label between the tensors
      remap w_m = detail::getNewLabels(w.order, v.order);

      // Do the contraction
      auto k = contract<Nv + Nw>(v, w.rename_dims(w_m), "");

      // The labels of the output tensor are the union of the input tensors labels
      std::string rorder = detail::union_dimensions(v.order, w.order);

      // The output tensor has the maximum size of the input tensors
      auto rdims = v_kvdim;
      for (const auto& it : w_kvdim)
	if (rdims.count(it.first) == 0 || rdims.at(it.first) == 1)
	  rdims[it.first] = it.second;

      // For the common labels, rename the singleton one
      remap k_m;
      for (const auto& it : w_kvdim)
      {
	if (v_kvdim.count(it.first) == 1 && it.second > 1)
	{
	  k_m[w_m.at(it.first)] = it.first;
	  k_m[it.first] = w_m.at(it.first);
	}
	else if (v_kvdim.count(it.first) == 0)
	{
	  k_m[w_m.at(it.first)] = it.first;
	}
      }

      // If the output tensor is not given create a new one
      if (!r)
      {
	auto r0 = v.template make_compatible<Nr>(rorder, rdims);
	k.rename_dims(k_m).copyTo(r0);
	return r0;
      }
      else
      {
	k.rename_dims(k_m).copyTo(r);
	return r;
      }
    }

    /// Compute the norm along some dimensions
    /// \param v: tensor
    /// \param order_t: labels not to contract (optional)
    /// \param order_rows: labels to contract (optional, either order_rows or order_t
    ///        should be provided)
    ///
    /// Example:
    ///
    ///   Tensor<2,Complex> t("cs", {{Nc,Ns}}), q("Ss", {{Ns,Ns}});
    ///   Tensor<2,Complex> r0 = contract<2>(t, q, "s"); // r0 dims are "cS"
    ///   Tensor<3,Complex> r1 = contract<3>(t, q, ""); // r1 dims are "csS"
    ///   Tensor<2,Complex> r2("cS", {{Nc,Ns}});
    ///   contract<2>(t, q, "s", CopyTo, r2); // r2 = q * s
    ///   Tensor<2,Complex> r3("cs", {{Nc,Ns}});
    ///   contract<2>(t, q, "s", CopyTo, r3, {{'s','S'}}); // r2 = q * s
    ///   contract<2>(t, q.rename_dims({{'s','S'},{'S','s'}}).conj(), "s", CopyTo, r3, {{'s','S'}}); // r2 = q * s^*

    template <std::size_t Nr, std::size_t Nv, typename T>
    Tensor<Nr, typename detail::real_type<T>::type> norm(const Tensor<Nv, T>& v,
							 Maybe<std::string> order_t = none,
							 Maybe<std::string> order_rows = none)
    {
      if (!order_t.hasSome() && !order_rows.hasSome())
	throw std::runtime_error(
	  "norm: invalid input, give at least either `order_t` or `order_rows`");

      // Compute the labels to contract
      std::string rorder = order_rows.hasSome()
			     ? order_rows.getSome()
			     : detail::remove_dimensions(v.order, order_t.getSome());
      std::string torder =
	order_t.hasSome() ? order_t.getSome() : detail::remove_dimensions(v.order, rorder);

      // Allocate the output on the host and spread the result to every process
      auto r = contract<Nr>(v.conj(), v, rorder, OnHost, OnEveryoneReplicated).reorder(torder);

      // Do the square root and return the result
      using Treal = typename detail::real_type<T>::type;
      return r.template transformWithCPUFun<Treal>(
	[](const typename detail::base_type<T>::type& t) { return std::sqrt(std::real(t)); });
    }

    /// Compute the Cholesky factor of `v' and contract its inverse with `w`
    /// \param v: tensor to compute the Cholesky factor
    /// \param order_rows: labels that are rows of the matrices to factor
    /// \param order_cols: labels that are columns of the matrices to factor
    /// \param w: the other tensor to contract
    /// \param labels_to_contract: labels dimensions to contract from `v` and `w`
    /// \param action: either to copy or add to the given output tensor if given (only `CopyTo' supported)
    /// \param r: optional given tensor where to put the resulting contraction

    template <std::size_t Nr, std::size_t Nv, std::size_t Nw, typename T>
    Tensor<Nr, T> cholInv(const Tensor<Nv, T>& v, const std::string& order_rows,
			  const std::string& order_cols, const Tensor<Nw, T>& w,
			  const std::string& labels_to_contract, Maybe<Action> action = none,
			  Tensor<Nr, T> r = {})
    {
      if (action.hasSome() != (bool)r)
	throw std::runtime_error("Invalid default value");

      // Compute the labels of the output tensor: v.order + w.order - labels_to_contract
      std::string rorder = detail::union_dimensions(v.order, w.order, labels_to_contract);
      if (Nr != rorder.size())
	throw std::runtime_error(
	  "cholInv: The dimension of the output tensor does not match the template argument");
      if (r && detail::union_dimensions(rorder, r.order) != rorder)
	throw std::runtime_error("cholInv: The given output tensor has an unexpected ordering");

      // If the output tensor is not given create a new one
      Tensor<Nr, T> r0;
      if (!r)
      {
	r0 = v.template like_this<Nr>(rorder, w.kvdim());
      }
      else
      {
	r0 = r;
      }

      // For now, only `CopyTo' action is supported
      if (action.hasSome() && action.getSome() != CopyTo)
	throw std::runtime_error("cholInv: unsupported action");

      // Do the contraction
      r0.cholInv(std::move(v), order_rows, order_cols, w);

      return r0;
    }

    /// Solve the linear systems within tensor `v' and right-hand-sides `w`
    /// \param v: tensor to compute the Cholesky factor
    /// \param order_rows: labels that are rows of the matrices to factor
    /// \param order_cols: labels that are columns of the matrices to factor
    /// \param w: the other tensor to contract
    /// \param labels_to_contract: labels dimensions to contract from `v` and `w`
    /// \param action: either to copy or add to the given output tensor if given (only `CopyTo' supported)
    /// \param r: optional given tensor where to put the resulting contraction

    template <std::size_t Nlabels, std::size_t Nr, std::size_t Nv, std::size_t Nw, typename T>
    Tensor<Nr, T> solve(const Tensor<Nv, T>& v, const std::string& order_rows,
			const std::string& order_cols, const Tensor<Nw, T>& w,
			const std::string& labels_to_contract, Maybe<Action> action = none,
			Tensor<Nr, T> r = {})
    {
      if (action.hasSome() != (bool)r)
	throw std::runtime_error("solve: Invalid default value");

      // Compute the labels of the output tensor: v.order + w.order - labels_to_contract
      std::string rorder = detail::union_dimensions(v.order, w.order, labels_to_contract);
      if (Nr != rorder.size())
	throw std::runtime_error(
	  "solve: The dimension of the output tensor does not match the template argument");
      if (r && detail::union_dimensions(rorder, r.order) != rorder)
	throw std::runtime_error("solve: The given output tensor has an unexpected ordering");
      if (Nlabels != labels_to_contract.size())
	throw std::runtime_error(
	  "solve: The length of `order_rows` does not match the template argument `Nrows`");
      if (order_rows.size() != order_cols.size())
	throw std::runtime_error("solve: unsupported ordering for the matrix");

      // If the output tensor is not given create a new one
      Tensor<Nr, T> r0;
      if (!r)
      {
	r0 = v.template like_this<Nr>(rorder, w.kvdim());
      }
      else
      {
	r0 = r;
      }

      // For now, only `CopyTo' action is supported
      if (action.hasSome() && action.getSome() != CopyTo)
	throw std::runtime_error("solve: unsupported action");

      // Compute the solution
      r0.solve(v, order_rows, order_cols, w);

      // Check the solution
      if (superbblas::getDebugLevel() > 0)
      {
	auto res = w.clone().scale(-1);
	remap m{};
	for (unsigned int i = 0; i < order_rows.size(); ++i)
	{
	  m[order_rows[i]] = order_cols[i];
	  m[order_cols[i]] = order_rows[i];
	}
	contract<Nw>(v, r0.rename_dims(m), labels_to_contract, AddTo, res.rename_dims(m));
	auto wnorms = norm<Nw - Nlabels>(w, none, order_cols);
	auto rnorms = norm<Nw - Nlabels>(res, none, order_cols);
	double err = 0;
	for (int i = 0, i1 = wnorms.volume(); i < i1; ++i)
	  err = std::max(err, (double)rnorms.data()[i] / wnorms.data()[i]);
	QDPIO::cout << "solve error: " << detail::tostr(err) << std::endl;
	auto eps = std::sqrt(std::numeric_limits<typename detail::real_type<T>::type>::epsilon());
	if (err > eps)
	  throw std::runtime_error(std::string("solve: too much error in dense solution, ") +
				   detail::tostr(err));
      }

      return r0;
    }

    /// Invert the matrices
    /// \param v: tensor to compute the inversion
    /// \param order_rows: labels that are rows of the matrices to factor
    /// \param order_cols: labels that are columns of the matrices to factor
    /// \param r: optional given tensor where to put the resulting contraction

    template <std::size_t N, typename T>
    Tensor<N, T> inv(const Tensor<N, T>& v, const std::string& order_rows,
		       const std::string& order_cols, Tensor<N, T> r = {})
    {
      if (r && detail::union_dimensions(v.order, r.order) != v.order)
	throw std::runtime_error("inv: The given output tensor has an unexpected ordering");
      if (order_rows.size() != order_cols.size())
	throw std::runtime_error("inv: unsupported ordering for the matrix");

      // If the output tensor is not given create a new one
      Tensor<N, T> r0;
      if (!r)
      {
	r0 = v.make_compatible();
      }
      else
      {
	r0 = r;
      }

      // Compute the solution
      r0.inv(v, order_rows, order_cols);

      // Check the solution
      if (superbblas::getDebugLevel() > 0)
      {
	remap m{};
	for (unsigned int i = 0; i < order_rows.size(); ++i)
	  m[order_rows[i]] = order_cols[i];
	auto dim = v.kvdim();
	char c = detail::get_free_label(v.order);
	dim[c] = 1;
	auto res = identity<N + 1, T>(dim, m).scale(-1);
	contract<N + 1>(v, r0.split_dimension(order_rows[0], std::string({c, order_rows[0]}), 1), m,
			AddTo, res);
	auto err = norm<1>(res, std::string(1, c)).get({0});
	QDPIO::cout << "inv error: " << detail::tostr(err) << std::endl;
	auto eps = std::sqrt(std::numeric_limits<typename detail::real_type<T>::type>::epsilon());
	if (err > eps)
	  throw std::runtime_error(std::string("inv: too much error in dense solution, ") +
				   detail::tostr(err));
      }

      return r0;
    }

    /// Elementwise division
    /// \param v: numerator
    /// \param w: denominator

    template <std::size_t N, typename T>
    Tensor<N, T> div(const Tensor<N, T>& v, const Tensor<N, T>& w)
    {
      auto r = v.make_compatible(none, {}, OnHost);
      v.copyTo(r);
      auto w0 = r.make_compatible();
      w.copyTo(w0);
      auto r_local = r.getLocal();
      auto w0_local = w0.getLocal();
      if (r_local)
      {
	auto rptr = r_local.data();
	auto w0ptr = w0_local.data();
	for (std::size_t i = 0, vol = r_local.volume(); i < vol; ++i) {
	  auto w0i =
	    detail::cond_conj(r_local.conjugate != w0_local.conjugate, w0ptr[i] * w0_local.scalar);
	  rptr[i] = std::norm(w0i) == 0 ? T{0} : rptr[i] / w0i;
	}
      }
      return r;
    }

    /// Compute the maximum for a small tensor
    /// \param v: tensor

    template <std::size_t N, typename T>
    typename detail::base_type<T>::type
    max(Tensor<N, T> v, typename detail::base_type<T>::type init =
			  std::numeric_limits<typename detail::base_type<T>::type>::min())
    {
      using value_type = typename detail::base_type<T>::type;
      if (v.dist != Glocal && v.dist != Local)
	v = v.make_sure(none, OnHost, OnEveryoneReplicated);
      else
	v = v.make_sure(none, OnHost);
      if (v.isSubtensor())
	v = v.clone();
      value_type r = init;
      v = v.getLocal();
      value_type* p = v.data();
      for (unsigned int i = 0, vol = v.volume(); i < vol; ++i)
	r = std::max(r, p[i]);
      return r;
    }

    /// Elementwise product
    /// \param v: numerator
    /// \param w: denominator

    template <std::size_t N, typename T>
    Tensor<N, T> mult(const Tensor<N, T>& v, const Tensor<N, T>& w)
    {
      auto r = v.make_compatible(none, {}, OnHost);
      v.copyTo(r);
      auto w0 = r.make_compatible();
      w.copyTo(w0);
      auto r_local = r.getLocal();
      auto w0_local = w0.getLocal();
      if (r_local)
      {
	auto rptr = r_local.data();
	auto w0ptr = w0_local.data();
	for (std::size_t i = 0, vol = r_local.volume(); i < vol; ++i)
	  rptr[i] = rptr[i] * detail::cond_conj(r_local.conjugate != w0_local.conjugate,
						w0ptr[i] * w0_local.scalar);
      }
      return r;
    }

    template <typename T>
    void* getQDPPtr(const T& t)
    {
#  if defined(QDP_IS_QDPJIT) && defined(SUPERBBLAS_USE_GPU)
      multi1d<QDPCache::ArgKey> v(1);
      v[0] = t.getId();
      void* r = QDP_get_global_cache().get_dev_ptrs(v)[0];
      assert(superbblas::detail::getPtrDevice(r) >= 0);
      return r;
#  else
      return t.getF();
#  endif
    }

    template <typename T>
    using LatticeColorVectorT = OLattice<PScalar<PColorVector<RComplex<T>, Nc>>>;

    template <typename T>
    Tensor<Nd + 2, std::complex<T>> asTensorView(const LatticeColorVectorT<T>& v)
    {
      using Complex = std::complex<T>;
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 2, Complex>("cxyztX", latticeSize<Nd + 2>("cxyztX"), OnHost,
				     OnEveryoneAsChroma, v_ptr);
    }

#  if !defined(QDP_IS_QDPJIT) || !defined(SUPERBBLAS_USE_GPU)
    inline Tensor<Nd + 3, Complex> asTensorView(const LatticeFermion& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 3, Complex>("csxyztX", latticeSize<Nd + 3>("csxyztX"), OnHost,
				     OnEveryoneAsChroma, v_ptr);
    }
#  else
    inline Tensor<Nd + 4, DIYComplex<REAL>> asTensorView(const LatticeFermion& v)
    {
      REAL* v_ptr = reinterpret_cast<REAL*>(getQDPPtr(v));
      return Tensor<Nd + 4, DIYComplex<REAL>>("xyztXsc.", latticeSize<Nd + 4>("xyztXsc."),
					      OnDefaultDevice, OnEveryoneAsChroma, v_ptr, '.');
    }
#  endif

#  if !defined(QDP_IS_QDPJIT) || !defined(SUPERBBLAS_USE_GPU)
    inline Tensor<Nd + 1, Complex> asTensorView(const LatticeComplex& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 1, Complex>("xyztX", latticeSize<Nd + 1>("xyztX"), OnHost,
				     OnEveryoneAsChroma, v_ptr);
    }
#  else
    inline Tensor<Nd + 2, DIYComplex<REAL>> asTensorView(const LatticeComplex& v)
    {
      REAL* v_ptr = reinterpret_cast<REAL*>(getQDPPtr(v));
      return Tensor<Nd + 2, DIYComplex<REAL>>("xyztX.", latticeSize<Nd + 2>("xyztX."),
					      OnDefaultDevice, OnEveryoneAsChroma, v_ptr, '.');
    }
#  endif

#  if !defined(QDP_IS_QDPJIT) || !defined(SUPERBBLAS_USE_GPU)
    inline Tensor<Nd + 3, Complex> asTensorView(const LatticeColorMatrix& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 3, Complex>("jixyztX",
				     latticeSize<Nd + 3>("jixyztX", {{'i', Nc}, {'j', Nc}}), OnHost,
				     OnEveryoneAsChroma, v_ptr);
    }
#  else
    inline Tensor<Nd + 4, DIYComplex<REAL>> asTensorView(const LatticeColorMatrix& v)
    {
      REAL* v_ptr = reinterpret_cast<REAL*>(getQDPPtr(v));
      return Tensor<Nd + 4, DIYComplex<REAL>>(
	"xyztXji.", latticeSize<Nd + 4>("xyztXji.", {{'i', Nc}, {'j', Nc}}), OnDefaultDevice,
	OnEveryoneAsChroma, v_ptr, '.');
    }
#  endif

    inline Tensor<Nd + 4, Complex> asTensorView(const LatticeColorVectorSpinMatrix& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 4, Complex>("cjixyztX",
				     latticeSize<Nd + 4>("cjixyztX", {{'i', Ns}, {'j', Ns}}),
				     OnHost, OnEveryoneAsChroma, v_ptr);
    }

    template <typename COMPLEX>
    Tensor<1, COMPLEX> asTensorView(std::vector<COMPLEX>& v,
				    Distribution dist = OnEveryoneReplicated)
    {
      return Tensor<1, COMPLEX>("i", Coor<1>{Index(v.size())}, OnHost, dist, v.data());
    }

    inline Tensor<2, Complex> asTensorView(SpinMatrix& smat)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(smat.getF());
      return Tensor<2, Complex>("ji", Coor<2>{Ns, Ns}, OnHost, OnEveryoneReplicated, v_ptr);
    }

    inline SpinMatrix SpinMatrixIdentity()
    {
      SpinMatrix one;
      // valgrind complains if all elements of SpinMatrix are not initialized!
      for (int i = 0; i < Ns; ++i)
	for (int j = 0; j < Ns; ++j)
	  pokeSpin(one, cmplx(Real(0), Real(0)), i, j);
      for (int i = 0; i < Ns; ++i)
	pokeSpin(one, cmplx(Real(1), Real(0)), i, i);
      return one;
    }

    template <typename COMPLEX = Complex>
    Tensor<2, COMPLEX> Gamma(int gamma, DeviceHost dev = OnDefaultDevice)
    {
      SpinMatrix g = QDP::Gamma(gamma) * SpinMatrixIdentity();
      Tensor<2, COMPLEX> r("ij", {Ns, Ns}, dev, OnEveryoneReplicated);
      asTensorView(g).copyTo(r);
      return r;
    }

    /// Broadcast a string from process zero
    inline std::string broadcast(const std::string& s)
    {
      // Broadcast the size of the string
      std::vector<float> size_orig(1, s.size()), size_dest(1, 0);
      asTensorView(size_orig, OnMaster).copyTo(asTensorView(size_dest));

      // Broadcast the content of the string
      std::vector<float> orig(s.begin(), s.end());
      orig.resize(size_dest[0]);
      std::vector<float> dest(size_dest[0]);
      asTensorView(orig, OnMaster).copyTo(asTensorView(dest));
      return std::string(dest.begin(), dest.end());
    }

    /// Broadcast a string from process zero
    inline int broadcast(int s)
    {
      std::vector<int> v(1, s), dest(1, 0);
      asTensorView(v, OnMaster).copyTo(asTensorView(dest));
      return dest[0];
    }

    /// Broadcast a string from process zero
    template <std::size_t N>
    Coor<N> broadcast(const Coor<N>& c)
    {
      std::vector<int> v(c.begin(), c.end()), dest(N);
      asTensorView(v, OnMaster).copyTo(asTensorView(dest));
      Coor<N> r;
      std::copy_n(dest.begin(), N, r.begin());
      return r;
    }

    /// Broadcast a string from process zero
    template <typename T>
    Maybe<T> broadcast(const Maybe<T>& c)
    {
      int has_something = broadcast(c.hasSome() ? 1 : 0);
      if (has_something == 1)
	return Maybe<T>(broadcast(c ? c.getSome() : T{}));
      return none;
    }

    /// Class for operating sparse tensors
    /// \tparam ND: number of domain dimensions
    /// \tparam NI: number of image dimensions
    /// \tparam T: datatype
    ///
    /// The class may support several variants of Column Sparse Row (CSR) format for representing
    /// sparse matrices, but for now only Block Sparse Row (BSR) with the same number of nonzeros
    /// on all rows is supported. Superbblas has some support for blocked ELL (BSR but with a negative
    /// column index for the unused blocks in a row), but most of the methods of this class aren't ready
    /// for that.
    ///
    /// Besides, this class implements an extension of the BSR in which the nonzero blocks are the result
    /// of the tensor product of two matrices one of them being constant among all edges on the same
    /// direction. This extension is referred as BSR Kronecker. When `kron_data` is given, the nonzeros
    /// should be ordered such that the nonzero blocks with the same `u` label are multiplied by the nonzero
    /// block in `kron_data` with that `u`.

    template <std::size_t ND, std::size_t NI, typename T>
    struct SpTensor {
      using value_type = typename detail::base_type<T>::type;
      static_assert(superbblas::supported_type<value_type>::value, "Not supported type");

    public:
      Tensor<ND, T> d;		   ///< Tensor example for the domain
      Tensor<NI, T> i;		   ///< Tensor example for the image
      Coor<ND> blkd;		   ///< blocking for the domain
      Coor<NI> blki;		   ///< blocking for the image
      Coor<ND> krond;		   ///< Kronecker blocking for the domain
      Coor<NI> kroni;		   ///< Kronecker blocking for the image
      Tensor<NI, int> ii;	   ///< Number of blocks in each row
      Tensor<NI + 2, int> jj;	   ///< Coordinate of the first element on each block
      Tensor<NI + ND + 1, T> data; ///< Nonzero values
      Tensor<NI + ND + 1, T> kron; ///< Nonzero values for the Kronecker values
      std::shared_ptr<superbblas::BSR_handle> handle; ///< suparbblas sparse tensor handle
      value_type scalar;			      ///< Scalar factor of the tensor
      bool isImgFastInBlock;			      ///< whether the BSR blocks are in row-major
      unsigned int nblockd;			      ///< Number of blocked domain dimensions
      unsigned int nblocki;			      ///< Number of blocked image dimensions
      unsigned int nkrond; ///< Number of Kronecker blocked domain dimensions
      unsigned int nkroni; ///< Number of Kronecker blocked image dimensions

      /// Low-level constructor with the Kronecker BSR extension
      SpTensor(Tensor<ND, T> d, Tensor<NI, T> i, Coor<ND> blkd, Coor<NI> blki, Coor<ND> krond,
	       Coor<NI> kroni, Tensor<NI, int> ii, Tensor<NI + 2, int> jj,
	       Tensor<NI + ND + 1, T> data, Tensor<NI + ND + 1, T> kron_data, value_type scalar,
	       bool isImgFastInBlock, unsigned int nblockd, unsigned int nblocki,
	       unsigned int nkrond, unsigned int nkroni)
	: d(d.make_eg()),
	  i(i.make_eg()),
	  blkd(blkd),
	  blki(blki),
	  krond(krond),
	  kroni(kroni),
	  ii(ii),
	  jj(jj),
	  data(data),
	  kron(kron_data),
	  scalar(scalar),
	  isImgFastInBlock(isImgFastInBlock),
	  nblockd(nblockd),
	  nblocki(nblocki),
	  nkrond(nkrond),
	  nkroni(nkroni)
      {
      }

      /// Low-level constructor without the Kronecker BSR extension
      SpTensor(Tensor<ND, T> d, Tensor<NI, T> i, Coor<ND> blkd, Coor<NI> blki, Tensor<NI, int> ii,
	       Tensor<NI + 2, int> jj, Tensor<NI + ND + 1, T> data, value_type scalar,
	       bool isImgFastInBlock, unsigned int nblockd, unsigned int nblocki)
	: SpTensor(d, i, blkd, blki, detail::ones<ND>(), detail::ones<NI>(), ii, jj, data,
		   Tensor<NI + ND + 1, T>(), scalar, isImgFastInBlock, nblockd, nblocki, 0, 0)
      {
      }

      /// Return a string describing the tensor
      /// \param ptr: pointer to the memory allocation
      /// \return: the string representing the tensor

      std::string repr() const
      {
	using namespace detail::repr;
	std::stringstream ss;
	ss << "SpTensor{";
	if (data.data())
	  ss << "data:" << data.data() << ", ";
	std::size_t sizemb = (ii.getLocal().volume() * sizeof(int) + //
			      jj.getLocal().volume() * sizeof(int) + //
			      data.getLocal().volume() * sizeof(value_type)) /
			     1024 / 1024;
	ss << "domain_order: " << d.order << ", domain_dim:" << d.dim << "image_order: " << i.order
	   << ", image_dim:" << i.dim << ", local_storage:" << sizemb << " MiB}";
	return ss.str();
      }

      /// Constructor
      /// \param d: example tensor for the domain
      /// \param i: example tensor for the image
      /// \param nblockd: the first `nblockd` domain labels will be blocked
      /// \param nblocki: the first `nblocki` image labels will blocked
      /// \param num_neighbors: number of nonzeros for each blocked row

      SpTensor(Tensor<ND, T> d, Tensor<NI, T> i, unsigned int nblockd, unsigned int nblocki,
	       unsigned int nkrond, unsigned int nkroni, unsigned int num_neighbors,
	       bool isImgFastInBlock = false)
	: d{d.make_eg()},
	  i{i.make_eg()},
	  scalar{value_type{1}},
	  isImgFastInBlock{isImgFastInBlock},
	  nblockd(nblockd),
	  nblocki(nblocki),
	  nkrond(nkrond),
	  nkroni(nkroni)
      {
	// Check that the examples are on the same device
	if (d.getDev() != i.getDev())
	  throw std::runtime_error("Please give example vectors on the same device");

	// Check that `d` and `i` are not subtensors
	if (this->d.isSubtensor() || this->i.isSubtensor())
	  throw std::runtime_error("unsupported subtensors for domain/image distributions");

	// Check that the domain and image labels are different and do not contain `u` or `~`
	detail::check_order<NI + ND + 2>(i.order + d.order + std::string("u~"));

	// Get the blocking and the Kronecker blocking
	krond = blkd = kvcoors<ND>(d.order, d.kvdim());
	for (unsigned int i = nblockd; i < ND; ++i)
	  blkd[i] = 1;
	kroni = blki = kvcoors<NI>(i.order, i.kvdim());
	for (unsigned int i = nblocki; i < NI; ++i)
	  blki[i] = 1;
	for (unsigned int i = 0; i < ND; ++i)
	  if (i < nblockd || i >= nblockd + nkrond)
	    krond[i] = 1;
	for (unsigned int i = 0; i < NI; ++i)
	  if (i < nblocki || i >= nblocki + nkroni)
	    kroni[i] = 1;

	// Create the tensor containing the number of neighbors for each blocking
	std::map<char, int> nonblki;
	for (unsigned int j = 0; j < NI; ++j)
	  nonblki[i.order[j]] = i.size[j] / blki[j] / kroni[j];
	ii = i.template make_compatible<NI, int>(none, nonblki);
	ii.set(num_neighbors);

	// Create the tensor containing the domain coordinates of the first nonzero in each block
	jj = ii.template make_compatible<NI + 2, int>(std::string("~u") + i.order,
						      {{'~', (int)ND}, {'u', (int)num_neighbors}});

	// Compute the data dimensions as
	//   image_blocked_dims + domain_dims + u + image_nonblocked_dims, for isImgFastInBlock
	//   domain_blocked_dims + image_blocked_dims + u + image_nonblockd_dims otherwise
	std::map<char, int> data_dims;
	for (unsigned int j = 0; j < NI; ++j)
	  data_dims[i.order[j]] = i.size[j] / kroni[j];
	for (unsigned int i = 0; i < ND; ++i)
	  data_dims[d.order[i]] = blkd[i];
	data_dims['u'] = num_neighbors;
	std::string data_order =
	  (isImgFastInBlock
	     ? std::string(i.order.begin(), i.order.begin() + nblocki + nkroni) + d.order
	     : std::string(d.order.begin(), d.order.begin() + nblockd + nkrond) +
		 std::string(i.order.begin(), i.order.begin() + nblocki + nkroni) +
		 std::string(d.order.begin() + nblockd + nkrond, d.order.end())) +
	  std::string("u") + std::string(i.order.begin() + nblocki + nkroni, i.order.end());
	data = ii.template make_compatible<NI + ND + 1, T>(data_order, data_dims);

	// Compute the Kronecker dimensions as `data`
	if (nkrond + nkroni > 0)
	{
	  std::map<char, int> kron_dims;
	  for (unsigned int i = 0; i < ND; ++i)
	    kron_dims[d.order[i]] = krond[i];
	  for (unsigned int j = 0; j < NI; ++j)
	    kron_dims[i.order[j]] = kroni[j];
	  kron_dims['u'] = num_neighbors;
	  std::string kron_order = data_order;
	  kron = data.like_this(kron_order, kron_dims, none, OnEveryoneReplicated);
	}

	std::string nonblock_img_labels(i.order.begin() + nblocki + nkroni, i.order.end());
	if (!ii.isDistributedAs(this->i, nonblock_img_labels) ||
	    !ii.isDistributedAs(jj, nonblock_img_labels) ||
	    !ii.isDistributedAs(data, nonblock_img_labels))
	  throw std::runtime_error("SpTensor: the dense tensors representing the sparse tensor "
				   "have incompatible distributions");
      }

      /// Empty constructor

      SpTensor()
	: blki{{}},
	  blkd{{}},
	  kroni{{}},
	  krond{{}},
	  scalar{0},
	  isImgFastInBlock{false},
	  nblockd{0},
	  nblocki{0},
	  nkrond{0},
	  nkroni{0}
      {
      }

      /// Return whether the tensor is not empty

      explicit operator bool() const noexcept
      {
	return (bool)d;
      }

      /// Return whether the sparse tensor has Kronecker form

      bool is_kronecker() const noexcept
      {
	return (bool)kron;
      }

      /// Construct the sparse operator
      void construct()
      {
	if ((ii.dist != OnEveryone && ii.dist != OnEveryoneAsChroma && ii.dist != Local) ||
	    ii.dist != jj.dist || ii.dist != data.dist ||
	    (kron && kron.dist != OnEveryoneReplicated))
	  throw std::runtime_error("SpTensor::construct: unexpected distribution of the data");

	// Superbblas needs the column coordinates to be local
	// Remove the local domain coordinates to jj
	const auto localFrom = d.p->localFrom();
	const auto domDim = d.dim;
	auto localjj =
	  jj.template transformWithCPUFunWithCoor<int>([&](const Coor<NI + 2>& c, const int& t) {
	    return (t - localFrom[c[0]] + domDim[c[0]]) % domDim[c[0]];
	  });

	std::string nonblock_img_labels(i.order.begin() + nblocki + nkroni, i.order.end());
	if (!ii.isDistributedAs(this->i, nonblock_img_labels) ||
	    !ii.isDistributedAs(localjj, nonblock_img_labels) ||
	    !ii.isDistributedAs(data, nonblock_img_labels))
	  throw std::runtime_error("SpTensor: the dense tensors representing the sparse tensor "
				   "have incompatible distributions");
	int* iiptr = ii.data();
	Coor<ND>* jjptr = (Coor<ND>*)localjj.data();
	// NOTE: despite jj being a vector of `int`, superbblas will use jj as a vector of Coor<NI>, so check that the alignment
	if (localjj.getLocal().volume() > 0 &&
	    superbblas::detail::align(alignof(Coor<NI>), sizeof(int), jjptr, sizeof(int)) ==
	      nullptr)
	  throw std::runtime_error("Ups! Look into this");
	const value_type* ptr = data.data();
	const value_type* kron_ptr = kron.data();
	MPI_Comm comm = (ii.dist == Local ? MPI_COMM_SELF : MPI_COMM_WORLD);
	superbblas::BSR_handle* bsr = nullptr;
	if (nkrond == 0 && nkroni == 0)
	{
	  superbblas::create_bsr<ND, NI, value_type>(
	    i.p->p.data(), i.dim, d.p->p.data(), d.dim, 1, blki, blkd, isImgFastInBlock, &iiptr,
	    &jjptr, &ptr, &data.ctx(), comm, superbblas::FastToSlow, &bsr);
	}
	else
	{
	  superbblas::create_kron_bsr<ND, NI, value_type>(
	    i.p->p.data(), i.dim, d.p->p.data(), d.dim, 1, blki, blkd, kroni, krond,
	    isImgFastInBlock, &iiptr, &jjptr, &ptr, &kron_ptr, &data.ctx(), comm,
	    superbblas::FastToSlow, &bsr);
	}
	handle = std::shared_ptr<superbblas::BSR_handle>(
	  bsr, [=](superbblas::BSR_handle* bsr) { destroy_bsr(bsr); });
      }

      /// Return a local support of the tensor

      SpTensor<ND, NI, T> getLocal() const
      {
	// Shortcut for empty and local tensors
	if (!*this || ii.dist == Local)
	  return *this;

	// Create the returning tensor
	SpTensor<ND, NI, T> r{d.getLocal(),
			      i.getLocal(),
			      nblockd,
			      nblocki,
			      nkrond,
			      nkroni,
			      (unsigned int)jj.kvdim().at('u'),
			      isImgFastInBlock};

	r.ii = ii.getLocal();
	const auto localFrom = d.p->localFrom();
	const auto domDim = d.dim;
	r.jj = jj.getLocal().template transformWithCPUFunWithCoor<int>(
	  [&](const Coor<NI + 2>& c, const int& t) {
	    return (t - localFrom[c[0]] + domDim[c[0]]) % domDim[c[0]];
	  });
	r.data = data.getLocal();
	r.kron = kron.getLocal();

	if (is_constructed())
	  r.construct();

	return r;
      }

      /// Return a local support of the tensor

      SpTensor<ND, NI, T> getGlocal() const
      {
	// Shortcut for empty and local tensors
	if (!*this || ii.dist == Local)
	  return *this;

	// Create the returning tensor
	SpTensor<ND, NI, T> r{d.getGlocal(),
			      i.getGlocal(),
			      nblockd,
			      nblocki,
			      nkrond,
			      nkroni,
			      (unsigned int)jj.kvdim().at('u'),
			      isImgFastInBlock};

	r.ii = ii.getGlocal();
	r.jj = jj.getGlocal();
	r.data = data.getGlocal();
	r.kron = kron.getGlocal();

	if (is_constructed())
	  r.construct();

	return r;
      }

      /// Split a dimension into another dimensions
      /// \param dom_dim_label: dominion dimension to split
      /// \param dom_new_labels: the labels of the new dominion dimensions
      /// \param dom_step: length of the first label in `dom_new_labels`
      /// \param img_dim_label: image dimension to split
      /// \param img_new_labels: the labels of the image new dimensions
      /// \param img_step: length of the first label in `img_new_labels`

      SpTensor<ND + 1, NI + 1, T>
      split_dimension(char dom_dim_label, const std::string& dom_new_labels, Index dom_step,
		      char img_dim_label, const std::string& img_new_labels, Index img_step) const
      {
	if (dom_new_labels.size() != 2)
	  throw std::runtime_error(
	    "split_dimension: invalid `dom_new_labels`, it should have size two");
	if (img_new_labels.size() != 2)
	  throw std::runtime_error(
	    "split_dimension: invalid `dom_new_labels`, it should have size two");
	if (d.kvdim().at(dom_dim_label) % dom_step != 0)
	  throw std::runtime_error(
	    "split_dimension: invalid `dom_step`, it should divide the dimension size");
	if (i.kvdim().at(img_dim_label) % img_step != 0)
	  throw std::runtime_error(
	    "split_dimension: invalid `img_step`, it should divide the dimension size");

	std::string::size_type d_pos = d.order.find(dom_dim_label);
	std::string::size_type i_pos = i.order.find(img_dim_label);

	if (blkd[d_pos] > 1 && (blkd[d_pos] % dom_step != 0 || krond[d_pos] % dom_step != 0))
	  throw std::runtime_error(
	    "split_dimension: invalid `dom_step`, it should divide the block size");
	if (blki[i_pos] > 1 && (blki[i_pos] % img_step != 0 || kroni[i_pos] % img_step != 0))
	  throw std::runtime_error(
	    "split_dimension: invalid `img_step`, it should divide the block size");

	// Transform the distribution of the domain and the image spaces
	// NOTE: blocking does not operate well with range intersection in the sense that
	//       intersection(range_a, range_b) != unblock(intersection(block(range_a), block(range_b))).
	//       A way to guarantee that is by enforcing that the first coordinate and the size of all ranges are
	//       divisible by the blocking. That's enforced by `coarse_support`.
	// FIXME: enforce that all contracted dense tensors with this sparse tensor have divisible partitions by
	//       the blocking.

	auto new_d = d.coarse_support({{dom_dim_label, dom_step}})
		       .split_dimension(dom_dim_label, dom_new_labels, dom_step);
	auto new_i = i.split_dimension(img_dim_label, img_new_labels, img_step);

	int new_blkd_pos = blkd[d_pos] == 1 ? 1 : dom_step;
	auto new_blkd = detail::insert_coor(blkd, d_pos, new_blkd_pos);
	new_blkd[d_pos + 1] /= new_blkd_pos;
	int new_blki_pos = blki[i_pos] == 1 ? 1 : img_step;
	auto new_blki = detail::insert_coor(blki, i_pos, new_blki_pos);
	new_blki[i_pos + 1] /= new_blki_pos;

	// Create the returning tensor
	SpTensor<ND + 1, NI + 1, T> r{
	  new_d,
	  new_i,
	  nblockd + (d_pos < nblockd ? 1 : 0),
	  nblocki + (i_pos < nblocki ? 1 : 0),
	  nkrond + (nblockd <= d_pos && d_pos < nblockd + nkrond ? 1 : 0),
	  nkroni + (nblocki <= i_pos && i_pos < nblocki + nkroni ? 1 : 0),
	  (unsigned int)jj.kvdim().at('u'),
	  isImgFastInBlock};

	ii.split_dimension(img_dim_label, img_new_labels, img_step).copyTo(r.ii);
	auto new_jj = r.jj.make_compatible(none, {}, OnHost);
	jj.split_dimension(img_dim_label, img_new_labels, img_step)
	  .copyTo(new_jj.kvslice_from_size({}, {{'~', (int)ND}}));
	{
	  auto local_new_jj = new_jj.getLocal();
	  int* p = local_new_jj.data();
	  auto new_dom_dim = detail::insert_coor(d.size, d_pos, dom_step);
	  new_dom_dim[d_pos + 1] /= dom_step;
	  std::size_t i1 = local_new_jj.volume() / (ND + 1);
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < i1; ++i)
	  {
	    Coor<ND> c;
	    std::copy_n(p + (ND + 1) * i, ND, c.begin());
	    Coor<ND + 1> new_c = detail::split_dimension(d_pos, c, new_dom_dim, detail::From);
	    std::copy_n(new_c.begin(), ND + 1, p + (ND + 1) * i);
	  }
	}
	new_jj.copyTo(r.jj);

	data.split_dimension(dom_dim_label, dom_new_labels, dom_step)
	  .split_dimension(img_dim_label, img_new_labels, img_step)
	  .copyTo(r.data);

	if (kron)
	{
	  kron.split_dimension(dom_dim_label, dom_new_labels, dom_step)
	    .split_dimension(img_dim_label, img_new_labels, img_step)
	    .copyTo(r.kron);
	}

	if (is_constructed())
	  r.construct();

	return r;
      }

      /// Return a slice of the tensor starting at coordinate `dom_kvfrom`, `img_kvfrom` and taking
      /// `dom_kvsize`, `img_kvsize` elements in each direction. The missing dimensions in `*_kvfrom`
      /// are set to zero and the missing directions in `*_kvsize` are set to the size of the tensor.
      ///
      /// \param dom_kvfrom: dictionary with the index of the first element in each domain direction
      /// \param dom_kvsize: dictionary with the number of elements in each domain direction
      /// \param img_kvfrom: dictionary with the index of the first element in each domain direction
      /// \param img_kvsize: dictionary with the number of elements in each domain direction
      /// \return: a copy of the tensor

      SpTensor<ND, NI, T> kvslice_from_size(const std::map<char, int>& dom_kvfrom = {},
					    const std::map<char, int>& dom_kvsize = {},
					    const std::map<char, int>& img_kvfrom = {},
					    const std::map<char, int>& img_kvsize = {}) const
      {
	// Check that we aren't slicing the blocking dimensions
	bool fail = false;
	std::string o_blk_d = std::string(d.order.begin(), d.order.begin() + nblockd);
	for (auto& it : dom_kvfrom)
	  if (detail::is_in(o_blk_d, it.first) && it.second != 0)
	    fail = true;
	auto dim_d = d.kvdim();
	for (auto& it : dom_kvsize)
	  if (detail::is_in(o_blk_d, it.first) && it.second != dim_d.at(it.first))
	    fail = true;

	std::string o_blk_i = std::string(i.order.begin(), i.order.begin() + nblocki);
	for (auto& it : img_kvfrom)
	  if (detail::is_in(o_blk_i, it.first) && it.second != 0)
	    fail = true;
	auto dim_i = i.kvdim();
	for (auto& it : img_kvsize)
	  if (detail::is_in(o_blk_i, it.first) && it.second != dim_i.at(it.first))
	    fail = true;

	if (fail)
	  throw std::runtime_error(
	    "SpTensor::kvslice_from_size: unsupported slicing on blocked dimensions");

	// We aren't free to redistribute `d` and `i`, because the support of the domain
	// in each process depends on the image support
	auto new_d = d.kvslice_from_size(dom_kvfrom, dom_kvsize).make_eg();
	auto new_i = i.kvslice_from_size(img_kvfrom, img_kvsize).make_eg();

	// Get the nonzeros in the slice
	auto ii_slice = ii.kvslice_from_size(img_kvfrom, img_kvsize).cloneOn(OnHost);
	auto new_ii = ii_slice.make_compatible(none, {}, OnHost);
	new_ii.set_zero();
	auto jj_slice = jj.kvslice_from_size(img_kvfrom, img_kvsize).cloneOn(OnHost);
	auto new_jj = jj_slice.make_compatible(none, {}, OnHost);
	auto new_jj_mask = jj_slice.template make_compatible<NI + 1, float>(
	  detail::remove_dimensions(jj_slice.order, "~"), {}, OnHost);
	new_jj_mask.set_zero();
	unsigned int num_neighbors = jj.kvdim().at('u');

	if (ii_slice.isSubtensor() || new_ii.isSubtensor() || jj_slice.isSubtensor() ||
	    new_jj.isSubtensor() || new_jj_mask.isSubtensor())
	{
	  throw std::runtime_error("This shouldn't happen");
	}
	if (!ii_slice.is_compatible(jj_slice) || !ii_slice.is_compatible(new_ii) ||
	    !ii_slice.is_compatible(new_jj) || !ii_slice.is_compatible(new_jj_mask))
	{
	  throw std::runtime_error("kvslice_from_size: hit corner case, sorry");
	}

	Tensor<1, float> dirs("u", {(int)num_neighbors}, OnHost, OnMaster);
	dirs.set_zero();
	auto dirs_local = dirs.getLocal();
	{
	  Coor<ND> from_dom = kvcoors<ND>(d.order, dom_kvfrom);
	  std::map<char, int> updated_dom_kvsize = d.kvdim();
	  for (const auto& it : dom_kvsize)
	    updated_dom_kvsize[it.first] = it.second;
	  Coor<ND> size_dom = kvcoors<ND>(d.order, updated_dom_kvsize);
	  auto ii_slice_local = ii_slice.getLocal();
	  int* ii_slice_ptr = ii_slice_local.data();
	  auto jj_slice_local = jj_slice.getLocal();
	  int* jj_slice_ptr = jj_slice_local.data();
	  auto new_ii_local = new_ii.getLocal();
	  int* new_ii_ptr = new_ii_local.data();
	  auto new_jj_local = new_jj.getLocal();
	  int* new_jj_ptr = new_jj_local.data();
	  auto new_jj_mask_local = new_jj_mask.getLocal();
	  float* new_jj_mask_ptr = new_jj_mask_local.data();
	  Coor<ND> size_nnz = d.size;
	  Tensor<1, float> dirs_global("u", {(int)num_neighbors}, OnHost, OnMaster);
	  dirs_global.set_zero();
	  auto dirs_local = dirs_global.getLocal();
	  for (unsigned int i = nblockd + nkrond; i < ND; ++i)
	    size_nnz[i] = 1;
	  for (std::size_t i = 0, i_acc = 0, i1 = new_ii_local.volume(); i < i1;
	       i_acc += ii_slice_ptr[i], ++i)
	  {
	    for (unsigned int j = i_acc, j1 = i_acc + ii_slice_ptr[i]; j < j1; ++j)
	    {
	      Coor<ND> from_nnz;
	      std::copy_n(jj_slice_ptr + j * ND, ND, from_nnz.begin());
	      Coor<ND> lfrom, lsize;
	      superbblas::detail::intersection(from_dom, size_dom, from_nnz, size_nnz, d.dim, lfrom,
					       lsize);
	      if (superbblas::detail::volume(lsize) == 0)
		continue;

	      using superbblas::detail::operator-;
	      Coor<ND> new_from_nnz = normalize_coor(from_nnz - from_dom, size_dom);
	      std::copy_n(new_from_nnz.begin(), ND, new_jj_ptr + (i_acc + new_ii_ptr[i]) * ND);
	      new_ii_ptr[i]++;
	      new_jj_mask_ptr[j] = 1;
	      if (dirs_local)
		dirs_local.data()[j - i_acc] = 1;
	    }
	    if (i > 0 && new_ii_ptr[i] != new_ii_ptr[0])
	      throw std::runtime_error("SpTensor::kvslice_from_size: unsupported slices ending up "
				       "in different number of nonzero values in each row");
	    if (i > 0 && kron)
	    {
	      for (unsigned int j = i_acc, j1 = i_acc + ii_slice_ptr[i]; j < j1; ++j)
		if (new_jj_mask_ptr[j] != new_jj_mask_ptr[j - i_acc])
		  throw std::runtime_error("SpTensor::kvslice_from_size unsupported slices ending "
					   "up in selecting different directions for each row");
	    }
	  }

	  // Make sure that all nodes with support have the same number of neighbors
	  if (Layout::nodeNumber() == 0 && new_ii_local.volume() == 0)
	    throw std::runtime_error("kvslice_from_size: unsupported distribution, master process "
				     "should have support on the origin tensor");
	  if (new_ii_local.volume() > 0)
	    num_neighbors = new_ii_ptr[0];
	  int global_num_neighbors = broadcast(num_neighbors);
	  if (new_ii_local.volume() > 0 && global_num_neighbors != num_neighbors)
	    throw std::runtime_error("SpTensor::kvslice_from_size: unsupported distribution");
	  num_neighbors = global_num_neighbors;

	  dirs = dirs_global.make_sure(none, OnDefaultDevice, OnEveryoneReplicated);
	}

	// Create the returning tensor
	SpTensor<ND, NI, T> r{new_d,  new_i,  nblockd,	     nblocki,
			      nkrond, nkroni, num_neighbors, isImgFastInBlock};
	new_ii.copyTo(r.ii);
	new_jj.kvslice_from_size({}, {{'u', num_neighbors}}).copyTo(r.jj);

	auto data_mask = data.create_mask();
	data_mask.set_zero();
	std::map<char, int> blk_m;
	auto data_dim = data.kvdim();
	for (unsigned int i = 0; i < ND; ++i)
	  blk_m[d.order[i]] = data_dim.at(d.order[i]);
	for (unsigned int i = 0; i < NI; ++i)
	  blk_m[this->i.order[i]] = (i < nblocki ? data_dim.at(this->i.order[i]) : 1);
	auto data_blk =
	  data.template like_this<NI + ND, float>("%", '%', "u", blk_m, none, OnEveryoneReplicated);
	data_blk.set(1);
	kronecker<NI + ND + 1>(new_jj_mask, data_blk)
	  .copyTo(data_mask.kvslice_from_size(img_kvfrom, img_kvsize));
	auto r_data_mask = r.data.create_mask();
	r_data_mask.set(1);
	r.data.set(detail::NaN<T>::get());
	data.kvslice_from_size(img_kvfrom, img_kvsize)
	  .copyToWithMask(r.data, data_mask.kvslice_from_size(img_kvfrom, img_kvsize), r_data_mask,
			  "u");

	if (kron)
	{
	  auto kron_mask = kron.create_mask();
	  kron_mask.set_zero();
	  auto blk_m = kron.kvdim();
	  blk_m.erase('u');
	  auto kron_blk = kron.template like_this<NI + ND, float>("%", '%', "u", blk_m, none,
								  OnEveryoneReplicated);
	  kron_blk.set(1);
	  kronecker<NI + ND + 1>(dirs, kron_blk)
	    .copyTo(kron_mask.kvslice_from_size(img_kvfrom, img_kvsize));
	  auto r_kron_mask = r.kron.create_mask();
	  r_kron_mask.set(1);
	  r.kron.set(detail::NaN<T>::get());
	  kron.kvslice_from_size(img_kvfrom, img_kvsize)
	    .copyToWithMask(r.kron, kron_mask.kvslice_from_size(img_kvfrom, img_kvsize),
			    r_kron_mask, "u");
	}

	if (is_constructed())
	  r.construct();

	// Do a test
	if (superbblas::getDebugLevel() > 0)
	{
	  auto x0 = d.template like_this<ND + 1>("%n", '%', "", {{'n', 2}});
	  x0.set_zero();
	  urand(x0.kvslice_from_size(dom_kvfrom, dom_kvsize), -1, 1);
	  auto y0 = i.template like_this<NI + 1>("%n", '%', "", {{'n', 2}});
	  contractWith(x0, {}, y0, {});
	  y0 = y0.kvslice_from_size(img_kvfrom, img_kvsize);

	  auto y = r.i.template like_this<NI + 1>("%n", '%', "", {{'n', 2}});
	  if (!is_constructed())
	    r.construct();
	  r.contractWith(x0.kvslice_from_size(dom_kvfrom, dom_kvsize), {}, y, {});

	  y0.scale(-1).addTo(y);
	  auto norm0 = norm<1>(y0, "n");
	  auto normdiff = norm<1>(y, "n");
	  double max_err = 0;
	  for (int i = 0, vol = normdiff.volume(); i < vol; ++i)
	    max_err = std::max(max_err, (double)normdiff.get({{i}}) / norm0.get({{i}}));
	  QDPIO::cout << "kvslice_from_size error: " << detail::tostr(max_err) << std::endl;
	}

	return r;
      }

      /// Reorder the domain and image orders
      /// \param new_dom_order: new ordering for the domain
      /// \param new_img_order: new ordering for the image
      /// \param remaining_char: if it isn't the null char, placeholder for the remaining dimensions

      SpTensor<ND, NI, T> reorder(const std::string& new_dom_order,
				  const std::string& new_img_order, char remaining_char = 0) const
      {
	if (remaining_char == '~' || remaining_char == 'u')
	  throw std::runtime_error("reorder: invalid remaining char, it shouldn't be `~` or `u`");

	std::string new_dom_order0 = d.get_order_for_reorder(new_dom_order, remaining_char);
	std::string new_img_order0 = i.get_order_for_reorder(new_img_order, remaining_char);
	auto new_d = d.reorder(new_dom_order0);
	auto new_i = i.reorder(new_img_order0);

	Coor<ND> d_perm = superbblas::detail::find_permutation(
	  detail::to_sb_order<ND>(d.order), detail::to_sb_order<ND>(new_dom_order0));
	Coor<NI> i_perm = superbblas::detail::find_permutation(
	  detail::to_sb_order<NI>(i.order), detail::to_sb_order<NI>(new_img_order0));
	auto new_blkd = superbblas::detail::reorder_coor(blkd, d_perm);
	auto new_blki = superbblas::detail::reorder_coor(blki, i_perm);
	auto new_krond = superbblas::detail::reorder_coor(krond, d_perm);
	auto new_kroni = superbblas::detail::reorder_coor(kroni, i_perm);

	// Check the blocking
	for (unsigned int i = 0; i < ND; ++i)
	  if ((i < nblockd && new_blkd[i] != new_d.size[i]) || (i >= nblockd && new_blkd[i] != 1))
	    throw std::runtime_error("reorder: invalid domain reordering, it is mixing blocking "
				     "and nonblocking dimensions");
	for (unsigned int i = 0; i < NI; ++i)
	  if ((i < nblocki && new_blki[i] != new_i.size[i]) || (i >= nblocki && new_blki[i] != 1))
	    throw std::runtime_error("reorder: invalid image reordering, it is mixing blocking "
				     "and nonblocking dimensions");
	for (unsigned int i = 0; i < ND; ++i)
	  if ((i >= nblockd && i < nblockd + nkrond && new_krond[i] != new_d.size[i]) ||
	      ((i < nblockd || i >= nblockd + nkrond) && new_krond[i] != 1))
	    throw std::runtime_error("reorder: invalid domain reordering, it is mixing blocking "
				     "and nonblocking dimensions");
	for (unsigned int i = 0; i < NI; ++i)
	  if ((i >= nblocki && i < nblocki + nkroni && new_kroni[i] != new_i.size[i]) ||
	      ((i < nblocki || i >= nblocki + nkroni) && new_kroni[i] != 1))
	    throw std::runtime_error("reorder: invalid image reordering, it is mixing blocking "
				     "and nonblocking dimensions");

	auto new_ii = ii.reorder(new_img_order0);
	auto new_jj = jj.reorder(std::string("~u") + new_img_order0);
	if (new_jj.order != jj.order)
	{
	  auto host_new_jj = new_jj.make_sure(none, OnHost);
	  auto local_new_jj = host_new_jj.getLocal();
	  int* new_p = local_new_jj.data();
	  for (std::size_t i = 0, i1 = local_new_jj.volume(); i < i1; i += ND)
	  {
	    Coor<ND> c;
	    std::copy_n(new_p + i, ND, c.begin());
	    Coor<ND> new_c = superbblas::detail::reorder_coor(c, d_perm);
	    std::copy_n(new_c.begin(), ND, new_p + i);
	  }
	  host_new_jj.copyTo(new_jj);
	}

	std::string data_order =
	  (isImgFastInBlock
	     ? std::string(new_i.order.begin(), new_i.order.begin() + nblocki + nkroni) +
		 new_d.order
	     : std::string(new_d.order.begin(), new_d.order.begin() + nblockd + nkrond) +
		 std::string(new_i.order.begin(), new_i.order.begin() + nblocki + nkroni) +
		 std::string(new_d.order.begin() + nblockd + nkrond, new_d.order.end())) +
	  std::string("u") + std::string(new_i.order.begin() + nblocki + nkrond, new_i.order.end());
	auto new_data = data.reorder(data_order);
	auto new_kron = kron ? kron.reorder(data_order) : kron;

	SpTensor<ND, NI, T> r(new_d, new_i, new_blkd, new_blki, new_krond, new_kroni, new_ii,
			      new_jj, new_data, new_kron, scalar, isImgFastInBlock, nblockd,
			      nblocki, nkrond, nkroni);
	if (is_constructed())
	  r.construct();

	return r;
      }

      /// Return a view of this tensor with an extra label for the real and the imaginary parts

      template <typename U = T>
      typename std::enable_if<!detail::is_diycomplex<U>::value && detail::is_complex<U>::value,
			      SpTensor<ND + 1, NI + 1, DIYComplex<typename U::value_type>>>::type
      toFakeReal() const
      {
	using newT = DIYComplex<typename U::value_type>;

	// Get the new domain and image
	char d_complexLabel = detail::get_free_label(d.order + i.order);
	auto new_d = d.toFakeReal(d_complexLabel);
	char i_complexLabel = detail::get_free_label(new_d.order + i.order);
	auto new_i = i.toFakeReal(i_complexLabel);

	// Create the returning tensor
	bool is_kron = (nkrond > 0 || nkroni > 0);
	SpTensor<ND + 1, NI + 1, newT> r{new_d,
					 new_i,
					 nblockd + 1,
					 nblocki + 1,
					 is_kron ? nkrond + 1 : 0,
					 is_kron ? nkroni + 1 : 0,
					 (unsigned int)jj.kvdim().at('u'),
					 isImgFastInBlock};

	// Copy the data to the new tensor
	// a) same number of nonzeros per row
	ii.copyTo(r.ii);
	// b) the nonzero blocks start at the same position
	r.jj.kvslice_from_size({}, {{'~', 1}}).set(0);
	jj.copyTo(r.jj.kvslice_from_size({{'~', 1}}, {{'~', ND}}));
	// c) each element in data xr+xi*i -> [xr  -xi; xi  xr]
	auto data0 = data.toFakeReal(d_complexLabel);
	data0.kvslice_from_size({}, {{d_complexLabel, 1}})
	  .copyTo(r.data.kvslice_from_size({}, {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	data0.kvslice_from_size({}, {{d_complexLabel, 1}})
	  .copyTo(r.data.kvslice_from_size({{d_complexLabel, 1}, {i_complexLabel, 1}},
					   {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	data0.kvslice_from_size({{d_complexLabel, 1}}, {{d_complexLabel, 1}})
	  .scale(-1)
	  .copyTo(r.data.kvslice_from_size({{d_complexLabel, 1}},
					   {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	data0.kvslice_from_size({{d_complexLabel, 1}}, {{d_complexLabel, 1}})
	  .copyTo(r.data.kvslice_from_size({{i_complexLabel, 1}},
					   {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	// d) each element in kron xr+xi*i -> [xr  -xi; xi  xr]
	if (kron)
	{
	  auto kron0 = kron.toFakeReal(d_complexLabel);
	  kron0.kvslice_from_size({}, {{d_complexLabel, 1}})
	    .copyTo(r.kron.kvslice_from_size({}, {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	  kron0.kvslice_from_size({}, {{d_complexLabel, 1}})
	    .copyTo(r.kron.kvslice_from_size({{d_complexLabel, 1}, {i_complexLabel, 1}},
					     {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	  kron0.kvslice_from_size({{d_complexLabel, 1}}, {{d_complexLabel, 1}})
	    .scale(-1)
	    .copyTo(r.kron.kvslice_from_size({{d_complexLabel, 1}},
					     {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	  kron0.kvslice_from_size({{d_complexLabel, 1}}, {{d_complexLabel, 1}})
	    .copyTo(r.kron.kvslice_from_size({{i_complexLabel, 1}},
					     {{d_complexLabel, 1}, {i_complexLabel, 1}}));
	}

	if (is_constructed())
	  r.construct();

	return r;
      }

      template <typename U = T>
      typename std::enable_if<detail::is_diycomplex<U>::value || !detail::is_complex<U>::value,
			      SpTensor<ND, NI, U>>::type
      toFakeReal() const
      {
	return *this;
      }

      /// Extend the support of each dimension by the given amount in each direction
      /// \param m: amount to extend the support for each process

      SpTensor<ND, NI, T> extend_support(const std::map<char, int>& m) const
      {
	std::map<char, int> md, mi;
	for (const auto& it : m)
	{
	  if (detail::is_in(d.order, it.first))
	    md[it.first] = it.second;
	  else if (detail::is_in(i.order, it.first))
	    mi[it.first] = it.second;
	  else
	    throw std::runtime_error("extend_support: unmatched label");
	}
	auto new_d = d.extend_support(md);
	auto new_i = i.extend_support(mi);

	// Create the returning tensor
	SpTensor<ND, NI, T> r{new_d,
			      new_i,
			      nblockd,
			      nblocki,
			      nkrond,
			      nkroni,
			      (unsigned int)jj.kvdim().at('u'),
			      isImgFastInBlock};

	// Populate the new tensor
	ii.copyTo(r.ii);
	jj.copyTo(r.jj);
	data.copyTo(r.data);
	if (kron)
	  kron.copyTo(r.kron);

	if (is_constructed())
	  r.construct();

	return r;
      }

      /// Return whether the sparse tensor has been constructed

      bool is_constructed() const
      {
	return (bool)handle;
      }

      /// Get where the tensor is stored

      DeviceHost getDev() const
      {
	return i.getDev();
      }

      // Contract the dimensions with the same label in this tensor and in `v` than do not appear on `w`.
      template <std::size_t Nv, std::size_t Nw, typename Tv, typename Tw,
		typename std::enable_if<
		  (!std::is_same<Tv, T>::value || !std::is_same<Tw, T>::value), bool>::type = true>
      void contractWith(Tensor<Nv, Tv> v, const remap& mv, const Tensor<Nw, Tw>& w,
			const remap& mw = {}, char power_label = 0) const
      {
	if (data.is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	if (!is_constructed())
	  throw std::runtime_error("invalid operation on an not constructed tensor");

	auto w0 = w.template cast_like<T>();
	contractWith(std::move(v).template cast<T>(), mv, w0);
	w0.copyTo(w);
      }

      template <std::size_t Nv, std::size_t Nw>
      void contractWith(Tensor<Nv, T> v, const remap& mv, Tensor<Nw, T> w, const remap& mw = {},
			char power_label = 0) const
      {
	if (data.is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	if (!is_constructed())
	  throw std::runtime_error("invalid operation on an not constructed tensor");

	// If either this tensor or v are on OnDevice, force both to be on the same device as this tensor.
	if (v.ctx().plat != data.ctx().plat)
	{
	  v = v.cloneOn(getDev());
	}

	if (getDev() != w.getDev())
	{
	  Tensor<Nw, T> aux = w.like_this(none, {}, getDev());
	  contractWith(v, mv, aux, mw, power_label);
	  aux.copyTo(w);
	  return;
	}

	// Check unsupported distributions for contraction
	if ((v.dist == Local) != (w.dist == Local) || (v.dist == Glocal) != (w.dist == Glocal) ||
	    (v.dist != Local && data.dist == Local))
	  throw std::runtime_error("contractWith: One of the contracted tensors or the output "
				   "tensor is local and others are not!");

	// We don't support conjugacy for now
	if (v.conjugate || w.conjugate)
	  throw std::runtime_error("contractWith: unsupported implicit conjugacy");

	// Check the power label
	if (power_label != 0 && detail::is_in(v.order, power_label) &&
	    v.kvdim().at(power_label) > 1)
	  throw std::runtime_error("contractWith: `power_label` for `v` does not have size one");
	if (power_label != 0 && !detail::is_in(w.order, power_label))
	  throw std::runtime_error("contractWith: `power_label` isn't in `w`");

	value_type* v_ptr = v.data();
	value_type* w_ptr = w.data_for_writing();
	std::string orderv = detail::update_order_and_check<Nv>(v.order, mv);
	std::string orderw = detail::update_order_and_check<Nw>(w.order, mw);
	superbblas::bsr_krylov<ND, NI, Nv, Nw, value_type>(
	  scalar * v.scalar / w.scalar, handle.get(), i.order.c_str(), d.order.c_str(),	       //
	  v.p->p.data(), 1, orderv.c_str(), v.from, v.size, v.dim, (const value_type**)&v_ptr, //
	  T{0}, w.p->p.data(), orderw.c_str(), w.from, w.size, w.dim, power_label,
	  (value_type**)&w_ptr, //
	  &data.ctx(), v.dist == Local ? MPI_COMM_SELF : MPI_COMM_WORLD, superbblas::FastToSlow,
	  nullptr, v.dist == Glocal);

	// Force synchronization in superbblas stream if the destination allocation isn't managed by superbblas
	if (!w.is_managed())
	  superbblas::sync(data.ctx());
      }

      void print(const std::string& name) const
      {
	std::stringstream ss;

	auto ii_host = ii.make_sure(none, OnHost, OnMaster).getLocal();
	auto jj_host = jj.make_sure(none, OnHost, OnMaster).getLocal();
	auto data_host = data.make_sure(none, OnHost, OnMaster).getLocal();
	auto kron_host = kron.make_sure(none, OnHost, OnMaster).getLocal();
	assert(!ii_host.isSubtensor() && !jj_host.isSubtensor() && !data_host.isSubtensor() &&
	       !kron_host.isSubtensor());

	std::size_t volblki = superbblas::detail::volume(blki);
	std::size_t volblkj = superbblas::detail::volume(blkd);
	std::size_t volkroni = superbblas::detail::volume(kroni);
	std::size_t volkronj = superbblas::detail::volume(krond);
	std::size_t volbi = volblki * volkroni;
	std::size_t volbj = volblkj * volkronj;

	auto ii_host_ptr = ii_host.data();
	auto jj_host_ptr = jj_host.data();
	auto data_host_ptr = data_host.data();
	auto kron_host_ptr = kron_host.data();

	// Print general tensor description in a matlab comment
	if (ii_host)
	  ss << "% " << repr() << std::endl;

	// If using the Kronecker format, reconstruct the nonzeros explicitly
	if (kron_host)
	{
	  // Print the Kronecker tensor (the spin-spin tensors)
	  int num_neighbors = data.kvdim().at('u');
	  ss << name << "_kron=reshape([";
	  for (std::size_t i = 0, kron_vol = kron_host.volume(); i < kron_vol; ++i)
	    detail::repr::operator<<(ss << " ", kron_host_ptr[i]);
	  ss << "], [" << (isImgFastInBlock ? volkroni : volkronj) << " "
	     << (isImgFastInBlock ? volkronj : volkroni) << " " << num_neighbors << "]);"
	     << std::endl;

	  // Print the data (the color-color tensors)
	  std::size_t numblks = data.volume() / volblki / volblkj;
	  ss << name << "_data0=reshape([";
	  for (std::size_t i = 0, data_vol = data_host.volume(); i < data_vol; ++i)
	    detail::repr::operator<<(ss << " ", data_host_ptr[i]);
	  ss << "], [" << (isImgFastInBlock ? volblki : volblkj) << " "
	     << (isImgFastInBlock ? volblkj : volblki) << " " << numblks << "]);" << std::endl;

	  // Preallocate and populate data in explicit format (no Kronecker format)
	  ss << name << "_data=zeros([" << (isImgFastInBlock ? volbi : volbj) << " "
	     << (isImgFastInBlock ? volbj : volbi) << " " << numblks << "]);" << std::endl;
	  ss << "for i=1:" << numblks << std::endl;
	  ss << "  " << name << "_data(:,:,i)=kron(squeeze(" << name << "_kron(:,:,mod(i-1,"
	     << num_neighbors << ")+1)), squeeze(" << name << "_data0(:,:,i)));" << std::endl;
	  ss << "end" << std::endl;
	}

	// Only master node prints
	if (ii_host)
	{
	  // Print for non-Kronecker variant
	  ss << name << "=sparse([";

	  // Print the row indices
	  for (std::size_t i = 0, iivol = ii_host.volume(); i < iivol; ++i)
	  {
	    for (unsigned int neighbor = 0, num_neighbors = ii_host_ptr[i];
		 neighbor < num_neighbors; ++neighbor)
	    {
	      if (isImgFastInBlock)
	      {
		for (unsigned int bj = 0; bj < volbj; ++bj)
		  for (unsigned int bi = 0; bi < volbi; ++bi)
		    ss << " " << i * volbi + bi + 1;
	      }
	      else
	      {
		for (unsigned int bi = 0; bi < volbi; ++bi)
		  for (unsigned int bj = 0; bj < volbj; ++bj)
		    ss << " " << i * volbi + bi + 1;
	      }
	    }
	  }
	  ss << "], [";

	  // Print the column indices
	  Stride<ND> dstrides =
	    superbblas::detail::get_strides<std::size_t>(d.size, superbblas::FastToSlow);
	  for (std::size_t j = 0, jjvol = jj_host.volume(); j < jjvol; j += ND)
	  {
	    Coor<ND> j_coor;
	    std::copy_n(jj_host_ptr + j, ND, j_coor.begin());
	    auto j_idx = superbblas::detail::coor2index(j_coor, d.size, dstrides);
	    if (isImgFastInBlock)
	    {
	      for (unsigned int bj = 0; bj < volbj; ++bj)
		for (unsigned int bi = 0; bi < volbi; ++bi)
		  ss << " " << j_idx + bj + 1;
	    }
	    else
	    {
	      for (unsigned int bi = 0; bi < volbi; ++bi)
		for (unsigned int bj = 0; bj < volbj; ++bj)
		  ss << " " << j_idx + bj + 1;
	    }
	  }
	  ss << "], ";

	  if (!kron)
	  {
	    ss << "[";

	    // Print the data values
	    for (std::size_t i = 0, data_vol = data_host.volume(); i < data_vol; ++i)
	      detail::repr::operator<<(ss << " ", data_host_ptr[i]);
	    ss << "]";
	  }
	  else
	  {
	    ss << name << "_data(:)";
	  }
	  ss << ");" << std::endl;
	}

	detail::log(1, ss.str());
      }

      /// Return a copy of the tensor in a different precision
      ///
      /// \tparam Q: new precision

      template <typename Q = T,
		typename std::enable_if<std::is_same<T, Q>::value, bool>::type = true>
      SpTensor<ND, NI, Q> cast() const
      {
	return *this;
      }

      template <typename Q = T,
		typename std::enable_if<!std::is_same<T, Q>::value, bool>::type = true>
      SpTensor<ND, NI, Q> cast() const
      {
	SpTensor<ND, NI, Q> r{d.template cast<Q>(),
			      i.template cast<Q>(),
			      blkd,
			      blki,
			      krond,
			      kroni,
			      ii,
			      jj,
			      data.template cast<Q>(),
			      kron.template cast<Q>(),
			      (Q)scalar,
			      isImgFastInBlock,
			      nblockd,
			      nblocki,
			      nkrond,
			      nkroni};
	if (is_constructed())
	  r.construct();
	return r;
      }
    };

    template <std::size_t N, typename T>
    struct StorageTensor {
      static_assert(superbblas::supported_type<T>::value, "Not supported type");

    public:
      std::string filename; ///< Storage file
      std::string metadata; ///< metadata
      std::string order;    ///< Labels of the tensor dimensions
      Coor<N> dim;	    ///< Length of the tensor dimensions
      Sparsity sparsity;    ///< Sparsity of the storage
      std::shared_ptr<superbblas::detail::Storage_context_abstract>
	ctx;	    ///< Superbblas storage handler
      Coor<N> from; ///< First active coordinate in the tensor
      Coor<N> size; ///< Number of active coordinates on each dimension
      T scalar;	    ///< Scalar factor of the tensor

      // Empty constructor
      StorageTensor()
	: filename{},
	  metadata{},
	  order(detail::getTrivialOrder(N)),
	  dim{{}},
	  sparsity(Dense),
	  ctx{},
	  from{{}},
	  size{{}},
	  scalar{0}
      {
      }

      // Create storage construct
      StorageTensor(const std::string& filename, const std::string& metadata,
		    const std::string& order, Coor<N> dim, Sparsity sparsity = Dense,
		    checksum_type checksum = checksum_type::NoChecksum)
	: filename(filename),
	  metadata(metadata),
	  order(order),
	  dim(dim),
	  sparsity(sparsity),
	  from{{}},
	  size{dim},
	  scalar{1}
      {
	checkOrder();
	superbblas::Storage_handle stoh;
	superbblas::create_storage<N, T>(dim, superbblas::FastToSlow, filename.c_str(),
					 metadata.c_str(), metadata.size(), checksum,
					 MPI_COMM_WORLD, &stoh);
	ctx = std::shared_ptr<superbblas::detail::Storage_context_abstract>(
	  stoh, [=](superbblas::detail::Storage_context_abstract* ptr) {
	    superbblas::close_storage<N, T>(ptr, MPI_COMM_WORLD);
	  });

	// If the tensor to store is dense, create the block here; otherwise, create the block on copy
	if (sparsity == Dense)
	{
	  superbblas::PartitionItem<N> p{Coor<N>{{}}, dim};
	  superbblas::append_blocks<N, T>(&p, 1, dim, stoh, MPI_COMM_WORLD, superbblas::FastToSlow);
	}
      }

      // Open storage construct
      StorageTensor(const std::string& filename, bool read_order = true,
		    const Maybe<std::string>& order_tag = none)
	: filename(filename), sparsity(Sparse), from{{}}, scalar{1}
      {
	// Read information from the storage
	superbblas::values_datatype values_dtype;
	std::vector<char> metadatav;
	std::vector<superbblas::IndexType> dimv;
	superbblas::read_storage_header(filename.c_str(), superbblas::FastToSlow, values_dtype,
					metadatav, dimv, MPI_COMM_WORLD);

	// Check that storage tensor dimension and value type match template arguments
	if (dimv.size() != N)
	  throw std::runtime_error(
	    "The storage tensor dimension does not match the template parameter N");
	if (superbblas::detail::get_values_datatype<T>() != values_dtype)
	  throw std::runtime_error("Storage type does not match template argument T");

	// Fill out the information of this class with storage header information
	std::copy(dimv.begin(), dimv.end(), dim.begin());
	size = dim;
	metadata = std::string(metadatav.begin(), metadatav.end());

	// Read the order
	if (read_order)
	{
	  std::istringstream is(metadata);
	  XMLReader xml_buf(is);
	  read(xml_buf, order_tag.getSome("order"), order);
	  checkOrder();
	}

	superbblas::Storage_handle stoh;
	superbblas::open_storage<N, T>(filename.c_str(), false /* don't allow writing */,
				       MPI_COMM_WORLD, &stoh);
	ctx = std::shared_ptr<superbblas::detail::Storage_context_abstract>(
	  stoh, [=](superbblas::detail::Storage_context_abstract* ptr) {
	    superbblas::close_storage<N, T>(ptr, MPI_COMM_WORLD);
	  });
      }

    protected:
      // Construct a slice/scale storage
      StorageTensor(const StorageTensor& t, const std::string& order, Coor<N> from, Coor<N> size,
		    T scalar)
	: filename(t.filename),
	  metadata(t.metadata),
	  order(order),
	  dim(t.dim),
	  ctx(t.ctx),
	  sparsity(t.sparsity),
	  from(normalize_coor(from, t.dim)),
	  size(size),
	  scalar{t.scalar}
      {
	checkOrder();
      }

    public:
      /// Return whether the tensor is not empty
      explicit operator bool() const noexcept
      {
	return superbblas::detail::volume(size) > 0;
      }

      // Return the dimensions of the tensor
      std::map<char, int> kvdim() const
      {
	std::map<char, int> d;
	for (unsigned int i = 0; i < N; ++i)
	  d[order[i]] = size[i];
	return d;
      }

      /// Rename dimensions
      StorageTensor<N, T> rename_dims(const SB::remap& m) const
      {
	return StorageTensor<N, T>(*this, detail::update_order_and_check<N>(order, m), this->from,
				   this->size);
      }

      // Return a slice of the tensor starting at coordinate `kvfrom` and taking `kvsize` elements in each direction.
      // The missing dimension in `kvfrom` are set to zero and the missing direction in `kvsize` are set to the active size of the tensor.
      StorageTensor<N, T> kvslice_from_size(const std::map<char, int>& kvfrom = {},
					    const std::map<char, int>& kvsize = {}) const
      {
	std::map<char, int> updated_kvsize = this->kvdim();
	for (const auto& it : kvsize)
	  updated_kvsize[it.first] = it.second;
	return slice_from_size(kvcoors<N>(order, kvfrom), kvcoors<N>(order, updated_kvsize));
      }

      // Return a slice of the tensor starting at coordinate `from` and taking `size` elements in each direction.
      StorageTensor<N, T> slice_from_size(Coor<N> from, Coor<N> size) const
      {
	for (unsigned int i = 0; i < N; ++i)
	{
	  if (size[i] > this->size[i])
	    throw std::runtime_error(
	      "The size of the slice cannot be larger than the original tensor");
	  if (normalize_coor(from[i], this->size[i]) + size[i] > this->size[i] &&
	      this->size[i] != this->dim[i])
	    throw std::runtime_error(
	      "Unsupported to make a view on a non-contiguous range on the tensor");
	}

	using superbblas::detail::operator+;
	return StorageTensor<N, T>(*this, order, this->from + from, size, scalar);
      }

      StorageTensor<N, T> scale(T s) const
      {
	return StorageTensor<N, T>(*this, order, from, scalar * s);
      }

      void release()
      {
	dim = {{}};
	ctx.reset();
	from = {{}};
	size = {{}};
	scalar = T{0};
	filename = "";
	metadata = "";
      }

      /// Check that the dimension labels are valid

      void checkOrder() const
      {
	// Check that all labels are different there are N
	detail::check_order<N>(order);

	for (auto s : size)
	  if (s < 0)
	    std::runtime_error("Invalid tensor size: it should be positive");
      }

      /// Preallocate space for the storage file
      /// \param size: expected final file size in bytes

      void preallocate(std::size_t size)
      {
	superbblas::preallocate_storage(ctx.get(), size);
      }

      /// Save content from the storage into the given tensor
      template <std::size_t Nw, typename Tw,
		typename std::enable_if<
		  detail::is_complex<T>::value == detail::is_complex<Tw>::value, bool>::type = true>
      void copyFrom(const Tensor<Nw, Tw>& w) const
      {
	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (wsize[i] > size[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	MPI_Comm comm = (w.dist == Local ? MPI_COMM_SELF : MPI_COMM_WORLD);

	// If the storage is sparse, add blocks for the new content
	if (sparsity == Sparse)
	{
	  superbblas::append_blocks<Nw, N, T>(w.p->p.data(), w.p->p.size(), w.order.c_str(), w.from,
					      w.size, w.dim, order.c_str(), from, ctx.get(), comm,
					      superbblas::FastToSlow);
	}

	Tw* w_ptr = w.data();
	superbblas::save<Nw, N, Tw, T>(detail::safe_div<Tw>(w.scalar, scalar), w.p->p.data(), 1,
				       w.order.c_str(), w.from, w.size, w.dim, (const Tw**)&w_ptr,
				       &w.ctx(), order.c_str(), from, ctx.get(), comm,
				       superbblas::FastToSlow);
      }

      /// Load content from the storage into the given tensor
      template <std::size_t Nw, typename Tw,
		typename std::enable_if<
		  detail::is_complex<T>::value == detail::is_complex<Tw>::value, bool>::type = true>
      void copyTo(const Tensor<Nw, Tw>& w) const
      {
	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > wsize[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	Tw* w_ptr = w.data_for_writing();
	MPI_Comm comm = (w.dist == Local ? MPI_COMM_SELF : MPI_COMM_WORLD);
	superbblas::load<N, Nw, T, Tw>(detail::safe_div<T>(scalar, w.scalar), ctx.get(),
				       order.c_str(), from, size, w.p->p.data(), 1, w.order.c_str(),
				       w.from, w.dim, &w_ptr, &w.ctx(), comm,
				       superbblas::FastToSlow, superbblas::Copy);
	if (!w.is_managed())
	  superbblas::sync(w.ctx());
      }
    };

    /// Return a tensor filled with the value of the function applied to each element
    /// \param order: dimension labels, they should start with "xyztX"
    /// \param from: coordinates of the first element
    /// \param size: length of each tensor dimension
    /// \param dim: length of each global dimension (which is usually equal to size)
    /// \param dev: either OnHost or OnDefaultDevice
    /// \param func: function (Coor<N-1>) -> COMPLEX
    /// \param zero_is_even: (optional) whether the first element (`from`) is an even site (usually true)

    template <std::size_t N, typename COMPLEX, typename Func>
    Tensor<N, COMPLEX> fillLatticeField(const std::string& order, const std::map<char, int>& from,
					const std::map<char, int>& size,
					const std::map<char, int>& dim, DeviceHost dev, Func func,
					bool zero_is_even = true)
    {
      using superbblas::detail::operator+;

      static_assert(N >= 5, "The minimum number of dimensions should be 5");
      if (order.size() < 5 || order.compare(0, 5, "xyztX") != 0)
	throw std::runtime_error("Wrong `order`, it should start with xyztX");

      // Get final object dimension
      Coor<N> dim_c = latticeSize<N>(order, dim);
      std::map<char, int> size0 = dim;
      for (const auto& it : size)
	size0[it.first] = it.second;
      Coor<N> size_c = latticeSize<N>(order, size0);
      Coor<N> from_c = kvcoors<N>(order, from);

      // Populate the tensor on CPU
      Tensor<N, COMPLEX> r(order, size_c, OnHost);
      Coor<N> local_latt_size = r.p->localSize(); // local dimensions for xyztX
      Stride<N> stride =
	superbblas::detail::get_strides<std::size_t>(local_latt_size, superbblas::FastToSlow);
      Coor<N> local_latt_from =
	r.p->localFrom(); // coordinates of first elements stored locally for xyztX
      local_latt_from = local_latt_from + from_c;
      std::size_t vol = superbblas::detail::volume(local_latt_size);
      Index nX = r.kvdim()['X'];
      COMPLEX* ptr = r.data();
      int d = (zero_is_even ? 0 : 1);

#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
      for (std::size_t i = 0; i < vol; ++i)
      {
	// Get the global coordinates
	Coor<N> c = normalize_coor(
	  superbblas::detail::index2coor(i, local_latt_size, stride) + local_latt_from, dim_c);

	// Translate even-odd coordinates to natural coordinates
	Coor<N - 1> coor;
	coor[0] = c[0] * nX + (c[1] + c[2] + c[3] + c[4] + d) % nX; // x
	coor[1] = c[1];						    // y
	coor[2] = c[2];						    // z
	coor[3] = c[3];						    // t
	std::copy_n(c.begin() + 5, N - 5, coor.begin() + 4);

	// Call the function
	ptr[i] = func(coor);
      }

      return r.make_sure(none, dev);
    }

    ///
    /// Operators
    ///

    /// Ordering of matrices

    enum ColOrdering {
      RowMajor,	   ///< row-major ordering, the fastest index is the column
      ColumnMajor, ///< row-major ordering, the fastest index is the row
    };

    /// Operator's layout
    enum OperatorLayout {
      NaturalLayout,	     ///< natural ordering
      XEvenOddLayout,	     ///< X:(x+y+z+t)%2, x:x/2, y:y, z:z, t:t
      XEvenOddLayoutZeroOdd, ///< X:(x+y+z+t+1)%2, x:x/2, y:y, z:z, t:t
      EvensOnlyLayout	     ///< x:x/2, y:y, z:z, t:t for all (x+y+z+t)%2==0
    };

    /// Return whether the layout is XEvenOddLayout or XEvenOddLayoutZeroOdd
    /// \param layout: layout to test

    namespace detail
    {
      inline bool isEvenOddLayout(OperatorLayout layout)
      {
	return layout == XEvenOddLayout || layout == XEvenOddLayoutZeroOdd;
      }
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

    /// Displacements of each site nonzero edge for every operator's site

    using NaturalNeighbors = std::vector<std::map<char, int>>;

    /// Matrix to contract on each direction

    template <typename T>
    using SpinMatrixDir = std::map<std::map<char, int>, Tensor<2, T>>;

    /// Special value to indicate that the operator is dense

    inline const NaturalNeighbors& DenseOperator()
    {
      static const NaturalNeighbors dense{{{(char)0, 0}}};
      return dense;
    }

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
      /// Whether the spin-color block nonzeros are the tensor product of spin-spin and color-color matrices
      bool kron;
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
	       ColOrdering preferred_col_ordering, bool kron, const std::string& id = "")
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
	  kron(kron),
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
	  kron(sp.is_kronecker()),
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
	  kron(op.is_kronecker()),
	  id(std::make_shared<std::string>(id))
      {
      }

      /// Return the local support of this tensor as a subset of the global tensor
      Operator<NOp, COMPLEX> getGlocal() const
      {
	return fop ? Operator<NOp, COMPLEX>{fop, d.getGlocal(), i.getGlocal(), fop_tconj, *this}
		   : Operator<NOp, COMPLEX>{
		       sp,	rd,	   max_power, d.getGlocal(), i.getGlocal(),
		       order_t, domLayout, imgLayout, neighbors,     preferred_col_ordering};
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
	  fop_tconj, i, d, fop, order_t, imgLayout, domLayout, neighbors, preferred_col_ordering,
	  kron};
      }

      /// Return whether the operator has transpose conjugate
      bool has_tconj() const
      {
	return !sp && fop_tconj;
      }

      /// Return whether the spin-color nonzero blocks are the tensor product of two matrices
      bool is_kronecker() const
      {
	return sp ? sp.is_kronecker() : kron;
      }

      /// Return compatible domain tensors
      /// \param col_order: order for the columns
      /// \param m: column dimension size

      template <std::size_t N, typename T = COMPLEX>
      Tensor<N, T> make_compatible_dom(const std::string& col_order,
				       const std::map<char, int>& m) const
      {
	return d.template make_compatible<N, T>(
	  preferred_col_ordering == ColumnMajor ? std::string("%") + col_order : col_order + "%",
	  '%', "", m);
      }

      /// Return compatible image tensors
      /// \param col_order: order for the columns
      /// \param m: column dimension size

      template <std::size_t N, typename T = COMPLEX>
      Tensor<N, T> make_compatible_img(const std::string& col_order,
				       const std::map<char, int>& m) const
      {
	return i.template make_compatible<N, T>(
	  preferred_col_ordering == ColumnMajor ? std::string("%") + col_order : col_order + "%",
	  '%', "", m);
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
	  auto y = make_compatible_img<N, T>(cols, t.kvdim());
	  if (t.dist == Glocal)
	    y = y.getGlocal();
	  sp.contractWith(t.rename_dims(mcols), rd, y.rename_dims(mcols), {});
	  return y;
	}
	else
	{
	  auto x =
	    t.template collapse_dimensions<NOp + 1>(cols, 'n', true).template make_sure<COMPLEX>();
	  auto y = make_compatible_img<NOp + 1>("n", {{'n', x.kvdim()['n']}});
	  if (t.dist == Glocal)
	    y = y.getGlocal();
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
	    order_t, domLayout, imgLayout, neighbors, preferred_col_ordering, is_kronecker());
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
	  (new_d.kvdim().at('X') != d.kvdim().at('X') && detail::isEvenOddLayout(domLayout)
	     ? EvensOnlyLayout
	     : domLayout);
	OperatorLayout new_imgLayout =
	  (new_i.kvdim().at('X') != i.kvdim().at('X') && detail::isEvenOddLayout(imgLayout)
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
	    order_t, new_domLayout, new_imgLayout, neighbors, preferred_col_ordering, is_kronecker());
	}
      }
    };

    namespace detail
    {
      enum BlockingAsSparseDimensions {
	ConsiderBlockingSparse, ///< Dimensions 0,1,2,3 will be sparse and part of the lattice
	ConsiderBlockingDense	///< Dimensions 0,1,2,3 will be dense and not lattice dimensions
      };

      /// Return the natural lattice dimensions
      /// \param dim: dimension for each label
      /// \param layout: operator's layout

      inline std::map<char, int>
      getNatLatticeDims(const std::map<char, int>& dim, OperatorLayout layout,
			BlockingAsSparseDimensions blockDims = ConsiderBlockingSparse)
      {
	int nX = (layout == EvensOnlyLayout ? 2 : (dim.count('X') == 1 ? dim.at('X') : 1));
	if (blockDims == ConsiderBlockingSparse)
	{
	  return std::map<char, int>{
	    {'x', dim.at('x') * (dim.count('0') == 1 ? dim.at('0') : 1) * nX},
	    {'y', dim.at('y') * (dim.count('1') == 1 ? dim.at('1') : 1)},
	    {'z', dim.at('z') * (dim.count('2') == 1 ? dim.at('2') : 1)},
	    {'t', dim.at('t') * (dim.count('3') == 1 ? dim.at('3') : 1)}};
	}
	else
	{
	  return std::map<char, int>{
	    {'x', dim.at('x') * nX}, {'y', dim.at('y')}, {'z', dim.at('z')}, {'t', dim.at('t')}};
	}
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
	Coor<Nd> nat_dims =
	  kvcoors<Nd>("xyzt", getNatLatticeDims(dim, layout, ConsiderBlockingDense));
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
	auto t = fillLatticeField<N, float>(
		   order, {}, real_dim, real_dim, OnDefaultDevice,
		   [&](Coor<N - 1> c) {
		     return (float)coloring.getColor({{c[0], c[1], c[2], c[3]}});
		   },
		   layout == XEvenOddLayout)
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
	else if (isEvenOddLayout(layout))
	{
	  if (dim.at('X') != 2)
	    throw std::runtime_error("getXOddityMask: invalid dimension size `X`");
	  Coor<4> c;
	  for (std::size_t i = 0, j = 0; i < N; ++i)
	    if (is_in("Xyzt", r.order[i]))
	      c[j++] = i;
	  if (layout == XEvenOddLayoutZeroOdd)
	    xoddity = (xoddity + 1) % 2;
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
      /// \param useKronFormat: whether to create a Kronecker BSR variant if the given operator is in that format
      /// \return: a pair of a sparse tensor and a remap; the sparse tensor has the same image
      ///          labels as the given operator and domain labels are indicated by the returned remap.
      ///
      /// NOTE: Encoding the Dirac-Wilson with the clover term into the Kronecker format gets convoluted.
      /// We treat differently the block-diagonal (the self direction) from the others (the x,y,z,t directions).
      /// The nonzeros of the latter directions are the addition of two matrices which are the result of
      /// the tensor product of a 4x4 (spin matrix) and a 3x3 (color matrix). One of the spin matrices is
      /// is the identity always and the other is the same for all nonzeros in a direction. The block-diagonal
      /// doesn't follow this pattern but it is block diagonal on the chirality, that is, there are nonzeros only
      /// for the combination of spin i,j such that floor(i/2) == floor(j/2).

      template <std::size_t NOp, typename COMPLEX,
		typename std::enable_if<(NOp >= Nd + 1), bool>::type = true>
      std::pair<SpTensor<NOp, NOp, COMPLEX>, remap>
      cloneUnblockedOperatorToSpTensor(const Operator<NOp, COMPLEX>& op, unsigned int power = 1,
				       ColOrdering coBlk = RowMajor, bool useKronFormat = true,
				       const std::string& prefix = "")
      {
	using value_type = typename detail::base_type<COMPLEX>::type;

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

	// The spin label if the spin-color matrices are the tensor product of a spin matrix
	// and a color matrix
	char kronecker_label = op.is_kronecker() && useKronFormat ? 's' : 0;

	// Create the ordering for the domain and the image where the dense dimensions indices run faster than the sparse dimensions.
	// If using Kronecker variant, the Kronecker label (the spin) runs the slowest of the dense labels.
	// NOTE: assuming that x,y,z,t are the only sparse dimensions; X remains sparse
	std::string sparse_labels("xyztX");
	std::string dense_labels = remove_dimensions(op.i.order, sparse_labels);
	if (kronecker_label)
	{
	  dense_labels = remove_dimensions(dense_labels, std::string(1, kronecker_label)) +
			 std::string(1, kronecker_label);
	}
	remap rd = getNewLabels(op.d.order, op.i.order + "u~0123");
	auto i = op.i.reorder(dense_labels + std::string("xyztX"))
		   .like_this(none, {}, OnDefaultDevice, OnEveryone)
		   .make_eg();

	// Get the blocking for the domain and the image
	std::map<char, int> blkd, blki;
	for (const auto& it : i.kvdim())
	  blki[it.first] = (is_in(dense_labels, it.first) ? it.second : 1);
	for (const auto& it : blki)
	  blkd[rd.at(it.first)] = it.second;

	// Check that the Kroneker label is a dense label if given
	if (kronecker_label != 0 && std::find(dense_labels.begin(), dense_labels.end(),
					      kronecker_label) == dense_labels.end())
	  throw std::runtime_error("The Kronecker label should be a blocking label");

	// Construct the probing vectors, which they have as the rows the domain labels and as
	// columns the domain blocking dimensions

	constexpr int Nblk = NOp - Nd - 1;
	std::map<char, int> blki_id;
	for (char c : dense_labels)
	  blki_id[c] = blki[c];
	remap rd_id;
	for (char c : dense_labels)
	  rd_id[c] = rd[c];
	auto t_blk =
	  identity<Nblk * 2, COMPLEX>(blki_id, rd_id).make_sure(none, none, OnEveryoneReplicated);

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
	else
	  throw std::runtime_error("Unsupported power");

	// Extend directions in case of using the Kronecker form
	if (kronecker_label != 0 && i.kvdim().at(kronecker_label) != 4)
	  throw std::runtime_error(
	    "Unsupported extraction of the Kronecker format from this operator");
	std::vector<Tensor<2, COMPLEX>> spin_matrix;
	if (kronecker_label != 0)
	{
	  std::vector<Coor<Nd>> new_neighbors;
	  int spin = i.kvdim().at(kronecker_label);
	  for (const auto& dir : neighbors)
	  {
	    if (dir == Coor<Nd>{{}})
	    {
	      // If self direction, create a single matrix on each combination of spins
	      // with the same chirality

	      for (int s = 0; s < spin * spin; ++s)
	      {
		int si = s % spin, sj = s / spin;
		Tensor<2, COMPLEX> mat(std::string(1, kronecker_label) +
					 std::string(1, rd.at(kronecker_label)),
				       {{spin, spin}}, OnHost, OnEveryoneReplicated);
		mat.set_zero();
		mat.set({{si, sj}}, 1);
		new_neighbors.push_back(dir);
		spin_matrix.push_back(mat);
	      }
	    }
	    else
	    {
	      // For the remaining directions, we capture the spin block diagonal on the first term
	      // and put an empty matrix on the second term so that will be guess further down
	      Tensor<2, COMPLEX> mat(std::string(1, kronecker_label) +
				       std::string(1, rd.at(kronecker_label)),
				     {{spin, spin}}, OnHost, OnEveryoneReplicated);
	      mat.set_zero();
	      for (int s = 0; s < spin; ++s)
		mat.set({{s, s}}, 1);
	      new_neighbors.push_back(dir);
	      spin_matrix.push_back(mat);
	      new_neighbors.push_back(dir);
	      spin_matrix.push_back(Tensor<2, COMPLEX>());
	    }
	  }
	  neighbors = new_neighbors;
	}

        // Chose a dense label that is not the spin
        const char color_label = 'c';

        // Extract the kronecker values with probing
	Tensor<3, COMPLEX> kron;
	std::vector<std::map<char, int>> nonzero_spins;
	if (kronecker_label != 0)
	{
	  unsigned int color = 0;

	  // Extracting the proving vector
	  auto t_l = colors.template transformWithCPUFun<COMPLEX>(
	    [&](float site_color) { return site_color == color ? value_type{1} : value_type{0}; });

	  // Skip empty masks
	  // NOTE: the call to split_dimension add a fake dimension that acts as columns
	  if (std::norm(norm<1>(t_l.split_dimension('X', "Xn", maxX), "n").get({{0}})) == 0)
	    throw std::runtime_error("Ups! We should do something more sophisticated here");

	  // Contracting the proving vector with the blocking components
	  auto site_size = blki;
	  for (char c : dense_labels)
	    site_size[rd[c]] = (c == kronecker_label ? blki[c] : 1);

	  // Compute the matvecs
	  auto mv =
	    op(contract<NOp + Nblk>(t_l, t_blk.kvslice_from_size({}, {{rd[color_label], 1}}), ""));

	  // Take a source
          auto source_coor =
              colors.find([&](float site_color) { return site_color == color; })
                  .getSome();
          std::map<char, int> source;
          for (std::size_t i = 0; i < colors.order.size(); ++i)
            if (is_in("xyztX", colors.order[i]))
              source[colors.order[i]] = source_coor[i];

          // Find the spin values for each direction
          std::string kron_order(3, 0);
          kron_order[0] = kronecker_label;
          kron_order[1] = rd[kronecker_label];
          kron_order[2] = 'u';
          kron = Tensor<3, COMPLEX>(kron_order,
                                    Coor<3>{blki[kronecker_label],
                                            blki[kronecker_label],
                                            (int)neighbors.size()},
                                    OnDefaultDevice, OnEveryoneReplicated);
          kron.set_zero();
          for (int mu = 0; mu < neighbors.size(); ++mu) {
	    // site = source + neighbors[mu], where the latter is in natural
	    // coordinates
	    const auto& coor_dir = neighbors[mu];
	    int sumdir = std::accumulate(coor_dir.begin(), coor_dir.end(), int{0});
	    std::map<char, int> site{{'X', source['X'] + sumdir},
				     {'x', (source['x'] * real_maxX + coor_dir[0]) / real_maxX},
				     {'y', source['y'] + coor_dir[1]},
				     {'z', source['z'] + coor_dir[2]},
				     {'t', source['t'] + coor_dir[3]}};
	    auto site_size = blki;
	    for (char c : dense_labels)
	      site_size[rd[c]] = (c == kronecker_label ? blki[c] : 1);

	    auto site_data =
	      mv.kvslice_from_size(site, site_size).make_sure(none, OnHost, OnEveryoneReplicated);

            if (!spin_matrix[mu]) {
              // Search for a nonzero element, we take the largest.
              // NOTE: don't take from spin block diagonal matrix, those
              // nonzeros are captured already
              auto spin_vals =
                  norm<2>(site_data, std::string(1, kronecker_label) +
                                         std::string(1, rd[kronecker_label]));
              int s_ref = 0;
              double val_ref = 0;
              for (int s = 0; s < blki[kronecker_label] * blki[kronecker_label];
                   ++s) {
                if (s % blki[kronecker_label] == s / blki[kronecker_label])
                  continue;
		double val = spin_vals.get(
		  kvcoors<2>(spin_vals.order, {{kronecker_label, s % blki[kronecker_label]},
					       {rd[kronecker_label], s / blki[kronecker_label]}}));
		if (val > val_ref) {
                  s_ref = s;
                  val_ref = val;
                }
              }

              // If the direction is empty, remove it!
              if (val_ref == 0) {
                neighbors.erase(neighbors.begin() + mu);
                spin_matrix.erase(spin_matrix.begin() + mu);
                mu--;
                continue;
              }

              nonzero_spins.push_back(
                  {{kronecker_label, s_ref % blki[kronecker_label]},
                   {rd[kronecker_label], s_ref / blki[kronecker_label]}});

              // Get the values
	      auto val0 = site_data.get(kvcoors<NOp + Nblk>(
		site_data.order, {{kronecker_label, s_ref % blki[kronecker_label]},
				  {rd[kronecker_label], s_ref / blki[kronecker_label]},
				  {color_label, 0},
				  {rd.at(color_label), 0}}));

	      for (int s = 0; s < blki[kronecker_label] * blki[kronecker_label];
                   ++s) {
                if (s % blki[kronecker_label] == s / blki[kronecker_label])
                  continue;
		kron.set(kvcoors<3>(kron.order, {{kronecker_label, s % blki[kronecker_label]},
						 {rd[kronecker_label], s / blki[kronecker_label]},
						 {'u', mu}}),
			 site_data.get(kvcoors<NOp + Nblk>(
			   site_data.order, {{kronecker_label, s % blki[kronecker_label]},
					     {rd[kronecker_label], s / blki[kronecker_label]},
					     {color_label, 0},
					     {rd.at(color_label), 0}})) /
			   val0);
	      }
            } else {
	      if (std::norm(
		    norm<1>(contract<NOp + Nblk + 1>(spin_matrix[mu].template reshape_dimensions<3>(
						       {{"s", "su"}}, {{'u', 1}}),
						     site_data, ""),
			    "u")
		      .get(Coor<1>{0})) == 0)
	      {
		neighbors.erase(neighbors.begin() + mu);
                spin_matrix.erase(spin_matrix.begin() + mu);
                mu--;
                continue;
	      }

	      // Copy the know spin matrix into sop.kron
              spin_matrix[mu].copyTo(
                  kron.kvslice_from_size({{'u', mu}}, {{'u', 1}}));

              // Set as the reference spin, the first nonzero
              for (int s = 0; s < blki[kronecker_label] * blki[kronecker_label];
                   ++s) {
                if (std::norm(spin_matrix[mu].get(
                        {{s % blki[kronecker_label],
                          s / blki[kronecker_label]}})) > 0) {
                  nonzero_spins.push_back(
                      {{kronecker_label, s % blki[kronecker_label]},
                       {rd[kronecker_label], s / blki[kronecker_label]}});
                  break;
                }
              }
            }
	  }
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
	const unsigned int Nkron = kronecker_label == 0 ? 0u : 1u;
	SpTensor<NOp, NOp, COMPLEX> sop{d_sop,
					i,
					Nblk - Nkron,
					Nblk - Nkron,
					Nkron,
					Nkron,
					(unsigned int)neighbors.size(),
					coBlk == ColumnMajor};

	// Copy the kronecker values
	if (kronecker_label != 0)
	{
	  kron.kvslice_from_size({}, {{'u', neighbors.size()}}).copyTo(sop.kron);
	}

	// Extract the nonzeros with probing
	sop.data.set_zero(); // all values may not be populated when using blocking
	for (unsigned int color = 0; color < num_colors; ++color)
	{
	  // Generate the proving vectors for the given color
	  auto t_l = colors.template transformWithCPUFun<COMPLEX>(
	    [&](float site_color) { return site_color == color ? value_type{1} : value_type{0}; });

	  // Skip empty masks
	  // NOTE: the call to split_dimension add a fake dimension that acts as columns
          if (std::norm(norm<1>(t_l.split_dimension('X', "Xn", maxX), "n")
                            .get({{0}})) == 0)
            continue;

          for (int color_idx = 0; color_idx < blki[color_label]; ++color_idx) {
            std::map<char, int> colorFrom{{rd[color_label], color_idx}};
            std::map<char, int> colorSize{{rd[color_label], 1}};

            // Contracting the proving vector with the blocking components
	    auto probs =
	      contract<NOp + Nblk>(t_l, t_blk.kvslice_from_size(colorFrom, colorSize), "");

	    // Compute the matvecs
            auto mv = op(std::move(probs));

            // Construct an indicator tensor where all blocking dimensions but
            // only the nodes colored `color` are copied
            auto color_mask = t_l.template transformWithCPUFun<float>(
                [](const value_type &t) { return (float)std::real(t); });
            auto sel_x_even = contract<NOp + Nblk>(
                contract<NOp>(color_mask, even_x_mask, ""),
                ones_blk.kvslice_from_size(colorFrom, colorSize), "");
            auto sel_x_odd = contract<NOp + Nblk>(
                contract<NOp>(color_mask, odd_x_mask, ""),
                ones_blk.kvslice_from_size(colorFrom, colorSize), "");

            // Populate the nonzeros by copying pieces from `mv` into sop.data.
            // We want to copy only the nonzeros in `mv`, which are `neighbors`
            // away from the nonzeros of `probs`.
            auto sop_data = sop.data.kvslice_from_size(colorFrom, colorSize);
            if (kronecker_label == 0) {
              latticeCopyToWithMask(mv, sop_data, 'u', neighbors,
                                    {{'X', real_maxX}}, sel_x_even, sel_x_odd);
            } else {
              std::map<char, int> single_spin{{kronecker_label, 1},
                                              {rd[kronecker_label], 1}};
              for (int dir = 0; dir < neighbors.size(); ++dir)
                latticeCopyToWithMask(
                    mv.kvslice_from_size(nonzero_spins[dir], single_spin),
                    sop_data.kvslice_from_size({{'u', dir}}, {{'u', 1}}), 'u',
                    std::vector<Coor<Nd>>(1, neighbors[dir]),
                    {{'X', real_maxX}},
                    sel_x_even.kvslice_from_size(nonzero_spins[dir],
                                                 single_spin),
                    sel_x_odd.kvslice_from_size(nonzero_spins[dir],
                                                single_spin));
            }
          }
        }

        // Populate the coordinate of the columns, that is, to give the domain coordinates of first nonzero in each
	// BSR nonzero block. Assume that we are processing nonzeros block for the image coordinate `c` on the
	// direction `dir`, that is, the domain coordinates will be (cx-dirx,cy-diry,cz-dirz,dt-dirt) in natural
	// coordinates. But we get the image coordinate `c` in even-odd coordinate, (cX,cx,cy,cz,ct), which has the
	// following natural coordinates (cx*2+(cX+cy+cz+ct)%2,cy,cz,ct). After subtracting the direction we get the
	// natural coordinates (cx*2+(cX+cy+cz+ct)%2-dirx,cy-diry,cz-dirz,ct-dirt), which corresponds to the following
	// even-odd coordinates ((cX-dirx-diry-dirz-dirt)%2,(cx*2+(cX+cy+cz+ct)%2-dirx)/2,cy-diry,cz-dirz,ct-dirt).

	Coor<Nd> real_dims = kvcoors<Nd>("xyzt", getNatLatticeDims(i.kvdim(), op.imgLayout));
	int d = op.imgLayout == XEvenOddLayoutZeroOdd ? 1 : 0;
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
	    int sumyzt = std::accumulate(c.begin() + 2 + Nblk + 1, c.end() - 1, d);
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
      /// \param useKronFormat: whether to create a Kronecker BSR variant if the given operator is in that format
      /// \return: a pair of a sparse tensor and a remap; the sparse tensor has the same image
      ///          labels as the given operator and domain labels are indicated by the returned remap.

      template <std::size_t NOp, typename COMPLEX,
		typename std::enable_if<(NOp > 9), bool>::type = true>
      std::pair<SpTensor<NOp, NOp, COMPLEX>, remap>
      cloneOperatorToSpTensor(const Operator<NOp, COMPLEX>& op, unsigned int power,
			      ColOrdering coBlk = RowMajor, bool useKronFormat = true,
			      const std::string& prefix = "")
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
	  coBlk,
	  op.is_kronecker()};

	// Get a sparse tensor representation of the operator
	unsigned int op_dist = getFurthestNeighborDistance(op);
	if (op_dist > 1 && power % op_dist != 0)
	  throw std::runtime_error("cloneOperatorToSpTensor: invalid power value, it isn't "
				   "divisible by the furthest neighbor distance");
	auto opdims = op.i.kvdim();
	auto t = cloneUnblockedOperatorToSpTensor(unblocked_op, std::min(power, op_dist), coBlk,
						  useKronFormat, prefix);
	remap rd = t.second;
	for (const auto& it :
	     detail::getNewLabels("0123", op.d.order + op.i.order + "0123" + t.first.d.order))
	  rd[it.first] = it.second;
	int max_op_power = (op_dist == 0 ? 0 : std::max(power / op_dist, 1u) - 1u);
	std::map<char, int> m_power{};
	for (char c : update_order("xyzt", rd))
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

      template <std::size_t NOp, typename COMPLEX,
		typename std::enable_if<(NOp <= 9), bool>::type = true>
      std::pair<SpTensor<NOp, NOp, COMPLEX>, remap>
      cloneOperatorToSpTensor(const Operator<NOp, COMPLEX>& op, unsigned int power,
			      ColOrdering coBlk = RowMajor, bool = true,
			      const std::string& prefix = "")
      {
	throw std::runtime_error("trying to clone an unblock operator with a blocking function");
      }

      /// Return an efficient operator application
      /// \param op: operator to extract the nonzeros from
      /// \param power: maximum distance to recover the nonzeros:
      ///               0, block diagonal; 1: near-neighbors...
      /// \param co: preferred ordering of dense input and output tensors
      /// \param coBlk: ordering of the nonzero blocks of the sparse operator

      template <std::size_t NOp, typename COMPLEX>
      Operator<NOp, COMPLEX>
      cloneOperator(const Operator<NOp, COMPLEX>& op, unsigned int power, ColOrdering co,
		    ColOrdering coBlk,
		    BlockingAsSparseDimensions blockingAsSparseDimensions = ConsiderBlockingSparse,
		    const std::string& prefix = "")
      {
	// If the operator is empty, just return itself
	if (op.d.volume() == 0 || op.i.volume() == 0)
	  return op;

	// Get a sparse tensor representation of the operator
	unsigned int op_dist = getFurthestNeighborDistance(op);
	if (op_dist > 1 && power % op_dist != 0)
	  throw std::runtime_error("cloneOperator: invalid power value");
	remap rd;
	SpTensor<NOp, NOp, COMPLEX> sop;
	if (blockingAsSparseDimensions == ConsiderBlockingSparse)
	{
	  auto t = cloneOperatorToSpTensor(op, power, coBlk, true /* use Kron format if possible */,
					   prefix);
	  sop = t.first;
	  rd = t.second;
	}
	else
	{
	  auto t = cloneUnblockedOperatorToSpTensor(op, power, coBlk,
						    true /* use Kron format if possible */, prefix);
	  sop = t.first;
	  rd = t.second;
	}

	// Construct the operator to return
	Operator<NOp, COMPLEX> rop{sop,	       rd,	     power,	   sop.i,	 sop.i,
				   op.order_t, op.domLayout, op.imgLayout, op.neighbors, co};

	// Skip tests if power < op_dist
	if (power < op_dist)
	  return rop;

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
      Operator<NOp, COMPLEX>
      cloneOperator(const Operator<NOp, COMPLEX>& op, ColOrdering co, ColOrdering coBlk,
		    BlockingAsSparseDimensions blockingAsSparseDimensions = ConsiderBlockingSparse,
		    const std::string& prefix = "")
      {
	return cloneOperator(op, getFurthestNeighborDistance(op), co, coBlk,
			     blockingAsSparseDimensions, prefix);
      }

    }

    namespace detail
    {
      inline std::mt19937_64& getSeed()
      {
	//  This is quick and dirty and nonreproducible if the lattice is distributed
	//  among the processes is different ways.
	static std::mt19937_64 twister_engine(10 + Layout::nodeNumber());
	return twister_engine;
      }
    }

    /// Modify with complex random uniformly distributed numbers with the real and the imaginary part between [a,b]
    /// \param t: tensor to fill with random numbers
    /// \param a: minimum random value
    /// \param b: maximum random value

    template <std::size_t N, typename T,
	      typename std::enable_if<detail::is_complex<T>::value, bool>::type = true>
    void urand(Tensor<N, T> t, typename T::value_type a = 0, typename T::value_type b = 1)
    {
      std::uniform_real_distribution<typename T::value_type> d(a, b);
      t.fillWithCPUFuncNoArgs(
	[&]() {
	  return T{d(detail::getSeed()), d(detail::getSeed())};
	},
	false);
    }

    /// Modify with random uniformly distributed numbers between [a,b]
    /// \param t: tensor to fill with random numbers
    /// \param a: minimum random value
    /// \param b: maximum random value

    template <std::size_t N, typename T,
	      typename std::enable_if<!detail::is_complex<T>::value, bool>::type = true>
    void urand(Tensor<N, T> t, typename detail::base_type<T>::type a = 0,
	       typename detail::base_type<T>::type b = 1)
    {
      std::uniform_real_distribution<typename detail::base_type<T>::type> d(a, b);
      t.fillWithCPUFuncNoArgs([&]() { return d(detail::getSeed()); }, false);
    }

    /// Modify with complex random normal distributed numbers
    /// \param t: tensor to fill with random numbers

    template <std::size_t N, typename T,
	      typename std::enable_if<detail::is_complex<T>::value, bool>::type = true>
    void nrand(Tensor<N, T> t)
    {
      std::normal_distribution<typename T::value_type> d{};
      t.fillWithCPUFuncNoArgs(
	[&]() {
	  return T{d(detail::getSeed()), d(detail::getSeed())};
	},
	false);
    }

    /// Modify with random normal distributed numbers
    /// \param t: tensor to fill with random numbers

    template <std::size_t N, typename T,
	      typename std::enable_if<!detail::is_complex<T>::value, bool>::type = true>
    void nrand(Tensor<N, T> t)
    {
      std::normal_distribution<typename detail::base_type<T>::type> d{};
      t.fillWithCPUFuncNoArgs([&]() { return d(detail::getSeed()); }, false);
    }

    /// Compute a shift of v onto the direction dir
    /// \param v: tensor to apply the displacement
    /// \param first_tslice: global index in the t direction of the first element
    /// \param len: step of the displacement
    /// \param dir: 0 is x; 1 is y...

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> shift(const Tensor<N, COMPLEX>& v, Index first_tslice, int len, int dir,
			     Maybe<Action> action = none, Tensor<N, COMPLEX> w = {})
    {
      if (dir < 0 || dir >= Nd - 1)
	throw std::runtime_error("Invalid direction");

      if (action.hasSome() != (bool)w)
	throw std::runtime_error("Invalid default value");

      // Address zero length case
      if (len == 0)
      {
	if (!w)
	  return v;
	v.doAction(action.getSome(), w);
	return w;
      }

      // NOTE: chroma uses the reverse convention for direction: shifting FORWARD moves the sites on the negative direction
      len = -len;

      const char dir_label[] = "xyz";
#  if QDP_USE_LEXICO_LAYOUT
      // If we are not using red-black ordering, return a view where the tensor is shifted on the given direction
      v = v.kvslice_from_size({{dir_label[dir], -len}});

      if (!w)
	return v;

      v.doAction(action, w);
      return w;

#  elif QDP_USE_CB2_LAYOUT
      // Assuming that v has support on the origin and destination lattice elements
      int dimX = v.kvdim()['X'];
      if (dimX != 2 && len % 2 != 0)
	throw std::runtime_error("Unsupported shift");

      if (dir != 0)
      {
	if (!w)
	  return v.kvslice_from_size({{'X', -len}, {dir_label[dir], -len}});
	v.doAction(action.getSome(), w.kvslice_from_size({{'X', len}, {dir_label[dir], len}}));
	return w;
      }
      else
      {
	auto dims = v.kvdim();
	int t = dims.at('t');
	if (t > 1 && t % 2 == 1)
	  throw std::runtime_error(
	    "The t dimension should be zero, one, or even when doing shifting on the X dimension");
	int maxX = dims.at('X');
	int maxY = std::min(2, dims.at('y'));
	int maxZ = std::min(2, dims.at('z'));
	int maxT = std::min(2, dims.at('t'));
	auto v_eo = v.split_dimension('y', "Yy", maxY)
		      .split_dimension('z', "Zz", maxZ)
		      .split_dimension('t', "Tt", maxT);
	Tensor<N, COMPLEX> r = w ? w : v.like_this();
	auto r_eo = r.split_dimension('y', "Yy", maxY)
		      .split_dimension('z', "Zz", maxZ)
		      .split_dimension('t', "Tt", maxT);
	//.make_writing_nonatomic();
	while (len < 0)
	  len += dims.at('x') * maxX;
	for (int T = 0; T < maxT; ++T)
	{
	  for (int Z = 0; Z < maxZ; ++Z)
	  {
	    for (int Y = 0; Y < maxY; ++Y)
	    {
	      for (int X = 0; X < maxX; ++X)
	      {
		auto v_eo_slice = v_eo.kvslice_from_size({{'X', X}, {'Y', Y}, {'Z', Z}, {'T', T}},
							 {{'X', 1}, {'Y', 1}, {'Z', 1}, {'T', 1}});
		auto r_eo_slice =
		  r_eo.kvslice_from_size({{'X', X + len},
					  {'x', (len + ((X + Y + Z + T + first_tslice) % 2)) / 2},
					  {'Y', Y},
					  {'Z', Z},
					  {'T', T}},
					 {{'Y', 1}, {'Z', 1}, {'T', 1}});
		v_eo_slice.doAction(action.getSome(CopyTo), r_eo_slice);
	      }
	    }
	  }
	}
	return r;
      }
#  else
      throw std::runtime_error("Unsupported layout");
#  endif
    }

    /// Compute a displacement of v onto the direction dir
    /// \param u: Gauge field
    /// \param v: tensor to apply the displacement
    /// \param first_tslice: global index in the t direction of the first element
    /// \param dir: 0: nothing; 1: forward x; -1: backward x; 2: forward y...

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> displace(const std::vector<Tensor<Nd + 3, COMPLEX>>& u, Tensor<N, COMPLEX> v,
				Index first_tslice, int dir, Maybe<Action> action = none,
				Tensor<N, COMPLEX> w = {})
    {
      if (std::abs(dir) > Nd)
	throw std::runtime_error("Invalid direction");

      if (action.hasSome() != (bool)w)
	throw std::runtime_error("Invalid default value");

      // Address the zero direction case
      if (dir == 0)
      {
	if (!w)
	  return v;
	v.doAction(action.getSome(), w);
	return w;
      }

      int d = std::abs(dir) - 1;    // space lattice direction, 0: x, 1: y, 2: z
      int len = (dir > 0 ? 1 : -1); // displacement unit direction
      assert(d < u.size());

      if (len > 0)
      {
	// Do u[d] * shift(x,d)
	Tensor<N, COMPLEX> r = w ? w : v.like_this();
	v = shift(std::move(v), first_tslice, len, d);
	r.contract(std::move(v), {}, NotConjugate, u[d], {{'j', 'c'}}, NotConjugate, {{'c', 'i'}},
		   action.getSome(CopyTo) == CopyTo ? 0.0 : 1.0);
	return r;
      }
      else
      {
	// Do shift(adj(u[d]) * x,d)
	Tensor<N, COMPLEX> r = v.like_this();
	r.contract(std::move(v), {}, NotConjugate, u[d], {{'i', 'c'}}, Conjugate, {{'c', 'j'}});
	return shift(std::move(r), first_tslice, len, d, action, w);
      }
    }

    /// Apply right nabla onto v on the direction dir
    /// \param u: Gauge field
    /// \param v: tensor to apply the derivative
    /// \param first_tslice: global index in the t direction of the first element
    /// \param dir: 0: nothing; 1: forward x; -1: backward x; 2: forward y...
    ///
    /// NOTE: the code returns U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> rightNabla(const std::vector<Tensor<Nd + 3, COMPLEX>>& u,
				  Tensor<N, COMPLEX> v, Index first_tslice, int dir)
    {
      auto r = displace(u, v, first_tslice, dir);
      displace(u, v, first_tslice, -dir).scale(-1).addTo(r);
      return r;
    }

    /// Compute a displacement of v onto the direction dir
    /// \param u: Gauge field
    /// \param v: tensor to apply the derivative
    /// \param first_tslice: global index in the t direction of the first element
    /// \param dir: 0: nothing; 1: forward x; -1: backward x; 2: forward y...
    /// \param moms: list of input momenta
    /// \param conjUnderAdd: if true, return a version, R(dir), so that
    ////       adj(R(dir)) * D(dir') == D(dir+dir'), where D(dir') is what this function returns
    ////       when conjUnderAdd is false.

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> leftRightNabla(const std::vector<Tensor<Nd + 3, COMPLEX>>& u,
				      Tensor<N, COMPLEX> v, Index first_tslice, int dir,
				      const std::vector<Coor<3>>& moms = {},
				      bool conjUnderAdd = false)
    {
      if (std::abs(dir) > Nd)
	throw std::runtime_error("Invalid direction");

      int d = std::abs(dir) - 1; // space lattice direction, 0: x, 1: y, 2: z

      // conj(phase)*displace(u, v, -dir) - phase*displace(u, v, dir)
      std::vector<COMPLEX> phases(moms.size());
      for (unsigned int i = 0; i < moms.size(); ++i)
      {

	typename COMPLEX::value_type angle = 2 * M_PI * moms[i][d] / Layout::lattSize()[d];
	phases[i] = COMPLEX{1} + COMPLEX{cos(angle), sin(angle)};
	if (conjUnderAdd)
	  phases[i] = std::sqrt(phases[i]);
      }

      // r = conj(phases) * displace(u, v, dir)
      Tensor<N, COMPLEX> r = v.like_this("c%xyzXtm", '%');
      r.contract(displace(u, v, first_tslice, -dir), {}, NotConjugate, asTensorView(phases),
		 {{'i', 'm'}}, Conjugate);

      // r = r - phases * displace(u, v, dir) if !ConjUnderAdd else r + phases * displace(u, v, dir)
      r.contract(displace(u, v, first_tslice, dir).scale(conjUnderAdd ? 1 : -1), {}, NotConjugate,
		 asTensorView(phases), {{'i', 'm'}}, NotConjugate, {}, 1.0);

      return r;
    }

    // template <std::size_t N, typename T>
    // class Transform : public Tensor<N,T> {
    // public:
    //   Transform(Tensor<N, T> t, remap input, remat output, std::string no_contract = std::string(),
    //     	bool conj = false)
    //     : Tensor<N, T>(t), input(input), output(output), no_contract(no_contract), conj(conj)
    //   {
    //   }

    //   Transform<N, T> rename_dims(remap m) const
    //   {
    //     remap new_input = input, new_output = output;
    //     std::string new_no_contract = new_constract;
    //     for (auto& it : new_input)
    //     {
    //       auto j = m.find(it->first);
    //       if (j != m.end())
    //         it->second = j->second;
    //     }
    //   	for (auto& it : new_output)
    //     {
    //       auto j = m.find(it->first);
    //       if (j != m.end())
    //         it->second = j->second;
    //     }
    //     for (auto& it : m) {
    //         auto n = new_no_contract.find(it->first);
    //         if (n != std::string::npos)
    //           new_no_contract[n] = it->second;
    //     }
    //     return Transform<N, T>(*this, new_input, new_output, new_no_contract, conj);
    //   }

    //   const remap input;  ///< rename dimensions before contraction
    //   const remap output; ///< rename dimensions after contraction
    //   const std::string
    //     no_contract; ///< dimension labels that should not be on the other contracted tensor
    //   const bool conj;

    //   virtual Transform<N, T> adj() const
    //   {
    //     return Transform<N, T>{*this, input, output, no_contract, !conj};
    //   }
    // };

    // template <std::size_t Nt>
    // Tensor<Nt, T> contractTo(Tensor<Nt, T> t) const
    // {
    // }

    // template <typename T>
    // class MatrixTransform : Transform<2, T>
    // {
    // public:
    //   MatrixTransform(Transform<2, T> t) : Transform<2, T>(t)
    //   {
    //   }

    //   MatrixTransform<N, T> adj() const override
    //   {
    //     remap m{{order[0], order[1]}, {order[1], order[0]}};
    //     return MatrixTransform<N, T>{
    //       Transform<N, T>{this->Tensor<N, T>::rename_dims(m), input, output, no_contract, !conj}};
    //   }
    // };

    //
    // Get colorvecs
    //

    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, Tensor<Nd, ComplexF>> MODS_t;
    typedef QDP::MapObjectDisk<KeyTimeSliceColorVec_t, Tensor<Nd, ComplexF>> MOD_t;

    // Represent either FILEDB or S3T handle for file containing distillation vectors

    struct ColorvecsStorage {
      std::shared_ptr<MODS_t> mod;	   // old storage
      StorageTensor<Nd + 2, ComplexD> s3t; // cxyztn
    };

    namespace ns_getColorvecs
    {
      /// Return the permutation from a natural layout to a red-black that is used by `applyPerm`
      /// \param t: index in the t-direction of the input elements
      /// NOTE: assuming the input layout is xyz and the output layout is xyzX for the given input
      ///       time-slice `t`

      inline std::vector<Index> getPermFromNatToRB(Index t)
      {
	if (Layout::nodeNumber() != 0)
	  return {};

	const Index x1 = Layout::lattSize()[0];
	const Index y1 = Layout::lattSize()[1];
	const Index z1 = Layout::lattSize()[2];
	std::vector<Index> perm(x1 * y1 * z1);

#  if QDP_USE_LEXICO_LAYOUT
	unsigned int n = x1 * y1 * z1;
#    ifdef _OPENMP
#      pragma omp parallel for schedule(static)
#    endif
	for (unsigned int i = 0; i < n; ++n)
	  perm[i] = i;

#  elif QDP_USE_CB2_LAYOUT
#    ifdef _OPENMP
#      pragma omp parallel for collapse(3) schedule(static)
#    endif
	for (unsigned int z = 0; z < z1; ++z)
	{
	  for (unsigned int y = 0; y < y1; ++y)
	  {
	    for (unsigned int x = 0; x < x1; ++x)
	    {
	      // index on natural ordering
	      Index i0 = x + y * x1 + z * x1 * y1;
	      // index in red-black
	      Index i1 = x / 2 + y * (x1 / 2) + z * (x1 * y1 / 2) +
			 ((x + y + z + t) % 2) * (x1 * y1 * z1 / 2);
	      perm[i1] = i0;
	    }
	  }
	}

#  else
	throw std::runtime_error("Unsupported layout");
#  endif

	return perm;
      }

      /// Apply a permutation generated by `getPermFromNatToRB`
      /// \param perm: permutation generated with getPermFromNatToRB
      /// \param tnat: input tensor with ordering cxyz
      /// \param trb: output tensor with ordering cxyzX

      template <typename T>
      void toRB(const std::vector<Index>& perm, Tensor<Nd, T> tnat, Tensor<Nd + 1, T> trb)
      {
	assert(tnat.order == "cxyz");
	assert(trb.order == "cxyzX");
	assert(tnat.p->localVolume() == perm.size() * Nc);

	unsigned int i1 = perm.size();
	const T* x = tnat.data();
	T* y = trb.data();

#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	for (unsigned int i = 0; i < i1; ++i)
	  for (unsigned int c = 0; c < Nc; ++c)
	    y[i * Nc + c] = x[perm[i] * Nc + c];
      }

      /// Apply a permutation generated by `getPermFromNatToRB`
      /// \param perm: permutation generated with getPermFromNatToRB
      /// \param tnat: input tensor with ordering cxyz
      /// \param trb: output tensor with ordering cxyzX

      template <typename T>
      void toNat(const std::vector<Index>& perm, Tensor<Nd + 1, T> trb, Tensor<Nd, T> tnat)
      {
	assert(tnat.order == "cxyz");
	assert(trb.order == "cxyzX");
	assert(tnat.p->localVolume() == perm.size() * Nc);

	unsigned int i1 = perm.size();
	T* x = tnat.data();
	const T* y = trb.data();

#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	for (unsigned int i = 0; i < i1; ++i)
	  for (unsigned int c = 0; c < Nc; ++c)
	    x[perm[i] * Nc + c] = y[i * Nc + c];
      }

      /// Return a lattice field with value exp(2*pi*(x./dim)'*phase) for each lattice site x
      /// \param phase: integer phase
      /// \param dev: device of the returned tensor

      template <typename T>
      Tensor<Nd + 1, T> getPhase(Coor<Nd - 1> phase, int tfrom, int tsize,
				 DeviceHost dev = OnDefaultDevice)
      {
	// Get spatial dimensions of the current lattice
	Coor<Nd> dim = latticeSize<Nd>("xyzX", {});
	dim[0] *= dim[3];
	return fillLatticeField<5, T>("xyztX", {{'t', tfrom}}, {{'t', tsize}}, {}, dev,
				      [=](Coor<Nd> c) {
					typename T::value_type phase_dot_coor = 0;
					for (int i = 0; i < Nd - 1; ++i)
					  phase_dot_coor += c[i] * 2 * M_PI * phase[i] / dim[i];

					return T{cos(phase_dot_coor), sin(phase_dot_coor)};
				      });
      }

#  if defined(BUILD_PRIMME)

      // Apply the laplacian operator on the spatial dimensions
      /// \param u: Gauge fields restricted to the same t-slice as chi and psi
      /// \param first_tslice: global t index of the zero t index
      /// \param chi: output vector
      /// \param psi: input vector

      inline void LaplacianOperator(const std::vector<Tensor<Nd + 3, ComplexD>>& u,
				    Index first_tslice, Tensor<Nd + 3, ComplexD> chi,
				    const Tensor<Nd + 3, ComplexD> psi)
      {
	int N = Nd - 1; // Only the spatial dimensions

	// chi = -2*N*psi
	psi.scale(-2 * N).copyTo(chi);

	for (int mu = 0; mu < N; ++mu)
	{
	  displace(u, psi, first_tslice, mu + 1, Action::AddTo, chi);
	  displace(u, psi, first_tslice, -(mu + 1), Action::AddTo, chi);
	}
      }

      // Auxiliary structure passed to PRIMME's matvec

      struct OperatorAux {
	const Operator<Nd + 2, ComplexD> op; // Operator in cxyztX
	const DeviceHost primme_dev;	     // where primme allocations are
      };

      // Wrapper for PRIMME of `LaplacianOperator`
      /// \param x: pointer to input vector
      /// \param ldx: leading dimension for `x`
      /// \param y: pointer to output vector
      /// \param ldy: leading dimension for `y`
      /// \param blockSize: number of input/output vectors
      /// \param ierr: output error state (zero means ok)

      extern "C" inline void primmeMatvec(void* x, PRIMME_INT* ldx, void* y, PRIMME_INT* ldy,
					  int* blockSize, primme_params* primme, int* ierr)
      {
	*ierr = -1;
	try
	{
	  // The implementation assumes that ldx and ldy is nLocal
	  if (*blockSize > 1 && (*ldx != primme->nLocal || *ldy != primme->nLocal))
	    throw std::runtime_error("We cannot play with the leading dimensions");

	  OperatorAux& opaux = *(OperatorAux*)primme->matrix;
	  const std::string order(opaux.op.d.order + std::string("n"));
	  Coor<Nd + 3> size = latticeSize<Nd + 3>(order, {{'n', *blockSize}, {'t', 1}});
	  Tensor<Nd + 3, ComplexD> tx(order, size, opaux.primme_dev, OnEveryone, (ComplexD*)x);
	  Tensor<Nd + 3, ComplexD> ty(order, size, opaux.primme_dev, OnEveryone, (ComplexD*)y);
	  assert(tx.getLocal().volume() == primme->nLocal * (*blockSize));
	  opaux.op(tx, ty);
	  assert(ty.allocation->pending_operations.size() == 0);
	  // Make sure cublas/hipblas handle operates on legacy stream for primme
	  superbblas::detail::gpuBlasCheck(SUPERBBLAS_GPUBLAS_SYMBOL(SetStream)(
	    *(superbblas::detail::GpuBlasHandle*)primme->queue, 0));
	  *ierr = 0;
	} catch (...)
	{
	}
      }

      /// Wrapper for PRIMME of a global sum for double
      /// \param sendBuf: pointer to input vector
      /// \param recvBuf: pointer to output vector
      /// \param count: number of elements in the input/output vector
      /// \param primme: pointer to the current primme_params
      /// \param ierr: output error state (zero means ok)

      extern "C" inline void primmeGlobalSum(void* sendBuf, void* recvBuf, int* count,
					     primme_params* primme, int* ierr)
      {
	if (sendBuf == recvBuf)
	{
	  *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPI_DOUBLE, MPI_SUM,
				MPI_COMM_WORLD) != MPI_SUCCESS;
	}
	else
	{
	  *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) !=
		  MPI_SUCCESS;
	}
      }

      /// Compute the eigenpairs of the laplacian operator on the spatial dimensions using PRIMME
      /// \param u: Gauge field
      /// \param from_tslice: index of the first t-slice to compute the eigenvectors from
      /// \param n_tslices: number of tslices to compute
      /// \param n_colorvecs: number of eigenpairs to compute
      /// \param order_: order of the output tensor for the eigenvectors
      /// \return: a pair of the eigenvectors and the eigenvalues

      inline std::pair<Tensor<Nd + 3, ComplexD>, std::vector<std::vector<double>>>
      computeColorvecs(const multi1d<LatticeColorMatrix>& u, int from_tslice, int n_tslices,
		       int n_colorvecs, const Maybe<const std::string>& order_ = none)
      {
	const std::string order = order_.getSome("cxyztXn");
	detail::check_order_contains(order, "cxyztXn");
	Tensor<Nd + 3, ComplexD> all_evecs(
	  order, latticeSize<Nd + 3>(order, {{'n', n_colorvecs}, {'t', n_tslices}}),
	  OnDefaultDevice, OnEveryone);
	std::vector<std::vector<double>> all_evals;

	for (Index t = 0; t < n_tslices; ++t)
	{
	  // Make a copy of the time-slicing of u[d] also supporting left and right
	  std::vector<Tensor<Nd + 3, ComplexD>> ut(Nd);
	  for (unsigned int d = 0; d < Nd - 1; d++)
	  {
	    ut[d] = asTensorView(u[d])
		      .kvslice_from_size({{'t', from_tslice + t}}, {{'t', 1}})
		      .toComplex()
		      .clone()
		      .template make_sure<ComplexD>("ijxyztX");
	  }

	  // If the 3D laplacian operator is big enough, run it on device
	  DeviceHost primme_dev = OnHost;
#    if defined(SUPERBBLAS_USE_GPU)
	  primme_dev = OnDefaultDevice;
#    endif

	  // Create an efficient representation of the laplacian operator
	  std::string order("cxyztX");
	  auto eg = Tensor<Nd + 2, ComplexD>(order, latticeSize<Nd + 2>(order, {{'t', 1}}),
					     primme_dev, OnEveryone)
		      .make_eg();
	  OperatorLayout op_layout =
	    ((from_tslice + t) % 2 == 0 ? XEvenOddLayout : XEvenOddLayoutZeroOdd);
	  auto laplacianOp = Chroma::SB::detail::cloneOperator(
	    Operator<Nd + 2, ComplexD>{
	      [&](Tensor<Nd + 3, ComplexD> x, Tensor<Nd + 3, ComplexD> y) {
		LaplacianOperator(ut, from_tslice + t, y, x);
	      },
	      eg, eg, nullptr, "", op_layout, op_layout,
	      detail::getNeighbors(eg.kvdim(), 1 /* near-neighbors links only */, op_layout),
	      ColumnMajor, false /* no kronecker op */},
	    ColumnMajor, RowMajor, Chroma::SB::detail::ConsiderBlockingDense, "laplacian");

	  // Create an auxiliary struct for the PRIMME's matvec
	  // NOTE: Please keep 'n' as the slowest index; the rows of vectors taken by PRIMME's matvec has dimensions 'cxyztX',
	  // and 'n' is the dimension for the columns.
	  OperatorAux opaux{laplacianOp, primme_dev};

	  // Make a bigger structure holding
	  primme_params primme;
	  primme_initialize(&primme);

	  // Get the global and local size of evec
	  std::size_t n = eg.volume();
	  std::size_t nLocal = eg.localVolume();

	  if (n_colorvecs > n)
	  {
	    std::cerr << "ERROR: the rank of the distillation basis shouldn't be larger than the "
			 "spatial dimensions"
		      << std::endl;
	    exit(1);
	  }

	  // Primme solver setup
	  primme.numEvals = n_colorvecs;
	  primme.printLevel = 0;
	  primme.n = n;
	  primme.eps = 1e-9;
	  primme.target = primme_largest;

	  // Set parallel settings
	  primme.nLocal = nLocal;
	  primme.numProcs = QDP::Layout::numNodes();
	  primme.procID = QDP::Layout::nodeNumber();
	  primme.globalSumReal = primmeGlobalSum;

	  // No preconditioner for my matrix
	  primme.matrixMatvec = primmeMatvec;
	  primme.matrix = &opaux;

	  // Set block size
	  if (n > 128)
	  {
	    primme.maxBasisSize = 64;
	    primme.maxBlockSize = 4;
	  }
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
	  const std::string evecs_order(eg.order + std::string("n"));
	  Tensor<Nd + 3, ComplexD> evecs(
	    evecs_order, latticeSize<Nd + 3>(evecs_order, {{'n', primme.numEvals}, {'t', 1}}),
	    primme_dev, OnEveryone);
	  assert(evecs.localVolume() == primme.nLocal * primme.numEvals);
#    if defined(SUPERBBLAS_USE_GPU)
	  superbblas::detail::GpuBlasHandle gpublas_handle =
	    superbblas::detail::getGpuBlasHandle(evecs.ctx().toGpu(0));
	  primme.queue = &gpublas_handle;
	  // Make sure cublas/hipblas handle operates on legacy stream for primme
	  superbblas::detail::gpuBlasCheck(SUPERBBLAS_GPUBLAS_SYMBOL(SetStream)(gpublas_handle, 0));
#    endif


	  // Call primme
	  int ret;
#    if defined(SUPERBBLAS_USE_GPU)
	  if (primme_dev == OnDefaultDevice)
	  {
	    ret = cublas_zprimme(evals.data(), evecs.data(), rnorms.data(), &primme);
	  }
	  else
#    endif
	  {
	    ret = zprimme(evals.data(), evecs.data(), rnorms.data(), &primme);
	  }

	  if (primme.procID == 0)
	  {
	    fprintf(stdout, " %d eigenpairs converged for tslice %d\n", primme.initSize,
		    from_tslice + t);
	    fprintf(stdout, "Tolerance : %-22.15E\n", primme.aNorm * primme.eps);
	    fprintf(stdout, "Iterations: %-d\n", (int)primme.stats.numOuterIterations);
	    fprintf(stdout, "Restarts  : %-d\n", (int)primme.stats.numRestarts);
	    fprintf(stdout, "Matvecs   : %-d\n", (int)primme.stats.numMatvecs);
	    fprintf(stdout, "Preconds  : %-d\n", (int)primme.stats.numPreconds);
	    fprintf(stdout, "T. ortho  : %g\n", primme.stats.timeOrtho);
	    fprintf(stdout, "T. matvec : %g\n", primme.stats.timeMatvec);
	    fprintf(stdout, "Total time: %g\n", primme.stats.elapsedTime);
	  }

	  if (ret != 0)
	  {
	    QDPIO::cerr << "Error: primme returned with nonzero exit status\n";
	    QDP_abort(1);
	  }

	  // Cleanup
	  primme_free(&primme);

	  // Check the residuals, |laplacian*v-lambda*v|_2<=|laplacian|*tol
	  if (evals.size() > 0)
	  {
	    auto r = evecs.like_this();
	    LaplacianOperator(ut, from_tslice + t, r, evecs);
	    std::vector<std::complex<double>> evals_cmpl(evals.begin(), evals.end());
	    contract<Nd + 3, Nd + 3, 1, ComplexD>(
	      evecs, asTensorView(evals_cmpl).rename_dims({{'i', 'n'}}).scale(-1), "", AddTo, r);
	    auto rnorm = norm<1>(r, "n");
	    for (int i = 0, vol = rnorm.volume(); i < vol; ++i)
	    {
	      if (rnorm.get({{i}}) > primme.stats.estimateLargestSVal * primme.eps * 10)
	      {
		QDPIO::cerr << "Error: primme returned eigenpairs with too much error\n";
		QDP_abort(1);
	      }
	    }
	  }

	  // Copy evecs into all_evecs
	  evecs.copyTo(all_evecs.kvslice_from_size({{'t', t}}, {{'t', 1}}));
	  all_evals.push_back(evals);
	}

	return {all_evecs, all_evals};
      }
#  else	 // BUILD_PRIMME
      inline std::pair<Tensor<Nd + 3, ComplexD>, std::vector<std::vector<double>>>
      computeColorvecs(const multi1d<LatticeColorMatrix>& u, int from_tslice, int n_tslices,
		       int n_colorvecs, const Maybe<const std::string>& order_ = none)
      {
	(void)u;
	(void)from_tslice;
	(void)n_tslices;
	(void)n_colorvecs;
	(void)order_;
	throw std::runtime_error("Functionality isn't available without compiling with PRIMME");
      }
#  endif // BUILD_PRIMME

      /// Read colorvecs from a FILEDB
      /// \param eigen_source: database handle
      /// \param decay_dir: something we assume is always three
      /// \param from_tslice: first tslice to read
      /// \param n_tslices: number of tslices to read
      /// \param n_colorvecs: number of eigenpairs to read
      /// \param order_: order of the output tensor for the eigenvectors
      /// \return: a tensor containing the eigenvectors

      template <typename COMPLEX = ComplexF>
      Tensor<Nd + 3, COMPLEX> getColorvecs(MODS_t& eigen_source, int decay_dir, int from_tslice,
					   int n_tslices, int n_colorvecs,
					   const Maybe<const std::string>& order_ = none,
					   DeviceHost dev = OnDefaultDevice)
      {
	const std::string order = order_.getSome("cxyztXn");
	detail::check_order_contains(order, "cxyztXn");

	from_tslice = normalize_coor(from_tslice, Layout::lattSize()[decay_dir]);

	// Allocate tensor to return
	Tensor<Nd + 3, COMPLEX> r(
	  order, latticeSize<Nd + 3>(order, {{'t', n_tslices}, {'n', n_colorvecs}}), dev);

	// Allocate a single time slice colorvec in natural ordering, as colorvec are stored
	Tensor<Nd, ComplexF> tnat("cxyz", latticeSize<Nd>("cxyz", {{'x', Layout::lattSize()[0]}}),
				  OnHost, OnMaster);

	// Allocate a single time slice colorvec in case of using RB ordering
	Tensor<Nd + 1, ComplexF> trb("cxyzX", latticeSize<Nd + 1>("cxyzX"), OnHost, OnMaster);

	// Allocate all colorvecs for the same time-slice
	Tensor<Nd + 2, ComplexF> t("cxyzXn", latticeSize<Nd + 2>("cxyzXn", {{'n', n_colorvecs}}),
				   OnHost, OnMaster);

	const int Nt = Layout::lattSize()[decay_dir];
	for (int t_slice = from_tslice, i_slice = 0; i_slice < n_tslices;
	     ++i_slice, t_slice = (t_slice + 1) % Nt)
	{
	  // Compute the permutation from natural ordering to red-black
	  std::vector<Index> perm = ns_getColorvecs::getPermFromNatToRB(t_slice);

	  for (int colorvec = 0; colorvec < n_colorvecs; ++colorvec)
	  {
	    // Read a single time-slice and colorvec
	    KeyTimeSliceColorVec_t key(t_slice, colorvec);
	    if (!eigen_source.exist(key))
	      throw std::runtime_error(
		"no colorvec exists with key t_slice= " + std::to_string(t_slice) +
		" colorvec= " + std::to_string(colorvec));
	    eigen_source.get(key, tnat);

	    // Correct ordering
	    ns_getColorvecs::toRB(perm, tnat, trb);

	    // t[n=colorvec] = trb
	    trb.copyTo(t.kvslice_from_size({{'n', colorvec}}, {{'n', 1}}));
	  }

	  // r[t=i_slice] = t, distribute the tensor from master to the rest of the nodes
	  t.copyTo(r.kvslice_from_size({{'t', i_slice}}));
	}

	return r;
      }

      /// Read colorvecs from a S3T file
      /// \param s3t: database handle
      /// \param u: gauge field
      /// \param decay_dir: something we assume is always three
      /// \param from_tslice: first tslice to read
      /// \param n_tslices: number of tslices to read
      /// \param n_colorvecs: number of eigenpairs to read
      /// \param order_: order of the output tensor for the eigenvectors
      /// \return: a tensor containing the eigenvectors

      template <typename COMPLEX = ComplexF>
      Tensor<Nd + 3, COMPLEX>
      getColorvecs(StorageTensor<Nd + 2, ComplexD> s3t, const multi1d<LatticeColorMatrix>& u,
		   int decay_dir, int from_tslice, int n_tslices, int n_colorvecs,
		   const Maybe<const std::string>& order_ = none, DeviceHost dev = OnDefaultDevice)
      {
	const std::string order = order_.getSome("cxyztXn");
	detail::check_order_contains(order, "cxyztXn");

	from_tslice = normalize_coor(from_tslice, Layout::lattSize()[decay_dir]);

	// Read the metadata and check that the file stores colorvecs and from a lattice of the same size
	std::istringstream is(s3t.metadata);
	XMLReader xml_buf(is);
	bool write_fingerprint = false;
	read(xml_buf, "/MODMetaData/fingerprint", write_fingerprint);
	GroupXML_t link_smear =
	  readXMLGroup(xml_buf, "/MODMetaData/LinkSmearing", "LinkSmearingType");

	// Smear the gauge field if needed
	multi1d<LatticeColorMatrix> u_smr = u;
	try
	{
	  std::istringstream xml_l(link_smear.xml);
	  XMLReader linktop(xml_l);
	  Handle<LinkSmearing> linkSmearing(TheLinkSmearingFactory::Instance().createObject(
	    link_smear.id, linktop, "/LinkSmearing"));
	  (*linkSmearing)(u_smr);
	} catch (const std::string& e)
	{
	  QDPIO::cerr << ": Caught Exception link smearing: " << e << std::endl;
	  QDP_abort(1);
	} catch (...)
	{
	  QDPIO::cerr << ": Caught unexpected exception" << std::endl;
	  QDP_abort(1);
	}

	// Allocate tensor with the content of s3t
	Tensor<Nd + 3, ComplexD> colorvecs_s3t(
	  order, latticeSize<Nd + 3>(order, {{'t', n_tslices}, {'n', n_colorvecs}}), dev);

	// Allocate a single time slice colorvec in natural ordering, as colorvec are stored
	Tensor<Nd, ComplexD> tnat("cxyz", latticeSize<Nd>("cxyz", {{'x', Layout::lattSize()[0]}}),
				  OnHost, OnMaster);

	// Allocate a single time slice colorvec in case of using RB ordering
	Tensor<Nd + 1, ComplexD> trb("cxyzX", latticeSize<Nd + 1>("cxyzX"), OnHost, OnMaster);

	const int Nt = Layout::lattSize()[decay_dir];
	for (int t_slice = from_tslice, i_slice = 0; i_slice < n_tslices;
	     ++i_slice, t_slice = (t_slice + 1) % Nt)
	{
	  // Compute the permutation from natural ordering to red-black
	  std::vector<Index> perm = ns_getColorvecs::getPermFromNatToRB(t_slice);

	  for (int colorvec = 0; colorvec < n_colorvecs; ++colorvec)
	  {
	    // Read a single time-slice and colorvec
	    tnat.set_zero();
	    s3t.kvslice_from_size({{'t', t_slice}, {'n', colorvec}}, {{'t', 1}, {'n', 1}})
	      .copyTo(tnat);

	    // Correct ordering
	    ns_getColorvecs::toRB(perm, tnat, trb);

	    // colorvecs_s3t[t=i_slice,n=colorvec] = trb
	    trb.copyTo(colorvecs_s3t.kvslice_from_size({{'t', i_slice}, {'n', colorvec}}));
	  }
	}

	// Compute the 2-norm of colorvecs_s3t and check that no vector is null
	auto colorvecs_s3t_norms = norm<2>(colorvecs_s3t, "nt");
	for (int t = 0; t < n_tslices; ++t)
	  for (int n = 0; n < n_colorvecs; ++n)
	    if (colorvecs_s3t_norms.get({n, t}) == 0)
	      throw std::runtime_error(
		"no colorvec exists with key t_slice= " + std::to_string(t + from_tslice) +
		" colorvec= " + std::to_string(n));

	if (write_fingerprint)
	{
	  // Compute the colorvecs
	  auto colorvecs =
	    ns_getColorvecs::computeColorvecs(u_smr, from_tslice, n_tslices, n_colorvecs, order_)
	      .first;

	  // We need to phase the individual eigenvectors so that the have the same phase as the
	  // s3t's colorvecs. That is, we need to apply a phase phi[i] to each eigenvector so that
	  //
	  //    colorvecs_s3t[i] = colorvecs[i] * phi[i].
	  //
	  // We have a subset of the s3t's colorvecs, so we restrict the above equation to that:
	  //
	  //    colorvecs_s3t[i]^\dagger * colorvecs_s3t[i] = colorvecs_s3t[i]^\dagger * colorvecs[i] * phi[i].
	  //
	  // Therefore, phi[i] = (colorvecs_s3t[i]^\dagger * colorvecs_s3t[i]) / (colorvecs_s3t[i]^\dagger * colorvecs[i])

	  auto ip = contract<2>(colorvecs_s3t.conj(), colorvecs,
				detail::remove_dimensions(colorvecs.order, "nt"), OnHost,
				OnEveryoneReplicated)
		      .reorder("nt");

	  auto phi = ip.like_this();
	  for (int t = 0; t < n_tslices; ++t)
	  {
	    for (int n = 0; n < n_colorvecs; ++n)
	    {
	      auto cv_norm = colorvecs_s3t_norms.get({n, t});
	      auto phi_i = cv_norm * cv_norm / ip.get({n, t});
	      if (std::fabs(std::fabs(phi_i) - 1) > 1e-4)
		throw std::runtime_error(
		  "The colorvec fingerprint does not correspond to current gates field");
	      phi.set({n, t}, phi_i);
	    }
	  }

	  // Apply the phase of the colorvecs in s3t to the computed colorvecs
	  colorvecs_s3t.contract(colorvecs, {}, NotConjugate, phi, {}, NotConjugate);
	}

	return colorvecs_s3t.make_sure<COMPLEX>();
      }
    }

    /// Read colorvecs from either a FILEDB or S3T file
    /// \param colorvec_files: filenames
    /// \return: a handle

    inline ColorvecsStorage openColorvecStorage(const std::vector<std::string>& colorvec_files)
    {
      ColorvecsStorage sto{}; // returned object

      std::string metadata; // the metadata content of the file

      // Try to open the file as a s3t database
      try
      {
	if (colorvec_files.size() == 1)
	  sto.s3t = StorageTensor<Nd + 2, ComplexD>(colorvec_files[0], true, "/MODMetaData/order");
	metadata = sto.s3t.metadata;
      } catch (...)
      {
      }

      // Try to open the files as a MOD database
      if (!sto.s3t)
      {
	sto.mod = std::make_shared<MODS_t>();
	sto.mod->setDebug(0);

	try
	{
	  // Open
	  sto.mod->open(colorvec_files);
	  sto.mod->getUserdata(metadata);
	} catch (std::bad_cast)
	{
	  QDPIO::cerr << ": caught dynamic cast error" << std::endl;
	  QDP_abort(1);
	} catch (const std::string& e)
	{
	  QDPIO::cerr << ": error extracting source_header: " << e << std::endl;
	  QDP_abort(1);
	} catch (const char* e)
	{
	  QDPIO::cerr << ": Caught some char* exception:" << std::endl;
	  QDPIO::cerr << e << std::endl;
	  QDP_abort(1);
	}
      }

      // Check that the file stores colorvecs and is from a lattice of the same size

      std::istringstream is(metadata);
      XMLReader xml_buf(is);

      std::string id;
      read(xml_buf, "/MODMetaData/id", id);
      if (id != "eigenVecsTimeSlice")
      {
	std::stringstream ss;
	ss << "The file `" << colorvec_files[0] << "' does not contain colorvecs";
	throw std::runtime_error(ss.str());
      }

      multi1d<int> spatialLayout(3);
      read(xml_buf, "/MODMetaData/lattSize", spatialLayout);
      if (spatialLayout[0] != Layout::lattSize()[0] || spatialLayout[1] != Layout::lattSize()[1] ||
	  spatialLayout[2] != Layout::lattSize()[2])
      {
	std::stringstream ss;
	ss << "The spatial dimensions of the colorvecs in `" << colorvec_files[0]
	   << "' do not much the current lattice";
	throw std::runtime_error(ss.str());
      }

      return sto;
    }

    /// Close a colorvec storage
    /// \param sto: colorvec storage handle

    inline void closeColorvecStorage(ColorvecsStorage& sto)
    {
      if (!sto.s3t)
      {
	sto.s3t.release();
      }
      else if (sto.mod)
      {
	sto.mod->close();
	sto.mod.reset();
      }
    }

    /// Phase colorvecs
    /// \param colorvecs: tensor with the colorvecs
    /// \param from_tslice: first tslice of the tensor
    /// \param phase: apply a phase to the eigenvectors
    /// \return: a tensor containing the eigenvectors phased

    template <typename COMPLEX>
    Tensor<Nd + 3, COMPLEX> phaseColorvecs(Tensor<Nd + 3, COMPLEX> colorvecs, int from_tslice,
					   Coor<Nd - 1> phase = {{}})
    {
      // Phase colorvecs if phase != (0,0,0)
      if (phase == Coor<Nd - 1>{{}})
	return colorvecs;

      Tensor<Nd + 1, COMPLEX> tphase = ns_getColorvecs::getPhase<COMPLEX>(
	phase, from_tslice, colorvecs.kvdim()['t'], colorvecs.getDev());
      Tensor<Nd + 3, COMPLEX> r = colorvecs.like_this();
      r.contract(colorvecs, {}, NotConjugate, tphase, {}, NotConjugate);
      return r;
    }

    /// Read colorvecs from a handle returned by `openColorvecStorage`
    /// \param sto: database handle
    /// \param u: gauge field
    /// \param decay_dir: something we assume is always three
    /// \param from_tslice: first tslice to read
    /// \param n_tslices: number of tslices to read
    /// \param n_colorvecs: number of eigenpairs to read
    /// \param order: order of the output tensor for the eigenvectors
    /// \param phase: apply a phase to the eigenvectors
    /// \return: a tensor containing the eigenvectors

    template <typename COMPLEX = ComplexF>
    Tensor<Nd + 3, COMPLEX>
    getColorvecs(const ColorvecsStorage& sto, const multi1d<LatticeColorMatrix>& u, int decay_dir,
		 int from_tslice, int n_tslices, int n_colorvecs,
		 const Maybe<const std::string>& order = none, Coor<Nd - 1> phase = {{}},
		 DeviceHost dev = OnDefaultDevice)
    {
      StopWatch sw;
      sw.reset();
      sw.start();

      if (decay_dir != 3)
	throw std::runtime_error("Only support for decay_dir being the temporal dimension");

      // Read the colorvecs with the proper function
      Tensor<Nd + 3, COMPLEX> r;
      if (sto.s3t)
	r = ns_getColorvecs::getColorvecs<COMPLEX>(sto.s3t, u, decay_dir, from_tslice, n_tslices,
						   n_colorvecs, order, dev);
      else if (sto.mod)
	r = ns_getColorvecs::getColorvecs<COMPLEX>(*sto.mod, decay_dir, from_tslice, n_tslices,
						   n_colorvecs, order, dev);

      // Phase colorvecs
      r = phaseColorvecs(r, from_tslice, phase);

      sw.stop();
      QDPIO::cout << "Time to read " << n_colorvecs << " colorvecs from " << n_tslices
		  << " time slices: " << sw.getTimeInSeconds() << " secs" << std::endl;

      return r;
    }

    /// Compute and store colorvecs
    /// \param colorvec_file: file to store the colorvecs
    /// \param link_smear: smearing gauge field options before building the laplacian
    /// \param u: gauge field
    /// \param from_tslice: first tslice to read
    /// \param n_tslices: number of tslices to read
    /// \param n_colorvecs: number of eigenpairs to read
    /// \param use_s3t_storage: if true S3T is used, otherwise FILEDB
    /// \param fingerprint: whether to store only a few sites of each colorvecs
    /// \param phase: apply a phase to the eigenvectors
    /// \param colorvec_file_src: if given, read the colorvecs from that file and if they
    ///        match the computed ones, they are the ones stored; this guarantee that the
    ///        that given smearing options were used to generate the colorvecs in `colorvec_file_src`

    inline void
    createColorvecStorage(const std::string& colorvec_file, GroupXML_t link_smear,
			  const multi1d<LatticeColorMatrix>& u, int from_tslice, int n_tslices,
			  int n_colorvecs, bool use_s3t_storage = false, bool fingerprint = false,
			  Coor<Nd - 1> phase = {{}},
			  const Maybe<std::vector<std::string>>& colorvec_file_src = none)
    {
      // Check input
      const int Nt = Layout::lattSize()[3];
      if (from_tslice < 0)
	throw std::runtime_error("The first t-slice to compute colorvecs is negative!");
      if (n_tslices < 0 || n_tslices > Nt)
	throw std::runtime_error(" The number of t-slices to compute colorvecs is negative or "
				 "greater than the t dimension of the lattice");

      // Smear the gauge field if needed
      multi1d<LatticeColorMatrix> u_smr = u;
      try
      {
	std::istringstream xml_l(link_smear.xml);
	XMLReader linktop(xml_l);
	Handle<LinkSmearing> linkSmearing(
	  TheLinkSmearingFactory::Instance().createObject(link_smear.id, linktop, link_smear.path));
	(*linkSmearing)(u_smr);
      } catch (const std::string& e)
      {
	QDPIO::cerr << ": Caught Exception link smearing: " << e << std::endl;
	QDP_abort(1);
      } catch (...)
      {
	QDPIO::cerr << ": Caught unexpected exception" << std::endl;
	QDP_abort(1);
      }

      // Some tasks read the eigenvalues from metadata but they not used; so we are going to give fake values
      multi1d<multi1d<double>> evals(n_colorvecs);
      for (int i = 0; i < n_colorvecs; ++i)
      {
	evals[i].resize(n_tslices);
	for (int t = 0; t < n_tslices; ++t)
	  evals[i][t] = 0;
      }

      // Open the DB and write metada
      MOD_t mod;
      StorageTensor<Nd + 2, ComplexD> sto;
      Coor<3> fingerprint_dim{{}};

      if (!use_s3t_storage)
      {
	XMLBufferWriter file_xml;

	push(file_xml, "MODMetaData");
	write(file_xml, "id", "eigenVecsTimeSlice");
	multi1d<int> spatialLayout(3);
	spatialLayout[0] = Layout::lattSize()[0];
	spatialLayout[1] = Layout::lattSize()[1];
	spatialLayout[2] = Layout::lattSize()[2];
	write(file_xml, "lattSize", spatialLayout);
	write(file_xml, "decay_dir", 3);
	write(file_xml, "num_vecs", n_colorvecs);
	write(file_xml, "Weights", evals);
	file_xml << link_smear.xml;
	pop(file_xml);

	mod.setDebug(0);

	mod.insertUserdata(file_xml.str());
	mod.open(colorvec_file, std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
      }
      else
      {
	std::string sto_order = "cxyztn"; // order for storing the colorvecs

	// If fingerprint, we store only the support of the colorvecs on a subset of the lattice;
	// compute the size of that subset
	for (int i = 0; i < 3; ++i)
	  fingerprint_dim[i] = std::min(4, Layout::lattSize()[i]);

	// Prepare metadata
	XMLBufferWriter file_xml;

	push(file_xml, "MODMetaData");
	write(file_xml, "id", "eigenVecsTimeSlice");
	multi1d<int> spatialLayout(3);
	spatialLayout[0] = Layout::lattSize()[0];
	spatialLayout[1] = Layout::lattSize()[1];
	spatialLayout[2] = Layout::lattSize()[2];
	write(file_xml, "lattSize", spatialLayout);
	write(file_xml, "decay_dir", 3);
	write(file_xml, "num_vecs", n_colorvecs);
	write(file_xml, "Weights", evals);
	write(file_xml, "order", sto_order);
	write(file_xml, "fingerprint", fingerprint);
	if (fingerprint)
	{
	  spatialLayout[0] = fingerprint_dim[0];
	  spatialLayout[1] = fingerprint_dim[1];
	  spatialLayout[2] = fingerprint_dim[2];
	  write(file_xml, "fingerprint_lattice", spatialLayout);
	}
	file_xml << link_smear.xml;
	pop(file_xml);

	// NOTE: file_xml has nonzero value only at the master node; so do a broadcast

	sto = StorageTensor<Nd + 2, ComplexD>(
	  colorvec_file, broadcast(file_xml.str()), sto_order,
	  latticeSize<Nd + 2>(sto_order, {{'n', n_colorvecs}, {'x', Layout::lattSize()[0]}}),
	  Sparse, superbblas::BlockChecksum);
      }

      // Open colorvec_file_src
      ColorvecsStorage colorvecsSto;
      if (colorvec_file_src.getSome({}).size() > 0)
	colorvecsSto = openColorvecStorage(colorvec_file_src.getSome());

      for (int i_tslice = 0; i_tslice < n_tslices; ++i_tslice, from_tslice = (from_tslice + 1) % Nt)
      {
	// Compute colorvecs
	std::string order = "cxyzXtn";
	auto colorvecs_and_evals =
	  ns_getColorvecs::computeColorvecs(u_smr, from_tslice, 1, n_colorvecs, order);
	auto colorvecs = colorvecs_and_evals.first;

	// Read the eigenvectors from another source if indicated
	if (colorvec_file_src.getSome({}).size() > 0)
	{
	  auto colorvecs_src =
	    getColorvecs<ComplexD>(colorvecsSto, u, 3, from_tslice, 1, n_colorvecs);

	  Tensor<2, ComplexD> ip("nt", Coor<2>{n_colorvecs, 1}, OnHost, OnEveryoneReplicated);
	  ip.contract(colorvecs, {}, Conjugate, colorvecs_src, {}, NotConjugate);
	  for (int n = 0; n < n_colorvecs; ++n)
	    if (std::fabs(std::fabs(ip.get({n, 0})) - 1) > 1e-4)
	      throw std::runtime_error(
		"The given colorvec does not correspond to current gates field and smearing");
	  colorvecs = colorvecs_src;
	}

	// Phase colorvecs
	colorvecs = phaseColorvecs(colorvecs, from_tslice, phase);

	// Compute the permutation from natural ordering to red-black
	std::vector<Index> perm = ns_getColorvecs::getPermFromNatToRB(from_tslice);

	// Store the colorvecs in natural order (not in red-black ordering)
	if (!use_s3t_storage)
	{
	  // Allocate a single time slice colorvec in natural ordering, as colorvec are stored
	  Tensor<Nd, ComplexF> tnat("cxyz", latticeSize<Nd>("cxyz", {{'x', Layout::lattSize()[0]}}),
				    OnHost, OnMaster);

	  // Allocate a single time slice colorvec in case of using RB ordering
	  Tensor<Nd + 1, ComplexF> trb("cxyzX", latticeSize<Nd + 1>("cxyzX"), OnHost, OnMaster);

	  for (int n = 0; n < n_colorvecs; ++n)
	  {
	    KeyTimeSliceColorVec_t time_key;
	    time_key.t_slice = from_tslice;
	    time_key.colorvec = n;
	    colorvecs.kvslice_from_size({{'t', 0}, {'n', n}}, {{'t', 1}, {'n', 1}}).copyTo(trb);
	    ns_getColorvecs::toNat(perm, trb, tnat);
	    mod.insert(time_key, tnat);
	  }
	}
	else
	{
	  // Allocate a single time slice colorvec in natural ordering, as colorvec are stored
	  Tensor<Nd, ComplexD> tnat("cxyz", latticeSize<Nd>("cxyz", {{'x', Layout::lattSize()[0]}}),
				    OnHost, OnMaster);

	  // Allocate a single time slice colorvec in case of using RB ordering
	  Tensor<Nd + 1, ComplexD> trb("cxyzX", latticeSize<Nd + 1>("cxyzX"), OnHost, OnMaster);

	  std::map<char, int> colorvec_size{};
	  if (fingerprint)
	    colorvec_size = std::map<char, int>{
	      {'x', fingerprint_dim[0]}, {'y', fingerprint_dim[1]}, {'z', fingerprint_dim[2]}};

	  for (int n = 0; n < n_colorvecs; ++n)
	  {
	    colorvecs.kvslice_from_size({{'t', 0}, {'n', n}}, {{'t', 1}, {'n', 1}}).copyTo(trb);
	    ns_getColorvecs::toNat(perm, trb, tnat);
	    sto.kvslice_from_size({{'t', from_tslice}, {'n', n}}, {{'t', 1}, {'n', 1}})
	      .copyFrom(tnat.kvslice_from_size({}, colorvec_size));
	  }
	}
      }

      if (!use_s3t_storage)
	mod.close();
    }

    //
    // High-level chroma operations
    //

    namespace detail
    {
      /// Path Node
      struct PathNode {
	std::map<int, PathNode> p; ///< following nodes
	int disp_index;		   ///< if >= 0, the index in the displacement list
      };

      /// Return the directions that are going to be use and the maximum number of displacements keep in memory
      inline void get_tree_mem_stats(const PathNode& disps, std::array<bool, Nd>& dirs,
				     unsigned int& max_rhs)
      {
	unsigned int max_rhs_sub = 0;
	for (const auto it : disps.p)
	{
	  unsigned int max_rhs_sub_it = 0;
	  get_tree_mem_stats(it.second, dirs, max_rhs_sub_it);
	  max_rhs_sub = std::max(max_rhs_sub, max_rhs_sub_it);

	  if (std::abs(it.first) <= Nd)
	    dirs[std::abs(it.first) - 1] = true;
	}

	if (disps.p.size() == 0)
	{
	  max_rhs = 0;
	}
	else if (disps.p.size() == 1)
	{
	  max_rhs = std::max(1u, max_rhs_sub);
	}
	else
	{
	  max_rhs = 1 + max_rhs_sub;
	}
      }

      const int path_separator = Nd + 1;

      /// Return the tree representing all paths
      /// \param paths: list of displacements
      /// \param allow_separator: allow a special direction Nd+1

      inline PathNode get_tree(const std::vector<std::vector<int>>& paths,
			       bool allow_separator = false)
      {
	PathNode r{{}, -1};
	int path_index = 0;
	for (const std::vector<int>& path : paths)
	{
	  PathNode* n = &r;
	  for (int d : path)
	  {
	    if (d == 0 || (std::abs(d) > Nd && (!allow_separator || d != path_separator)))
	      throw std::runtime_error("Invalid direction: " + std::to_string(d));

	    auto it = n->p.find(d);
	    if (it != n->p.end())
	      n = &it->second;
	    else
	    {
	      n->p[d] = PathNode{{}, -1};
	      n = &n->p[d];
	    }
	  }
	  if (n->disp_index < 0)
	    n->disp_index = path_index++;
	}
	return r;
      }

    }

    namespace ns_doMomGammaDisp_contractions
    {
      using namespace detail;

      /// Contract two LatticeFermion with different momenta, gammas, and displacements.
      /// \param leftconj: left lattice fermion tensor, cSxyzXN
      /// \param right: right lattice fermion tensor, csxyzXn
      /// \param disps: tree of displacements/derivatives
      /// \param deriv: if true, do left-right nabla derivatives
      /// \param gammas: tensor with spins, QSg
      /// \param moms: list of momenta
      /// \param max_rhs: maximum number of vectors hold in memory
      /// \param r tensor holding the contractions, sqnNmgd where
      ///        q and N (s and n) are the spin and vector from left (right) vectors, m is the momentum
      ///        index, g is the gamma index, and d is the displacement index
      /// \param: disp_indices: dictionary that map each `d` index in r displacement index.

      template <typename COMPLEX, std::size_t Nleft, std::size_t Nright, std::size_t Nout>
      void doMomGammaDisp_contractions(const std::vector<Tensor<Nd + 3, Complex>>& u,
				       const Tensor<Nleft, COMPLEX> leftconj,
				       Tensor<Nright, COMPLEX> right, Index first_tslice,
				       const PathNode& disps, bool deriv, Tensor<3, COMPLEX> gammas,
				       const std::vector<Coor<Nd - 1>>& moms, int max_rhs,
				       Tensor<Nout, COMPLEX> r, std::vector<int>& disp_indices)
      {
	max_rhs = std::max(1, max_rhs);

	if (disps.disp_index >= 0)
	{
	  detail::log(1, "contracting for disp_index=" + std::to_string(disps.disp_index));
	  // Contract the spatial components and the color of the leftconj and right tensors
	  Tensor<Nout, COMPLEX> aux =
	    r.template like_this<Nout, COMPLEX>("mNQqnSst%", '%', "gd", {{'S', Ns}, {'Q', Ns}});
	  aux.contract(leftconj, {}, Conjugate, right, {}, NotConjugate, {});

	  // Contract the spin components S and Q with the gammas, and put the result on r[d=disp_indices.size()]
	  Tensor<Nout - 1, COMPLEX> aux0 =
	    r.template like_this<Nout - 1, COMPLEX>("gmNqnst%", '%', "d");
	  aux0.contract(gammas, {}, NotConjugate, aux, {}, NotConjugate);
	  aux0.copyTo(r.kvslice_from_size({{'d', disp_indices.size()}}, {{'d', 1}}));

	  // Annotate on disp_indices the displacement being computed for the current `d`
	  disp_indices.push_back(disps.disp_index);
	}

	// Apply displacements on the right and call recursively
	const int num_vecs = right.kvdim()['n'];
	unsigned int node_disp = 0;
	for (const auto it : disps.p)
	{
	  detail::log(1, "push on direction " + std::to_string(it.first));
	  // Apply displacement on the right vectors
	  // NOTE: avoid that the memory requirements grow linearly with the number of displacements
	  //       by killing the reference to `right` as soon as possible
	  Tensor<Nright, COMPLEX> right_disp =
	    !deriv ? displace(u, right, first_tslice, it.first)
		   : leftRightNabla(u, right, first_tslice, it.first, moms);
	  if (node_disp == disps.p.size() - 1)
	    right.release();
	  doMomGammaDisp_contractions(u, leftconj, std::move(right_disp), first_tslice, it.second,
				      deriv, gammas, moms, max_rhs - num_vecs, r, disp_indices);
	  node_disp++;
	  detail::log(1, "pop direction");
	}
      }
    }

    using CoorMoms = std::vector<Coor<3>>;

    template <typename COMPLEX = Complex>
    using Moms = std::pair<Tensor<Nd + 2, COMPLEX>, std::vector<Coor<3>>>;

    /// Contract two LatticeFermion with different momenta, gammas, and displacements.
    /// \param leftconj: left lattice fermion tensor, cxyzXNQqt
    /// \param right: right lattice fermion tensor, cxyzXnSst
    /// \param first_tslice: first time-slice in leftconj and right
    /// \param moms: momenta to apply
    /// \param moms_first: index of the first momenta to apply
    /// \param num_moms: number of momenta to apply (if none, apply all of them)
    /// \param gammas: list of gamma matrices to apply
    /// \param disps: list of displacements/derivatives
    /// \param deriv: if true, do left-right nabla derivatives
    /// \param max_rhs: maximum number of vectors hold in memory
    /// \param order_out: coordinate order of the output tensor, a permutation of nNSQmgd where
    ///        q and N (s and n) are the spin and vector from left (right) vectors, m is the momentum
    ///        index, g is the gamma index, and d is the displacement index
    /// \return: a pair made of a tensor sqnNmgd and a vector that is a dictionary that map each `d`
    ///        index in the tensor with an input displacement index.

    template <std::size_t Nout, std::size_t Nleft, std::size_t Nright, typename COMPLEX>
    std::pair<Tensor<Nout, COMPLEX>, std::vector<int>> doMomGammaDisp_contractions(
      const multi1d<LatticeColorMatrix>& u, Tensor<Nleft, COMPLEX> leftconj,
      Tensor<Nright, COMPLEX> right, Index first_tslice, const CoorMoms& moms, int first_mom,
      Maybe<int> num_moms, const std::vector<Tensor<2, COMPLEX>>& gammas,
      const std::vector<std::vector<int>>& disps, bool deriv,
      const std::string& order_out = "gmNndsqt", Maybe<int> max_active_tslices = none,
      DeviceHost dev = OnDefaultDevice)
    {
      detail::check_order_contains(order_out, "gmNndsqt");
      detail::check_order_contains(leftconj.order, "cxyzXNQqt");
      detail::check_order_contains(right.order, "cxyzXnSst");

      if (right.kvdim()['t'] != leftconj.kvdim()['t'])
	throw std::runtime_error("The t component of `right' and `left' does not match");
      int Nt = right.kvdim()['t'];

      int max_t = max_active_tslices.getSome(Nt);
      if (max_t <= 0)
	max_t = Nt;

      // Form a tree with the displacement paths
      detail::PathNode tree_disps = ns_doMomGammaDisp_contractions::get_tree(disps);

      // Get what directions are going to be used and the maximum number of displacements in memory
      std::array<bool, Nd> active_dirs{{}};
      unsigned int max_active_disps = 0;
      detail::get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

      // Number of moments to apply
      int numMom = num_moms.getSome(moms.size());
      if (first_mom + numMom > moms.size())
	throw std::runtime_error("Invalid range of momenta");

      // Allocate output tensor
      std::map<char, int> r_size = {{'t', Nt},
				    {'n', right.kvdim()['n']},
				    {'s', right.kvdim()['s']},
				    {'N', leftconj.kvdim()['N']},
				    {'q', leftconj.kvdim()['q']},
				    {'m', numMom},
				    {'g', gammas.size()},
				    {'d', disps.size()}};
      for (char c : detail::remove_dimensions(order_out, "gmNndsqt"))
	r_size[c] = leftconj.kvdim()[c];
      Tensor<Nout, COMPLEX> r(order_out, kvcoors<Nout>(order_out, r_size));

      // Create mom_list
      std::vector<Coor<Nd - 1>> mom_list(moms.begin() + first_mom,
					 moms.begin() + first_mom + numMom);

      // Copy all gammas into a single tensor
      Tensor<3, COMPLEX> gammast("gQS", {(Index)gammas.size(), Ns, Ns}, dev, OnEveryoneReplicated);
      for (unsigned int g = 0; g < gammas.size(); g++)
      {
	gammas[g]
	  .rename_dims({{'i', 'Q'}, {'j', 'S'}})
	  .copyTo(gammast.kvslice_from_size({{'g', g}}, {{'g', 1}}));
      }

      // Iterate over time-slices
      std::vector<int> disp_indices;

      for (int tfrom = 0, tsize = std::min(max_t, Nt); tfrom < Nt;
	   tfrom += tsize, tsize = std::min(max_t, Nt - tfrom))
      {
	// Make tsize one or even
	if (tsize > 1 && tsize % 2 != 0)
	  --tsize;

	detail::log(1, "contracting " + std::to_string(tsize) +
			 " tslices from tslice= " + std::to_string(tfrom));

	disp_indices.resize(0);

	// Copy moms into a single tensor
	std::string momst_order = "mxyzXt";
	Tensor<Nd + 2, COMPLEX> momst(
	  momst_order, latticeSize<Nd + 2>(momst_order, {{'t', tsize}, {'m', numMom}}), dev);
	for (int m = 0; m < numMom; ++m)
	{
	  ns_getColorvecs::getPhase<COMPLEX>(moms[first_mom + m], first_tslice + tfrom, tsize,
					     momst.getDev())
	    .copyTo(momst.kvslice_from_size({{'m', m}}, {{'m', 1}}));
	}

	// Apply momenta conjugated to the left tensor and rename the spin components s and Q to q and Q,
	// and the colorvector component n to N
	Tensor<Nleft + 1, COMPLEX> moms_left = leftconj.template like_this<Nleft + 1>(
	  "mQNqc%xyzXt", '%', "", {{'m', numMom}, {'t', tsize}});
	moms_left.contract(std::move(momst), {}, Conjugate,
			   leftconj.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}}), {},
			   NotConjugate);
	if (tfrom + tsize >= Nt)
	  leftconj.release();

	// Make a copy of the time-slicing of u[d] also supporting left and right
	std::vector<Tensor<Nd + 3, Complex>> ut(Nd);
	for (unsigned int d = 0; d < Nd - 1; d++)
	{
	  if (!active_dirs[d])
	    continue;

	  // NOTE: This is going to create a tensor with the same distribution of the t-dimension as leftconj and right
	  ut[d] = asTensorView(u[d])
		    .kvslice_from_size({{'t', first_tslice + tfrom}}, {{'t', tsize}})
		    .toComplex()
		    .make_sure(none, dev);
	}

	// Do the thing
	auto this_right = right.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});
	if (tfrom + tsize >= Nt)
	  right.release();
	auto this_r = r.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});
	if (!deriv)
	{
	  ns_doMomGammaDisp_contractions::doMomGammaDisp_contractions(
	    ut, std::move(moms_left), std::move(this_right), first_tslice + tfrom, tree_disps,
	    deriv, gammast, mom_list, 0, this_r, disp_indices);
	}
	else
	{
	  throw std::runtime_error("Derivatives are not implemented! Sorry!");
	  // std::vector<COMPLEX> ones(moms.numMom(), COMPLEX(1));
	  // std::string right_moms_order = std::string(right.order.begin(), right.order.size()) + "m";
	  // Tensor<Nright + 1, COMPLEX> right_moms =
	  //   right.like_this<Nright + 1>(right_moms_order.c_str());
	  // right_moms.contract(asTensorView(ones), {{'i', 'm'}}, NotConjugate, std::move(right), {},
	  // 		    NotConjugate);
	  // doMomGammaDisp_contractions(u, gammast_moms_left, right_moms, tree_disps, deriv, mom_list,
	  // 			    max_rhs, r, disp_indices);
	}
      }

      return {r, disp_indices};
    }

    /// Callback function for each displacement/derivate, and chunk of time-slices and momenta
    /// Arguments of the callback:
    /// \param tensor: output tensor with order ijkmt
    /// \param disp: index of the displacement/derivative
    /// \param first_timeslice: index of the first time-slice in the tensor
    /// \param first_mom: index of the first momentum in the tensor

    template <typename COMPLEX = Complex>
    using ColorContractionFn = std::function<void(Tensor<5, COMPLEX>, int, int, int)>;

    namespace ns_doMomDisp_colorContractions
    {
      using namespace detail;

      /// Return the tree representing all paths
      inline PathNode get_tree(const std::vector<std::array<std::vector<int>, 3>>& paths)
      {
	// Concatenate the three paths in each displacement
	std::vector<std::vector<int>> paths_out;
	for (const std::array<std::vector<int>, 3>& tripletpath : paths)
	{
	  std::vector<int> p;
	  for (int i = 0; i < 3; ++i)
	  {
	    p.insert(p.end(), tripletpath[i].begin(), tripletpath[i].end());
	    if (i < 2)
	      p.push_back(path_separator);
	  }
	  paths_out.push_back(p);
	}
	return detail::get_tree(paths_out, true);
      }

      /// Contract three LatticeColorvec with different momenta and displacements.
      /// Auxiliary function traversing the tree for disps2.
      /// \param colorvecs: lattice color tensor on several t_slices, ctxyzXn
      /// \param disps: tree of displacements/derivatives for colorvecs
      /// \param deriv: if true, do right nabla derivatives
      /// \param moms: momenta tensor on several t_slices, mtxyzX
      /// \param first_mom: index of the first momentum being computed
      /// \param max_cols: maximum number from colorvecs[0] to be contracted at once
      /// \param order_out: coordinate order of the output tensor, a permutation of ijkmt
      /// \param call: function to call for each combination of disps0, disps1,
      ///        and disps2.

      template <typename COMPLEX, std::size_t Nin>
      void doMomDisp_colorContractions(const std::vector<Tensor<Nd + 3, COMPLEX>>& u,
				       std::array<Tensor<Nin, COMPLEX>, 3> colorvecs,
				       Index first_tslice, const PathNode& disps, bool deriv,
				       int current_colorvec, const Moms<COMPLEX> moms,
				       int first_mom, int max_cols, const std::string& order_out,
				       DeviceHost dev, Distribution dist,
				       const ColorContractionFn<COMPLEX>& call)
      {
	if (disps.disp_index >= 0)
	{
	  detail::log(1, "contracting for disp_index=" + std::to_string(disps.disp_index));

	  // Create the output tensor
	  Tensor<5, COMPLEX> colorvec012m =
	    colorvecs[0].template like_this<5, COMPLEX>(order_out,
							{{'i', colorvecs[0].kvdim()['n']},
							 {'j', colorvecs[1].kvdim()['n']},
							 {'k', colorvecs[2].kvdim()['n']},
							 {'m', moms.first.kvdim()['m']}},
							dev, dist);

	  // Contract colorvec2 and moms
	  Tensor<Nd + 4, COMPLEX> colorvec2m =
	    colorvecs[2]
	      .template like_this<Nd + 4, COMPLEX>("ncm%xyzXt", '%', "",
						   {{'m', moms.first.kvdim()['m']}})
	      .rename_dims({{'n', 'k'}});
	  colorvec2m.contract(colorvecs[2].rename_dims({{'n', 'k'}}), {}, NotConjugate, moms.first,
			      {}, NotConjugate);

	  int imax = max_cols;
	  if (imax <= 0)
	    imax = colorvecs[0].kvdim()['n'];

	  for (int i0 = 0, i1 = colorvecs[0].kvdim()['n'], isize = std::min(imax, i1); i0 < i1;
	       i0 += isize, isize = std::min(imax, i1 - i0))
	  {
	    // Color-contract colorvec0 and colorvec1
	    Tensor<Nin + 1, COMPLEX> colorvec01 =
	      colorvecs[0]
		.template like_this<Nin + 1, COMPLEX>(
		  "njcxyzXt%", '%', "", {{'n', isize}, {'j', colorvecs[1].kvdim()['n']}})
		.rename_dims({{'n', 'i'}});
	    auto colorvec0 =
	      colorvecs[0].rename_dims({{'n', 'i'}}).kvslice_from_size({{'i', i0}}, {{'i', isize}});
	    auto colorvec1 = colorvecs[1].rename_dims({{'n', 'j'}});
	    colorvec01.contract(colorvec0.kvslice_from_size({{'c', 2}}), {}, NotConjugate,
				colorvec1.kvslice_from_size({{'c', 1}}), {}, NotConjugate);
	    colorvec01.contract(colorvec0.kvslice_from_size({{'c', 1}}), {}, NotConjugate,
				colorvec1.kvslice_from_size({{'c', 2}}), {}, NotConjugate, {}, -1);
	    colorvec0.release();
	    colorvec1.release();

	    // Contract colorvec01 and colorvec2m
	    colorvec012m.kvslice_from_size({{'i', i0}}, {{'i', isize}})
	      .contract(std::move(colorvec01), {}, NotConjugate, colorvec2m, {}, NotConjugate);
	  }

	  // Do whatever
	  call(std::move(colorvec012m), disps.disp_index, first_tslice, first_mom);
	}

	// Apply displacements on colorvec2 and call recursively
	unsigned int node_disp = 0;
	for (const auto it : disps.p)
	{
	  detail::log(1, "for disps, push on direction " + std::to_string(it.first));
	  // Apply displacement on the current colorvec
	  // NOTE: avoid that the memory requirements grow linearly with the number of displacements
	  //       by killing the reference to `colorvec2` as soon as possible
	  std::array<Tensor<Nin, COMPLEX>, 3> colorvecs_disp = colorvecs;
	  int this_current_colorvec = current_colorvec;
	  if (abs(it.first) <= Nd)
	    colorvecs_disp[current_colorvec] =
	      !deriv ? displace(u, colorvecs[current_colorvec], first_tslice, it.first)
		     : rightNabla(u, colorvecs[current_colorvec], first_tslice, it.first);
	  else if (it.first == path_separator)
	    ++this_current_colorvec;
	  else
	    throw std::runtime_error("Invalid direction");
	  if (node_disp == disps.p.size() - 1)
	    for (auto& i : colorvecs)
	      i.release();
	  doMomDisp_colorContractions(u, std::move(colorvecs_disp), first_tslice, it.second, deriv,
				      this_current_colorvec, moms, first_mom, max_cols, order_out,
				      dev, dist, call);
	  node_disp++;
	  detail::log(1, "for disps, pop direction");
	}
      }
    }

    /// Contract three LatticeColorvec with different momenta and displacements.
    /// \param colorvecs: lattice color tensor on several t_slices, ctxyzXn
    /// \param moms: momenta tensor on several t_slices, mtxyzX
    /// \param disps: list of displacements/derivatives
    /// \param deriv: if true, do right nabla derivatives
    /// \param call: function to call for each combination of disps0, disps1, and disps2
    /// \param order_out: coordinate order of the output tensor, a permutation of ijkmt where
    ///        i, j, and k are the n index in colorvecs; and m is the momentum index

    template <std::size_t Nin, typename COMPLEX>
    void doMomDisp_colorContractions(
      const multi1d<LatticeColorMatrix>& u, Tensor<Nin, COMPLEX> colorvec, const CoorMoms& moms,
      Index first_tslice, const std::vector<std::array<std::vector<int>, 3>>& disps, bool deriv,
      const ColorContractionFn<COMPLEX>& call, Maybe<int> max_active_tslices = none,
      Maybe<int> max_active_momenta = none, Maybe<int> max_cols = none,
      const Maybe<std::string>& order_out = none, Maybe<DeviceHost> dev = none,
      Maybe<Distribution> dist = none)
    {
      const std::string order_out_str = order_out.getSome("ijkmt");
      detail::check_order_contains(order_out_str, "ijkmt");
      detail::check_order_contains(colorvec.order, "cxyzXtn");

      // Form a tree with the displacement paths
      detail::PathNode tree_disps = ns_doMomDisp_colorContractions::get_tree(disps);

      // Get what directions are going to be used and the maximum number of displacements in memory
      std::array<bool, Nd> active_dirs{{}};
      unsigned int max_active_disps = 0;
      detail::get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

      // Check that all tensors have the same number of time
      int Nt = colorvec.kvdim()['t'];

      int max_t = max_active_tslices.getSome(Nt);
      if (max_t <= 0)
	max_t = Nt;

      int Nmom = moms.size();
      int max_active_moms = max_active_momenta.getSome(Nmom);
      if (max_active_moms <= 0)
	max_active_moms = Nmom;

      // Iterate over time-slices
      for (int tfrom = 0, tsize = std::min(max_t, Nt); tfrom < Nt;
	   tfrom += tsize, tsize = std::min(max_t, Nt - tfrom))
      {
	detail::log(1, "color contracting " + std::to_string(tsize) +
			 " tslices from tslice= " + std::to_string(first_tslice + tfrom));

	// Make a copy of the time-slicing of u[d] also supporting left and right
	std::vector<Tensor<Nd + 3, COMPLEX>> ut(Nd);
	for (unsigned int d = 0; d < Nd - 1; d++)
	{
	  if (!active_dirs[d])
	    continue;

	  // NOTE: This is going to create a tensor with the same distribution of the t-dimension as colorvec and moms
	  ut[d] = asTensorView(u[d])
		    .kvslice_from_size({{'t', first_tslice + tfrom}}, {{'t', tsize}})
		    .toComplex();
	}

	// Get the time-slice for colorvec
	auto this_colorvec = colorvec.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});

	// Loop over the momenta
	for (int mfrom = 0, msize = std::min(max_active_moms, Nmom); mfrom < Nmom;
	     mfrom += msize, msize = std::min(max_active_moms, Nmom - mfrom))
	{
	  auto this_moms =
	    this_colorvec.template like_this<Nd + 2, COMPLEX>("xyzXtm", {{'m', msize}});
	  for (int m = 0; m < msize; ++m)
	    ns_getColorvecs::getPhase<COMPLEX>(moms[mfrom + m], first_tslice + tfrom, tsize,
					       this_moms.getDev())
	      .copyTo(this_moms.kvslice_from_size({{'m', m}}, {{'m', 1}}));

	  if (tfrom + tsize >= Nt && mfrom + msize >= Nmom)
	  {
	    colorvec.release();
	  }
	  std::vector<Coor<3>> moms_list(moms.begin() + mfrom, moms.begin() + mfrom + msize);
	  if (!deriv)
	  {
	    ns_doMomDisp_colorContractions::doMomDisp_colorContractions<COMPLEX, Nin>(
	      ut, {this_colorvec, this_colorvec, this_colorvec}, first_tslice + tfrom, tree_disps,
	      deriv, 0, {this_moms, moms_list}, mfrom, max_cols.getSome(0), order_out_str,
	      dev.getSome(OnDefaultDevice), dist.getSome(OnEveryoneReplicated), call);
	  }
	  else
	  {
	    // When using derivatives, each momenta has a different effect
	    std::vector<COMPLEX> ones(msize, COMPLEX(1));
	    Tensor<Nin + 1, COMPLEX> this_colorvec_m =
	      this_colorvec.template like_this<Nin + 1>("%m", '%', "", {{'m', msize}});
	    this_colorvec_m.contract(this_colorvec, {}, NotConjugate, asTensorView(ones),
				     {{'i', 'm'}}, NotConjugate);
	    ns_doMomDisp_colorContractions::doMomDisp_colorContractions<COMPLEX, Nin + 1>(
	      ut, {this_colorvec_m, this_colorvec_m, this_colorvec_m}, first_tslice + tfrom,
	      tree_disps, deriv, 0, {this_moms, moms_list}, mfrom, max_cols.getSome(0),
	      order_out_str, dev.getSome(OnDefaultDevice), dist.getSome(OnEveryoneReplicated),
	      call);
	  }
	}
      }
    }

    /// Callback function for each displacement/derivate, and chunk of time-slices and momenta
    /// Arguments of the callback:
    /// \param tensor: output tensor with order ijmt, where i and j are the right and left colorvec indices,
    ///        m is the momentum index, and t is the t-slice
    /// \param disp: index of the displacement/derivative
    /// \param first_timeslice: index of the first time-slice in the tensor
    /// \param first_mom: index of the first momentum in the tensor

    template <typename COMPLEX = Complex>
    using ContractionFn = std::function<void(Tensor<4, COMPLEX>, int, int, int)>;

    namespace ns_doMomDisp_contractions
    {
      using namespace detail;

      /// Contract two LatticeColorvec with different momenta and displacements.
      /// Auxiliary function traversing the tree for disps2.
      /// \param colorvecs: lattice color tensor on several t_slices, ctxyzXn
      /// \param disps: tree of displacements/derivatives for colorvecs
      /// \param moms: momenta tensor on several t_slices, mtxyzX
      /// \param first_mom: index of the first momentum being computed
      /// \param order_out: coordinate order of the output tensor, a permutation of ijkmt
      /// \param call: function to call for each combination of displacement, t-slice, and momentum

      template <typename COMPLEX, std::size_t Nleft, std::size_t Nright>
      void doMomDisp_contractions(const std::vector<Tensor<Nd + 3, COMPLEX>>& u,
				  Tensor<Nleft, COMPLEX> left, Tensor<Nright, COMPLEX> right,
				  Index first_tslice, const PathNode& disps, bool deriv,
				  const std::vector<Coor<Nd - 1>>& moms, int first_mom,
				  const std::string& order_out, DeviceHost dev, Distribution dist,
				  const ContractionFn<COMPLEX>& call)
      {
	if (disps.disp_index >= 0)
	{
	  detail::log(1, "contracting for disp_index=" + std::to_string(disps.disp_index));

	  // Contract left and right
	  auto this_right = right.rename_dims({{'n', 'i'}});
	  auto this_left = left.rename_dims({{'n', 'j'}});
	  Tensor<4, COMPLEX> r = this_left.template like_this<4, COMPLEX>(
	    "jimt", {{'i', this_right.kvdim()['i']}}, dev, dist);
	  r.contract(std::move(this_left), {}, Conjugate, std::move(this_right), {}, NotConjugate);

	  // Do whatever
	  call(std::move(r), disps.disp_index, first_tslice, first_mom);
	}

	// Apply displacements on right and call recursively
	unsigned int node_disp = 0;
	for (const auto it : disps.p)
	{
	  detail::log(1, "for disps, push on direction " + std::to_string(it.first));
	  // Apply displacement on the right colorvec
	  // NOTE: avoid that the memory requirements grow linearly with the number of displacements
	  //       by killing the reference to `right` as soon as possible
	  Tensor<Nright, COMPLEX> right_disp =
	    !deriv ? displace(u, right, first_tslice, it.first)
		   : leftRightNabla(u, right, first_tslice, it.first, moms);
	  if (node_disp == disps.p.size() - 1)
	    right.release();
	  doMomDisp_contractions(u, left, std::move(right_disp), first_tslice, it.second, deriv,
				 moms, first_mom, order_out, dev, dist, call);
	  node_disp++;
	  detail::log(1, "for disps, pop direction");
	}
      }
    }

    /// Contract three LatticeColorvec with different momenta and displacements.
    /// It computes
    ///    \eta_j^\dagger exp(- i left_phase \cdot x) \Gamma \eta_k exp(i right_phase \cdot x)
    /// where \eta_i is the ith colorvec, x is a lattice site, and \Gamma is a combination of
    /// derivatives and momenta.
    /// \param colorvecs: lattice color tensor on several t_slices, ctxyzXn
    /// \param left_phase: phase to the left colorvecs
    /// \param right_phase: phase to the right colorvecs
    /// \param moms: momenta tensor on several t_slices, mtxyzX
    /// \param disps: list of displacements/derivatives
    /// \param deriv: if true, do left-right nabla derivatives
    /// \param call: function to call for each combination of disps0, disps1, and disps2
    /// \param order_out: coordinate order of the output tensor, a permutation of ijmt where
    ///        i and j are the n index in the right and left colorvec respectively; and
    ///        m is the momentum index

    template <std::size_t Nin, typename COMPLEX>
    void doMomDisp_contractions(const multi1d<LatticeColorMatrix>& u, Tensor<Nin, COMPLEX> colorvec,
				Coor<3> left_phase, Coor<3> right_phase, const CoorMoms& moms,
				Index first_tslice, const std::vector<std::vector<int>>& disps,
				bool deriv, const ContractionFn<COMPLEX>& call,
				const Maybe<std::string>& order_out = none,
				Maybe<DeviceHost> dev = none, Maybe<Distribution> dist = none,
				int max_tslices_in_contraction = 0, int max_moms_in_contraction = 0)
    {
      const std::string order_out_str = order_out.getSome("ijmt");
      detail::check_order_contains(order_out_str, "ijmt");
      detail::check_order_contains(colorvec.order, "cxyzXtn");

      // Form a tree with the displacement paths
      detail::PathNode tree_disps = detail::get_tree(disps);

      // Get what directions are going to be used and the maximum number of displacements in memory
      std::array<bool, Nd> active_dirs{{}};
      unsigned int max_active_disps = 0;
      detail::get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

      // Check that all tensors have the same number of time
      int Nt = colorvec.kvdim()['t'];

      if (max_tslices_in_contraction <= 0)
	max_tslices_in_contraction = Nt;
      int Nmom = moms.size();
      if (max_moms_in_contraction <= 0)
	max_moms_in_contraction = Nmom;

      // Iterate over time-slices
      for (int tfrom = 0, tsize = std::min(Nt, max_tslices_in_contraction); tfrom < Nt;
	   tfrom += tsize, tsize = std::min(max_tslices_in_contraction, Nt - tfrom))
      {
	// Make tsize one or even
	if (tsize > 1 && tsize % 2 != 0)
	  --tsize;

	detail::log(1, "contracting " + std::to_string(tsize) +
			 " tslices from tslice= " + std::to_string(tfrom));

	// Make a copy of the time-slicing of u[d] also supporting left and right
	std::vector<Tensor<Nd + 3, COMPLEX>> ut(Nd);
	for (unsigned int d = 0; d < Nd - 1; d++)
	{
	  if (!active_dirs[d])
	    continue;

	  // NOTE: This is going to create a tensor with the same distribution of the t-dimension as colorvec and moms
	  ut[d] = asTensorView(u[d])
		    .kvslice_from_size({{'t', first_tslice + tfrom}}, {{'t', tsize}})
		    .toComplex();
	}

	// Get the time-slice for colorvec
	auto this_colorvec = colorvec.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});

	// Apply the phases
	auto this_colorvec_phase_right =
	  phaseColorvecs(this_colorvec, first_tslice + tfrom, right_phase);
	auto this_colorvec_phase_left =
	  phaseColorvecs(this_colorvec, first_tslice + tfrom, left_phase);

	// Loop over the momenta
	for (int mfrom = 0, msize = std::min(max_moms_in_contraction, Nmom); mfrom < Nmom;
	     mfrom += msize, msize = std::min(max_moms_in_contraction, Nmom - mfrom))
	{

	  auto this_moms =
	    this_colorvec_phase_left.template like_this<Nd + 2, COMPLEX>("xyzXtm", {{'m', msize}});
	  for (int m = 0; m < msize; ++m)
	    ns_getColorvecs::getPhase<COMPLEX>(moms[mfrom + m], first_tslice + tfrom, tsize,
					       this_moms.getDev())
	      .copyTo(this_moms.kvslice_from_size({{'m', m}}, {{'m', 1}}));

	  // Apply left phase and momenta conjugated to the left tensor
	  // NOTE: look for the minus sign on left_phase in the doc of this function
	  Tensor<Nin + 1, COMPLEX> moms_left =
	    this_colorvec.template like_this<Nin + 1>("mc%xyzXt", '%', "", {{'m', msize}});
	  moms_left.contract(std::move(this_moms), {}, Conjugate, this_colorvec_phase_left, {},
			     NotConjugate);

	  if (tfrom + tsize >= Nt && mfrom + msize >= Nmom)
	  {
	    colorvec.release();
	  }

	  auto this_moms_coors =
	    std::vector<Coor<Nd - 1>>(moms.begin() + mfrom, moms.begin() + mfrom + msize);
	  if (!deriv)
	  {
	    ns_doMomDisp_contractions::doMomDisp_contractions<COMPLEX>(
	      ut, std::move(moms_left), this_colorvec_phase_right, first_tslice + tfrom, tree_disps,
	      deriv, this_moms_coors, mfrom, order_out_str, dev.getSome(OnDefaultDevice),
	      dist.getSome(OnEveryoneReplicated), call);
	  }
	  else
	  {
	    // When using derivatives, each momenta has a different effect
	    std::vector<COMPLEX> ones(msize, COMPLEX(1));
	    Tensor<Nin + 1, COMPLEX> this_colorvec_m =
	      this_colorvec.template like_this<Nin + 1>("%m", '%', "", {{'m', msize}});
	    this_colorvec_m.contract(this_colorvec_phase_right, {}, NotConjugate,
				     asTensorView(ones), {{'i', 'm'}}, NotConjugate);
	    ns_doMomDisp_contractions::doMomDisp_contractions<COMPLEX>(
	      ut, std::move(moms_left), std::move(this_colorvec_m), first_tslice + tfrom,
	      tree_disps, deriv, this_moms_coors, mfrom, order_out_str,
	      dev.getSome(OnDefaultDevice), dist.getSome(OnEveryoneReplicated), call);
	  }
	}
      }
    }

    /// Call the destroy list and clean superbblas cache

    inline void finish()
    {
      // Show performance reports
      // NOTE: QDPIO::cout doesn't support iomanip format declarations, so use std::cout on master node instead
      if (Layout::nodeNumber() == 0)
      {
	superbblas::reportTimings(std::cout);
	superbblas::reportCacheUsage(std::cout);
      }

      // Clear internal superbblas caches
      superbblas::clearCaches();

      // Call the destroy list
      for (const auto& f : detail::getDestroyList())
	f();
      detail::getDestroyList().clear();
    }

    /// Return the smallest interval containing the union of two intervals
    /// \param from0: first element of the first interval
    /// \param size0: length of the first interval
    /// \param from1: first element of the second interval
    /// \param size1: length of the second interval
    /// \param dim: dimension length
    /// \param fromr: (output) first element of the union of the two intervals
    /// \param sizer: (output) length of the union of the two intervals

    inline void union_interval(Index from0, Index size0, Index from1, Index size1, Index dim,
			       Index& fromr, Index& sizer)
    {
      // Check inputs
      if (size0 > dim || size1 > dim)
	throw std::runtime_error(
	  "Invalid interval to union! Some of input intervals exceeds the lattice dimension");

      // Normalize from and take as from0 the leftmost interval of the two input intervals
      from0 = normalize_coor(from0, dim);
      from1 = normalize_coor(from1, dim);
      if (from0 > from1)
      {
	std::swap(from0, from1);
	std::swap(size0, size1);
      }

      // If some interval is empty, return the other
      if (size0 == 0)
      {
	fromr = from1;
	sizer = size1;
      }
      else if (size1 == 0)
      {
	fromr = from0;
	sizer = size0;
      }
      else
      {
	// Return the shortest interval resulting from the leftmost point of the
	// first interval and the rightmost point of both intervals, and the
	// leftmost point of the second interval and the rightmost point of both
	// intervals

	Index fromra = from0;
	Index sizera = std::max(from0 + size0, from1 + size1) - from0;
	Index fromrb = from1;
	Index sizerb = std::max(from0 + dim + size0, from1 + size1) - from1;
	fromr = (sizera <= sizerb ? fromra : fromrb);
	sizer = (sizera <= sizerb ? sizera : sizerb);
      }

      // Normalize the output if the resulting interval is the whole dimension
      if (sizer >= dim)
      {
	fromr = 0;
	sizer = dim;
      }
    }

    /// Return a multi1d from std::vector
    /// \param v: vector to convert

    template <typename T>
    multi1d<T> tomulti1d(const std::vector<T>& v)
    {
      multi1d<T> r(v.size());
      for (int i = 0; i < v.size(); ++i)
	r[i] = v[i];
      return r;
    }

    /// Return a multi1d from std::array
    /// \param v: array to convert

    template <typename T, std::size_t N>
    multi1d<T> tomulti1d(const std::array<T, N>& v)
    {
      multi1d<T> r(v.size());
      for (int i = 0; i < v.size(); ++i)
	r[i] = v[i];
      return r;
    }

    /// Return all momenta with magnitude squared within a range
    /// \param min_mom2: all returned momenta should have this magnitude squared at least
    /// \param max_mom2: all returned momenta should have up to this magnitude squared

    inline CoorMoms getMomenta(int min_mom2, int max_mom2)
    {
      static_assert(Nd == 4);
      int max_component = (int)std::sqrt((float)max_mom2) + 1;
      CoorMoms r;
      for (int i = -max_component; i <= max_component; ++i)
      {
	for (int j = -max_component; j <= max_component; ++j)
	{
	  for (int k = -max_component; k <= max_component; ++k)
	  {
	    int mom_magnitude2 = i * i + j * j + k * k;
	    if (min_mom2 <= mom_magnitude2 && mom_magnitude2 <= max_mom2)
	      r.push_back(Coor<3>{i, j, k});
	  }
	}
      }
      return r;
    }

    /// Return a list of momenta as std::vector<Coor<3>> from std::vector<std::vector<int>>
    /// \param v: list of momenta to transform

    inline CoorMoms getMomenta(const std::vector<std::vector<int>>& v)
    {
      static_assert(Nd == 4);
      CoorMoms r;
      for (const auto vi : v)
      {
	Coor<3> c;
	std::copy_n(vi.begin(), 3, c.begin());
	r.push_back(c);
      }
      return r;
    }

  }
}

namespace QDP
{
  //! Binary input
  template <std::size_t N, typename T>
  void read(BinaryReader& bin, Chroma::SB::Tensor<N, T> t)
  {
    t.binaryRead(bin);
  }

  //! Binary output
  template <std::size_t N, typename T>
  inline void write(BinaryWriter& bin, const Chroma::SB::Tensor<N, T> t)
  {
    t.binaryWrite(bin);
  }
}

#endif // BUILD_SB
#endif // __INCLUDE_SUPERB_CONTRACTIONS__
