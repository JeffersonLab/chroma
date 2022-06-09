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
#  include "meas/smear/link_smearing_factory.h"
#  include "qdp.h"
#  include "qdp_map_obj_disk_multiple.h"
#  include "superbblas.h"
#  include "util/ferm/key_timeslice_colorvec.h"
#  include "util/ft/sftmom.h"
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

#  if defined(BUILD_MAGMA)
#    include "magma_v2.h"
#  endif

namespace Chroma
{

  namespace SB
  {

    using Index = superbblas::IndexType;
    using Complex = std::complex<REAL>;
    using ComplexD = std::complex<REAL64>;
    using ComplexF = std::complex<REAL32>;
    template <std::size_t N>
    using Coor = superbblas::Coor<N>;
    using checksum_type = superbblas::checksum_type;

    /// Where to store the tensor (see class Tensor)
    enum DeviceHost {
      OnHost,	      ///< on cpu memory
      OnDefaultDevice ///< on GPU memory if possible
    };

    /// How to distribute the tensor (see class Tensor)
    enum Distribution {
      OnMaster,		    ///< Fully supported on node with index zero
      OnEveryone,	    ///< Distribute the lattice dimensions (x, y, z, t) as chroma does
      OnEveryoneReplicated, ///< All nodes have a copy of the tensor
      Local		    ///< Non-collective
    };

    /// Whether complex conjugate the elements before contraction (see Tensor::contract)
    enum Conjugation { NotConjugate, Conjugate };

    /// Whether the tensor is dense or sparse (see StorageTensor)
    enum Sparsity { Dense, Sparse };

    /// Whether to copy or add the values into the destination tensor (see Tensor::doAction)
    enum Action { CopyTo, AddTo };

    /// Auxiliary class for initialize Maybe<T> with no value
    struct None {
    };

    /// Class for optional values
    template <typename T>
    struct Maybe {
      /// opt_val.first is whether a value was set, and opt_val.second has the value if that's the case
      std::pair<bool, T> opt_val;

      /// Constructor without a value
      Maybe() : opt_val{false, {}}
      {
      }

      /// Constructor without a value
      Maybe(None) : Maybe()
      {
      }

      /// Constructor with a value
      template <typename Q,
		typename std::enable_if<std::is_convertible<Q, T>::value, bool>::type = true>
      Maybe(const Q& t) : opt_val{true, T(t)}
      {
      }

      /// Return whether it has been initialized with a value
      bool hasSome() const
      {
	return opt_val.first;
      }

      /// Return whether it has been initialized with a value
      explicit operator bool() const noexcept
      {
	return hasSome();
      }

      /// Return the value if it has been initialized with some
      T getSome() const
      {
	if (opt_val.first)
	  return opt_val.second;
	throw std::runtime_error("W!");
      }

      /// Return the value if it has been initialized with some; otherwise return `def`
      T getSome(T def) const
      {
	if (opt_val.first)
	  return opt_val.second;
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
	: superbblas::detail::tracker<superbblas::detail::Cpu>(funcName, superbblas::detail::Cpu{})
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

    inline int normalize_coor(int v, int dim)
    {
      return (v + dim * (v < 0 ? -v / dim + 1 : 0)) % dim;
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

      /// Reshape several dimensions into another dimension
      /// \param ncollapse: number of dimensions to collapse starting from each old dimension
      /// \param nsplit: number of dimensions to split starting from each old dimension
      /// \param old_dim: dimensions of the old tensor
      /// \param new_dim: maximum dimension size for the new tensor
      /// \param t: either `From` (first element) or `Size` (number of elements in each dimension)
      /// \param c: coordinate to transform

      template <std::size_t Nout, std::size_t N>
      Coor<Nout> reshape_dimensions(const Coor<N>& ncollapse, const Coor<N>& nsplit,
				    const Coor<N>& old_dim, const Coor<Nout>& new_dim, CoorType t,
				    const Coor<N>& c)
      {
	Coor<Nout> r;
	unsigned int ri = 0;
	for (unsigned int ci = 0; ci < N; ++ci)
	{
	  if (ncollapse[ci] == 1 && nsplit[ci] == 1) {
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
		    throw std::runtime_error(
		      "reshape_dimensions: unsupported to collapse a range with holes");
		  odd_dim_watched = true;
		}
		if (t == From)
		  idx += c[ci + k] * stride;
		else
		  idx *= c[ci + k];
		stride *= old_dim[ci + k];
	      }
	    }
	    ci += ncollapse[ci] - 1;

	    // Split the new dimension into it[2] new dimensions
	    {
	      Index stride = 1;
	      for (unsigned int k = 0; k < nsplit[ci]; ++k)
	      {
		if (!(t == From || idx < stride || idx % stride == 0))
		  throw std::runtime_error("reshape_dimension: Not supporting for this partition");
		if (t == From)
		  r[ri + k] = (idx / stride) % new_dim[ri + k];
		else
		  r[ri + k] = std::min((idx + stride - 1) / stride, new_dim[ri + k]);
		stride *= new_dim[ri + k];
	      }
	    }
	    ri += nsplit[ci];
	  }
	}

	// Return the new coordinates
	return r;
      }


#  if defined(BUILD_MAGMA)
      // Return a MAGMA context
      inline std::shared_ptr<magma_queue_t> getMagmaContext(Maybe<int> device = none)
      {
	static std::shared_ptr<magma_queue_t> queue;
	if (!queue)
	{
	  // Start MAGMA and create a queue
	  int dev = device.getSome(-1);
	  if (dev < 0)
	  {
#    ifdef SUPERBBLAS_USE_CUDA
	    superbblas::detail::cudaCheck(cudaGetDevice(&dev));
#    elif defined(SUPERBBLAS_USE_HIP)
	    superbblas::detail::hipCheck(hipGetDevice(&dev));
#    else
#      error superbblas was not build with support for GPUs
#    endif
	  }
	  magma_init();
	  magma_queue_t q;
	  magma_queue_create(dev, &q);
	  queue = std::make_shared<magma_queue_t>(q);
	}
	return queue;
      }
#  endif

      // Return a context on either the host or the device
      inline std::shared_ptr<superbblas::Context> getContext(DeviceHost dev)
      {
	// Creating GPU context can be expensive; so do it once
	static std::shared_ptr<superbblas::Context> cudactx;
	static std::shared_ptr<superbblas::Context> cpuctx;
	if (!cpuctx)
	  cpuctx = std::make_shared<superbblas::Context>(superbblas::createCpuContext());

	switch (dev)
	{
	case OnHost: return cpuctx;
	case OnDefaultDevice:
#  ifdef SUPERBBLAS_USE_GPU
	  if (!cudactx)
	  {

	    int dev = -1;
#    ifdef SUPERBBLAS_USE_CUDA
	    superbblas::detail::cudaCheck(cudaGetDevice(&dev));
#    elif defined(SUPERBBLAS_USE_HIP)
	    superbblas::detail::hipCheck(hipGetDevice(&dev));
#    else
#      error unsupported GPU platform
#    endif

#    if defined(BUILD_MAGMA)
	    // Force the creation of the queue before creating a superbblas context (otherwise cublas complains)
	    getMagmaContext();
#    endif

	    // Workaround on a potential issue in qdp-jit: avoid passing through the pool allocator
#    if defined(QDP_IS_QDPJIT)
	    if (jit_config_get_max_allocation() != 0)
	    {
	      cudactx = std::make_shared<superbblas::Context>(superbblas::createGpuContext(
		dev,

		// Make superbblas use the same memory allocator for gpu as any other qdp-jit lattice object
		[](std::size_t size, superbblas::platform plat) -> void* {
		  if (size == 0)
		    return nullptr;
		  if (plat == superbblas::CPU)
		    return malloc(size);
		  void* ptr = nullptr;
		  QDP_get_global_cache().addDeviceStatic(&ptr, size, true);
		  assert(superbblas::detail::getPtrDevice(ptr) >= 0);
		  return ptr;
		},

		// The corresponding deallocator
		[](void* ptr, superbblas::platform plat) {
		  if (ptr == nullptr)
		    return;
		  if (plat == superbblas::CPU)
		    free(ptr);
		  else
		    QDP_get_global_cache().signoffViaPtr(ptr);
		}));
	    }
	    else
#    endif // defined(QDP_IS_QDPJIT)
	    {
	      cudactx = std::make_shared<superbblas::Context>(superbblas::createGpuContext(dev));
	    }
	  }
	  return cudactx;
#  else	 // SUPERBBLAS_USE_GPU
	  return cpuctx;
#  endif // SUPERBBLAS_USE_GPU
	}
	throw std::runtime_error("Unsupported `DeviceHost`");
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
	  switch (dist)
	  {
	  case OnMaster: p = all_tensor_on_master(dim); break;
	  case OnEveryone: p = partitioning_chroma_compatible(order, dim); break;
	  case OnEveryoneReplicated: p = all_tensor_replicated(dim); break;
	  case Local:
	    p = local(dim);
	    isLocal = true;
	    break;
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
	TensorPartition<Nout> reshape_dimensions(const Coor<N>& ncollapse, const Coor<N>& nsplit,
						 const Coor<Nout>& new_dim) const
	{
	  typename TensorPartition<Nout>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back(
	      {detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, new_dim, From, i[0]),
	       detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, new_dim, Size, i[1])});
	  return TensorPartition<Nout>{
	    detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, new_dim, Size, dim), r,
	    isLocal};
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
	    r.push_back({lfrom, lsize});
	  }
	  return TensorPartition<N>{size, r, isLocal};
	}

	/// Return a partition with the local portion of the tensor

	TensorPartition<N> get_local_partition() const
	{
	  return TensorPartition<N>{
	    localSize(), PartitionStored(1, superbblas::PartitionItem<N>{{{}, localSize()}}), true};
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
	  }
	  return fs;
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

	template <typename Ostream>
	Ostream& operator<<(Ostream& s, Distribution dist)
	{
	  switch (dist)
	  {
	  case OnMaster: s << "OnMaster"; break;
	  case OnEveryone: s << "OnEveryone"; break;
	  case OnEveryoneReplicated: s << "OnEveryoneReplicated"; break;
	  case Local: s << "Local"; break;
	  }
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

      /// is_complex<T>::value is true if `T` is complex

      template <typename T>
      struct is_complex : std::false_type {
      };

      template <typename T>
      struct is_complex<std::complex<T>> : std::true_type {
      };

      /// real_type<T>::type is T::value_type if T is complex; otherwise it is T

      template <typename T>
      struct real_type {
	using type = T;
      };

      template <typename T>
      struct real_type<std::complex<T>> {
	using type = T;
      };

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
    }

    /// Class for operating dense tensors

    template <std::size_t N, typename T>
    struct Tensor {
      static_assert(superbblas::supported_type<T>::value, "Not supported type");

    public:
      std::string order;			///< Labels of the tensor dimensions
      Coor<N> dim;				///< Length of the tensor dimensions
      std::shared_ptr<superbblas::Context> ctx; ///< Tensor storage information (device/host)
      std::shared_ptr<T> data;			///< Pointer to the tensor storage
      std::shared_ptr<detail::TensorPartition<N>>
	p;		 ///< Distribution of the tensor among the processes
      Distribution dist; ///< Whether the tensor is stored on the cpu or a device
      Coor<N> from;	 ///< First active coordinate in the tensor
      Coor<N> size;	 ///< Number of active coordinates on each dimension
      Coor<N> strides;	 ///< Displacement for the next element along every direction
      T scalar;		 ///< Scalar factor of the tensor
      bool conjugate;	 ///< Whether the values are implicitly conjugated
      bool eg;		 ///< Whether this tensor is an example

      /// Return a string describing the tensor
      /// \param ptr: pointer to the memory allocation
      /// \return: the string representing the tensor

      std::string repr(T* ptr = nullptr) const
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
	     Distribution dist = OnEveryone)
	: Tensor(order, dim, dev, dist,
		 std::make_shared<detail::TensorPartition<N>>(
		   detail::TensorPartition<N>(order, dim, dist)))
      {
      }

      /// Empty constructor

      Tensor()
	: order(detail::getTrivialOrder(N)),
	  dim{},
	  ctx(detail::getContext(OnHost)),
	  p(std::make_shared<detail::TensorPartition<N>>(
	    detail::TensorPartition<N>(detail::getTrivialOrder(N), {}, OnEveryoneReplicated))),
	  dist(OnEveryoneReplicated),
	  from{},
	  size{},
	  strides{},
	  scalar{0},
	  conjugate{false},
	  eg{false}
      {
      }

      /// Constructor for bringing the memory allocation (see `asTensorView`)
      /// \param order: dimension labels of the tensor
      /// \param dim: size for each dimension
      /// \param dev: where to allocate the content on the GPU if available (`OnDefaultDevice`)
      ///        or on CPU always (`OnHost`)
      /// \param dist: how to distribute the tensor, see `Distribution`

      Tensor(const std::string& order, Coor<N> dim, DeviceHost dev, Distribution dist,
	     std::shared_ptr<T> data)
	: order(order),
	  dim(dim),
	  ctx(detail::getContext(dev)),
	  data(data),
	  dist(dist),
	  from{},
	  size(dim),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow)),
	  scalar{1},
	  conjugate{false},
	  eg{false}
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
      /// \param ctx: superbblas context
      /// \param data: memory allocation
      /// \param p: partition of the tensor among the processes
      /// \param dist: how to distribute the tensor, see `Distribution`
      /// \param from: coordinate of the first element in this view
      /// \param size: elements in each direction in this view
      /// \param scalar: scalar factor of this view
      /// \param conjugate: whether the elements are implicitly conjugated

      Tensor(const std::string& order, Coor<N> dim, std::shared_ptr<superbblas::Context> ctx,
	     std::shared_ptr<T> data, std::shared_ptr<detail::TensorPartition<N>> p,
	     Distribution dist, Coor<N> from, Coor<N> size, T scalar, bool conjugate, bool eg)
	: order(order),
	  dim(dim),
	  ctx(ctx),
	  data(data),
	  p(p),
	  dist(dist),
	  from(normalize_coor(from, dim)),
	  size(size),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow)),
	  scalar(scalar),
	  conjugate(conjugate),
	  eg(eg)
      {
	checkOrder();
      }

      /// Internal constructor, used by `make_suitable_for_contraction`
      /// \param order: dimension labels of the tensor
      /// \param dim: size for each dimension
      /// \param dist: how to distribute the tensor, see `Distribution`
      /// \param p: partition of the tensor among the processes

      Tensor(const std::string& order, Coor<N> dim, DeviceHost dev, Distribution dist,
	     std::shared_ptr<detail::TensorPartition<N>> p)
	: order(order),
	  dim(dim),
	  ctx(detail::getContext(dev)),
	  p(p),
	  dist(dist),
	  from{},
	  size(dim),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow)),
	  scalar{1},
	  conjugate{false},
	  eg{false}
      {
	checkOrder();
	superbblas::Context ctx0 = *ctx;
	std::string s = repr();
	detail::log(1, "allocating " + s);
	T* ptr = superbblas::allocate<T>(p->localVolume(), *ctx);
	detail::log_mem();
	data = std::shared_ptr<T>(ptr, [=](const T* ptr) {
	  superbblas::deallocate(ptr, ctx0);
	  detail::log(1, "deallocated " + s);
	  detail::log_mem();
	});

	// TEMP!!!
	///set(detail::NaN<T>::get());
      }

    protected:
      /// Internal constructor, used by functions making slices, eg. `kvslice_from_size`
      /// \param order: dimension labels of the tensor
      /// \param from: coordinate of the first element in this view
      /// \param size: elements in each direction in this view

      Tensor(const Tensor& t, const std::string& order, Coor<N> from, Coor<N> size)
	: order(order),
	  dim(t.dim),
	  ctx(t.ctx),
	  data(t.data),
	  p(t.p),
	  dist(t.dist),
	  from(normalize_coor(from, t.dim)),
	  size(size),
	  strides(t.strides),
	  scalar{t.scalar},
	  conjugate{t.conjugate},
	  eg{t.eg}
      {
	checkOrder();
      }

      /// Internal constructor, used by `scale` and `conj`
      /// \param scalar: scalar factor of this view
      /// \param conjugate: whether the elements are implicitly conjugated

      Tensor(const Tensor& t, T scalar, bool conjugate)
	: order(t.order),
	  dim(t.dim),
	  ctx(t.ctx),
	  data(t.data),
	  p(t.p),
	  dist(t.dist),
	  from(t.from),
	  size(t.size),
	  strides(t.strides),
	  scalar{scalar},
	  conjugate{conjugate},
	  eg{false}
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
	return (from != Coor<N>{} || size != dim);
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

      T get(Coor<N> coor) const
      {
	if (ctx->plat != superbblas::CPU)
	  throw std::runtime_error(
	    "Unsupported to `get` elements from tensors not stored on the host");
	if (dist == OnEveryone)
	  throw std::runtime_error(
	    "Unsupported to `get` elements on a distributed tensor; change the distribution to "
	    "be supported on master, replicated among all nodes, or local");
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// coor[i] = coor[i] + from[i]
	for (unsigned int i = 0; i < N; ++i)
	  coor[i] = normalize_coor(normalize_coor(coor[i], size[i]) + from[i], dim[i]);

	return detail::cond_conj(conjugate,
				 data.get()[detail::coor2index<N>(coor, dim, strides)] * scalar);
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

      void set(Coor<N> coor, T v)
      {
	if (ctx->plat != superbblas::CPU)
	  throw std::runtime_error(
	    "Unsupported to `get` elements from tensors not stored on the host");
	if (dist == OnEveryone)
	  throw std::runtime_error(
	    "Unsupported to `set` elements on a distributed tensor; change the distribution to "
	    "be supported on master, replicated among all nodes, or local");
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// coor[i] = coor[i] + from[i]
	for (unsigned int i = 0; i < N; ++i)
	  coor[i] = normalize_coor(normalize_coor(coor[i], size[i]) + from[i], dim[i]);

	data.get()[detail::coor2index<N>(coor, dim, strides)] =
	  detail::cond_conj(conjugate, v) / scalar;
      }

      /// Modify the content this tensor with the result of a function on each element
      /// \param func: function () -> COMPLEX
      /// \param threaded: whether to run threaded

      template <typename Func>
      void fillWithCPUFuncNoArgs(Func func, bool threaded = true)
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	std::size_t vol = t.getLocal().volume();
	T* ptr = t.data.get();

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
      /// \param func: function (Coor<N>) -> T
      /// \param threaded: whether to run threaded

      template <typename Func>
      void fillCpuFunCoor(Func func, bool threaded = true)
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	using superbblas::detail::operator+;

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	std::size_t vol = t.getLocal().volume();
	T* ptr = t.data.get();
	/// Number of elements in each direction for the local part
	Coor<N> local_size = t.getLocal().size;
	/// Stride for the local volume
	Coor<N> local_stride = superbblas::detail::get_strides(local_size, superbblas::FastToSlow);
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
      /// \param func: function T -> Tr
      /// \param threaded: whether to run threaded

      template <typename Tr, typename Func>
      Tensor<N, Tr> transformWithCPUFun(Func func, bool threaded = true) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	auto r = t.template like_this<N, Tr>();
	assert(!r.isSubtensor() && !t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	T* tptr = t.data.get();
	Tr* rptr = r.data.get();

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
      /// \param func: function (Coor<N>, T) -> Tr
      /// \param threaded: whether to run threaded

      template <typename Tr, typename FuncWithCoor>
      Tensor<N, Tr> transformWithCPUFunWithCoor(FuncWithCoor func, bool threaded = true) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	using superbblas::detail::operator+;

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	auto r = t.template like_this<N, Tr>();
	assert(!r.isSubtensor() && !t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	T* tptr = t.data.get();
	Tr* rptr = r.data.get();
	/// Number of elements in each direction for the local part
	Coor<N> local_size = t.getLocal().size;
	/// Stride for the local volume
	Coor<N> local_stride = superbblas::detail::get_strides(local_size, superbblas::FastToSlow);
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

      /// Apply the function to each tensor element
      /// \param func: function T -> void
      /// \param threaded: whether to run threaded

      template <typename Func>
      void foreachWithCPUFun(Func func, bool threaded = true) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	auto t = isSubtensor() ? cloneOn(OnHost) : make_sure(none, OnHost);
	assert(!t.isSubtensor());
	std::size_t vol = t.getLocal().volume();
	T* tptr = t.data.get();

	if (threaded)
	{
#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
	  for (std::size_t i = 0; i < vol; ++i)
	    func(tptr[i]);
	}
	else
	{
	  for (std::size_t i = 0; i < vol; ++i)
	    func(tptr[i]);
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

      void set(T v)
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
	return Tensor<Nn, Tn>(new_order_, kvcoors<Nn>(new_order_, new_kvdim, 0, ThrowOnMissing),
			      new_dev.getSome(getDev()), new_dist.getSome(dist));
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
	std::map<char, int> new_kvdim = kvdim();
	for (const auto& it : kvsize)
	  new_kvdim[it.first] = it.second;
	std::string new_order_ =
	  detail::remove_dimensions(get_order_for_reorder(new_order, remaining_char), remove_dims);
	return Tensor<Nn, Tn>(new_order_, kvcoors<Nn>(new_order_, new_kvdim, 0, ThrowOnMissing),
			      new_dev.getSome(getDev()), new_dist.getSome(dist));
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

	Tensor<N, Tn> r = like_this<N, Tn>(none, {}, new_dev);
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
	  order, size, ctx, std::shared_ptr<T>(),
	  std::make_shared<detail::TensorPartition<N>>(p->get_subpartition(from, size)), dist, {},
	  size, T{1}, false /* not conjugate */, true /* is eg */);
      }

    private:
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

    public:
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
	Tensor<N, T> r = like_this(new_order1);
	if (is_eg())
	  r = r.make_eg();
	else
	  copyTo(r);
	return r;
      }

      /// Return whether the tensor has complex components although being stored with a non-complex type `T`

      bool isFakeReal() const
      {
	return order.find('.') != std::string::npos;
      }

      /// Check that the dimension labels are valid

      void checkOrder() const
      {
	// Check that all labels are different there are N
	detail::check_order<N>(order);

	/// Throw exception if this a fake real tensor but with a complex type `T`
	if (isFakeReal() && detail::is_complex<T>::value)
	  throw std::runtime_error("Invalid tensor: it is fake real and complex!");

	for (auto s : size)
	  if (s < 0)
	    std::runtime_error("Invalid tensor size: it should be positive");
      }

      /// Return a fake real view of this tensor

      template <typename U = T,
		typename std::enable_if<detail::is_complex<U>::value, bool>::type = true>
      Tensor<N + 1, typename U::value_type> toFakeReal() const
      {
	assert(!isFakeReal());

	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	std::string new_order = "." + order;
	Coor<N + 1> new_from = {0};
	std::copy_n(from.begin(), N, new_from.begin() + 1);
	Coor<N + 1> new_size = {2};
	std::copy_n(size.begin(), N, new_size.begin() + 1);
	Coor<N + 1> new_dim = {2};
	std::copy_n(dim.begin(), N, new_dim.begin() + 1);
	if (std::fabs(std::imag(scalar)) != 0)
	  throw std::runtime_error(
	    "Unsupported conversion to fake real tensors with an implicit complex scale");
	using new_T = typename T::value_type;
	new_T new_scalar = std::real(scalar);
	auto this_data = data;
	auto new_data =
	  std::shared_ptr<new_T>((new_T*)data.get(), [=](const new_T* ptr) { (void)this_data; });
	auto new_p = std::make_shared<detail::TensorPartition<N + 1>>(p->insert_dimension(0, 2));

	return Tensor<N + 1, new_T>(new_order, new_dim, ctx, new_data, new_p, dist, new_from,
				    new_size, new_scalar, conjugate, eg);
      }

      template <typename U = T,
		typename std::enable_if<!detail::is_complex<U>::value, bool>::type = true>
      Tensor<N - 1, std::complex<U>> toComplex(bool allow_cloning = true) const
      {
	assert(isFakeReal() && kvdim()['.'] == 2);

	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	std::size_t dot_pos = order.find('.');
	std::string new_order = detail::remove_coor(order, dot_pos);

	if (dot_pos != 0)
	{
	  if (allow_cloning)
	    return reorder("." + new_order).toComplex(false);
	  else
	    throw std::runtime_error("Not allow to create a new tensor in `toComplex`");
	}

	Coor<N - 1> new_from = detail::remove_coor(from, dot_pos);
	Coor<N - 1> new_size = detail::remove_coor(size, dot_pos);
	Coor<N - 1> new_dim = detail::remove_coor(dim, dot_pos);
	using new_T = std::complex<T>;
	new_T new_scalar = new_T{scalar};
	auto this_data = data;
	auto new_data =
	  std::shared_ptr<new_T>((new_T*)data.get(), [=](const new_T* ptr) { (void)this_data; });
	auto new_p = std::make_shared<detail::TensorPartition<N - 1>>(p->remove_dimension(dot_pos));

	return Tensor<N - 1, new_T>(new_order, new_dim, ctx, new_data, new_p, dist, new_from,
				    new_size, new_scalar, conjugate, eg);
      }

      template <typename U = T,
		typename std::enable_if<!detail::is_complex<U>::value, bool>::type = true>
      Tensor<N, U> toFakeReal() const
      {
	assert(isFakeReal());
	return *this;
      }

      template <typename U = T,
		typename std::enable_if<detail::is_complex<U>::value, bool>::type = true>
      Tensor<N, U> toComplex(bool allow_cloning = true) const
      {
	(void)allow_cloning;
	assert(!isFakeReal());
	return *this;
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

	return Tensor<Nout, T>(new_order, new_p->dim, ctx, data, new_p, dist,
			       detail::split_dimension(pos, from, d, From),
			       detail::split_dimension(pos, size, d, Size), scalar, conjugate, eg);
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

	return Tensor<Nout, T>(new_order, new_p->dim, ctx, data, new_p, dist,
			       detail::collapse_dimensions<Nout>(pos, from, dim, From),
			       detail::collapse_dimensions<Nout>(pos, size, dim, Size), scalar,
			       conjugate, eg);
      }

      /// Rearrange several dimensions into new ones
      /// \param m: maps from suborder of the current tensor to new orders
      /// \param allow_copy: whether to allow to return a reordered copy of the current tensor

      template <std::size_t Nout, typename std::enable_if<(N > 0), bool>::type = true>
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
	  return reorder(old_order).template reshape_dimensions<Nout>(m, new_dim, false);

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
	Coor<Nout> d = detail::reshape_dimensions<Nout>(
	  ncollapse, nsplit, dim, kvcoors<Nout>(new_order, new_dim0, 0, ThrowOnMissing), Size, dim);

	// Transform the partition
	auto new_p = std::make_shared<detail::TensorPartition<Nout>>(
	  p->template reshape_dimensions<Nout>(ncollapse, nsplit, d));

	return Tensor<Nout, T>(
	  new_order, new_p->dim, ctx, data, new_p, dist,
	  detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, d, From, from),
	  detail::reshape_dimensions<Nout>(ncollapse, nsplit, dim, d, Size, size), scalar,
	  conjugate, eg);
      }

      /// Copy/add this tensor into the given one
      /// NOTE: if this tensor or the given tensor is fake real, force both to be fake real

      template <std::size_t Nw, typename Tw,
		typename std::enable_if<
		  detail::is_complex<T>::value != detail::is_complex<Tw>::value, bool>::type = true>
      void doAction(Action action, Tensor<Nw, Tw> w) const
      {
	toFakeReal().doAction(action, w.toFakeReal());
      }

      /// Return the local support of this tensor
      Tensor<N, T> getLocal() const
      {
	// Compute the size of the intersection of the current view and the local support
	Coor<N> lfrom, lsize;
	superbblas::detail::intersection(p->localFrom(), p->localSize(), from, size, dim, lfrom,
					 lsize);

	// If the current process has no support, return the empty tensor
	if (superbblas::detail::volume(lsize) == 0)
	  return Tensor<N, T>{};

	using superbblas::detail::operator-;
	return Tensor<N, T>(order, p->localSize(), ctx, data,
			    std::make_shared<detail::TensorPartition<N>>(p->get_local_partition()),
			    Local, normalize_coor(from - p->localFrom(), dim), lsize, scalar,
			    conjugate, eg);
      }

      /// Set zero
      void set_zero()
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	T* ptr = this->data.get();
	MPI_Comm comm = (dist == OnMaster || dist == Local ? MPI_COMM_SELF : MPI_COMM_WORLD);
	if (dist != OnMaster || Layout::nodeNumber() == 0)
	  superbblas::copy<N, N>(T{0}, p->p.data(), 1, order.c_str(), from, size, dim,
				 (const T**)&ptr, nullptr, &*ctx, p->p.data(), 1, order.c_str(),
				 from, dim, &ptr, nullptr, &*ctx, comm, superbblas::FastToSlow,
				 superbblas::Copy);
      }

      /// Return whether the given tensor has the same distribution as this one
      template <typename Tv>
      bool is_distributed_like(Tensor<N, Tv> v) const
      {
	return from == v.from && size == v.size && dim == v.dim && p->p == v.p->p;
      }

      /// Return whether the given tensor has the same distribution as this one
      template <std::size_t Nv, typename Tv, typename std::enable_if<N != Nv, bool>::type = true>
      bool is_distributed_like(Tensor<Nv, Tv>) const
      {
	return false;
      }

      /// Return whether the given tensor has the same distribution as this one
      Tensor<N, float> create_mask() const
      {
	Tensor<N, float> m{order, dim, getDev(), dist, p};
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
      void doAction(Action action, Tensor<Nw, Tw> w, Maybe<Tensor<Nm, Tm>> m = none,
		    Maybe<Tensor<Nwm, Twm>> wm = none) const
      {
	if (is_eg() || w.is_eg() || (m.hasSome() && m.getSome().is_eg()) ||
	    (wm.hasSome() && wm.getSome().is_eg()))
	  throw std::runtime_error("Invalid operation from an example tensor");

	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > wsize[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	if (action == AddTo && w.scalar != Tw{1})
	  throw std::runtime_error(
	    "Not allowed to add to a tensor whose implicit scalar factor is not one");

	if (conjugate != w.conjugate)
	  throw std::runtime_error(
	    "Not allowed to copy or add tensor with different implicit conjugacy");

	if ((dist == Local && w.dist != Local) || (dist != Local && w.dist == Local))
	{
	  getLocal().doAction(action, w.getLocal(), m);
	  return;
	}

	// Compute masks
	float *m0ptr = nullptr, *m1ptr = nullptr;
	Tensor<N, float> m0;
	Tensor<Nw, float> m1;
	if (m.hasSome() || wm.hasSome())
	{
	  if (m.hasSome())
	  {
	    if (is_distributed_like(m.getSome()))
	    {
	      m0 = m.getSome();
	    }
	    else
	    {
	      m0 = create_mask();
	      m.getSome().copyTo(m0);
	    }
	  }

	  if (wm.hasSome())
	  {
	    if (w.is_distributed_like(wm.getSome()))
	    {
	      m1 = wm.getSome();
	    }
	    else
	    {
	      m1 = w.create_mask();
	      wm.getSome().copyTo(m1);
	    }
	  }

	  if (m.hasSome() && !wm.hasSome())
	    m0.copyTo(m1);
	  if (!m.hasSome() && wm.hasSome())
	    m1.copyTo(m0);

	  m0ptr = m0.data.get();
	  m1ptr = m1.data.get();
	}

	T* ptr = this->data.get();
	Tw* w_ptr = w.data.get();
	MPI_Comm comm =
	  ((dist == OnMaster && w.dist == OnMaster) || dist == Local ? MPI_COMM_SELF
								     : MPI_COMM_WORLD);
	if (dist != OnMaster || w.dist != OnMaster || Layout::nodeNumber() == 0)
	{
	  superbblas::copy<N, Nw>(
	    detail::safe_div<T>(scalar, w.scalar), p->p.data(), 1, order.c_str(), from, size, dim,
	    (const T**)&ptr, (const float**)&m0ptr, &*ctx, w.p->p.data(), 1, w.order.c_str(),
	    w.from, w.dim, &w_ptr, (const float**)&m1ptr, &*w.ctx, comm, superbblas::FastToSlow,
	    action == CopyTo ? superbblas::Copy : superbblas::Add);
	}
      }

      /// Copy this tensor into the given one
      template <std::size_t Nw, typename Tw>
      void copyTo(Tensor<Nw, Tw> w) const
      {
	doAction(CopyTo, w);
      }

      /// Copy this tensor into the given one but only the elements where the mask is nonzero
      template <std::size_t Nw, typename Tw, std::size_t Nm, typename Tm, std::size_t Nwm,
		typename Twm>
      void copyToWithMask(Tensor<Nw, Tw> w, Tensor<Nm, Tm> m, Tensor<Nwm, Twm> wm) const
      {
	doAction<Nw, Tw, Nm, Tm, Nwm, Twm>(CopyTo, w, m, wm);
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

	Tensor<N, T> r(order, dim, getDev(), OnEveryone, new_p);
	copyTo(r);
	return r;
      }

      // Contract the dimensions with the same label in `v` and `w` than do not appear on `this` tensor.
      template <std::size_t Nv, std::size_t Nw>
      void contract(Tensor<Nv, T> v, const remap& mv, Conjugation conjv, Tensor<Nw, T> w,
		    const remap& mw, Conjugation conjw, const remap& mr = {}, T beta = T{0})
      {
	if (is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	// If either v or w is on OnDevice, force both to be on device
	if (v.ctx->plat != w.ctx->plat)
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
	  Tensor<N, T> aux =
	    std::norm(beta) == 0 ? like_this(none, {}, v.getDev()) : cloneOn(v.getDev());
	  aux.contract(v, mv, conjv, w, mw, conjw, mr, beta);
	  aux.copyTo(*this);
	  return;
	}

	if ((v.dist == Local) != (w.dist == Local) || (w.dist == Local) != (dist == Local))
	  throw std::runtime_error(
	    "One of the contracted tensors or the output tensor is local and others are not!");

	if ((v.dist == OnMaster && w.dist == OnEveryone) ||
	    (v.dist == OnEveryone && w.dist == OnMaster))
	  throw std::runtime_error("Incompatible layout for contractions: one of the tensors is on "
				   "the master node and the other is distributed");

	// Avoid replicating the same computation in all nodes
	// TODO: check whether superbblas does this already
	if ((v.dist == OnMaster || v.dist == OnEveryoneReplicated) &&
	    (w.dist == OnMaster || w.dist == OnEveryoneReplicated))
	{
	  v = v.make_sure(none, none, OnMaster);
	  w = w.make_sure(none, none, OnMaster);
	}

	if (v.dist == OnEveryone && w.dist == OnEveryoneReplicated)
	  w = w.make_suitable_for_contraction(v);

	if (v.dist == OnEveryoneReplicated && w.dist == OnEveryone)
	  v = v.make_suitable_for_contraction(w);

	T* v_ptr = v.data.get();
	T* w_ptr = w.data.get();
	T* ptr = this->data.get();
	std::string orderv_ = detail::update_order_and_check<Nv>(v.order, mv);
	std::string orderw_ = detail::update_order_and_check<Nw>(w.order, mw);
	std::string order_ = detail::update_order_and_check<N>(order, mr);
	bool conjv_ = (((conjv == Conjugate) xor v.conjugate) xor conjugate);
	bool conjw_ = (((conjw == Conjugate) xor w.conjugate) xor conjugate);
	superbblas::contraction<Nv, Nw, N>(
	  detail::cond_conj(conjv_, v.scalar) * detail::cond_conj(conjw_, w.scalar) / scalar, //
	  v.p->p.data(), v.dim, 1, orderv_.c_str(), conjv_, (const T**)&v_ptr, &*v.ctx,	      //
	  w.p->p.data(), w.dim, 1, orderw_.c_str(), conjw_, (const T**)&w_ptr, &*w.ctx,	      //
	  detail::cond_conj(conjugate, beta), p->p.data(), dim, 1, order_.c_str(), &ptr, &*ctx,
	  MPI_COMM_WORLD, superbblas::FastToSlow);
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
	if (v.ctx->plat != w.ctx->plat)
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
	  Tensor<N, T> aux = like_this(none, {}, v.getDev());
	  aux.cholInv(v, order_rows, order_cols, w);
	  aux.copyTo(*this);
	  return;
	}

	// v is going to be modified and is reference, make a clone
	if (v.data.use_count() > 1)
	  v = v.clone();

	if ((v.dist == Local) != (w.dist == Local) || (w.dist == Local) != (dist == Local))
	  throw std::runtime_error(
	    "One of the contracted tensors or the output tensor is local and others are not!");

	if (v.dist == OnEveryone && w.dist == OnEveryoneReplicated)
	  w = w.make_suitable_for_contraction(v);

	if (v.dist == OnEveryoneReplicated && w.dist == OnEveryone)
	  v = v.make_suitable_for_contraction(w);

	if (std::fabs(std::imag(v.scalar)) != 0 || std::real(v.scalar) < 0)
	  throw std::runtime_error("cholInv: unsupported a negative or imaginary scale");

	T* v_ptr = v.data.get();
	T* w_ptr = w.data.get();
	T* ptr = this->data.get();
	superbblas::cholesky<Nv>(v.p->p.data(), v.dim, 1, v.order.c_str(), &v_ptr,
				 order_rows.c_str(), order_cols.c_str(), &*v.ctx, MPI_COMM_WORLD,
				 superbblas::FastToSlow);
	superbblas::trsm<Nv, Nw, N>(
	  w.scalar / std::sqrt(v.scalar) / scalar, //
	  v.p->p.data(), v.dim, 1, v.order.c_str(), (const T**)&v_ptr, order_rows.c_str(),
	  order_cols.c_str(),
	  &*v.ctx,								//
	  w.p->p.data(), w.dim, 1, w.order.c_str(), (const T**)&w_ptr, &*w.ctx, //
	  p->p.data(), dim, 1, order.c_str(), &ptr, &*ctx, MPI_COMM_WORLD, superbblas::FastToSlow);
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
	if (v.ctx->plat != w.ctx->plat)
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
	  Tensor<N, T> aux = like_this(none, {}, v.getDev());
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

	if (v.dist == OnEveryone && w.dist == OnEveryoneReplicated)
	  w = w.make_suitable_for_contraction(v);

	if (v.dist == OnEveryoneReplicated && w.dist == OnEveryone)
	  v = v.make_suitable_for_contraction(w);

	T* v_ptr = v.data.get();
	T* w_ptr = w.data.get();
	T* ptr = this->data.get();
	superbblas::gesm<Nv, Nw, N>(
	  w.scalar / v.scalar / scalar, //
	  v.p->p.data(), v.dim, 1, v.order.c_str(), (const T**)&v_ptr, order_rows.c_str(),
	  order_cols.c_str(),
	  &*v.ctx,								//
	  w.p->p.data(), w.dim, 1, w.order.c_str(), (const T**)&w_ptr, &*w.ctx, //
	  p->p.data(), dim, 1, order.c_str(), &ptr, &*ctx, MPI_COMM_WORLD, superbblas::FastToSlow);
      }

      /// Return a view of this tensor where the elements are scaled by the given argument
      /// \param s: scaling factor
      /// \return: a new view (it doesn't create a copy of the tensor)

      Tensor<N, T> scale(T s) const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");
	return Tensor<N, T>(*this, scalar * detail::cond_conj(conjugate, s), conjugate);
      }

      /// Return a view of this tensor where the elements are conjuated
      /// \return: a new view (it doesn't create a copy of the tensor)

      Tensor<N, T> conj() const
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");
	return Tensor<N, T>(*this, scalar, !conjugate);
      }

      void release()
      {
	dim = {};
	data.reset();
	p.reset();
	ctx.reset();
	from = {};
	size = {};
	strides = {};
	scalar = T{0};
	conjugate = false;
	eg = false;
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
	  Tensor<N, Tn> r = like_this(new_order, {}, new_dev, new_dist);
	  if (is_eg())
	  {
	    r = r.make_eg();
	  }
	  else
	  {
	    r.conjugate = conjugate;
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
	  copyTo(r);
	}
	return r;
      }

      /// Extend the support of each dimension by the given amount in each direction
      /// \param m: amount to extend the support for each process
      /// \return a new tensor with the extension

      Tensor<N, T> extend_support(const std::map<char, int>& m) const
      {
	Tensor<N, T> r{
	  order, dim, getDev(), dist,
	  std::make_shared<detail::TensorPartition<N>>(p->extend_support(kvcoors<N>(order, m, 0)))};
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
	return (ctx->plat != superbblas::CPU ? OnDefaultDevice : OnHost);
#  else
	return OnDefaultDevice;
#  endif
      }

      void binaryRead(BinaryReader& bin)
      {
	if (is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");
	if (ctx->plat != superbblas::CPU)
	  throw std::runtime_error("Only supported to read on `OnHost` tensors");
	if (dist != OnMaster)
	  throw std::runtime_error("Only supported to read on `OnMaster` tensors");
	if (!isContiguous())
	  throw std::runtime_error("Only supported contiguous views in memory");
	if (scalar != T{1} || conjugate)
	  throw std::runtime_error(
	    "Not allowed for tensor with a scale not being one or implicitly conjugated");

	// Only on primary node read the data
	std::size_t vol = volume();
	std::size_t disp = detail::coor2index<N>(from, dim, strides);
	std::size_t word_size = sizeof(typename detail::WordType<T>::type);
	bin.readArrayPrimaryNode((char*)&data.get()[disp], word_size, sizeof(T) / word_size * vol);
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
	if (scalar != T{1} || isSubtensor() || ctx->plat != superbblas::CPU)
	{
	  cloneOn(OnHost).binaryWrite(bin);
	  return;
	}

	// Write the local data
	std::size_t vol = p->localVolume();
	std::size_t word_size = sizeof(typename detail::WordType<T>::type);
	bin.writeArrayPrimaryNode((char*)data.get(), word_size, sizeof(T) / word_size * vol);
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
	  ss << "% " << repr(data.get()) << std::endl;
	  ss << "% dist=" << p->p << std::endl;
	  ss << name << "=reshape([";
	  assert(!t_host.isSubtensor());
	  std::size_t vol = volume();
	  for (std::size_t i = 0; i < vol; ++i)
	  {
	    //using detail::repr::operator<<;
	    ss << " ";
	    detail::repr::operator<<(ss, t_host.data.get()[i]);
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
    void latticeCopyToWithMask(Tensor<Nv, Tv> v, Tensor<Nw, Tw> w, char label_mu,
			       const std::vector<Coor<Nd>>& disps,
			       const std::map<char, int>& real_dims, Tensor<Nm, Tm> mask_even,
			       Tensor<Nm, Tm> mask_odd)
    {
      // Shortcuts
      if (disps.size() == 0)
	return;

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

    /// Contract some dimension of the given tensors
    /// \param v: one tensor to contract
    /// \param w: the other tensor to contract
    /// \param labels_to_contract: labels dimensions to contract from `v` and `w`
    /// \param action: either to copy or add to the given output tensor if given
    /// \param r: optional given tensor where to put the resulting contraction
    /// \param mr: map from the given `r` to the labels of the contraction
    /// \param beta: scale on `r` if the `action` in `AddTo`
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
    Tensor<Nr, T> contract(Tensor<Nv, T> v, Tensor<Nw, T> w, const std::string& labels_to_contract,
			   Maybe<Action> action = none, Tensor<Nr, T> r = Tensor<Nr, T>{},
			   const remap& mr = {}, T beta = T{1})
    {
      if (action.hasSome() != (bool)r)
	throw std::runtime_error("Invalid default value");

      // Compute the labels of the output tensor: v.order + w.order - labels_to_contract
      std::string rorder = detail::union_dimensions(v.order, w.order, labels_to_contract);
      if (Nr != rorder.size())
	throw std::runtime_error(
	  "contract: The dimension of the output tensor does not match the template argument");
      if ((bool)r && detail::union_dimensions(rorder, r.order) != rorder)
	throw std::runtime_error("contract: The given output tensor has an unexpected ordering");

      // If the output tensor is not given create a new one
      Tensor<Nr, T> r0;
      if (!r)
      {
	r0 = v.template like_this<Nr>(rorder, w.kvdim());
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
    Tensor<Nr, typename detail::real_type<T>::type>
    norm(Tensor<Nv, T> v, Maybe<std::string> order_t = none, Maybe<std::string> order_rows = none)
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
      auto r = contract<Nr>(v.conj(), v, rorder).make_sure(torder, OnHost, OnEveryoneReplicated);

      // Do the square root and return the result
      using Treal = typename detail::real_type<T>::type;
      return r.template transformWithCPUFun<Treal>(
	[](const T& t) { return std::sqrt(std::real(t)); });
    }

    /// Compute the maximum for a small tensor
    /// \param v: tensor

    template <std::size_t N, typename T>
    T max(Tensor<N, T> v)
    {
      v = v.make_sure(none, OnHost, OnEveryoneReplicated);
      if (v.isSubtensor())
	v = v.clone();
      T r = T{0};
      T* p = v.data.get();
      for (unsigned int i = 0, vol = v.volume(); i < vol; ++i)
	r = std::max(r, p[i]);
      return r;
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
    Tensor<Nr, T> cholInv(Tensor<Nv, T> v, const std::string& order_rows,
			  const std::string& order_cols, Tensor<Nw, T> w,
			  const std::string& labels_to_contract, Maybe<Action> action = none,
			  Maybe<Tensor<Nr, T>> r = none)
    {
      if (action.hasSome() != r.hasSome())
	throw std::runtime_error("Invalid default value");

      // Compute the labels of the output tensor: v.order + w.order - labels_to_contract
      std::string rorder = detail::union_dimensions(v.order, w.order, labels_to_contract);
      if (Nr != rorder.size())
	throw std::runtime_error(
	  "cholInv: The dimension of the output tensor does not match the template argument");
      if (r && detail::union_dimensions(rorder, r.getSome().order) != rorder)
	throw std::runtime_error("cholInv: The given output tensor has an unexpected ordering");

      // If the output tensor is not given create a new one
      Tensor<Nr, T> r0;
      if (!r)
      {
	r0 = v.template like_this<Nr>(rorder, w.kvdim());
      }
      else
      {
	r0 = r.getSome();
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
    Tensor<Nr, T> solve(Tensor<Nv, T> v, const std::string& order_rows,
			const std::string& order_cols, Tensor<Nw, T> w,
			const std::string& labels_to_contract, Maybe<Action> action = none,
			Maybe<Tensor<Nr, T>> r = none)
    {
      if (action.hasSome() != r.hasSome())
	throw std::runtime_error("solve: Invalid default value");

      // Compute the labels of the output tensor: v.order + w.order - labels_to_contract
      std::string rorder = detail::union_dimensions(v.order, w.order, labels_to_contract);
      if (Nr != rorder.size())
	throw std::runtime_error(
	  "solve: The dimension of the output tensor does not match the template argument");
      if (r && detail::union_dimensions(rorder, r.getSome().order) != rorder)
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
	r0 = r.getSome();
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
	  err = std::max(err, (double)rnorms.data.get()[i] / wnorms.data.get()[i]);
	QDPIO::cout << "solve error: " << detail::tostr(err) << std::endl;
	auto eps = std::sqrt(std::numeric_limits<typename detail::real_type<T>::type>::epsilon());
	if (err > eps)
	  throw std::runtime_error(std::string("solve: too much error in dense solution, ") +
				   detail::tostr(err));
      }

      return r0;
    }

    template <typename T>
    void* getQDPPtr(const T& t)
    {
#  if defined(QDP_IS_QDPJIT) && defined(SUPERBBLAS_USE_GPU)
      std::vector<QDPCache::ArgKey> v(1, t.getId());
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
      return Tensor<Nd + 2, Complex>("cxyztX", latticeSize<Nd + 2>("cxyztX"), OnHost, OnEveryone,
				     std::shared_ptr<Complex>(v_ptr, [](Complex*) {}));
    }

#  if !defined(QDP_IS_QDPJIT) || !defined(SUPERBBLAS_USE_GPU)
    inline Tensor<Nd + 3, Complex> asTensorView(const LatticeFermion& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 3, Complex>("csxyztX", latticeSize<Nd + 3>("csxyztX"), OnHost, OnEveryone,
				     std::shared_ptr<Complex>(v_ptr, [](Complex*) {}));
    }
#  else
    inline Tensor<Nd + 4, REAL> asTensorView(const LatticeFermion& v)
    {
      REAL* v_ptr = reinterpret_cast<REAL*>(getQDPPtr(v));
      return Tensor<Nd + 4, REAL>("xyztXsc.", latticeSize<Nd + 4>("xyztXsc."), OnDefaultDevice,
				  OnEveryone, std::shared_ptr<REAL>(v_ptr, [](REAL*) {}));
    }
#  endif

#  if !defined(QDP_IS_QDPJIT) || !defined(SUPERBBLAS_USE_GPU)
    inline Tensor<Nd + 1, Complex> asTensorView(const LatticeComplex& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 1, Complex>("xyztX", latticeSize<Nd + 1>("xyztX"), OnHost, OnEveryone,
				     std::shared_ptr<Complex>(v_ptr, [](Complex*) {}));
    }
#  else
    inline Tensor<Nd + 2, REAL> asTensorView(const LatticeComplex& v)
    {
      REAL* v_ptr = reinterpret_cast<REAL*>(getQDPPtr(v));
      return Tensor<Nd + 2, REAL>("xyztX.", latticeSize<Nd + 2>("xyztX."), OnDefaultDevice,
				  OnEveryone, std::shared_ptr<REAL>(v_ptr, [](REAL*) {}));
    }
#  endif

#  if !defined(QDP_IS_QDPJIT) || !defined(SUPERBBLAS_USE_GPU)
    inline Tensor<Nd + 3, Complex> asTensorView(const LatticeColorMatrix& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 3, Complex>("jixyztX",
				     latticeSize<Nd + 3>("jixyztX", {{'i', Nc}, {'j', Nc}}), OnHost,
				     OnEveryone, std::shared_ptr<Complex>(v_ptr, [](Complex*) {}));
    }
#  else
    inline Tensor<Nd + 4, REAL> asTensorView(const LatticeColorMatrix& v)
    {
      REAL* v_ptr = reinterpret_cast<REAL*>(getQDPPtr(v));
      return Tensor<Nd + 4, REAL>(
	"xyztXji.", latticeSize<Nd + 4>("xyztXji.", {{'i', Nc}, {'j', Nc}}), OnDefaultDevice,
	OnEveryone, std::shared_ptr<REAL>(v_ptr, [](REAL*) {}));
    }
#  endif

    inline Tensor<Nd + 4, Complex> asTensorView(const LatticeColorVectorSpinMatrix& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 4, Complex>(
	"cjixyztX", latticeSize<Nd + 4>("cjixyztX", {{'i', Ns}, {'j', Ns}}), OnHost, OnEveryone,
	std::shared_ptr<Complex>(v_ptr, [](Complex*) {}));
    }

    template <typename COMPLEX>
    Tensor<1, COMPLEX> asTensorView(std::vector<COMPLEX>& v,
				    Distribution dist = OnEveryoneReplicated)
    {
      return Tensor<1, COMPLEX>("i", Coor<1>{Index(v.size())}, OnHost, dist,
				std::shared_ptr<COMPLEX>(v.data(), [](COMPLEX*) {}));
    }

    inline Tensor<2, Complex> asTensorView(SpinMatrix& smat)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(smat.getF());
      return Tensor<2, Complex>("ji", Coor<2>{Ns, Ns}, OnHost, OnEveryoneReplicated,
				std::shared_ptr<Complex>(v_ptr, [](Complex*) {}));
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

    /// Class for operating sparse tensors
    /// \tparam ND: number of domain dimensions
    /// \tparam NI: number of image dimensions
    /// \tparam T: datatype

    template <std::size_t ND, std::size_t NI, typename T>
    struct SpTensor {
      static_assert(superbblas::supported_type<T>::value, "Not supported type");

    public:
      Tensor<ND, T> d;		   ///< Tensor example for the domain
      Tensor<NI, T> i;		   ///< Tensor example for the image
      Coor<ND> blkd;		   ///< blocking for the domain
      Coor<NI> blki;		   ///< blocking for the image
      Tensor<NI, int> ii;	   ///< Number of blocks in each row
      Tensor<NI + 2, int> jj;	   ///< Coordinate of the first element on each block
      Tensor<NI + ND + 1, T> data; ///< Nonzero values
      std::shared_ptr<superbblas::BSR_handle> handle; ///< suparbblas sparse tensor handle
      T scalar;					      ///< Scalar factor of the tensor
      const bool isImgFastInBlock;		      ///< whether the BSR blocks are in row-major
      const unsigned int nblockd;		      ///< Number of blocked domain dimensions
      const unsigned int nblocki;		      ///< Number of blocked image dimensions

      /// Low-level constructor
      SpTensor(Tensor<ND, T> d, Tensor<NI, T> i, Coor<ND> blkd, Coor<NI> blki, Tensor<NI, int> ii,
	       Tensor<NI + 2, int> jj, Tensor<NI + ND + 1, T> data, T scalar, bool isImgFastInBlock,
	       unsigned int nblockd, unsigned int nblocki)
	: d(d.make_eg()),
	  i(i.make_eg()),
	  blkd(blkd),
	  blki(blki),
	  ii(ii),
	  jj(jj),
	  data(data),
	  scalar(scalar),
	  isImgFastInBlock(isImgFastInBlock),
	  nblockd(nblockd),
	  nblocki(nblocki)
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
	if (data.data)
	  ss << "data:" << data.data << ", ";
	std::size_t sizemb = (ii.getLocal().volume() * sizeof(int) + //
			      jj.getLocal().volume() * sizeof(int) + //
			      data.getLocal().volume() * sizeof(T)) /
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
	       unsigned int num_neighbors, bool isImgFastInBlock = false)
	: d{d.make_eg()},
	  i{i.make_eg()},
	  scalar{T{1}},
	  isImgFastInBlock{isImgFastInBlock},
	  nblockd(nblockd),
	  nblocki(nblocki)
      {
	// Check that the examples are on the same device
	if (d.getDev() != i.getDev())
	  throw std::runtime_error("Please give example vectors on the same device");

	// Check that `d` and `i` are not subtensors
	if (this->d.isSubtensor() || this->i.isSubtensor())
	  throw std::runtime_error("unsupported subtensors for domain/image distributions");

	// Check that the domain and image labels are different and do not contain `u` or `~`
	detail::check_order<NI + ND + 2>(i.order + d.order + std::string("u~"));

	// Check the blocking
	blkd = kvcoors<ND>(d.order, d.kvdim());
	for (unsigned int i = nblockd; i < ND; ++i)
	  blkd[i] = 1;
	blki = kvcoors<NI>(i.order, i.kvdim());
	for (unsigned int i = nblocki; i < NI; ++i)
	  blki[i] = 1;

	// Create the tensor containing the number of neighbors for each blocking
	std::map<char, int> nonblki;
	for (unsigned int j = 0; j < NI; ++j)
	  nonblki[i.order[j]] = i.size[j] / blki[j];
	ii = i.template like_this<NI, int>(none, nonblki);
	ii.set(num_neighbors);

	// Create the tensor containing the domain coordinates of the first nonzero in each block
	jj = ii.template like_this<NI + 2, int>(
	  "~u%", '%', "", std::map<char, int>{{'~', (int)ND}, {'u', (int)num_neighbors}});

	// Compute the data dimensions as
	//   image_blocked_dims + domain_dims + u + image_nonblocked_dims, for isImgFastInBlock
	//   domain_blocked_dims + image_blocked_dims + u + image_nonblockd_dims otherwise
	std::map<char, int> data_dims = i.kvdim();
	for (unsigned int i = 0; i < ND; ++i)
	  data_dims[d.order[i]] = blkd[i];
	data_dims['u'] = num_neighbors;
	std::string data_order =
	  (isImgFastInBlock ? std::string(i.order.begin(), i.order.begin() + nblocki) + d.order
			    : std::string(d.order.begin(), d.order.begin() + nblockd) +
				std::string(i.order.begin(), i.order.begin() + nblocki) +
				std::string(d.order.begin() + nblockd, d.order.end())) +
	  std::string("u") + std::string(i.order.begin() + nblocki, i.order.end());
	data = i.template like_this<NI + ND + 1, T>(data_order, data_dims);
      }

      /// Empty constructor

      SpTensor() : blki{}, blkd{}, scalar{0}, isImgFastInBlock{false}, nblockd{0}, nblocki{0}
      {
      }

      /// Construct the sparse operator
      void construct()
      {
	// Superbblas needs the column coordinates to be local
	// Remove the local domain coordinates to jj
	const auto localFrom = d.p->localFrom();
	const auto domDim = d.dim;
	auto localjj =
	  jj.template transformWithCPUFunWithCoor<int>([&](const Coor<NI + 2>& c, const int& t) {
	    return (t - localFrom[c[0]] + domDim[c[0]]) % domDim[c[0]];
	  });

	int* iiptr = ii.data.get();
	Coor<ND>* jjptr = (Coor<ND>*)localjj.data.get();
	// NOTE: despite jj being a vector of `int`, superbblas will use jj as a vector of Coor<NI>, so check that the alignment
	if (localjj.getLocal().volume() > 0 &&
	    superbblas::detail::align(alignof(Coor<NI>), sizeof(int), jjptr, sizeof(int)) ==
	      nullptr)
	  throw std::runtime_error("Ups! Look into this");
	const T* ptr = data.data.get();
	superbblas::BSR_handle* bsr = nullptr;
	superbblas::create_bsr<ND, NI, T>(i.p->p.data(), i.dim, d.p->p.data(), d.dim, 1, blki, blkd,
					  isImgFastInBlock, &iiptr, &jjptr, &ptr, &*data.ctx,
					  MPI_COMM_WORLD, superbblas::FastToSlow, &bsr);
	handle = std::shared_ptr<superbblas::BSR_handle>(
	  bsr, [=](superbblas::BSR_handle* bsr) { destroy_bsr(bsr); });
      }

      /// Split a dimension into another dimensions
      /// \param dom_dim_label: dominion dimension to split
      /// \param dom_new_labels: the labels of the new dominion dimensions
      /// \param dom_step: length of the first label in `dom_new_labels`
      /// \param img_dim_label: image dimension to split
      /// \param img_new_labels: the labels of the image new dimensions
      /// \param img_step: length of the first label in `img_new_labels`

      SpTensor<ND + 1, NI + 1, T> split_dimension(char dom_dim_label, const std::string& dom_new_labels,
						  Index dom_step, char img_dim_label,
						  const std::string& img_new_labels,
						  Index img_step) const
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

	if (blkd[d_pos] > 1 && blkd[d_pos] % dom_step != 0)
	  throw std::runtime_error(
	    "split_dimension: invalid `dom_step`, it should divide the block size");
	if (blki[i_pos] > 1 && blki[i_pos] % img_step != 0)
	  throw std::runtime_error(
	    "split_dimension: invalid `img_step`, it should divide the block size");

	auto new_d = d.split_dimension(dom_dim_label, dom_new_labels, dom_step);
	auto new_i = i.split_dimension(img_dim_label, img_new_labels, img_step);

	int new_blkd_pos = blkd[d_pos] == 1 ? 1 : dom_step;
	auto new_blkd = detail::insert_coor(blkd, d_pos, new_blkd_pos);
	new_blkd[d_pos + 1] /= new_blkd_pos;
	int new_blki_pos = blki[i_pos] == 1 ? 1 : img_step;
	auto new_blki = detail::insert_coor(blki, i_pos, new_blki_pos);
	new_blki[i_pos + 1] /= new_blki_pos;

	auto new_ii = ii.split_dimension(img_dim_label, img_new_labels, img_step);
	auto new_jj_dim = new_ii.kvdim();
	new_jj_dim['~'] = (int)ND + 1;
	auto new_jj = jj.template like_this<NI + 3>(std::string("~u") + new_ii.order, new_jj_dim);
	{
	  auto local_jj = jj.make_sure(none, OnHost).getLocal();
	  auto local_new_jj = new_jj.make_sure(none, OnHost).getLocal();
	  assert(local_jj.volume() * (ND + 1) == local_new_jj.volume() * ND);
	  const int *p = local_jj.data.get();
	  int* new_p = local_new_jj.data.get();
	  auto new_dom_dim = detail::insert_coor(d.size, d_pos, dom_step);
	  for (std::size_t i = 0, i1 = local_jj.volume() / ND; i < i1; ++i)
	  {
	    Coor<ND> c;
	    std::copy_n(p + ND * i, ND, c.begin());
	    Coor<ND + 1> new_c = detail::split_dimension(d_pos, c, new_dom_dim, detail::From);
	    std::copy_n(new_c.begin(), ND + 1, new_p + (ND + 1) * i);
	  }
	  local_new_jj.copyTo(new_jj.getLocal());
	}

	auto new_data = data.split_dimension(dom_dim_label, dom_new_labels, dom_step)
			  .split_dimension(img_dim_label, img_new_labels, img_step);

	SpTensor<ND + 1, NI + 1, T> r(new_d, new_i, new_blkd, new_blki, new_ii, new_jj, new_data,
				      scalar, isImgFastInBlock, nblockd + (d_pos < nblockd ? 1 : 0),
				      nblocki + (i_pos < nblocki ? 1 : 0));
	if (is_constructed())
	  r.construct();

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

	auto new_d = d.reorder(new_dom_order, remaining_char);
	auto new_i = i.reorder(new_img_order, remaining_char);

	Coor<ND> d_perm = superbblas::detail::find_permutation(d.order, new_dom_order);
	Coor<NI> i_perm = superbblas::detail::find_permutation(i.order, new_img_order);
	auto new_blkd = superbblas::detail::reorder_coor(blkd, d_perm);
	auto new_blki = superbblas::detail::reorder_coor(blki, i_perm);

	// Check the blocking
	for (unsigned int i = 0; i < ND; ++i)
	  if (i >= nblockd && new_blkd[i] != new_d.size[i])
	    throw std::runtime_error("reorder: invalid domain reordering, it is mixing blocking "
				     "and nonblocking dimensions");
	for (unsigned int i = 0; i < NI; ++i)
	  if (i >= nblocki && new_blki[i] != new_i.size[i])
	    throw std::runtime_error("reorder: invalid image reordering, it is mixing blocking "
				     "and nonblocking dimensions");

	auto new_ii = ii.reorder(new_img_order, remaining_char);
	auto new_jj = jj.reorder(std::string("~u") + new_img_order, remaining_char);
	if (new_jj.order != jj.order)
	{
	  auto local_new_jj = new_jj.make_sure(none, OnHost).getLocal();
	  int* new_p = local_new_jj.data.get();
	  for (std::size_t i = 0, i1 = local_new_jj.volume() / ND; i < i1; ++i)
	  {
	    Coor<ND> c;
	    std::copy_n(new_p + ND * i, ND, c.begin());
	    Coor<ND> new_c = superbblas::detail::reorder_coor(c, d_perm);
	    std::copy_n(new_c.begin(), new_p + ND * i, ND);
	  }
	  local_new_jj.copyTo(new_jj.getLocal());
	}

	std::string data_order =
	  (isImgFastInBlock ? std::string(i.order.begin(), i.order.begin() + nblocki) + d.order
			    : std::string(d.order.begin(), d.order.begin() + nblockd) +
				std::string(i.order.begin(), i.order.begin() + nblocki) +
				std::string(d.order.begin() + nblockd, d.order.end())) +
	  std::string("u") + std::string(i.order.begin() + nblocki, i.order.end());
	auto new_data = data.reorder(data_order);

	SpTensor<ND + 1, NI + 1, T> r(new_d, new_i, new_blkd, new_blki, new_ii, new_jj, new_data,
				      scalar, isImgFastInBlock, nblockd, nblocki);
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
      template <std::size_t Nv, std::size_t Nw>
      void contractWith(Tensor<Nv, T> v, const remap& mv, Tensor<Nw, T> w,
			const remap& mw = {}) const
      {
	if (data.is_eg() || v.is_eg() || w.is_eg())
	  throw std::runtime_error("Invalid operation from an example tensor");

	if (!is_constructed())
	  throw std::runtime_error("invalid operation on an not constructed tensor");

	// If either this tensor or v are on OnDevice, force both to be on the same device as this tensor.
	if (v.ctx->plat != data.ctx->plat)
	{
	  v = v.cloneOn(getDev());
	}

	if (getDev() != w.getDev())
	{
	  Tensor<Nw, T> aux = w.like_this(none, {}, getDev());
	  contractWith(v, mv, aux, mw);
	  aux.copyTo(w);
	  return;
	}

	// Check unsupported distributions for contraction
	if ((v.dist == Local) != (w.dist == Local) || (w.dist == Local) != (data.dist == Local))
	  throw std::runtime_error(
	    "One of the contracted tensors or the output tensor is local and others are not!");

	if ((v.dist == OnMaster && w.dist == OnEveryone) ||
	    (v.dist == OnEveryone && w.dist == OnMaster))
	  throw std::runtime_error("Incompatible layout for contractions: one of the tensors is on "
				   "the master node and the other is distributed");

	// We don't support conjugacy for now
	if (v.conjugate || w.conjugate)
	  throw std::runtime_error("contractWith: unsupported implicit conjugacy");

	T* v_ptr = v.data.get();
	T* w_ptr = w.data.get();
	std::string orderv = detail::update_order_and_check<Nv>(v.order, mv);
	std::string orderw = detail::update_order_and_check<Nw>(w.order, mw);
	superbblas::bsr_krylov<ND, NI, Nv, Nw, T>(
	  scalar * v.scalar / w.scalar, handle.get(), i.order.c_str(), d.order.c_str(), //
	  v.p->p.data(), 1, orderv.c_str(), v.from, v.size, v.dim, (const T**)&v_ptr,	//
	  w.p->p.data(), orderw.c_str(), w.from, w.size, w.dim,
	  0 /* krylov power dim (none for now) */,
	  (T**)&w_ptr, //
	  &*data.ctx, MPI_COMM_WORLD, superbblas::FastToSlow);
      }
    };

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
	  dim{},
	  sparsity(Dense),
	  ctx{},
	  from{},
	  size{},
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
	  from{},
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
	  superbblas::PartitionItem<N> p{Coor<N>{}, dim};
	  superbblas::append_blocks<N, T>(&p, 1, dim, stoh, MPI_COMM_WORLD, superbblas::FastToSlow);
	}
      }

      // Open storage construct
      StorageTensor(const std::string& filename, bool read_order = true,
		    const Maybe<std::string>& order_tag = none)
	: filename(filename), sparsity(Sparse), from{}, scalar{1}
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
	superbblas::open_storage<N, T>(filename.c_str(), MPI_COMM_WORLD, &stoh);
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
	dim = {};
	ctx.reset();
	from = {};
	size = {};
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
      void copyFrom(Tensor<Nw, Tw> w) const
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

	Tw* w_ptr = w.data.get();
	superbblas::save<Nw, N, Tw, T>(detail::safe_div<Tw>(w.scalar, scalar), w.p->p.data(), 1,
				       w.order.c_str(), w.from, w.size, w.dim, (const Tw**)&w_ptr,
				       &*w.ctx, order.c_str(), from, ctx.get(), comm,
				       superbblas::FastToSlow);
      }

      /// Load content from the storage into the given tensor
      template <std::size_t Nw, typename Tw,
		typename std::enable_if<
		  detail::is_complex<T>::value == detail::is_complex<Tw>::value, bool>::type = true>
      void copyTo(Tensor<Nw, Tw> w) const
      {
	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > wsize[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	Tw* w_ptr = w.data.get();
	MPI_Comm comm = (w.dist == Local ? MPI_COMM_SELF : MPI_COMM_WORLD);
	superbblas::load<N, Nw, T, Tw>(detail::safe_div<T>(scalar, w.scalar), ctx.get(),
				       order.c_str(), from, size, w.p->p.data(), 1, w.order.c_str(),
				       w.from, w.dim, &w_ptr, &*w.ctx, comm, superbblas::FastToSlow,
				       superbblas::Copy);
      }
    };

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
    void urand(Tensor<N, T>& t, typename T::value_type a = 0, typename T::value_type b = 1)
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
    void urand(Tensor<N, T>& t, T a = 0, T b = 1)
    {
      std::uniform_real_distribution<T> d(a, b);
      fillWithCPUFuncNoArgs(
	t, [&]() { return d(detail::getSeed()); }, false);
    }

    /// Return an identity matrix
    /// \param dim: length for each of the row dimensions
    /// \param m: labels map from the row to the column dimensions
    /// \param order: order of the rows for the returned tensor
    /// \return: tensor with ordering `order`+`m[order]`

    template <std::size_t N, typename T>
    Tensor<N * 2, T> identity(const std::map<char, int>& dim, const remap& m,
			      Maybe<std::string> order)
    {
      // Get the order for the rows
      std::string orows;
      if (order)
      {
	orows = order.getSome();
      }
      else
      {
	for (const auto& it : dim)
	  orows.push_back(it.first);
      }

      // Get the order for the columns
      std::string ocols = detail::update_order(orows, m);

      // Get the dimensions of the returned tensor
      std::map<char, int> tdim = dim;
      for (const auto& it : dim)
	tdim[m.at(it.first)] = it.second;

      // Create the identity tensor
      Tensor<N * 2, T> t{orows + ocols, kvcoors<N * 2>(orows + ocols, tdim, 0, ThrowOnMissing),
			 OnHost, OnMaster};
      t.set_zero();
      if (t.getLocal())
      {
	T* p = t.data.get();
	for (unsigned int i = 0, vol = detail::volume(dim, orows); i < vol; ++i)
	  p[vol * i + i] = T{1};
      }

      return t;
    }

    /// Return a tensor filled with the value of the function applied to each element
    /// \param order: dimension labels, they should start with "xyztX"
    /// \param size: length of each dimension
    /// \param dev: either OnHost or OnDefaultDevice
    /// \param func: function (Coor<N-1>) -> COMPLEX

    template <std::size_t N, typename COMPLEX, typename Func>
    Tensor<N, COMPLEX> fillLatticeField(const std::string& order, const std::map<char, int>& size,
					DeviceHost dev, Func func)
    {
      static_assert(N >= 5, "The minimum number of dimensions should be 5");
      if (order.size() < 5 || order.compare(0, 5, "xyztX") != 0)
	throw std::runtime_error("Wrong `order`, it should start with xyztX");

      // Get final object dimension
      Coor<N> dim = latticeSize<N>(order, size);

      // Populate the tensor on CPU
      Tensor<N, COMPLEX> r(order, dim, OnHost);
      Coor<N> local_latt_size = r.p->localSize(); // local dimensions for xyztX
      Coor<N> stride = superbblas::detail::get_strides(local_latt_size, superbblas::FastToSlow);
      Coor<N> local_latt_from =
	r.p->localFrom(); // coordinates of first elements stored locally for xyztX
      std::size_t vol = superbblas::detail::volume(local_latt_size);
      Index nX = r.kvdim()['X'];
      COMPLEX* ptr = r.data.get();

#  ifdef _OPENMP
#    pragma omp parallel for schedule(static)
#  endif
      for (std::size_t i = 0; i < vol; ++i)
      {
	// Get the global coordinates
	using superbblas::detail::operator+;
	Coor<N> c = normalize_coor(
	  superbblas::detail::index2coor(i, local_latt_size, stride) + local_latt_from, dim);

	// Translate even-odd coordinates to natural coordinates
	Coor<N - 1> coor;
	coor[0] = c[0] * nX + (c[1] + c[2] + c[3] + c[4]) % nX; // x
	coor[1] = c[1];						// y
	coor[2] = c[2];						// z
	coor[3] = c[3];						// t
	std::copy_n(c.begin() + 5, N - 5, coor.begin() + 4);

	// Call the function
	ptr[i] = func(coor);
      }

      return r.make_sure(none, dev);
    }

    /// Compute a shift of v onto the direction dir
    /// \param v: tensor to apply the displacement
    /// \param first_tslice: global index in the t direction of the first element
    /// \param len: step of the displacement
    /// \param dir: 0 is x; 1 is y...

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> shift(const Tensor<N, COMPLEX> v, Index first_tslice, int len, int dir,
			     Maybe<Action> action = none, Maybe<Tensor<N, COMPLEX>> w = none)
    {
      if (dir < 0 || dir >= Nd - 1)
	throw std::runtime_error("Invalid direction");

      if (action.hasSome() != w.hasSome())
	throw std::runtime_error("Invalid default value");

      // Address zero length case
      if (len == 0)
      {
	if (!w.hasSome())
	  return v;
	v.doAction(action.getSome(), w.getSome());
	return w.getSome();
      }

      // NOTE: chroma uses the reverse convention for direction: shifting FORWARD moves the sites on the negative direction
      len = -len;

      const char dir_label[] = "xyz";
#  if QDP_USE_LEXICO_LAYOUT
      // If we are not using red-black ordering, return a view where the tensor is shifted on the given direction
      v = v.kvslice_from_size({{dir_label[dir], -len}});

      if (!w.hasSome())
	return v;

      v.doAction(action, w.getSome());
      return w.getSome();

#  elif QDP_USE_CB2_LAYOUT
      // Assuming that v has support on the origin and destination lattice elements
      int dimX = v.kvdim()['X'];
      if (dimX != 2 && len % 2 != 0)
	throw std::runtime_error("Unsupported shift");

      if (dir != 0)
      {
	if (!w.hasSome())
	  return v.kvslice_from_size({{'X', -len}, {dir_label[dir], -len}});
	v.doAction(action.getSome(),
		   w.getSome().kvslice_from_size({{'X', len}, {dir_label[dir], len}}));
	return w.getSome();
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
	Tensor<N, COMPLEX> r = w.hasSome() ? w.getSome() : v.like_this();
	auto r_eo = r.split_dimension('y', "Yy", maxY)
		      .split_dimension('z', "Zz", maxZ)
		      .split_dimension('t', "Tt", maxT);
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
				Maybe<Tensor<N, COMPLEX>> w = none)
    {
      if (std::abs(dir) > Nd)
	throw std::runtime_error("Invalid direction");

      if (action.hasSome() != w.hasSome())
	throw std::runtime_error("Invalid default value");

      // Address the zero direction case
      if (dir == 0)
      {
	if (!w.hasSome())
	  return v;
	v.doAction(action.getSome(), w.getSome());
	return w.getSome();
      }

      int d = std::abs(dir) - 1;    // space lattice direction, 0: x, 1: y, 2: z
      int len = (dir > 0 ? 1 : -1); // displacement unit direction
      assert(d < u.size());

      if (len > 0)
      {
	// Do u[d] * shift(x,d)
	Tensor<N, COMPLEX> r = w.hasSome() ? w.getSome() : v.like_this();
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
	const T* x = tnat.data.get();
	T* y = trb.data.get();

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
	T* x = tnat.data.get();
	const T* y = trb.data.get();

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
      Tensor<Nd + 1, T> getPhase(Coor<Nd - 1> phase, DeviceHost dev = OnDefaultDevice)
      {
	// Get spatial dimensions of the current lattice
	Coor<Nd> dim = latticeSize<Nd>("xyzX", {});
	dim[0] *= dim[3];
	return fillLatticeField<5, T>("xyztX", {}, dev, [=](Coor<Nd> c) {
	  typename T::value_type phase_dot_coor = 0;
	  for (int i = 0; i < Nd - 1; ++i)
	    phase_dot_coor += c[i] * 2 * M_PI * phase[i] / dim[i];

	  return T{cos(phase_dot_coor), sin(phase_dot_coor)};
	});
      }

      // NOTE: for now, the GPU version requires MAGMA
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

	// I have no idea how to do this....
	using MaybeTensor = Maybe<Tensor<Nd + 3, ComplexD>>;

	for (int mu = 0; mu < N; ++mu)
	{
	  displace(u, psi, first_tslice, mu + 1, Action::AddTo, MaybeTensor(chi));
	  displace(u, psi, first_tslice, -(mu + 1), Action::AddTo, MaybeTensor(chi));
	}
      }

      // Auxiliary structure passed to PRIMME's matvec

      struct OperatorAux {
	const std::vector<Tensor<Nd + 3, ComplexD>> u; // Gauge fields
	const Index first_tslice;		       // global t index
	const std::string order;		       // Laplacian input/output tensor's order
	const DeviceHost primme_dev;		       // where primme allocations are
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
	  Coor<Nd + 3> size = latticeSize<Nd + 3>(opaux.order, {{'n', *blockSize}, {'t', 1}});
	  Tensor<Nd + 3, ComplexD> tx(opaux.order, size, opaux.primme_dev, OnEveryone,
				      std::shared_ptr<ComplexD>((ComplexD*)x, [](ComplexD*) {}));
	  Tensor<Nd + 3, ComplexD> ty(opaux.order, size, opaux.primme_dev, OnEveryone,
				      std::shared_ptr<ComplexD>((ComplexD*)y, [](ComplexD*) {}));
	  LaplacianOperator(opaux.u, opaux.first_tslice, ty, tx);
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
		      .template make_sure<ComplexD>("ijxyztX");
	  }

#    if defined(SUPERBBLAS_USE_CUDA) && defined(BUILD_MAGMA)
	  DeviceHost primme_dev = OnDefaultDevice;
#    else
	  DeviceHost primme_dev = OnHost;
#    endif

	  // Create an auxiliary struct for the PRIMME's matvec
	  // NOTE: Please keep 'n' as the slowest index; the rows of vectors taken by PRIMME's matvec has dimensions 'cxyztX',
	  // and 'n' is the dimension for the columns.
	  OperatorAux opaux{ut, from_tslice + t, "cxyztXn", primme_dev};

	  // Make a bigger structure holding
	  primme_params primme;
	  primme_initialize(&primme);

	  // Get the global and local size of evec
	  std::size_t n, nLocal;
	  {
	    Tensor<Nd + 3, Complex> aux_tensor(
	      opaux.order, latticeSize<Nd + 3>(opaux.order, {{'n', 1}, {'t', 1}}), OnDefaultDevice,
	      OnEveryone);
	    n = aux_tensor.volume();
	    nLocal = aux_tensor.getLocal().volume();
	  }

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
	  Tensor<Nd + 3, ComplexD> evecs(
	    opaux.order, latticeSize<Nd + 3>(opaux.order, {{'n', primme.numEvals}, {'t', 1}}),
	    primme_dev, OnEveryone);
#    if defined(SUPERBBLAS_USE_CUDA) && defined(BUILD_MAGMA)
	  primme.queue = &*detail::getMagmaContext();
#    endif

	  // Call primme
#    if defined(SUPERBBLAS_USE_CUDA) && defined(BUILD_MAGMA)
	  int ret = magma_zprimme(evals.data(), evecs.data.get(), rnorms.data(), &primme);
#    else
	  int ret = zprimme(evals.data(), evecs.data.get(), rnorms.data(), &primme);
#    endif

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
	    LaplacianOperator(opaux.u, opaux.first_tslice, r, evecs);
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
				detail::remove_dimensions(colorvecs.order, "nt"))
		      .make_sure("nt", OnHost, OnEveryoneReplicated);

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
					   Coor<Nd - 1> phase = {})
    {
      // Phase colorvecs if phase != (0,0,0)
      if (phase == Coor<Nd - 1>{})
	return colorvecs;

      Tensor<Nd + 1, COMPLEX> tphase =
	ns_getColorvecs::getPhase<COMPLEX>(phase, colorvecs.getDev())
	  .kvslice_from_size({{'t', from_tslice}}, {{'t', colorvecs.kvdim()['t']}});
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
    Tensor<Nd + 3, COMPLEX> getColorvecs(const ColorvecsStorage& sto,
					 const multi1d<LatticeColorMatrix>& u, int decay_dir,
					 int from_tslice, int n_tslices, int n_colorvecs,
					 const Maybe<const std::string>& order = none,
					 Coor<Nd - 1> phase = {}, DeviceHost dev = OnDefaultDevice)
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
			  Coor<Nd - 1> phase = {},
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
      Coor<3> fingerprint_dim{};

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

    template <typename COMPLEX = Complex>
    using Moms = std::pair<Tensor<Nd + 2, COMPLEX>, std::vector<Coor<3>>>;

    /// Copy several momenta into a single tensor
    /// \param decay_dir: something that should be three
    /// \param moms: momenta to apply
    /// \param first_mom: first momentum to extract
    /// \param num_moms: number of momenta to extract
    /// \param first_tslice: first time-slice to extract
    /// \param num_tslice: number of time-slices to extract
    /// \param order_out: coordinate order of the output tensor, a permutation of mxyzXt
    /// \return: the tensor with the momenta

    template <typename COMPLEX = Complex>
    Moms<COMPLEX> getMoms(int decay_dir, const SftMom& moms, Maybe<int> first_mom = none,
			  Maybe<int> num_moms = none, Maybe<Index> first_tslice = none,
			  Maybe<int> num_tslices = none, const std::string& order_out = "mxyzXt")
    {
      // Copy moms into a single tensor
      const int Nt = Layout::lattSize()[decay_dir];
      int tfrom = first_tslice.getSome(0);	   // first tslice to extract
      int tsize = num_tslices.getSome(Nt);	   // number of tslices to extract
      int mfrom = first_mom.getSome(0);		   // first momentum to extract
      int msize = num_moms.getSome(moms.numMom()); // number of momenta to extract

      Tensor<Nd + 2, COMPLEX> momst(order_out,
				    latticeSize<Nd + 2>(order_out, {{'t', tsize}, {'m', msize}}));
      for (unsigned int mom = 0; mom < msize; ++mom)
      {
	asTensorView(moms[mfrom + mom])
	  .kvslice_from_size({{'t', tfrom}}, {{'t', tsize}})
	  .copyTo(momst.kvslice_from_size({{'m', mom}}, {{'m', 1}}));
      }

      // Create mom_list
      std::vector<Coor<Nd - 1>> mom_list(msize);
      for (unsigned int mom = 0; mom < msize; ++mom)
      {
	for (unsigned int i = 0; i < 3; ++i)
	  mom_list[mom][i] = moms.numToMom(mfrom + mom)[i];
      }

      return {momst, mom_list};
    }

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
      Tensor<Nright, COMPLEX> right, Index first_tslice, const SftMom& moms, int first_mom,
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
      std::array<bool, Nd> active_dirs{};
      unsigned int max_active_disps = 0;
      detail::get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

      // Number of moments to apply
      int numMom = num_moms.getSome(moms.numMom());
      if (first_mom + numMom > moms.numMom())
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
      std::vector<Coor<Nd - 1>> mom_list(numMom);
      for (unsigned int mom = 0; mom < numMom; ++mom)
      {
	for (unsigned int i = 0; i < Nd - 1; ++i)
	  mom_list[mom][i] = moms.numToMom(first_mom + mom)[i];
      }

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
	for (unsigned int mom = 0; mom < numMom; ++mom)
	{
	  asTensorView(moms[first_mom + mom])
	    .kvslice_from_size({{'t', tfrom + first_tslice}}, {{'t', tsize}})
	    .copyTo(momst.kvslice_from_size({{'m', mom}}, {{'m', 1}}));
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
      const multi1d<LatticeColorMatrix>& u, Tensor<Nin, COMPLEX> colorvec, Moms<COMPLEX> moms,
      Index first_tslice, const std::vector<std::array<std::vector<int>, 3>>& disps, bool deriv,
      const ColorContractionFn<COMPLEX>& call, Maybe<int> max_active_tslices = none,
      Maybe<int> max_active_momenta = none, Maybe<int> max_cols = none,
      const Maybe<std::string>& order_out = none, Maybe<DeviceHost> dev = none,
      Maybe<Distribution> dist = none)
    {
      const std::string order_out_str = order_out.getSome("ijkmt");
      detail::check_order_contains(order_out_str, "ijkmt");
      detail::check_order_contains(colorvec.order, "cxyzXtn");
      detail::check_order_contains(moms.first.order, "xyzXtm");

      // Form a tree with the displacement paths
      detail::PathNode tree_disps = ns_doMomDisp_colorContractions::get_tree(disps);

      // Get what directions are going to be used and the maximum number of displacements in memory
      std::array<bool, Nd> active_dirs{};
      unsigned int max_active_disps = 0;
      detail::get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

      // Check that all tensors have the same number of time
      int Nt = colorvec.kvdim()['t'];
      if (Nt != moms.first.kvdim()['t'])
	throw std::runtime_error("The t component of `colorvec' and `moms' does not match");

      int max_t = max_active_tslices.getSome(Nt);
      if (max_t <= 0)
	max_t = Nt;

      int Nmom = moms.first.kvdim()['m'];
      int max_active_moms = max_active_momenta.getSome(Nmom);
      if (max_active_moms <= 0)
	max_active_moms = Nmom;

      // Iterate over time-slices
      for (int tfrom = 0, tsize = std::min(max_t, Nt); tfrom < Nt;
	   tfrom += tsize, tsize = std::min(max_t, Nt - tfrom))
      {
	// Make tsize one or even
	if (tsize > 1 && tsize % 2 != 0)
	  --tsize;

	detail::log(1, "color contracting " + std::to_string(tsize) +
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

	// Loop over the momenta
	for (int mfrom = 0, msize = std::min(max_active_moms, Nmom); mfrom < Nmom;
	     mfrom += msize, msize = std::min(max_active_moms, Nmom - mfrom))
	{
	  auto this_moms = moms.first.kvslice_from_size({{'t', tfrom}, {'m', mfrom}},
							{{'t', tsize}, {'m', msize}});
	  if (tfrom + tsize >= Nt && mfrom + msize >= Nmom)
	  {
	    colorvec.release();
	    moms.first.release();
	  }
	  std::vector<Coor<3>> moms_list(moms.second.begin() + mfrom,
					 moms.second.begin() + mfrom + msize);
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
				Coor<3> left_phase, Coor<3> right_phase, Moms<COMPLEX> moms,
				Index first_tslice, const std::vector<std::vector<int>>& disps,
				bool deriv, const ContractionFn<COMPLEX>& call,
				const Maybe<std::string>& order_out = none,
				Maybe<DeviceHost> dev = none, Maybe<Distribution> dist = none)
    {
      const std::string order_out_str = order_out.getSome("ijmt");
      detail::check_order_contains(order_out_str, "ijmt");
      detail::check_order_contains(colorvec.order, "cxyzXtn");
      detail::check_order_contains(moms.first.order, "xyzXtm");

      // Form a tree with the displacement paths
      detail::PathNode tree_disps = detail::get_tree(disps);

      // Get what directions are going to be used and the maximum number of displacements in memory
      std::array<bool, Nd> active_dirs{};
      unsigned int max_active_disps = 0;
      detail::get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

      // Check that all tensors have the same number of time
      int Nt = colorvec.kvdim()['t'];
      if (Nt != moms.first.kvdim()['t'])
	throw std::runtime_error("The t component of `colorvec' and `moms' does not match");

      // Iterate over time-slices
      for (int tfrom = 0, tsize = Nt; tfrom < Nt; tfrom += tsize, tsize = Nt - tfrom)
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
	auto this_moms = moms.first.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});

	// Apply left phase and momenta conjugated to the left tensor
	// NOTE: look for the minus sign on left_phase in the doc of this function
	int Nmom = moms.first.kvdim()['m'];
	Tensor<Nin + 1, COMPLEX> moms_left =
	  colorvec.template like_this<Nin + 1>("mc%xyzXt", '%', "", {{'m', Nmom}});
	moms_left.contract(std::move(this_moms), {}, Conjugate,
			   phaseColorvecs(this_colorvec, first_tslice + tfrom, left_phase), {},
			   NotConjugate);

	// Apply the right phase
	this_colorvec = phaseColorvecs(this_colorvec, first_tslice + tfrom, right_phase);

	if (tfrom + tsize >= Nt)
	{
	  colorvec.release();
	  moms.first.release();
	}

	if (!deriv)
	{
	  ns_doMomDisp_contractions::doMomDisp_contractions<COMPLEX>(
	    ut, std::move(moms_left), this_colorvec, first_tslice + tfrom, tree_disps, deriv,
	    moms.second, 0, order_out_str, dev.getSome(OnDefaultDevice),
	    dist.getSome(OnEveryoneReplicated), call);
	}
	else
	{
	  // When using derivatives, each momenta has a different effect
	  std::vector<COMPLEX> ones(Nmom, COMPLEX(1));
	  Tensor<Nin + 1, COMPLEX> this_colorvec_m =
	    this_colorvec.template like_this<Nin + 1>("%m", '%', "", {{'m', Nmom}});
	  this_colorvec_m.contract(std::move(this_colorvec), {}, NotConjugate, asTensorView(ones),
				   {{'i', 'm'}}, NotConjugate);
	  ns_doMomDisp_contractions::doMomDisp_contractions<COMPLEX>(
	    ut, std::move(moms_left), this_colorvec_m, first_tslice + tfrom, tree_disps, deriv,
	    moms.second, 0, order_out_str, dev.getSome(OnDefaultDevice),
	    dist.getSome(OnEveryoneReplicated), call);
	}
      }
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
