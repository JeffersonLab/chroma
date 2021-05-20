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

#  include "actions/ferm/fermacts/fermact_factory_w.h"
#  include "actions/ferm/fermacts/fermacts_aggregate_w.h"
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
#  include <mutex>
#  include <thread>
#  include <set>
#  include <sstream>
#  include <stdexcept>
#  include <string>
#  include <type_traits>

#  ifndef M_PI
#    define M_PI                                                                                   \
      3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068L
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

    /// Where to store the tensor (see class Tensor)
    enum DeviceHost {
      OnHost,	      ///< on cpu memory
      OnDefaultDevice ///< on GPU memory if possible
    };

    /// How to distribute the tensor (see class Tensor)
    enum Distribution {
      OnMaster,		    ///< Fully supported on node with index zero
      OnEveryone,	    ///< Distributed the lattice dimensions (x, y, z, t) as chroma does
      OnEveryoneReplicated, ///< All nodes have a copy of the tensor
      Local		    ///< Non-collective
    };

    /// Whether complex conjugate the elements before contraction (see Tensor::contract)
    enum Conjugation { NotConjugate, Conjugate };

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
      template <typename Q>
      Maybe(Q t) : opt_val{true, T(t)}
      {
      }

      /// Return whether it has been initialized with a value
      bool hasSome() const
      {
	return opt_val.first;
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

	std::set<char> s;
	for (unsigned int i = 0; i < N; ++i)
	{
	  if (!s.insert(order[i]).second)
	  {
	    std::stringstream ss;
	    ss << "Invalid order: some label names are repeated `" << order << "`";
	    throw std::runtime_error(ss.str());
	  }
	}
      }
    }

    enum Throw_kvcoors { NoThrow, ThrowOnUnmatchLabel, ThrowOnMissing };

    template <std::size_t N>
    Coor<N> kvcoors(const std::string& order, std::map<char, int> m, Index missing = 0,
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
    Coor<N> latticeSize(const std::string& order, std::map<char, int> m = {})
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

      // Throw an error if `order` does not contain a label in `should_contain`
      inline std::string remove_dimensions(const std::string& order, const std::string& remove_dims)
      {
	std::string out;
	out.reserve(order.size());
	for (char c : order)
	  if (remove_dims.find(c) == std::string::npos)
	    out.push_back(c);
	return out;
      }

      template <std::size_t N>
      std::string update_order(std::string order, remap m)
      {
	for (std::size_t i = 0; i < N; ++i)
	{
	  auto it = m.find(order[i]);
	  if (it != m.end())
	    order[i] = it->second;
	}
	check_order<N>(order);
	return order;
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
#  ifdef QDP_IS_QDPJIT
	  if (!cudactx)
	  {
	    int dev = -1;
	    superbblas::detail::cudaCheck(cudaGetDevice(&dev));
	    cudactx = std::make_shared<superbblas::Context>(superbblas::createCudaContext(
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
	  return cudactx;
#  else
	  return cpuctx;
#  endif
	}
	throw std::runtime_error("Unsupported `DeviceHost`");
      }

      /// Return if two devices are the same

      inline bool is_same(DeviceHost a, DeviceHost b)
      {
#  ifdef QDP_IS_QDPJIT
	return a == b;
#  else
	// Without gpus, OnHost and OnDefaultDevice means on cpu.
	return true;
#  endif
      }

      /// Return an ordering with labels 0, 1, ...
      inline std::string getTrivialOrder(std::size_t N)
      {
	std::string r(N, 0);
	for (std::size_t i = 0; i < N; ++i)
	  r[i] = i % 128;
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

	TensorPartition(Coor<N> dim, const PartitionStored& p, bool isLocal = false)
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
	  return TensorPartition<N + 1>{insert_coor(dim, pos, dim_size), r};
	}

	/// Remove a non-distributed dimension

	TensorPartition<N - 1> remove_dimension(std::size_t pos) const
	{
	  typename TensorPartition<N - 1>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	    r.push_back({remove_coor(i[0], pos), remove_coor(i[1], pos)});
	  return TensorPartition<N - 1>{remove_coor(dim, pos), r};
	}

	/// Split a dimension into a non-distributed dimension and another dimension

	TensorPartition<N + 1> split_dimension(std::size_t pos, Index step) const
	{
	  typename TensorPartition<N + 1>::PartitionStored r;
	  r.reserve(p.size());
	  for (const auto& i : p)
	  {
	    if (i[1][pos] % step != 0 && i[1][pos] > step)
	      throw std::runtime_error("Unsupported splitting a dimension with an uneven lattice "
				       "portions in all processes");
	    r.push_back(
	      {insert_coor(replace_coor(i[0], pos, i[0][pos] % step), pos + 1, i[0][pos] / step),
	       insert_coor(replace_coor(i[1], pos, std::min(i[1][pos], step)), pos + 1,
			   (i[1][pos] + step - 1) / step)});
	  }
	  return TensorPartition<N + 1>{
	    insert_coor(replace_coor(dim, pos, std::min(dim[pos], step)), pos + 1, dim[pos] / step),
	    r};
	}

	/// Return a partition with the local portion of the tensor

	TensorPartition<N> get_local_partition() const
	{
	  return TensorPartition<N>{
	    localSize(), PartitionStored(1, superbblas::PartitionItem<N>{{{}, localSize()}}), true};
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

      using Sessions = std::unordered_map<std::thread::id, superbblas::Session>;

      inline Sessions& getSessions()
      {
	static Sessions sessions(16);
	return sessions;
      }

      inline std::mutex& getSessionsMutex()
      {
	static std::mutex m;
	return m;
      }

      inline superbblas::Session getSession()
      {
	std::lock_guard<std::mutex> g(getSessionsMutex());
	auto it = getSessions().find(std::this_thread::get_id());
	return (it != getSessions().end() ? it->second : 0);
      }

      inline void setSession(superbblas::Session session)
      {
	std::lock_guard<std::mutex> g(getSessionsMutex());
	getSessions()[std::this_thread::get_id()] = session;
      }

      inline void finishSession()
      {
	std::lock_guard<std::mutex> g(getSessionsMutex());
	getSessions().erase(std::this_thread::get_id());
      }

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

      inline void log(int level, std::string s)
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
	   << superbblas::getCpuMemUsed() / 1024 / 1024
	   << " MiB   GPU: " << superbblas::getGpuMemUsed() / 1024 / 1024 << " MiB";
	log(1, ss.str());
      }

      /// is_complex<T>::value is true if `T` is complex

      template <typename T>
      struct is_complex : std::false_type {
      };

      template <typename T>
      struct is_complex<std::complex<T>> : std::true_type {
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
    }

    /// RAII class for setting a secondary session for doing works in the background.
    /// NOTE: only one thread should be using the secondary session at the same time

    class AuxiliarySession
    {
    public:
      AuxiliarySession()
      {
	detail::setSession(1);
      }

      ~AuxiliarySession()
      {
	detail::finishSession();
      }
    };

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

      // Return a string describing the tensor
      std::string repr(T* ptr = nullptr) const
      {
	using namespace detail::repr;
	std::stringstream ss;
	ss << "Tensor{";
	if (ptr)
	  ss << "data:" << ptr << ", ";
	std::size_t sizemb = p->localVolume() * sizeof(T) / 1024 / 1024;
	ss << "order:" << order << ", dim:" << dim << ", dist:" << dist
	   << ", local_storage:" << sizemb << " MiB}";
	return ss.str();
      }

      // Construct used by non-Chroma tensors
      Tensor(const std::string& order, Coor<N> dim, DeviceHost dev = OnDefaultDevice,
	     Distribution dist = OnEveryone)
	: order(order),
	  dim(dim),
	  ctx(detail::getContext(dev)),
	  dist(dist),
	  from{},
	  size(dim),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow)),
	  scalar{1}
      {
	checkOrder();
	superbblas::Context ctx0 = *ctx;
	p = std::make_shared<detail::TensorPartition<N>>(
	  detail::TensorPartition<N>(order, dim, dist));
	std::string s = repr();
	detail::log(1, std::string("allocating ") + s);
	T* ptr = superbblas::allocate<T>(p->localVolume(), *ctx);
	detail::log_mem();
	data = std::shared_ptr<T>(ptr, [=](const T* ptr) {
	  superbblas::deallocate(ptr, ctx0);
	  detail::log(1, std::string("deallocated ") + s);
	  detail::log_mem();
	});
      }

      // Empty constructor
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
	  scalar{0}
      {
      }

      // Construct used by Chroma tensors (see `asTensorView`)
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
	  scalar{1}
      {
	checkOrder();

	// For now, TensorPartition creates the same distribution as chroma for tensor with
	// dimensions divisible by chroma logical dimensions
	p = std::make_shared<detail::TensorPartition<N>>(
	  detail::TensorPartition<N>(order, dim, dist));
      }

      // Construct for toFakeReal
      Tensor(const std::string& order, Coor<N> dim, std::shared_ptr<superbblas::Context> ctx,
	     std::shared_ptr<T> data, std::shared_ptr<detail::TensorPartition<N>> p,
	     Distribution dist, Coor<N> from, Coor<N> size, T scalar)
	: order(order),
	  dim(dim),
	  ctx(ctx),
	  data(data),
	  p(p),
	  dist(dist),
	  from(normalize_coor(from, dim)),
	  size(size),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow)),
	  scalar(scalar)
      {
	checkOrder();
      }

    protected:
      // Construct a slice of a tensor
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
	  scalar{t.scalar}
      {
	checkOrder();
      }

      // Construct a scaled tensor
      Tensor(const Tensor& t, T scalar)
	: order(t.order),
	  dim(t.dim),
	  ctx(t.ctx),
	  data(t.data),
	  p(t.p),
	  dist(t.dist),
	  from(t.from),
	  size(t.size),
	  strides(t.strides),
	  scalar{scalar}
      {
	checkOrder();
      }

    public:
      /// Return whether the tensor is not empty
      explicit operator bool() const noexcept
      {
	return superbblas::detail::volume(size) > 0;
      }

      // Return whether from != 0 or size != dim
      bool isSubtensor() const
      {
	return (from != Coor<N>{} || size != dim);
      }

      // Return the dimensions of the tensor
      std::map<char, int> kvdim() const
      {
	std::map<char, int> d;
	for (unsigned int i = 0; i < N; ++i)
	  d[order[i]] = size[i];
	return d;
      }

      // Get an element of the tensor
      T get(Coor<N> coor) const
      {
	if (ctx->plat != superbblas::CPU)
	  throw std::runtime_error(
	    "Unsupported to `get` elements from tensors not stored on the host");
	if (dist != OnMaster && dist != Local)
	  throw std::runtime_error("Unsupported to `get` elements on tensor that are not local or "
				   "not being fully supported on the master node");

	// coor[i] = coor[i] + from[i]
	for (unsigned int i = 0; i < N; ++i)
	  coor[i] = normalize_coor(normalize_coor(coor[i], size[i]) + from[i], dim[i]);

	return data.get()[detail::coor2index<N>(coor, dim, strides)] * scalar;
      }

      /// Rename dimensions
      Tensor<N, T> rename_dims(SB::remap m) const
      {
	return Tensor<N, T>(*this, detail::update_order<N>(order, m), this->from, this->size);
      }

      // Return a slice of the tensor starting at coordinate `kvfrom` and taking `kvsize` elements in each direction.
      // The missing dimension in `kvfrom` are set to zero and the missing direction in `kvsize` are set to the active size of the tensor.
      Tensor<N, T> kvslice_from_size(std::map<char, int> kvfrom = {},
				     std::map<char, int> kvsize = {}) const
      {
	std::map<char, int> updated_kvsize = this->kvdim();
	for (const auto& it : kvsize)
	  updated_kvsize[it.first] = it.second;
	return slice_from_size(kvcoors<N>(order, kvfrom), kvcoors<N>(order, updated_kvsize));
      }

      // Return a slice of the tensor starting at coordinate `from` and taking `size` elements in each direction.
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

      template <std::size_t Nn = N, typename Tn = T>
      Tensor<Nn, Tn> like_this(Maybe<std::string> new_order = none, std::map<char, int> kvsize = {},
			       Maybe<DeviceHost> new_dev = none,
			       Maybe<Distribution> new_dist = none) const
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

      template <std::size_t Nn = N, typename Tn = T>
      Tensor<Nn, Tn> like_this(std::string new_order, char remaining_char,
			       std::string remove_dims = "", std::map<char, int> kvsize = {},
			       Maybe<DeviceHost> new_dev = none,
			       Maybe<Distribution> new_dist = none) const
      {
	std::map<char, int> new_kvdim = kvdim();
	for (const auto& it : kvsize)
	  new_kvdim[it.first] = it.second;
	std::string new_order_ =
	  detail::remove_dimensions(get_order_for_reorder(new_order, remaining_char), remove_dims);
	return Tensor<Nn, Tn>(new_order_, kvcoors<Nn>(new_order_, new_kvdim, 0, ThrowOnMissing),
			      new_dev.getSome(getDev()), new_dist.getSome(dist));
      }

      /// Return a copy of this tensor, possible with a new precision `nT`

      template <typename Tn = T>
      Tensor<N, Tn> clone() const
      {
	return cloneOn<Tn>(getDev());
      }

      /// Return a copy of this tensor on device `new_dev`, possible with a new precision `nT`
      /// \param new_dev: device that will hold the new tensor

      template <typename Tn = T>
      Tensor<N, Tn> cloneOn(DeviceHost new_dev) const
      {
	Tensor<N, Tn> r = like_this<N, Tn>(none, {}, new_dev);
	copyTo(r);
	return r;
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
				    new_size, new_scalar);
      }

      template <typename U = T,
		typename std::enable_if<!detail::is_complex<U>::value, bool>::type = true>
      Tensor<N - 1, std::complex<U>> toComplex(bool allow_cloning = true) const
      {
	assert(isFakeReal() && kvdim()['.'] == 2);

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
				    new_size, new_scalar);
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

      /// Return a fake real view of this tensor

      Tensor<N + 1, T> split_dimension(char dim_label, std::string new_labels, Index step) const
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

	// Check the other arguments
	if (new_labels.size() != 2)
	  throw std::runtime_error("`new_labels` should have two labels!");

	if (step < 1)
	  throw std::runtime_error("`step` cannot be zero or negative");

	if (size[pos] % step != 0 && size[pos] > step)
	  throw std::runtime_error("Not supporting `split_dimension` for this lattice dimensions");

	// Set the new characteristics of the tensor
	std::string new_order =
	  insert_coor(replace_coor(order, pos, new_labels[0]), pos + 1, new_labels[1]);
	Coor<N + 1> new_from =
	  insert_coor(replace_coor(from, pos, from[pos] % step), pos + 1, from[pos] / step);
	Coor<N + 1> new_size = insert_coor(replace_coor(size, pos, std::min(size[pos], step)),
					   pos + 1, (size[pos] + step - 1) / step);
	Coor<N + 1> new_dim = insert_coor(replace_coor(dim, pos, std::min(dim[pos], step)), pos + 1,
					  (dim[pos] + step - 1) / step);

	auto new_p =
	  std::make_shared<detail::TensorPartition<N + 1>>(p->split_dimension(pos, step));

	return Tensor<N + 1, T>(new_order, new_dim, ctx, data, new_p, dist, new_from, new_size,
				scalar);
      }

      /// Copy this tensor into the given one
      template <std::size_t Nw, typename Tw,
		typename std::enable_if<
		  detail::is_complex<T>::value != detail::is_complex<Tw>::value, bool>::type = true>
      void copyTo(Tensor<Nw, Tw> w) const
      {
	toFakeReal().copyTo(w.toFakeReal());
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
			    Local, normalize_coor(from - p->localFrom(), dim), lsize, scalar);
      }

      /// Copy this tensor into the given one
      template <std::size_t Nw, typename Tw,
		typename std::enable_if<
		  detail::is_complex<T>::value == detail::is_complex<Tw>::value, bool>::type = true>
      void copyTo(Tensor<Nw, Tw> w) const
      {
	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > wsize[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	if ((dist == Local && w.dist != Local) || (dist != Local && w.dist == Local))
	{
	  getLocal().copyTo(w.getLocal());
	  return;
	}

	T* ptr = this->data.get();
	Tw* w_ptr = w.data.get();
	MPI_Comm comm =
	  ((dist == OnMaster && w.dist == OnMaster) || dist == Local ? MPI_COMM_SELF
								     : MPI_COMM_WORLD);
	if (dist != OnMaster || w.dist != OnMaster || Layout::nodeNumber() == 0)
	{
	  superbblas::copy<N, Nw>(detail::safe_div<T>(scalar, w.scalar), p->p.data(), 1,
				  order.c_str(), from, size, (const T**)&ptr, &*ctx, w.p->p.data(),
				  1, w.order.c_str(), w.from, &w_ptr, &*w.ctx, comm,
				  superbblas::FastToSlow, superbblas::Copy, detail::getSession());
	}
      }

      // Add `this` tensor into the given one
      template <std::size_t Nw, typename Tw>
      void addTo(Tensor<Nw, Tw> w) const
      {
	Coor<N> wsize = kvcoors<N>(order, w.kvdim(), 1, NoThrow);
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > wsize[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	if (w.scalar != T{1})
	  throw std::runtime_error("Not allowed to addTo to tensor with a scalar not being one");

	T* ptr = this->data.get();
	Tw* w_ptr = w.data.get();
	superbblas::copy<N, Nw>(scalar, p->p.data(), 1, order.c_str(), from, size, (const T**)&ptr,
				&*ctx, w.p->p.data(), 1, w.order.c_str(), w.from, &w_ptr, &*w.ctx,
				MPI_COMM_WORLD, superbblas::FastToSlow, superbblas::Add,
				detail::getSession());
      }

      // Contract the dimensions with the same label in `v` and `w` than do not appear on `this` tensor.
      template <std::size_t Nv, std::size_t Nw>
      void contract(Tensor<Nv, T> v, remap mv, Conjugation conjv, Tensor<Nw, T> w, remap mw,
		    Conjugation conjw, remap mr = {}, T beta = T{0})
      {
	// If either v or w is on OnDevice, force both to be on device
	if (v.ctx->plat != w.ctx->plat)
	{
	  if (v.getDev() != OnDefaultDevice)
	    v = v.cloneOn(OnDefaultDevice);
	  if (w.getDev() != OnDefaultDevice)
	    w = w.cloneOn(OnDefaultDevice);
	}

	// Superbblas tensor contraction is shit and those not deal with subtensors (for now)
	if (v.isSubtensor())
	  v = v.clone();
	if (w.isSubtensor())
	  w = w.clone();
	if (isSubtensor())
	{
	  Tensor<N, T> aux = like_this();
	  aux.contract(v, mv, conjv, w, mw, conjw, mr);
	  aux.copyTo(*this);
	  return;
	}

	T* v_ptr = v.data.get();
	T* w_ptr = w.data.get();
	T* ptr = this->data.get();
	std::string orderv_ = detail::update_order<Nv>(v.order, mv);
	std::string orderw_ = detail::update_order<Nw>(w.order, mw);
	std::string order_ = detail::update_order<N>(order, mr);
	superbblas::contraction<Nv, Nw, N>(
	  v.scalar * w.scalar / scalar, v.p->p.data(), 1, orderv_.c_str(), conjv == Conjugate,
	  (const T**)&v_ptr, &*v.ctx, w.p->p.data(), 1, orderw_.c_str(), conjw == Conjugate,
	  (const T**)&w_ptr, &*w.ctx, beta, p->p.data(), 1, order_.c_str(), &ptr, &*ctx,
	  MPI_COMM_WORLD, superbblas::FastToSlow, detail::getSession());
      }

      Tensor<N, T> scale(T s) const
      {
	return Tensor<N, T>(*this, scalar * s);
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
      }

      // Return whether the current view is contiguous in memory
      bool isContiguous() const
      {
	// Meaningless for tensors not been fully supported on a single node
	if (dist != OnMaster && dist != Local)
	  return false;

	if (superbblas::detail::volume(size) > 0 && N > 1)
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
      Tensor<N, Tn> make_sure(Maybe<std::string> new_order = none, Maybe<DeviceHost> new_dev = none,
			      Maybe<Distribution> new_dist = none) const
      {
	if (new_order.getSome(order) != order ||
	    !detail::is_same(new_dev.getSome(getDev()), getDev()) || new_dist.getSome(dist) != dist)
	{
	  Tensor<N, Tn> r = like_this(new_order, {}, new_dev, new_dist);
	  copyTo(r);
	  return r;
	}
	else
	{
	  return *this;
	}
      }

      template <typename Tn = T,
		typename std::enable_if<!std::is_same<T, Tn>::value, bool>::type = true>
      Tensor<N, Tn> make_sure(Maybe<std::string> new_order = none, Maybe<DeviceHost> new_dev = none,
			      Maybe<Distribution> new_dist = none) const
      {
	Tensor<N, Tn> r = like_this<Tn>(new_order, {}, new_dev, new_dist);
	copyTo(r);
	return r;
      }

      /// Get where the tensor is stored

      DeviceHost getDev() const
      {
#  ifdef QDP_IS_QDPJIT
	return (ctx->plat != superbblas::CPU ? OnDefaultDevice : OnHost);
#  else
	return OnDefaultDevice;
#  endif
      }

      void binaryRead(BinaryReader& bin)
      {
	if (ctx->plat != superbblas::CPU)
	  throw std::runtime_error("Only supported to read on `OnHost` tensors");
	if (dist != OnMaster)
	  throw std::runtime_error("Only supported to read on `OnMaster` tensors");
	if (!isContiguous())
	  throw std::runtime_error("Only supported contiguous views in memory");
	if (scalar != T{1})
	  throw std::runtime_error("Not allowed for tensor with a scale not being one");

	// Only on primary node read the data
	std::size_t vol = superbblas::detail::volume(size);
	std::size_t disp = detail::coor2index<N>(from, dim, strides);
	std::size_t word_size = sizeof(typename detail::WordType<T>::type);
	bin.readArrayPrimaryNode((char*)&data.get()[disp], word_size, sizeof(T) / word_size * vol);
      }

      void binaryWrite(BinaryWriter& bin) const
      {
	// If the writing is collective, the root process needs to hold the whole tensor
	if (!bin.isLocal() && dist != OnMaster)
	  throw std::runtime_error("For collective writing, the tensor should be `OnMaster`");

	// If the writing is non-collective, the tensor should be local
	if (bin.isLocal() && dist != Local)
	  throw std::runtime_error("For non-collective writing, the tensor should be `Local`");

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

      void print(std::string name) const
      {
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
	  std::size_t vol = superbblas::detail::volume(size);
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

    inline void* getQDPPtrFromId(int id)
    {
#  ifdef QDP_IS_QDPJIT
      std::vector<QDPCache::ArgKey> v(id, 1);
      return QDP_get_global_cache().get_kernel_args(v, false)[0];
#  else
      return nullptr;
#  endif
    }

    template <typename T>
    void* getQDPPtr(const T& t)
    {
#  ifdef QDP_IS_QDPJIT
      std::vector<QDPCache::ArgKey> v(1, t.getId());
      void* r = QDP_get_global_cache().get_kernel_args(v, false)[0];
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

#  ifndef QDP_IS_QDPJIT
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

#  ifndef QDP_IS_QDPJIT
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

#  ifndef QDP_IS_QDPJIT
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
    Tensor<1, COMPLEX> asTensorView(std::vector<COMPLEX>& v)
    {
      return Tensor<1, COMPLEX>("i", Coor<1>{Index(v.size())}, OnHost, OnEveryoneReplicated,
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

    /// Return a tensor filled with the value of the function applied to each element
    /// \param order: dimension labels, they should start with "xyztX"
    /// \param size: length of each dimension
    /// \param dev: either OnHost or OnDefaultDevice
    /// \param func: function (Coor<N-1>) -> COMPLEX

    template <std::size_t N, typename COMPLEX, typename Func>
    Tensor<N, COMPLEX> fillLatticeField(std::string order, std::map<char, int> size, DeviceHost dev,
					Func func)
    {
      static_assert(N >= Nd + 1, "The minimum number of dimensions should be Nd+1");
      if (order.size() < Nd + 1 || order.compare(0, 5, "xyztX") != 0)
	throw std::runtime_error("Wrong `order`, it should start with xyztX");

      // Get final object dimension
      Coor<N> dim = latticeSize<N>(order, size);

      // Populate the tensor on CPU
      Tensor<N, COMPLEX> r(order, dim, OnHost);
      Coor<N> local_latt_size = r.p->localSize(); // local dimensions for xyztX
      Coor<N> stride =
	superbblas::detail::get_strides<Nd + 1>(local_latt_size, superbblas::FastToSlow);
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
	coor[0] = c[0] * 2 + (c[1] + c[2] + c[3] + c[4]) % nX; // x
	coor[1] = c[1];					       // y
	coor[2] = c[2];					       // z
	coor[3] = c[3];					       // t
	std::copy_n(c.begin() + 5, N - (Nd - 1), coor.begin() + 4);

	// Call the function
	ptr[i] = func(coor);
      }

      return r.make_sure(none, dev);
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
      void applyPerm(const std::vector<Index>& perm, Tensor<Nd, T> tnat, Tensor<Nd + 1, T> trb)
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
    }

    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, Tensor<Nd, ComplexF>> MODS_t;

    /// Get colorvecs(t,n) for t=from_slice..(from_slice+n_tslices-1) and n=0..(n_colorvecs-1)

    template <typename COMPLEX = ComplexF>
    Tensor<Nd + 3, COMPLEX>
    getColorvecs(MODS_t& eigen_source, int decay_dir, int from_tslice, int n_tslices,
		 int n_colorvecs, Maybe<const std::string> order_ = none, Coor<Nd - 1> phase = {})
    {
      using namespace ns_getColorvecs;

      StopWatch sw;
      sw.reset();
      sw.start();

      if (decay_dir != 3)
	throw std::runtime_error("Only support for decay_dir being the temporal dimension");
      const std::string order = order_.getSome("cxyztXn");
      detail::check_order_contains(order, "cxyztXn");

      // Phase colorvecs if phase != (0,0,0)
      bool phasing = (phase != Coor<Nd - 1>{});

      from_tslice = normalize_coor(from_tslice, Layout::lattSize()[decay_dir]);

      // Allocate tensor to return
      std::string r_order = phasing ? "cnxyztX" : order;
      Tensor<Nd + 3, COMPLEX> r(
	r_order, latticeSize<Nd + 3>(r_order, {{'t', n_tslices}, {'n', n_colorvecs}}));

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
	std::vector<Index> perm = getPermFromNatToRB(t_slice);

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
	  applyPerm(perm, tnat, trb);

	  // t[n=colorvec] = trb
	  trb.copyTo(t.kvslice_from_size({{'n', colorvec}}, {{'n', 1}}));
	}

	// r[t=i_slice] = t, distribute the tensor from master to the rest of the nodes
	t.copyTo(r.kvslice_from_size({{'t', i_slice}}));
      }

      // Apply phase
      if (phasing)
      {
	Tensor<Nd + 1, COMPLEX> tphase =
	  getPhase<COMPLEX>(phase)
	    .kvslice_from_size({{'t', from_tslice}}, {{'t', n_tslices}})
	    .reorder("xyztX");
	Tensor<Nd + 3, COMPLEX> rp = r.like_this("cnxyztX");
	rp.contract(r, {}, NotConjugate, tphase, {}, NotConjugate);
	r = rp.reorder(order);
      }

      sw.stop();
      QDPIO::cout << "Time to read " << n_colorvecs << " colorvecs from " << n_tslices
		  << " time slices: " << sw.getTimeInSeconds() << " secs" << std::endl;

      return r;
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
    doInversion(const SystemSolver<LatticeFermion>& PP, const Tensor<Nd + 3, COMPLEX_CHI> chi,
		int t_source, int first_tslice_out, int n_tslice_out, std::vector<int> spin_sources,
		int max_rhs, const std::string& order_out = "cSxyztXns")
    {
      detail::check_order_contains(order_out, "cSxyztXns");
      if (chi.kvdim()['t'] != 1)
	throw std::runtime_error("Expected one time-slice");
      const int num_vecs = chi.kvdim()['n'];

      if (n_tslice_out > Layout::lattSize()[3])
	throw std::runtime_error("Too many tslices");

      Tensor<Nd + 5, COMPLEX_OUT> psi(
	order_out,
	latticeSize<Nd + 5>(
	  order_out,
	  {{'t', n_tslice_out}, {'S', Ns}, {'s', spin_sources.size()}, {'n', num_vecs}}));

      int max_step = std::max(num_vecs, max_rhs);
      std::vector<std::shared_ptr<LatticeFermion>> chis(max_step), quark_solns(max_step);
      for (int col = 0; col < max_step; col++)
	chis[col].reset(new LatticeFermion);
      for (int col = 0; col < max_step; col++)
	quark_solns[col].reset(new LatticeFermion);

      StopWatch snarss1;
      snarss1.reset();
      snarss1.start();

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

      snarss1.stop();
      QDPIO::cout << "Time to compute inversions for " << spin_sources.size()
		  << " spin sources and " << num_vecs
		  << " colorvecs : " << snarss1.getTimeInSeconds() << " secs" << std::endl;

      return psi;
    }

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> shift(const Tensor<N, COMPLEX> v, Index first_tslice, int len, int dir)
    {
      if (dir < 0 || dir >= Nd - 1)
	throw std::runtime_error("Invalid direction");

      if (len == 0)
	return v;

      // NOTE: chroma uses the reverse convention for direction: shifting FORWARD moves the sites on the negative direction
      len = -len;

      const char dir_label[] = "xyz";
#  if QDP_USE_LEXICO_LAYOUT
      // If we are not using red-black ordering, return a view where the tensor is shifted on the given direction
      return v.kvslice_from_size({{dir_label[dir], -len}});

#  elif QDP_USE_CB2_LAYOUT
      // Assuming that v has support on the origin and destination lattice elements
      if (v.kvdim()['X'] != 2 && len % 2 != 0)
	throw std::runtime_error("Unsupported shift");

      Tensor<N, COMPLEX> r = v.like_this();
      if (dir != 0)
      {
	v.copyTo(r.kvslice_from_size({{'X', len}, {dir_label[dir], len}}));
      }
      else
      {
	auto v_eo = v.split_dimension('y', "Yy", 2)
		      .split_dimension('z', "Zz", 2)
		      .split_dimension('t', "Tt", 2);
	auto r_eo = r.split_dimension('y', "Yy", 2)
		      .split_dimension('z', "Zz", 2)
		      .split_dimension('t', "Tt", 2);
	while (len < 0)
	  len += v.kvdim()['x'] * 2;
	for (int T = 0; T < 2; ++T)
	  for (int Z = 0; Z < 2; ++Z)
	    for (int Y = 0; Y < 2; ++Y)
	      for (int X = 0; X < 2; ++X)
		v_eo
		  .kvslice_from_size({{'X', X}, {'Y', Y}, {'Z', Z}, {'T', T}},
				     {{'X', 1}, {'Y', 1}, {'Z', 1}, {'T', 1}})
		  .copyTo(
		    r_eo.kvslice_from_size({{'X', X + len},
					    {'x', (len + ((X + Y + Z + T + first_tslice) % 2)) / 2},
					    {'Y', Y},
					    {'Z', Z},
					    {'T', T}},
					   {{'Y', 1}, {'Z', 1}, {'T', 1}}));
      }
      return r;
#  else
      throw std::runtime_error("Unsupported layout");
#  endif
    }

    /// Compute displace/derivate of v in the direction dir
    /// \param u: Gauge field
    /// \param v: tensor to apply the derivative
    /// \param dir: 0: nothing; 1: forward x; -1: backward x; 2: forward y...
    /// \patam deriv: whether use derivative
    /// \param moms: list of input momenta
    /// \param conjUnderAdd: if true, return a version, R(dir), so that
    ////       adj(R(dir)) * D(dir') == D(dir+dir'), where D(dir') is what this function returns
    ////       when conjUnderAdd is false.

    template <typename COMPLEX, std::size_t N>
    Tensor<N, COMPLEX> displace(const std::vector<Tensor<Nd + 3, Complex>>& u, Tensor<N, COMPLEX> v,
				Index first_tslice, int dir, bool deriv = false,
				std::vector<Coor<Nd - 1>> moms = {}, bool conjUnderAdd = false)
    {
      if (std::abs(dir) > Nd)
	throw std::runtime_error("Invalid direction");

      if (dir == 0)
	return v;

      int d = std::abs(dir) - 1;    // space lattice direction, 0: x, 1: y, 2: z
      int len = (dir > 0 ? 1 : -1); // displacement unit direction

      Tensor<N, COMPLEX> r = v.like_this("cnSs%xyzXt", '%');
      if (!deriv)
      {
	assert(d < u.size());

	if (conjUnderAdd)
	  len *= -1;

	v = v.reorder("cnSs%xyzXt", '%');
	if (len > 0)
	{
	  // Do u[d] * shift(x,d)
	  v = shift(std::move(v), first_tslice, len, d);
	  r.contract(u[d], {{'j', 'c'}}, NotConjugate, std::move(v), {}, NotConjugate,
		     {{'c', 'i'}});
	}
	else
	{
	  // Do shift(adj(u[d]) * x,d)
	  r.contract(u[d], {{'i', 'c'}}, Conjugate, std::move(v), {}, NotConjugate, {{'c', 'j'}});
	  r = shift(std::move(r), first_tslice, len, d);
	}
      }
      else
      {
	// conj(phase)*displace(u, v, -length, d) - phase*displace(u, v, length, d)
	std::vector<COMPLEX> phases(moms.size());
	for (unsigned int i = 0; i < moms.size(); ++i)
	{

	  double angle = 2 * M_PI * moms[i][d] / Layout::lattSize()[d];
	  phases[i] = COMPLEX{1} + COMPLEX{cos(angle), sin(angle)};
	  if (conjUnderAdd)
	    phases[i] = std::sqrt(phases[i]);
	}

	// r = conj(phases) * displace(u, v, dir)
	r.contract(displace(u, v, first_tslice, -dir), {}, NotConjugate, asTensorView(phases),
		   {{'i', 'm'}}, Conjugate);
	// r = r - phases * displace(u, v, dir) if !ConjUnderAdd else r + phases * displace(u, v, dir)
	r.contract(displace(u, v, first_tslice, dir, false).scale(conjUnderAdd ? 1 : -1), {},
		   NotConjugate, asTensorView(phases), {{'i', 'm'}}, NotConjugate, {}, 1.0);
      }
      return r;
    }

    namespace ns_doMomGammaDisp_contractions
    {
      /// Path Node
      struct PathNode {
	std::map<char, PathNode> p; ///< following nodes
	int disp_index;		    ///< if >= 0, the index in the displacement list
      };

      /// Return the tree representing all paths
      inline PathNode get_tree(const std::vector<std::vector<int>>& paths)
      {
	PathNode r{{}, -1};
	int path_index = 0;
	for (const std::vector<int>& path : paths)
	{
	  PathNode* n = &r;
	  for (char d : path)
	  {
	    if (d == 0 || std::abs(d) > Nd)
	      throw std::runtime_error("Invalid direction: " + std::to_string((int)d));

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

      /// Contract two LatticeFermion with different momenta, gammas, and displacements.
      /// \param leftconj: left lattice fermion tensor, cSxyzXN
      /// \param right: right lattice fermion tensor, csxyzXn
      /// \param disps: tree of displacements/derivatives
      /// \param deriv: whether use derivatives instead of displacements
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
	  detail::log(1, std::string("contracting for disp_index=") +
			   std::to_string(disps.disp_index));
	  // Contract the spatial components and the color of the leftconj and right tensors
	  Tensor<Nout, COMPLEX> aux =
	    r.template like_this<Nout, COMPLEX>("mNQqnSst%", '%', "gd", {{'S', Ns}, {'Q', Ns}});
	  aux.contract(leftconj, {}, Conjugate, right.reorder("cxyzXnSst%", '%'), {}, NotConjugate,
		       {});

	  // Contract the spin components S and Q with the gammas, and put the result on r[d=disp_indices.size()]
	  aux = aux.reorder("QSmNqnst");
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
	  detail::log(1, std::string("push on direction ") + std::to_string(it.first));
	  // Apply displacement on the right vectors
	  // NOTE: avoid that the memory requirements grow linearly with the number of displacements
	  //       by killing the reference to `right` as soon as possible
	  Tensor<Nright, COMPLEX> right_disp =
	    (node_disp < disps.p.size() - 1)
	      ? displace(u, right, first_tslice, it.first, deriv, moms)
	      : displace(u, std::move(right), first_tslice, it.first, deriv, moms);
	  doMomGammaDisp_contractions(u, leftconj, std::move(right_disp), first_tslice, it.second,
				      deriv, gammas, moms, max_rhs - num_vecs, r, disp_indices);
	  node_disp++;
	  detail::log(1, std::string("pop direction"));
	}
      }
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
    /// \param deriv: whether use derivatives instead of displacements
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
      std::vector<std::vector<int>> disps, bool deriv, const std::string& order_out = "gmNndsqt",
      Maybe<int> max_active_tslices = none)
    {
      using namespace ns_doMomGammaDisp_contractions;

      detail::check_order_contains(order_out, "gmNndsqt");
      detail::check_order_contains(leftconj.order, "cxyzXNQqt");
      detail::check_order_contains(right.order, "cxyzXnSst");

      assert(right.kvdim()['t'] == leftconj.kvdim()['t']);
      int Nt = right.kvdim()['t'];
      if (Nt % 2 != 0)
	throw std::runtime_error("The number of time-slices should be even");

      int max_t = max_active_tslices.getSome(Nt);
      if (max_t <= 0)
	max_t = Nt;
      max_t = max_t + (max_t % 2);
      max_t = std::min(Nt, max_t);

      // Form a tree with the displacement paths
      PathNode tree_disps = get_tree(disps);

      // Get what directions are going to be used and the maximum number of displacements in memory
      std::array<bool, Nd> active_dirs{};
      unsigned int max_active_disps = 0;
      get_tree_mem_stats(tree_disps, active_dirs, max_active_disps);

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
      Tensor<3, COMPLEX> gammast("gQS", {(Index)gammas.size(), Ns, Ns}, OnDefaultDevice,
				 OnEveryoneReplicated);
      for (unsigned int g = 0; g < gammas.size(); g++)
      {
	gammas[g]
	  .rename_dims({{'i', 'Q'}, {'j', 'S'}})
	  .copyTo(gammast.kvslice_from_size({{'g', g}}, {{'g', 1}}));
      }

      // Iterate over time-slices
      std::vector<int> disp_indices;
      leftconj = leftconj.reorder("QNqc%xyzXt", '%');

      for (int tfrom = 0, tsize = std::min(max_t, Nt); tfrom < Nt;
	   tfrom += tsize, tsize = std::min(max_t, Nt - tfrom))
      {
	detail::log(1,
		    std::string("contracting " + std::to_string(tsize) + " tslices from tslice= ") +
		      std::to_string(tfrom));

	disp_indices.resize(0);

	// Copy moms into a single tensor
	std::string momst_order = "mxyzXt";
	Tensor<Nd + 2, COMPLEX> momst(
	  momst_order, latticeSize<Nd + 2>(momst_order, {{'t', tsize}, {'m', numMom}}));
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
	moms_left = moms_left.reorder("cxyzXmNQqt%", '%');

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
		    .reorder("ijxyzXt");
	}

	// Do the thing
	auto this_right = right.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});
	if (tfrom + tsize >= Nt)
	  right.release();
	auto this_r = r.kvslice_from_size({{'t', tfrom}}, {{'t', tsize}});
	if (!deriv)
	{
	  doMomGammaDisp_contractions(ut, std::move(moms_left), std::move(this_right),
				      first_tslice + tfrom, tree_disps, deriv, gammast, mom_list, 0,
				      this_r, disp_indices);
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
