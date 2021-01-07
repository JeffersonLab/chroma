// -*- C++ -*-
/*! \file                                                                    
 * \brief Copy, permuting, contracting tensors with superbblas
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_SUPERB_CONTRACTIONS__
#define __INCLUDE_SUPERB_CONTRACTIONS__

// Activate the MPI support in Superbblas
#define SUPERBBLAS_USE_MPI

#include "chromabase.h"
#include "qdp.h"
#include "qdp_map_obj_disk_multiple.h"
#include "superbblas.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <map>
#include <memory>
#include <stdexcept>

namespace Chroma
{

  namespace SB
  {

    using Index = superbblas::IndexType;
    using Complex = std::complex<REAL>;
    using ComplexF = std::complex<REAL32>;
    template <std::size_t N>
    using Coor = superbblas::Coor<N>;
    template <std::size_t N>
    using Order = superbblas::Order<N>;

    enum DeviceHost { OnHost, OnDefaultDevice };
    enum Distribution { OnMaster, OnEveryone };
    enum Conjugation { NotConjugate, Conjugate };

    template <std::size_t N>
    Coor<N> kvcoors(const char* order, std::map<char, int> m, Index missing = 0)
    {
      if (std::strlen(order) != N)
	throw std::runtime_error(
	  "The length of the string `order` should match the template argument `N`");
      Coor<N> r;
      for (std::size_t i = 0; i < N; ++i)
      {
	auto it = m.find(order[i]);
	if (it != m.end())
	  r[i] = it->second;
	else
	  r[i] = missing;
      }
      return r;
    }

    template <std::size_t N>
    Coor<N> kvcoors(Order<N> order, std::map<char, int> m, Index missing = 0)
    {
      Coor<N> r;
      for (std::size_t i = 0; i < N; ++i)
      {
	auto it = m.find(order[i]);
	if (it != m.end())
	  r[i] = it->second;
	else
	  r[i] = missing;
      }
      return r;
    }

    template <std::size_t N>
    Coor<N> latticeSize(const char* order, std::map<char, int> m = {})
    {
      if (std::strlen(order) != N)
	throw std::runtime_error(
	  "The length of the string `order` should match the template argument `N`");

      std::map<char, int> m0 = {{'x', Layout::lattSize()[0]},
				{'y', Layout::lattSize()[1]},
				{'z', Layout::lattSize()[2]},
				{'t', Layout::lattSize()[3]},
				{'s', Ns},
				{'c', Nc}};
      m0.insert(m.begin(), m.end());
      return kvcoors<N>(order, m0);
    }

    // Replace a label by another label
    using remap = std::map<char, char>;

    namespace detail
    {
      using namespace superbblas::detail;

      template <std::size_t N>
      Order<N> toOrder(const char* order)
      {
	if (std::strlen(order) != N)
	  throw std::runtime_error(
	    "The length of the string `order` should match the template argument `N`");
	Order<N> o;
	std::copy_n(order, N, o.begin());
	return o;
      }

      template <std::size_t N>
      Order<N + 1> toOrderStr(const Order<N>& order)
      {
	Order<N + 1> o;
	std::copy_n(order.begin(), N, o.begin());
	o[N] = 0;
	return o;
      }

      // Throw an error if `order` does not contain a label in `should_contain`
      inline void check_order_contains(const char* order, const char* should_contain)
      {
	bool ok = false;
	if (std::strlen(order) == std::strlen(should_contain))
	{
	  int n = std::strlen(order);
	  unsigned int i;
	  for (i = 0; i < n; ++i)
	  {
	    unsigned int j;
	    for (j = 0; j < n && order[i] != should_contain[j]; ++j)
	      ;
	    if (j >= n)
	      break;
	  }
	  if (i >= n)
	    ok = true;
	}
	if (!ok)
	{
	  std::stringstream ss;
	  ss << "The input order `" << order
	     << "` is missing one of this labels: " << should_contain;
	  throw std::runtime_error(ss.str());
	}
      }

      // Return the equivalent value of the coordinate `v` in the interval [0, dim[ for a periodic
      // dimension with length `dim`.

      int normalize_coor(int v, int dim)
      {
	return (v + dim * (v < 0 ? -v / dim + 1 : 0)) % dim;
      }

      template <std::size_t N>
      Order<N> update_order(Order<N> order, remap m)
      {
	for (std::size_t i = 0; i < N; ++i)
	{
	  auto it = m.find(order[i]);
	  if (it != m.end())
	    order[i] = it->second;
	}
	return order;
      }

      template <std::size_t N>
      class TensorPartition
      {
      public:
	Coor<N> dim;
	superbblas::From_size<N> p;
	TensorPartition(Order<N> order, Coor<N> dim, Distribution dist) : dim(dim)
	{
	  switch (dist)
	  {
	  case OnMaster: p = all_tensor_on_master(dim); break;
	  case OnEveryone: p = partitioning_chroma_compatible(order, dim); break;
	  }
	}

	TensorPartition(Coor<N> dim, superbblas::From_size<N> p) : dim(dim), p(p)
	{
	}

      private:
	/// Return a partitioning where the root node has support for the whole tensor
	/// \param dim: dimension size for the tensor

	static superbblas::From_size<N> all_tensor_on_master(Coor<N> dim)
	{
	  int nprocs = Layout::numNodes();
	  superbblas::From_size<N> fs(nprocs);
	  if (1 <= nprocs)
	    fs[0][1] = dim;
	  return fs;
	}

	/// Return a partitioning for a tensor of `dim` dimension onto a grid of processes
	/// \param order: dimension labels (use x, y, z, t for lattice dimensions)
	/// \param dim: dimension size for the tensor

	static superbblas::From_size<N> partitioning_chroma_compatible(Order<N> order, Coor<N> dim)
	{
	  // Get the number of procs use in each dimension; for know we put as many as chroma
	  // put onto the lattice dimensions
	  multi1d<int> procs_ = Layout::logicalSize();
	  Coor<N> procs = kvcoors(
	    order, {{'x', procs_[0]}, {'y', procs_[1]}, {'z', procs_[2]}, {'t', procs_[3]}}, 1);

	  // For each proc, get its coordinate in procs (logical coordinate) and compute the
	  // fair range of the tensor supported on the proc
	  int num_procs = Layout::numNodes();
	  superbblas::From_size<N> fs(num_procs);
	  for (int rank = 0; rank < num_procs; ++rank)
	  {
	    multi1d<int> cproc_ = Layout::getLogicalCoordFrom(rank);
	    Coor<N> cproc = kvcoors(
	      order, {{'x', cproc_[0]}, {'y', cproc_[1]}, {'z', cproc_[2]}, {'t', cproc_[3]}});
	    for (unsigned int i = 0; i < N; ++i)
	    {
	      // First coordinate in process with rank 'rank' on dimension 'i'
	      fs[rank][0][i] = dim[i] / procs[i] * cproc[i] + std::min(cproc[i], dim[i] % procs[i]);
	      // Number of elements in process with rank 'cproc[i]' on dimension 'i'
	      fs[rank][1][i] = dim[i] / procs[i] + (dim[i] % procs[i] > cproc[i] ? 1 : 0);
	    }
	  }
	  return fs;
	}
      };

    }

    template <std::size_t N, typename T>
    class Tensor
    {
      static_assert(superbblas::supported_type<T>::value, "Not supported type");

    public:
      // Construct used by non-Chroma tensors
      Tensor(const char* order_, Coor<N> dim, DeviceHost dev = OnDefaultDevice,
	     Distribution dist = OnEveryone)
	: order(detail::toOrder<N>(order_)),
	  dim(dim),
	  ctx(getContext(dev)),
	  dist(dist),
	  from{},
	  size(dim),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow))
      {

	data = std::shared_ptr<T[]>(superbblas::allocate<T>(detail::volume<N>(dim), ctx),
				    [&](const T* ptr) { superbblas::deallocate(ptr, ctx); });
	p = std::make_shared<detail::TensorPartition<N>>(
	  detail::TensorPartition<N>(order, dim, dist));
      }

      // Construct used by Chroma tensors (see `asTensorView`)
      Tensor(const char* order_, Coor<N> dim, DeviceHost dev, Distribution dist,
	     std::shared_ptr<T[]> data)
	: order(detail::toOrder<N>(order_)),
	  dim(dim),
	  ctx(getContext(dev)),
	  data(data),
	  dist(dist),
	  from{},
	  size(dim),
	  strides(detail::get_strides<N>(dim, superbblas::FastToSlow))
      {

	// For now, TensorPartition creates the same distribution as chroma for tensor with
	// dimensions divisible by chroma logical dimensions
	p = std::make_shared<detail::TensorPartition<N>>(
	  detail::TensorPartition<N>(order, dim, dist));
      }

    protected:
      // Construct a slice of a tensor
      Tensor(const Tensor& t, Order<N> order, Coor<N> from, Coor<N> size)
	: order(order),
	  dim(t.dim),
	  ctx(t.ctx),
	  data(t.data),
	  p(t.p),
	  dist(t.dist),
	  from(from),
	  size(size),
	  strides(t.strides)
      {
      }

    public:
      // Get an element of the tensor
      T& get(Coor<N> coor) const
      {
	if (ctx.plat != superbblas::CPU)
	  throw std::runtime_error(
	    "Unsupported to `get` elements from tensors not stored on the host");
	if (dist != OnMaster)
	  throw std::runtime_error(
	    "Unsupported to `get` elements on tensor not being fully on the master node");

	using superbblas::detail::operator+;
	return data.get()[detail::coor2index<N>(coor + from, dim, strides)];
      }

      // Return a slice of the tensor starting at coordinate `kvfrom` and taking `kvsize` elements in each direction.
      // The missing dimension in `kvfrom` are set to zero and the missing direction in `kvsize` are set to the active size of the tensor.
      Tensor<N, T> kvslice_from_size(std::map<char, int> kvfrom = {},
				     std::map<char, int> kvsize = {}) const
      {
	std::map<char, int> updated_kvsize;
	for (unsigned int i = 0; i < N; ++i)
	  updated_kvsize[order[i]] = size[i];
	updated_kvsize.insert(kvsize.begin(), kvsize.end());
	return slice_from_size(kvcoors(order, kvfrom), kvcoors(order, updated_kvsize));
      }

      // Return a slice of the tensor starting at coordinate `from` and taking `size` elements in each direction.
      Tensor<N, T> slice_from_size(Coor<N> from, Coor<N> size) const
      {
	for (unsigned int i = 0; i < N; ++i)
	{
	  if (size[i] > this->size[i])
	    throw std::runtime_error(
	      "The size of the slice cannot be larger than the original tensor");
	  if (from[i] + size[i] >= this->size[i] && this->size[i] != this->dim[i])
	    throw std::runtime_error(
	      "Unsupported to make a view on a non-contiguous range on the tensor");
	}
	return Tensor<N, T>(*this, order, from, size);
      }

      // Copy `this` tensor into the given one
      template <std::size_t Nw, typename Tw>
      void copyTo(Tensor<Nw, Tw> w) const
      {
	for (unsigned int i = 0; i < N; ++i)
	  if (size[i] > w.size[i])
	    throw std::runtime_error("The destination tensor is smaller than the source tensor");

	T* ptr = this->data.get();
	Tw* w_ptr = w.data.get();
	Order<N + 1> order_ = detail::toOrderStr(order);
	Order<Nw + 1> orderw_ = detail::toOrderStr(w.order);
	superbblas::copy<N, Nw>(p->p, 1, &order_[0], from, size, (const T**)&ptr, &ctx, w.p->p, 1,
				&orderw_[0], w.from, &w_ptr, &w.ctx, MPI_COMM_WORLD,
				superbblas::FastToSlow);
      }

      // Contract the dimensions with the same label in `v` and `w` than do not appear on `this` tensor.
      template <std::size_t Nv, std::size_t Nw>
      void contract(Tensor<Nv, T> v, remap mv, Conjugation conjv, Tensor<Nw, T> w, remap mw,
		    Conjugation conjw, remap mr = {})
      {
	T* v_ptr = v.data.get();
	T* w_ptr = w.data.get();
	T* ptr = this->data.get();
	Order<Nv + 1> orderv_ = detail::toOrderStr(detail::update_order(v.order, mv));
	Order<Nw + 1> orderw_ = detail::toOrderStr(detail::update_order(w.order, mw));
	Order<N + 1> order_ = detail::toOrderStr(detail::update_order(order, mr));
	superbblas::contraction<Nv, Nw, N>(
	  v.p->p, 1, &orderv_[0], conjv == Conjugate, (const T**)&v_ptr, &v.ctx, w.p->p, 1,
	  &orderw_[0], conjw == Conjugate, (const T**)&w_ptr, &w.ctx, p->p, 1, &order_[0], &ptr,
	  &ctx, MPI_COMM_WORLD, superbblas::FastToSlow);
      }

      const Order<N> order;		///< Labels of the tensor dimensions
      const Coor<N> dim;		///< Length of the tensor dimensions
      const superbblas::Context ctx;	///< Tensor storage information (device/host)
      std::shared_ptr<T[]> data;	///< Pointer to the tensor storage
      std::shared_ptr<detail::TensorPartition<N>>
	p;				///< Distribution of the tensor among the processes
      const Distribution dist;		///< Whether the tensor is stored on the cpu or a device
      const Coor<N> from;		///< First active coordinate in the tensor
      const Coor<N> size;		///< Number of active coordinates on each dimension
      const Coor<N> strides;		///< Displacement for the next element along every direction

      // Return whether the current view is contiguous in memory
      bool isContiguous() const
      {
	if (superbblas::detail::volume<N>(size) > 0)
	  for (unsigned int i = 1; i < N; ++i)
	    if (from[i] != 0 || size[i] != dim[i])
	      return false;
	return true;
      }

      // Return a context on either the host or the device
      static superbblas::Context getContext(DeviceHost dev)
      {
	// Creating GPU context can be expensive; so do it once
	static std::shared_ptr<superbblas::Context> cudactx;

	switch (dev)
	{
	case OnHost: return superbblas::createCpuContext();
#ifdef QDP_IS_QDPJIT
	case OnDefaultDevice:
	  if (!cudactx)
	    cudactx = std::make_shared<superbblas::Context>(superbblas::createCudaContext());
	  return *cudactx;
#else
	case OnDefaultDevice: return superbblas::createCpuContext();
#endif
	}
	throw std::runtime_error("Unsupported `DeviceHost`");
      }

      void binaryRead(BinaryReader& bin)
      {
	if (dist != OnMaster)
	  throw std::runtime_error("Only supported to read on `OnMaster` tensors");
	if (!isContiguous())
	  throw std::runtime_error("Only supported contiguous views in memory");

	// Only on primary node read the data
	std::size_t vol = superbblas::detail::volume<N>(size);
	std::size_t disp = detail::coor2index<N>(from, dim, strides);
	bin.readArrayPrimaryNode((char*)&data.get()[disp], sizeof(T), vol);
      }

      void binaryWrite(BinaryWriter& bin) const
      {
	if (dist != OnMaster)
	  throw std::runtime_error("Only supported to write from `OnMaster` tensors");
	if (!isContiguous())
	  throw std::runtime_error("Only supported contiguous views in memory");

	// Only on primary node write the data
	std::size_t vol = superbblas::detail::volume<N>(size);
	std::size_t disp = detail::coor2index<N>(from, dim, strides);
	bin.writeArrayPrimaryNode((char*)&data.get()[disp], sizeof(T), vol);
      }
    };

    Tensor<Nd + 2, Complex> asTensorView(LatticeFermion& v)
    {
      Complex* v_ptr = reinterpret_cast<Complex*>(v.getF());
      return Tensor<Nd + 2, Complex>("csxyzt", latticeSize<Nd + 2>("csxyzt"), OnDefaultDevice,
				     OnEveryone,
				     std::shared_ptr<Complex[]>(v_ptr, [](Complex*) {}));
    }

    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, Tensor<Nd + 1, ComplexF>> MODS_t;

    Tensor<Nd + 2, ComplexF> getColorvecs(MODS_t& eigen_source, int decay_dir, int from_tslice,
					  int n_tslices, int n_colorvecs,
					  const char order[] = "cxyztn")
    {
      if (decay_dir != 3)
	throw std::runtime_error("Only support for decay_dir being the temporal dimension");
      detail::check_order_contains(order, "cxyztn");

      from_tslice = detail::normalize_coor(from_tslice, Layout::lattSize()[decay_dir]);

      Tensor<Nd + 2, ComplexF> r(order, latticeSize<Nd + 2>(order, {{'n', n_colorvecs}}));
      Tensor<Nd + 1, ComplexF> t("cxyzn", latticeSize<Nd + 1>("cxyzn", {{'n', n_colorvecs}}),
				 OnHost, OnMaster);
      //Tensor<Nd, ComplexF> v("cxyz", latticeSize<Nd>("cxyz"), OnHost, OnMaster);
      for (int t_slice = from_tslice, i_slice = 0, n_tslices = Layout::lattSize()[decay_dir];
	   i_slice < n_tslices; ++i_slice, t_slice = (t_slice + 1) % n_tslices)
      {
	for (int colorvec = 0; colorvec < n_colorvecs; ++colorvec)
	{
	  KeyTimeSliceColorVec_t key(t_slice, colorvec);
	  Tensor<Nd + 1, ComplexF> v = t.kvslice_from_size({{'n', colorvec}}, {{'n', 1}});
	  eigen_source.get(key, v);
	  // v.copyTo(t.kvslice_from_size({{'n', colorvec}})); // t[n=colorvec] = v
	}
	t.copyTo(r.kvslice_from_size({{'t', i_slice}})); // r[t=i_slice] = t
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

#endif // __INCLUDE_SUPERB_CONTRACTIONS__
