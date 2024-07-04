// -*- C++ -*-
/*! \file                                                                    
 * \brief Greedy coloring for Lattice probing
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#include "meas/hadron/greedy_coloring.h"
#include "chromabase.h"
#include "qdp_layout.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace Chroma
{

  namespace
  {
    // Number of dimensions
    constexpr unsigned int D = 4;
    // Array index type
    typedef unsigned int IndexType;
    // Coordinate type
    typedef std::array<IndexType, D> CoorType;
    // Vector of indices
    typedef std::vector<IndexType> Indices;
    // Vector of coordinates
    typedef std::vector<CoorType> Coors;

    // Return the permutation [0,1,...]

// Avoid intel compiler explosion
#ifdef __INTEL_COMPILER
#  pragma intel optimization_level 0
#endif
    CoorType naturalOrder()
    {
      CoorType c;
      for (unsigned int i = 0; i < c.size(); i++)
	c[i] = i;
      return c;
    }

    // Return the permuation [N-1,N-2,...]

// Avoid intel compiler explosion
#ifdef __INTEL_COMPILER
#  pragma intel optimization_level 0
#endif
    CoorType antiNaturalOrder()
    {
      CoorType c;
      for (unsigned int i = 0; i < c.size(); i++)
	c[i] = c.size() - i - 1;
      return c;
    }

    // Return the indices associated to each coordinate
    // \param coors: input coordinates
    // \param dim: lattice dimension
    // \param order: ordering of the dimensions, anti-natural ordering by default
    //
    // Return a vector with indices of the passed coordinates in
    // the same order.
    //
    // NOTE: we used anti-natural order, the last coordinate moves the fastest

    Indices coor2index(const Coors& coors, const CoorType& dim,
		       const CoorType& order = antiNaturalOrder())
    {
      // Quick exit
      if (dim.size() <= 0)
	return Indices();

      // Output array
      Indices indices(coors.size());

      // p(order(i)) = prod(dim(order(1:i)))
      CoorType p;
      p[order[0]] = 1;
      for (unsigned int i = 1; i < dim.size(); i++)
	p[order[i]] = p[order[i - 1]] * dim[order[i - 1]];

// indices(i) = inner product of coor(i) and p
// NOTE: every coordinate value is normalize to modulus the lattice dimension
#ifdef _OPENMP
#  pragma omp parallel for schedule(static)
#endif
      for (unsigned int i = 0; i < coors.size(); i++)
      {
	IndexType r = 0;
	for (unsigned int j = 0; j < dim.size(); j++)
	  r += (coors[i][j] % dim[j]) * p[j];
	indices[i] = r;
      }

      return indices;
    }

    // Return the coordinates associated to each index
    // \param indices: input vertex indices
    // \param dim: lattice dimension
    // \param order: ordering of the dimensions, anti-natural ordering by default
    //
    // Return a vector with the coordinates of the passed indices in
    // the same order.
    //
    // NOTE: we used anti-natural order, the last coordinate moves the fastest

    Coors index2coor(const Indices& indices, const CoorType& dim,
		     const CoorType& order = antiNaturalOrder())
    {
      // Quick exit
      if (dim.size() <= 0)
	return Coors();

      // Output array
      Coors coors(indices.size());

      // p(i) = prod(dim(end:-1:i))
      CoorType p;
      p[order[0]] = 1;
      for (unsigned int i = 1; i < dim.size(); i++)
	p[order[i]] = p[order[i - 1]] * dim[order[i - 1]];

      // coors(i,j) = indices(i) / p(i)
      // NOTE: every coordinate value is normalize to modulus the lattice dimension
      for (unsigned int i = 0; i < indices.size(); i++)
	for (unsigned int j = 0; j < dim.size(); j++)
	  coors[i][j] = (indices[i] / p[j]) % dim[j];

      return coors;
    }

    // Return the neighbors of the given vertices
    // \param coors: input vertices
    // \param dim: lattice dimensions
    //
    // Return the neighbors' coordinates of the vertices 'coors'.

    Coors neighbors(const Coors coors, const CoorType& dim)
    {
      if (dim.size() <= 0)
	return Coors();

      // Number of neighbors
      int k = 0;
      for (unsigned int i = 0; i < dim.size(); i++)
	k += (dim[i] == 2 ? 1 : (dim[i] > 2 ? 2 : 0));
      Coors b(coors.size() * k);
      unsigned int j = 0;
      for (CoorType coor : coors)
      {
	for (unsigned int i = 0; i < dim.size(); i++)
	{
	  if (dim[i] >= 2)
	  {
	    b[j] = coor;
	    b[j][i] = (b[j][i] + 1) % dim[i];
	    j++;
	  }
	  if (dim[i] > 2)
	  {
	    b[j] = coor;
	    b[j][i] = (dim[i] + b[j][i] - 1) % dim[i];
	    j++;
	  }
	}
      }
      assert(j == b.size());
      return b;
    }

    // Return the number of vertices in a lattice
    // \param dim: lattice dimensions

    std::size_t volume(const CoorType& dim)
    {
      if (dim.size() <= 0)
	return 0;

      std::size_t vol = dim[0];
      for (unsigned int i = 1; i < dim.size(); ++i)
	vol *= dim[i];
      return vol;
    }

    // Sort a vector and remove duplicate values
    // \param v: input vector

    template <typename T>
    std::vector<T> unique_and_sort(const std::vector<T>& v)
    {
      std::vector<T> r(v);
      std::sort(r.begin(), r.end());
      r.erase(std::unique(r.begin(), r.end()), r.end());
      return r;
    }

    // Add c to all passed coordinates
    // \param c: single coordinate
    // \param coors: various coordinate differences
    // \param dim: lattice dimension
    //
    // This function is used to get all the neighbors of vertex 'c' if coors vector
    // comes from 'neighbors_upto_distance'.
    //
    // NOTE: every returned coordinate value is normalize to modulus the lattice dimension

    Coors plus(const CoorType& c, const Coors& coors, const CoorType& dim)
    {
      Coors r(coors.size());
#ifdef _OPENMP
#  pragma omp parallel for schedule(static)
#endif
      for (unsigned int i = 0; i < coors.size(); ++i)
	for (unsigned int j = 0; j < dim.size(); ++j)
	  r[i][j] = (dim[j] + coors[i][j] + c[j]) % dim[j];
      return r;
    }

    IndexType normalize_coor(IndexType coor, IndexType dim)
    {
      return (dim == 0 ? 0 : coor % dim);
    }

// Avoid intel compiler explosion
#ifdef __INTEL_COMPILER
#  pragma intel optimization_level 0
#endif
    CoorType normalize_coor(const CoorType& coor, const CoorType& dim)
    {
      CoorType r;
      for (unsigned int i = 0; i < coor.size(); ++i)
	r[i] = normalize_coor(coor[i], dim[i]);
      return r;
    }

    IndexType euclidian_dist_squared(const CoorType& a, const CoorType& b, const CoorType& dim)
    {
      CoorType d;
      for (unsigned int i = 0; i < d.size(); ++i)
	d[i] = normalize_coor(dim[i] + a[i] - b[i], dim[i]);
      for (unsigned int i = 0; i < d.size(); ++i)
	d[i] = std::min(d[i], dim[i] - d[i]);
      IndexType dist = 0;
      for (unsigned int i = 0; i < d.size(); ++i)
	dist += d[i] * d[i];
      return dist;
    }

// Avoid intel compiler explosion
#ifdef __INTEL_COMPILER
#  pragma intel optimization_level 0
#endif
    CoorType minus(const CoorType& a, const CoorType& b)
    {
      CoorType r;
      for (unsigned int i = 0; i < a.size(); ++i)
	r[i] = a[i] - b[i];
      return r;
    }

    // Return all neighbors up to a given distance WITH THE SAME PARITY
    // \param dists: shifts to consider
    // \param power: all neighbors up to this distance
    //
    // Return a vector of coordinate differences to all neighbors up to
    // the given distance (power) and the same parity as the neighbors at the
    // given distance. The returned vector is used as input to the
    // function 'plus' to get all neighbors to a given vertex.
    //
    // NOTE: no performance requirements for this function; the function
    //       'plus' is doing the heavy lifting.

    Coors neighbors_upto_distance(const Coors& dists, unsigned int power, const CoorType& dim)
    {
      // Find all neighbors of the vertex at origin up to the
      // given distances; a regular code should do something like:
      //   for i=1:dist, vertices=union(vertices, neighbors(vertices))
      // BUT get_motive ONLY CARES ABOUT EVEN VERTICES
      Coors centers;
      for (const auto& dist : dists)
      {
	centers.push_back(dist);
	centers.push_back(normalize_coor(minus(dim, dist), dim));
      }
      Coors neighbors_pattern = centers, prev;
      for (unsigned int i = 0; i < power; i++)
      {
	prev = neighbors_pattern;
	neighbors_pattern = neighbors(neighbors_pattern, dim);
	neighbors_pattern = index2coor(unique_and_sort(coor2index(neighbors_pattern, dim)), dim);
      }
      neighbors_pattern.insert(neighbors_pattern.end(), prev.begin(), prev.end());

      // Filter our neighbors further than power distance in euclidian metric from the centers
      Coors filter;
      for (const CoorType& c : neighbors_pattern)
	if (coor2index(Coors(1, c), dim)[0] != 0 &&
	    (euclidian_dist_squared(c, centers[0], dim) <= power * power ||
	     euclidian_dist_squared(c, centers[1], dim) <= power * power))
	  filter.push_back(c);

      // Remove duplicate vertices
      return index2coor(unique_and_sort(coor2index(filter, dim)), dim);
    }

    // Return k-distance coloring
    // \param dist: distance, k
    // \param dim: lattice dimensions
    // \param num_colors: return the number of colors
    //
    // Return a vector with color value for every vertex following
    // the ordering of coor2index.
    //
    // Strategy:
    // * First color even vertices with greedy coloring
    // * Then color odd vertices copying the coloring of the even vertices

    Indices get_colors(const Coors& dists, unsigned int power, const CoorType& dim,
		       unsigned int& num_colors)
    {
      Coors neighbors_rel = neighbors_upto_distance(dists, power, dim);

      const unsigned int vol = volume(dim);
      std::vector<unsigned int> color(vol);

      // For some reason, natural ordering works better for the tested lattices;
      // also, it seems beneficial to swap the first and the last dimensions if
      // the latter is the smallest
      CoorType perm = naturalOrder();
      if (dim[0] > dim[3])
	std::swap(perm[0], perm[3]);

      // First color even vertices with greedy coloring
      num_colors = 0;
      for (int oddity = 0; oddity < 2; oddity++)
      {
	for (unsigned int i = 0; i < vol; i++)
	{
	  // Get coordinates of vertex with index i
	  CoorType c = index2coor(Indices(1, i), dim, perm)[0];

	  // Only process even vertices: mod(sum(c),2) == 0
	  if (std::accumulate(c.begin(), c.end(), 0) % 2 != oddity)
	    continue;

	  // Get neighbors of c
	  Indices c_neighbors = coor2index(plus(c, neighbors_rel, dim), dim);

	  // Get the colors of the c's neighbors
	  std::vector<bool> used_colors(num_colors + 1);
	  for (unsigned int j = 0; j < c_neighbors.size(); ++j)
	    used_colors[color[c_neighbors[j]]] = true;

	  // Find the first free color
	  // NOTE: color 0 means no color assigned
	  unsigned int j;
	  for (j = 1; j <= num_colors && used_colors[j]; j++)
	    ;

	  // Assign the free color
	  unsigned int idx = coor2index(Coors(1, c), dim)[0];
	  color[idx] = j;

	  // Track number of used colors
	  num_colors = std::max(num_colors, color[idx]);
	}
      }

      // Make zero the first color index
      for (unsigned int i = 0; i < vol; i++)
	color[i]--;

      //assert(check_coloring(color, dist, power, dim));
      return color;
    }

    bool check_coloring(const Indices& color, const Coors& dists, unsigned int power,
			const CoorType& dim)
    {
      Coors neighbors_rel = neighbors_upto_distance(dists, power, dim);

      const unsigned int vol = volume(dim);
      for (unsigned int i = 0; i < vol; i++)
      {
	// Process element i
	CoorType c = index2coor(Indices(1, i), dim)[0];

	// Get neighbors of c
	Indices c_neighbors = coor2index(plus(c, neighbors_rel, dim), dim);

	// Check that c has different color than its neighbors
	for (unsigned int j = 0; j < c_neighbors.size(); ++j)
	  if (c_neighbors[j] != i && color[i] == color[c_neighbors[j]])
	    return false;
      }

      return true;
    }
  }

  // Construct a k-distance coloring
  void Coloring::construct(const std::vector<std::array<int, 4>>& distances, unsigned int power,
			   const CoorType& latt_size, bool build_local)
  {
    // Get the absolute value of the distances
    Coors abs_distances;
    for (const auto& dist : distances)
      abs_distances.push_back(
	CoorType{(unsigned int)std::abs(dist[0]), (unsigned int)std::abs(dist[1]),
		 (unsigned int)std::abs(dist[2]), (unsigned int)std::abs(dist[3])});

    // Compute the maximum shift/distance requested
    CoorType max_distance{{}};
    for (const auto& dist : abs_distances)
      for (unsigned int i = 0; i < dist.size(); ++i)
	max_distance[i] = std::max(max_distance[i], dist[i]);

    // Compute the tile size; the tile size should be divisible by the lattice size and
    // greater or equal than 2*(dist+power)
    for (unsigned int i = 0; i < latt_size.size(); i++)
    {
      tile_size[i] = std::min(2 * (max_distance[i] + power), latt_size[i]);
      while (latt_size[i] % tile_size[i] != 0)
	tile_size[i]++;
    }

    // Get colors for all nodes
    colors = get_colors(abs_distances, power, tile_size, num_colors);

    if (build_local)
    {
      // Store the colors of the local nodes
      int this_node = Layout::nodeNumber();
      local_colors.resize(Layout::sitesOnNode());
      for (unsigned int i = 0; i < Layout::sitesOnNode(); i++)
      {
	// Local coordinates of node i
	multi1d<int> x = Layout::siteCoords(this_node, i);

	CoorType c;
	for (unsigned int j = 0; j < c.size(); j++)
	  c[j] = x[j];

	local_colors[i] = colors[coor2index(Coors(1, c), tile_size)[0]];
      }
    }
  }

  // Construct a k-distance coloring
  Coloring::Coloring(const std::vector<std::array<int, 4>>& distances, unsigned int power)
  {
    // Get lattice dimensions
    CoorType latt_size;
    for (unsigned int i = 0; i < latt_size.size(); i++)
      latt_size[i] = Layout::lattSize()[i];

    construct(distances, power, latt_size, true);
  }

  // Read the coloring from a file
  Coloring::Coloring(const std::string& filename)
  {
    // Get lattice dimensions
    CoorType latt_size;
    for (unsigned int i = 0; i < latt_size.size(); i++)
      latt_size[i] = Layout::lattSize()[i];

    // Read colors for all nodes from 'filename' which are
    // stored in natural order
    std::size_t vol = volume(latt_size);
    Indices colors(vol);
    {
      std::ifstream f(filename);
      IndexType idx;
      std::size_t i = 0;
      for (i = 0; f >> idx; ++i)
      {
	if (i >= vol)
	{
	  throw std::runtime_error("The coloring file has too many rows!");
	}
	colors[coor2index(index2coor(Indices(1, i), latt_size, naturalOrder()), latt_size)[0]] =
	  idx - 1;
      }
      if (i < vol)
      {
	throw std::runtime_error("Missing rows from the coloring file!");
      }
    }

    num_colors = *std::max_element(colors.begin(), colors.end()) + 1;

    // Store the colors of the local nodes
    int this_node = Layout::nodeNumber();
    local_colors.resize(Layout::sitesOnNode());
    for (unsigned int i = 0; i < Layout::sitesOnNode(); i++)
    {
      // Local coordinates of node i
      multi1d<int> x = Layout::siteCoords(this_node, i);

      CoorType c;
      for (unsigned int j = 0; j < c.size(); j++)
	c[j] = x[j];

      local_colors[i] = colors[coor2index(Coors(1, c), latt_size)[0]];
    }
  }

  // Return a probing vector for the given color
  void Coloring::getVec(LatticeInteger& vec, unsigned int color) const
  {
    if (color >= num_colors)
      throw std::runtime_error("Invalid color value");
    if (local_colors.size() == 0)
      throw std::runtime_error("Invalid function");

    int node = Layout::nodeNumber();
    for (int s(0); s < Layout::sitesOnNode(); s++)
    {
      multi1d<int> x = Layout::siteCoords(node, s);
      Integer v = local_colors[s] == color ? 1 : 0;
      pokeSite(vec, v, x);
    }
  }

  // Return the color for a site
  unsigned int Coloring::getColor(const std::array<int, 4>&  coor) const
  {
    CoorType coor0{(unsigned int)coor[0], (unsigned int)coor[1], (unsigned int)coor[2],
		   (unsigned int)coor[3]};
    return colors[coor2index(Coors(1, coor0), tile_size)[0]];
  }

  // Return all neighbors
  std::vector<std::array<int, 4>> Coloring::all_neighbors(unsigned int farthest_neighbor,
							  const std::array<int, 4>& dim)
  {
    std::array<unsigned int, 4> dimu{(unsigned int)dim[0], (unsigned int)dim[1],
				     (unsigned int)dim[2], (unsigned int)dim[3]};
    Coors neighbors_coors(1);
    for (unsigned int i = 0; i < farthest_neighbor; i++)
    {
      auto new_neighbors_coors = neighbors(neighbors_coors, dimu);
      neighbors_coors.insert(neighbors_coors.end(), new_neighbors_coors.begin(),
			     new_neighbors_coors.end());
      neighbors_coors = index2coor(unique_and_sort(coor2index(neighbors_coors, dimu)), dimu);
    }
    std::vector<std::array<int, 4>> r;
    for (const auto& i : neighbors_coors)
      r.push_back(std::array<int, 4>{(int)i[0], (int)i[1], (int)i[2], (int)i[3]});
    return r;
  }
}
