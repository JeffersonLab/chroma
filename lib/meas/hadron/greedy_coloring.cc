// -*- C++ -*-                                                                   
/*! \file                                                                    
 * \brief Greedy coloring for Lattice probing
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#include "chromabase.h"
#include "meas/hadron/greedy_coloring.h"
#include "qdp_layout.h"
#include <stdexcept>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <iterator>

namespace Chroma {

namespace {
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

	// Return the indices associated to each coordinate
	// \param coors: input coordinates
	// \param dim: lattice dimension
	//
	// Return a vector with indices of the passed coordinates in
	// the same order.
	//
	// NOTE: we used anti-natural order, the last coordinate moves the fastest

	Indices coor2index(const Coors &coors, const CoorType dim) {
		// Quick exit
		if (dim.size() <= 0) return Indices();

		// Output array
		Indices indices(coors.size());

		// p(i) = prod(dim(end:-1:i))
		CoorType p;
		p.back() = 1;
		for (int i=p.size()-1; i>=1; i--) p[i-1] = p[i]*dim[i];

		// indices(i) = inner product of coor(i) and p
		// NOTE: every coordinate value is normalize to modulus the lattice dimension
		for (unsigned int i=0; i<coors.size(); i++) {
			IndexType r = 0;
			for (unsigned int j=0; j<dim.size(); j++) r += (coors[i][j] % dim[j]) * p[j];
			indices[i] = r;
		}

		return indices;
	}

	// Return the coordinates associated to each index
	// \param indices: input vertex indices
	// \param dim: lattice dimension
	//
	// Return a vector with the coordinates of the passed indices in
	// the same order.
	//
	// NOTE: we used anti-natural order, the last coordinate moves the fastest

	Coors index2coor(const Indices &indices, const CoorType dim) {
		// Quick exit
		if (dim.size() <= 0) return Coors();

		// Output array
		Coors coors(indices.size());

		// p(i) = prod(dim(end:-1:i))
		CoorType p;
		p.back() = 1;
		for (int i=p.size()-1; i>=1; i--) p[i-1] = p[i]*dim[i];

		// coors(i,j) = indices(i) / p(i)
		// NOTE: every coordinate value is normalize to modulus the lattice dimension
		for (unsigned int i=0; i<indices.size(); i++)
			for (unsigned int j=0; j<dim.size(); j++)
				coors[i][j] = (indices[i] / p[j]) % dim[j];

		return coors;
	}

	// Return the neighbors of the given vertices
	// \param coors: input vertices
	// \param dim: lattice dimensions
	//
	// Return the neighbors' coordinates of the vertices 'coors'.

	Coors neighbors(const Coors coors, const CoorType dim) {
		if (dim.size() <= 0) return Coors();

		// Number of neighbors
		int k = 0;
		for (unsigned int i=0; i<dim.size(); i++) k += (dim[i] == 2 ? 1 : (dim[i] > 2 ? 2 : 0));
		Coors b(coors.size() * k);
		unsigned int j = 0;
		for (CoorType coor: coors) {
			for (unsigned int i=0; i<dim.size(); i++) {
				if (dim[i] >= 2) {
					b[j] = coor;
					b[j][i] = (b[j][i] + 1) % dim[i];
					j++;
				}
				if (dim[i] > 2) {
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

	std::size_t volume(const CoorType dim) {
		if (dim.size() <= 0) return 0;

		std::size_t vol = dim[0];
		for (unsigned int i=1; i<dim.size(); ++i)
			vol *= dim[i];
		return vol;
	}

	// Sort a vector and remove duplicate values
	// \param v: input vector

	template<typename T> std::vector<T> unique_and_sort(const std::vector<T>& v) {
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

	Coors plus(const CoorType& c, const Coors& coors, const CoorType& dim) {
		Coors r(coors.size());
		for (unsigned int i=0; i<coors.size(); ++i)
			for (unsigned int j=0; j<dim.size(); ++j)
				r[i][j] = (dim[j] + coors[i][j] + c[j]) % dim[j];
		return r;
	}

	// Return all neighbors up to a given distance WITH THE SAME PARITY
	// \param dist: distance
	//
	// Return a vector of coordinate differences to all neighbors up to
	// the given distance and the same parity as the neighbors at the
	// given distance. The returned vector is used as input to the
	// function 'plus' to get all neighbors to a given vertex.
	//
	// NOTE: no performance requirements for this function; the function
	//       'plus' is doing the heavy lifting.

	Coors neighbors_upto_distance(unsigned int dist, const CoorType& dim) {
		// Find all neighbors of the vertex at origin up to the
		// given distance; a regular code should do something like:
		//   for i=1:dist, vertices=union(vertices, neighbors(vertices))
		// BUT get_motive ONLY CARES ABOUT EVEN VERTICES
		Coors neighbors_pattern(1);
		for (unsigned int i=0; i<dist; i++)
			neighbors_pattern = neighbors(neighbors_pattern, dim);

		// Remove duplicate vertices
		neighbors_pattern = index2coor(unique_and_sort(coor2index(neighbors_pattern, dim)), dim);

		return neighbors_pattern;
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

	Indices get_motive(unsigned int dist, const CoorType& dim, unsigned int& num_colors) {
		dist = dist/2*2;
		Coors neighbors_rel = neighbors_upto_distance(dist, dim);

		const unsigned int vol = volume(dim);
		std::vector<unsigned int> color(vol);

		// First color even vertices with greedy coloring
		num_colors = 0;
		for (unsigned int i=0; i<vol; i++) {
			// Get coordinates of vertex with index i
			CoorType c = index2coor(Indices(1, i), dim)[0]; 

			// Only process even vertices: mod(sum(c),2) == 0
			if (std::accumulate(c.begin(), c.end(), 0) % 2 != 0) continue;

			// Get neighbors of c
			Indices c_neighbors = coor2index(plus(c, neighbors_rel, dim), dim);

			// Get the colors of the c's neighbors
			std::vector<bool> used_colors(num_colors+1);
			for (unsigned int j=0; j<c_neighbors.size(); ++j)
				used_colors[color[c_neighbors[j]]] = true;

			// Find the first free color
			// NOTE: color 0 means no color assigned
			unsigned int j;
			for (j=1; j<=num_colors && used_colors[j]; j++);

			// Assign the free color
			color[i] = j;

			// Track number of used colors
			num_colors = std::max(num_colors, color[i]);
		}

		// Color odd vertices as the color of the forward neighbor on
		// the first direction
		for (unsigned int i=0; i<vol; i++) {
			// Get coordinates of vertex i
			CoorType c = index2coor(Indices(1, i), dim)[0];

			// Only process odd vertices
			if (std::accumulate(c.begin(), c.end(), 0) % 2 != 1) continue;

			// Get the color of the forward neighbor at the first direction
			c[0]++;
			color[i] = color[coor2index(Coors(1, c), dim)[0]] + num_colors;
		}
 		num_colors *= 2;

		// Make zero the first color index
		for (unsigned int i=0; i<vol; i++) color[i]--;

		//assert(check_coloring(color, dist, dim));
		return color;
	}

	// Return k-distance coloring
	// \param dist: distance, k
	// \param dim: lattice dimensions
	//
	// Return a vector with color value for every vertex following
	// the anti-natural ordering.
	//
	// Strategy:
	// * First color a tile of dimension 2^dist
	// * Then color the lattice following the tile

	Indices get_colors(unsigned int dist, const CoorType& dim, unsigned int& num_colors) {
		// The tile has dimension 2^dist if 2^dist is divisible by the original lattice
		// dimension; otherwise the tile is as large as the lattice dimension 
		CoorType dim_motive;
		unsigned int d = (1 << dist);
		for (unsigned int i=0; i<dim.size(); i++)
			dim_motive[i] = (dim[i] % d == 0 ? d : dim[i]);

		// dist-coloring the tile
		Indices color_motive = get_motive(dist, dim_motive, num_colors);

		// Color the original lattice following the coloring of the tile
		const unsigned int vol = volume(dim);
		std::vector<unsigned int> color(vol);
		for (unsigned int i=0; i<vol; i++) {
			color[i] = color_motive[coor2index(index2coor(Indices(1, i), dim), dim_motive)[0]];
		}

		return color;
	}

	bool check_coloring(const Indices& color, unsigned int dist, const CoorType& dim) {
		Coors neighbors_rel = neighbors_upto_distance(dist, dim);

		const unsigned int vol = volume(dim);
		for (unsigned int i=0; i<vol; i++) {
			// Process element i
			CoorType c = index2coor(Indices(1, i), dim)[0];

			// Get neighbors of c
			Indices c_neighbors = coor2index(plus(c, neighbors_rel, dim), dim);

			// Check that c has different color than its neighbors
			for (unsigned int j=0; j<c_neighbors.size(); ++j)
				if (c_neighbors[j] != i && color[i] == color[c_neighbors[j]])
					return false;
		}

		return true;
	}
}

// Construct a k-distance coloring
Coloring::Coloring(unsigned int distance) {
	// Get lattice dimensions
	CoorType latt_size;
	for (unsigned int i=0; i<latt_size.size(); i++)
		latt_size[i] = Layout::lattSize()[i];
	
	// Get colors for all nodes
	Indices colors = get_colors(distance, latt_size, num_colors);

	// Store the colors of the local nodes
	int this_node = Layout::nodeNumber();
	local_colors.resize(Layout::sitesOnNode());
	for (unsigned int i=0; i<Layout::sitesOnNode(); i++) {
		// Local coordinates of node i
		multi1d<int> x = Layout::siteCoords(this_node, i);

		CoorType c;
		for (unsigned int j=0; j<c.size(); j++) c[j] = x[j];

		local_colors[i] = colors[coor2index(Coors(1, c), latt_size)[0]];
	}
}

// Return a probing vector for the given color
void Coloring::getVec(LatticeInteger& vec, unsigned int color) const {
	if (color >= num_colors)
		throw std::runtime_error("Invalid color value");

	int node = Layout::nodeNumber();
	for( int s(0); s<Layout::sitesOnNode(); s++) {
		multi1d<int> x = Layout::siteCoords(node,s) ;
		Integer v = local_colors[s] == color ? 1 : 0;
		pokeSite(vec,v,x);
	}
}

}
