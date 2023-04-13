// -*- C++ -*-
/*! \file                                                                    
 * \brief Greedy coloring for Lattice probing
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_GREEDY_COLORING__
#define __INCLUDE_GREEDY_COLORING__

#include "chromabase.h"

#include <array>
#include <vector>

namespace Chroma
{

  // Interface for computes distance-k coloring for toroidal lattices
  struct Coloring {
    // Construct a k-distance coloring
    Coloring(unsigned int distance, unsigned int power);
    // Construct a k-distance coloring
    Coloring(unsigned int distance, unsigned int power, std::array<int, 4> dim)
    {
      std::array<unsigned int, 4> dimu{(unsigned int)dim[0], (unsigned int)dim[1],
				       (unsigned int)dim[2], (unsigned int)dim[3]};
      construct(distance, power, dimu, false);
    }

    // Reading the coloring from a file
    Coloring(const std::string& filename);

    // Return a probing vector for the given color
    void getVec(LatticeInteger& vec, unsigned int color) const;

    // Return the color for each node
    unsigned int getColor(std::array<int, 4> dim) const;

    // Return the number of colors
    unsigned int numColors() const
    {
      return num_colors;
    }

    static std::vector<std::array<int, 4>> all_neighbors(unsigned int farthest_neighbor,
							 std::array<int, 4> dim);

  private:
    std::vector<unsigned int> colors;
    std::vector<unsigned int> local_colors;
    std::array<unsigned int, 4> tile_size;
    unsigned int num_colors;
    void construct(unsigned int distance, unsigned int power, std::array<unsigned int, 4> latt_size,
		   bool build_local);
  };

}
#endif // __INCLUDE_GREEDY_COLORING__
