// -*- C++ -*-                                                                   
/*! \file                                                                    
 * \brief Greedy coloring for Lattice probing
 *                                                                             
 * Hadron spectrum calculations utilities
 */

#ifndef __INCLUDE_GREEDY_COLORING__
#define __INCLUDE_GREEDY_COLORING__

#include "chromabase.h"

namespace Chroma {

// Interface for computes distance-k coloring for toroidal lattices
struct Coloring {
	// Construct a k-distance coloring
	Coloring(unsigned int distance);
	// Reading the coloring from a file
	Coloring(const std::string& filename);

	// Return a probing vector for the given color
	void getVec(LatticeInteger& vec, unsigned int color) const;

	// Return the number of colors
	unsigned int numColors() const { return num_colors; }

	private:
	multi1d<unsigned int> local_colors;
	unsigned int num_colors;
};

}
#endif // __INCLUDE_GREEDY_COLORING__
