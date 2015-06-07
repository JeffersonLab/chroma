
/*! \file
 *  \brief Read in a configuration written in CERN format.
 */

#ifndef __readcern_h__
#define __readcern_h__

#include "io/cern_io.h"

namespace Chroma {

    //  Read gauge field in CERN format into "u".
    //  Details about CERN format:
    //     - little endian, ints are 4 bytes, doubles are 8 bytes
    //     - 4 ints NT,NX,NY,NZ
    //     - 1 double for average plaquette
    //     - links as SU3 matrices (row major, double precision complex)
    //         in the following order:
    //            - the 8 links in directions +0,-0,...,+3,-3 at the first 
    //              odd point, the second odd point, and so on. 
    //            - The order of the point (x0,x1,x2,x3) with
    //               Cartesian coordinates in the range 
    //                 0<=x0<N0,...,0<=x3<N3 is determined by the index
    //                   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0,
    //               where N0,N1,N2,N3 are the global lattice sizes


void readCERN(multi1d<LatticeColorMatrix>& u, const std::string& cfg_file);


}  // end namespace Chroma

#endif
