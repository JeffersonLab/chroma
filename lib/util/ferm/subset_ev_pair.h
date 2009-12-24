// -*- C++ -*-
/*! \file
 *  \brief Holds of vectors and eigenvalues
 */

#ifndef __EV_PAIR_H__
#define __EV_PAIR_H__


#include "chromabase.h"

namespace Chroma 
{

  //! Weights for subset of vectors
  /*! \ingroup ferm */
  struct SubsetVectorWeight_t
  {
    multi1d<Real> weights;
  };

  //! A Pair type.
  template<typename T>
  struct EVPair {
    SubsetVectorWeight_t eigenValue;
    T                    eigenVector;
  };

  template<typename T>
  void read(BinaryReader& bin_in, EVPair<T>& evpair)
  {
    read(bin_in, evpair.eigenVector);
    read(bin_in, evpair.eigenValue.weights);
  }

  template<typename T>
  void write(BinaryWriter& bin_out, const EVPair<T>& evpair)
  {
    write(bin_out, evpair.eigenVector);
    write(bin_out, evpair.eigenValue.weights);
  }

};


#endif

