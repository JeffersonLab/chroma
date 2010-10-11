// -*- C++ -*-
/*! \file
 *  \brief Holds of vectors and weights
 */

#ifndef __subset_vectors_h__
#define __subset_vectors_h__

#include "chromabase.h"
#include "util/ferm/subset_ev_pair.h"
#include "qdp_map_obj.h"


namespace Chroma
{
 
  //! Reader
  /*! \ingroup ferm */
  void read(XMLReader& xml, const std::string& path, SubsetVectorWeight_t& param);
 
  //! Writer
  /*! \ingroup ferm */
  void write(XMLWriter& xml, const std::string& path, const SubsetVectorWeight_t& param);


  //! Extract eigenvalues from a map, and arrange them in a different format
  multi1d<SubsetVectorWeight_t> getEigenValues(const QDP::MapObject<int,EVPair<LatticeColorVector> >& eigen_source, int num_vecs);
  
  
}


#endif
