/*! \file
 *  \brief Holds of vectors and weights
 */

#include "chromabase.h"
#include "util/ferm/subset_vectors.h"

using namespace QDP;

namespace Chroma 
{

  // Reader
  void read(XMLReader& xml, const std::string& path, SubsetVectorWeight_t& param)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, path, param.weights);
  }
 
  // Writer
  void write(XMLWriter& xml, const std::string& path, const SubsetVectorWeight_t& param)
  {
    //push(xml, path);
    write(xml, path, param.weights);
    //pop(xml);
  }


  // Extract eigenvalues from a map, and arrange them in a different format
  multi1d<SubsetVectorWeight_t> getEigenValues(const MapObject<int,EVPair<LatticeColorVector> >& eigen_source, int num_vecs)
  {
    multi1d<SubsetVectorWeight_t> evals(num_vecs);

    // Basic sanity check
    if (num_vecs > eigen_source.size())
    {
      QDPIO::cerr << __func__ << ": requesting more vectors than are in the eigenvector map\n";
      QDP_abort(1);
    }

    // A straight lookup is fine - will fail if not right number of eigenvalues
    for(int i=0; i < num_vecs; ++i)
    {
      EVPair<LatticeColorVector> tmpvec; eigen_source.get(i, tmpvec);

      evals[i] = tmpvec.eigenValue;
    }

    return evals;
  }
  
}  // end namespace Chroma
