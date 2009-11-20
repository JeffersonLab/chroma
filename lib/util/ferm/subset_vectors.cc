// $Id: subset_vectors.cc,v 1.2 2008-09-26 01:19:41 jbulava Exp $
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

}  // end namespace Chroma
