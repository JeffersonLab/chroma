#ifndef __QUDA_EIGENVECTOR_PARAMS_H__
#define __QUDA_EIGENVECTOR_PARAMS_H__

#include "chromabase.h"
#include <quda.h>
#include <enum_quda.h>
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"
#include "actions/ferm/invert/quda_solvers/xml_array_utils.h"



namespace Chroma
{



  // Need maps to read QudaEigType and QudaSpectrumType
  void read(XMLReader& xml_in, const std::string& path, QudaEigParam& p);
  void write(XMLWriter& xml_in, const std::string& path, const QudaEigParam& p);

}
#endif
