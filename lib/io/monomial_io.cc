#include "chromabase.h"
#include "io/monomial_io.h"
#include "update/molecdyn/monomial_factory.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma { 

void read(XMLReader& xml, 
	  const std::string& path, 
	  Handle< ExactMonomial< multi1d<LatticeColorMatrix>, 
	                         multi1d<LatticeColorMatrix> > >& mon_handle )
{
  XMLReader paramtop(xml, path);
  std::string monomial_name;
  try { 
    read( paramtop, "./Name", monomial_name);
    mon_handle = TheExactMonomialFactory::Instance().createObject( monomial_name, xml, path );
  }
  catch( const std::string& e ) { 
    QDPIO::cerr << "Error Reading Exact Monommial: " << e << endl;
    QDP_abort(1);
  }
}

};
