#include "chromabase.h"
#include "io/monomial_io.h"
#include "update/molecdyn/monomial/monomial_factory.h"


namespace Chroma { 

void read(XMLReader& xml, 
	  const std::string& path, 
	  Handle< Monomial< multi1d<LatticeColorMatrix>, 
	                    multi1d<LatticeColorMatrix> > >& mon_handle )
{
  XMLReader paramtop(xml, path);
  std::string monomial_name;
  try { 
    read( paramtop, "./Name", monomial_name);
    mon_handle = TheMonomialFactory::Instance().createObject( monomial_name, xml, path );
  }
  catch( const std::string& e ) { 
    QDPIO::cerr << "Error Reading Monommial: " << e << endl;
    QDP_abort(1);
  }
}

void read(XMLReader& xml, 
	  const std::string& path, 
	  Handle< ExactMonomial< multi1d<LatticeColorMatrix>, 
	                         multi1d<LatticeColorMatrix> > >& mon_handle )
{

  XMLReader paramtop(xml, path);
  std::string monomial_name;
  Monomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* m;

  try { 
    read( paramtop, "./Name", monomial_name);
    m = TheMonomialFactory::Instance().createObject( monomial_name, xml, path );
  }
  catch( const std::string& e ) { 
    QDPIO::cerr << "Error Reading Monommial: " << e << endl;
    QDP_abort(1);
  }

  ExactMonomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* em;
  em = dynamic_cast< ExactMonomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* >(m);
  if( em == 0 ) { 
    QDPIO::cerr << "Failed to downcast monomial to exact monomial " << endl;
    QDP_abort(1);
  }

  mon_handle = em;

}

};
