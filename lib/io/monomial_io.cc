#include "chromabase.h"
#include "io/monomial_io.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "meas/inline/io/named_objmap.h"

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

//! Read a named monomial from an XML reader, usa factory to create and assign the pointer to a handle in the named object map of monomial handles.
/*! \ingroup io */
void readNamedMonomial(XMLReader& xml,
		       const std::string& path,
		       std::string& monomial_id)

{

  typedef multi1d<LatticeColorMatrix> LCM;

  try { 
    XMLReader paramtop(xml, path);
    
    try { 
      XMLReader named_object_xml(paramtop, "./NamedObject");
      read(named_object_xml, "./monomial_id", monomial_id);
    }
    catch(const std::string& e) { 
      QDPIO::cerr << "Failed to find NamedObject tag or monomial ID in readNamedMonomial" << endl << flush;
      QDP_abort(1);
    }
    
    // Create the ID in the named Object space:
    TheNamedObjMap::Instance().create< Handle<Monomial<LCM,LCM> > >(monomial_id);
    XMLBufferWriter file_xml;
    push(file_xml, "DummyFileXML");
    pop(file_xml);

    XMLBufferWriter record_xml;
    push(record_xml, "Monomial");
    record_xml << paramtop; // Everything in current scope int
    pop(record_xml);

    TheNamedObjMap::Instance().get(monomial_id).setFileXML(file_xml);
    TheNamedObjMap::Instance().get(monomial_id).setRecordXML(record_xml);

    std::string monomial_name;
    read( paramtop, "./Name", monomial_name);

    
    Handle< Monomial<LCM,LCM> > mon_handle( TheMonomialFactory::Instance().createObject( monomial_name, xml, path ) );

    // Here I add the handle to the named object map...
    TheNamedObjMap::Instance().getData< Handle< Monomial<LCM,LCM> > >(monomial_id) = mon_handle;
  }
  catch( const std::string &e) { 
    QDPIO::cerr << "Caught exception with message: " << e << endl;
    QDP_abort(1);
  }
}

void readNamedMonomialArray(XMLReader& xml, 
			    const std::string& path)
{
  try { 
    std::string monomial_id;
    XMLReader arraytop(xml, path);
    int n_items = arraytop.count("./elem");
    for(int i=1; i <= n_items; i++) {
      std::ostringstream os;
      os << "./elem["<< i << "]";
      readNamedMonomial(arraytop, os.str(), monomial_id);
      QDPIO::cout << "Read Monomial with monomial id: " << monomial_id <<endl;

    }
  }
  catch( const std::string& e) {
    QDPIO::cout << "Caught Exception with message: " << e << endl << flush;
    QDP_abort(1);
  }
}

};
