
#include "chromabase.h"
#include "io/bc_io.h"

using namespace QDP;
using namespace std;

void read(XMLReader& xml, const string& path, SimpleBCType_t& t)
{
  START_CODE();

  string token;

  try { 
    read(xml, path, token);
  }
  catch(const string& e ) { 
    QDPIO::cerr << "Caught exception while reading XML: " << e << endl;
    QDP_abort(1);
  }

  if( token == "ANTIPERIODIC" ) { 
    t = BC_ANTIPERIODIC;
  }
  else if( token == "DIRICHLET" ) { 
    t = BC_DIRICHLET;
  }
  else if( token == "PERIODIC" ) { 
    t = BC_PERIODIC;
  }
  else { 
    QDPIO::cerr << token << " is not a supported simple BC type " << endl;
    QDP_abort(1);
  }

  END_CODE();
}

void write(XMLWriter& xml, const string& path, const SimpleBCType_t& t)
{
  START_CODE();

  string token;
  switch( t ) { 
  case BC_ANTIPERIODIC: 
    token = "ANTIPERIODIC";
    break;
  case BC_DIRICHLET:
    token = "DIRICHLET";
    break;
  case BC_PERIODIC:
    token = "PERIODIC";
    break;
  default:
    QDPIO::cerr << "Unsupported simple BC type " << endl;
    QDP_abort(1);
    break;
  }

  try {
    write(xml, path, token);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception writing XML: " << e << endl;
    QDP_abort(1);
  }
  
 END_CODE();
}
