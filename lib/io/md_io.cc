#include "chromabase.h"
#include "io/md_io.h"

using namespace QDP;
using namespace std;

void read(XMLReader& xml, const string& path, MDIntegratorType_t& t)
{
  string token;

  try { 
    read(xml, path, token);
  }
  catch(const string& e ) {
    QDPIO::cout << "Caught exception while reading XML " << e << endl;
    QDP_abort(1);
  }

  if ( token == "PQP_LEAPFROG" ) { 
    t = MD_PQP_LEAPFROG;
  }
  else if ( token == "QPQ_LEAPFROG" ) { 
    t = MD_QPQ_LEAPFROG;
  }
  else { 
    QDPIO::cerr << " Unsupported MDIntegratorType_t: " << token << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const MDIntegratorType_t& t) 
{
  string token; 
  
  switch( t ) { 
  case MD_PQP_LEAPFROG: 
    token = "PQP_LEAPFROG";
    break;
  case MD_QPQ_LEAPFROG:
    token = "QPQ_LEAPFROG";
    break;
  default:
    QDPIO::cerr<< "Unsupported MDIntegratorType_t : " << t << endl;
    QDP_abort(1);
  }

  try { 
    write(xml, path, token);
  }
  catch( const string& e ) { 
    QDPIO::cerr << "Caught exception while writing XML: " << e << endl;
    QDP_abort(1);
  }
}

LeapfrogParams::LeapfrogParams(XMLReader& top) 
{
  try { 
    read(top, "./IntegratorType", my_type);
    
    switch( my_type ) { 
    case MD_PQP_LEAPFROG: 
      break;
    case MD_QPQ_LEAPFROG:
      break;
    default:
      QDPIO::cerr << "Internal error: Trying to read non leapfrog params with Leapfrog params" << endl;
      QDP_abort(1);
      break;
    }

    read(top, "./dt", dt);
    read(top, "./tau", tau);
  }
  catch( const string& e ) { 
    QDPIO::cerr << "Caught exception while reading XML: "<< e << endl;
    QDP_abort(1);
  }
}
  
void write(XMLWriter& xml, const string& path, const LeapfrogParams& p) 
{
  try { 
    push(xml, path);
    write(xml, "IntegratorType", p.getType());
    write(xml, "dt", p.getStepSize());
    write(xml, "tau", p.getTrajLength());
    pop(xml);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception while writing XML: " << e << endl;
    QDP_abort(1);
  }
}

MDIntegratorParamsBase* readMDIntegratorParams(XMLReader& xml, const string& path)
{
  MDIntegratorType_t my_type;
  try {
    XMLReader top(xml, path);
  
    read(top, "./IntegratorType", my_type);

    switch(my_type) { 
    case MD_PQP_LEAPFROG:
      return new LeapfrogParams(top);
      break;
    case MD_QPQ_LEAPFROG:
      return new LeapfrogParams(top);
      break;
    default:
      QDPIO::cerr << "Unsupported MDIntegratorType_t " << my_type << endl;
      QDP_abort(1);
    }
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception while reading XML " << e << endl;
    QDP_abort(1);
  }

  return 0;
}

void write(XMLWriter& xml, const string& path, const MDIntegratorParamsBase& p) 
{
  switch(p.getType()) { 
  case MD_PQP_LEAPFROG:
    // Intentional fall through
  case MD_QPQ_LEAPFROG: 
    {
      const LeapfrogParams& lp=dynamic_cast<const LeapfrogParams&>(p);
      write(xml, path, lp);
    }
    break;
  default:
    QDPIO::cerr << "Unknown MDIntegratorType_t: " << p.getType() << endl;
    QDP_abort(1);
  }
}

  
