#include "chromabase.h"
#include "io/enum_io/enum_io.h"
#include "io/gaugebc_io.h"
#include <string>

using namespace QDP;
using namespace std;
using namespace Chroma;
    
GaugeBCPeriodicParams::GaugeBCPeriodicParams(XMLReader& xml)
{
  START_CODE();

  GaugeBCType my_type;
  try { 
    read(xml, "BCType", my_type);
  }
  catch(const string& e) { 
    QDPIO::cerr << "CaughtException reading BCType: " << e << endl;
    QDP_abort(1);
  }

  if( my_type != GAUGEBC_ALL_PERIODIC ) { 
    QDPIO::cerr << "Internal error: GaugeBCPeriodicParams trying to read different kind of BC: " << my_type << endl;
    QDP_abort(1);
  }

  END_CODE();
}

void write(XMLWriter& xml, const string& path, const GaugeBCPeriodicParams& p) 
{

  START_CODE();

  try { 
    push(xml, path); 
    write(xml, "BCType", p.getType());
    pop(xml);
  }
  catch(const string& e ) { 
    QDPIO::cerr << "Caught exception while writing XML" << e << endl;
  }

  END_CODE();

}

GaugeBCSimpleParams::GaugeBCSimpleParams(XMLReader& xml)
{
  START_CODE();

  GaugeBCType my_type;
  try { 
    read(xml, "BCType", my_type);
  }
  catch(const string& e) { 
    QDPIO::cerr << "CaughtException reading BCType: " << e << endl;
    QDP_abort(1);
  }

  if( my_type != GAUGEBC_SIMPLE ) { 
    QDPIO::cerr << "Internal error: GaugeBCSimpleParams trying to read different kind of BC: " << my_type << endl;
    QDP_abort(1);
  }

  try { 
    read(xml, "boundary", boundary);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception while reading XML: " << e << endl;
    QDP_abort(1);
  }

  END_CODE();
    
}

void write(XMLWriter& xml, const string& path, const GaugeBCSimpleParams& p) 
{

  START_CODE();

  try { 
    push(xml, path); 
    write(xml, "BCType", p.getType());
    write(xml, "boundary", p.getBoundary());
    pop(xml);
  }
  catch(const string& e ) { 
    QDPIO::cerr << "Caught exception while writing XML" << e << endl;
  }

  END_CODE();

}

GaugeBCSchrParams::GaugeBCSchrParams(XMLReader& xml)
{
  START_CODE();

  GaugeBCType my_type;
  try { 
    read(xml, "BCType", my_type);
  }
  catch(const string& e) { 
    QDPIO::cerr << "CaughtException reading BCType: " << e << endl;
    QDP_abort(1);
  }

  switch( my_type ) { 
  case GAUGEBC_SCHROEDINGER_1LINK:
    type_t = my_type;
    break;
  case GAUGEBC_SCHROEDINGER_2LINK:
    type_t = my_type;
    break;
  default:
    QDPIO::cerr << "Internal error: GaugeBCSchrParams trying to read non Schroedinger BC: " << my_type << endl;
    QDP_abort(1);
  }

  try { 
    read(xml, "SchrFunType", SchrFun);
    read(xml, "SchrPhiMult", SchrPhiMult);
  }
  catch( const string& e) { 
    QDPIO::cerr << "Caught exception while reading XML: " << e << endl;
    QDP_abort(1);
  }

  END_CODE();
}

void write(XMLWriter& xml, const string& path, const GaugeBCSchrParams& p) 
{

  START_CODE();

  try { 
    push(xml, path); 
    write(xml, "BCType", p.getType());
    write(xml, "SchrFunType", p.getSchrFun());
    write(xml, "SchrPhiMult", p.getSchrPhiMult());

    pop(xml);
  }
  catch(const string& e ) { 
    QDPIO::cerr << "Caught exception while writing XML" << e << endl;
  }

  END_CODE();

}


// toplevel write
void write(XMLWriter& xml, const string& path, const GaugeBCParamsBase& p) 
{
  START_CODE();

  GaugeBCType my_type;

  my_type = p.getType();

  switch( my_type ) {
  case GAUGEBC_ALL_PERIODIC: 
    {
      const GaugeBCPeriodicParams& pp = 
	dynamic_cast<const GaugeBCPeriodicParams&>(p);

      write(xml, path, pp);
    }
    break;
  case GAUGEBC_SCHROEDINGER_1LINK:
    // DELIBERATE FALL THROUGH
  case GAUGEBC_SCHROEDINGER_2LINK:
    {
      const GaugeBCSchrParams& pp = 
	dynamic_cast<const GaugeBCSchrParams&>(p);
      write(xml, path, pp);
    }
    break;
  case GAUGEBC_SIMPLE:
    { 
      const GaugeBCSimpleParams& pp = 
	dynamic_cast<const GaugeBCSimpleParams&>(p);
      
      write(xml, path, pp);
    }
    break;
  default:
    QDPIO::cerr << "Unsupported GaugeBCParams " << endl;
    QDP_abort(1);
  }

  END_CODE();

}


// Toplevel reader
GaugeBCParamsBase* readGaugeBCParams(XMLReader& xml, const string& path)
{
  START_CODE();

  XMLReader top(xml, path);

  enum GaugeBCType bc_type;

  try { 
    read( top, "BCType", bc_type );
  }
  catch( const string& e ) { 
    QDPIO::cerr << "Unable to read FermAct: " << e << endl;
    QDP_abort(1);
  }

  switch(bc_type) { 
  case GAUGEBC_ALL_PERIODIC:
    return new GaugeBCPeriodicParams(top);
    break;
  case GAUGEBC_SCHROEDINGER_1LINK:
    // DELIBERATE FALL THROUGH
  case GAUGEBC_SCHROEDINGER_2LINK:   
    return new GaugeBCSchrParams(top);
    break;
  case GAUGEBC_SIMPLE:
    return new GaugeBCSimpleParams(top);
    break;
  default:
    QDPIO::cerr << "Unsupported GaugeBCParams " << endl;
    QDP_abort(1);
  }

  END_CODE();

  return 0;
}
