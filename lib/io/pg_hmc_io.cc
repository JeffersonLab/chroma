#include "chromabase.h"
#include "io/pg_hmc_io.h"

#include <string>

using namespace QDP;
using namespace std;

void read(XMLReader& xml, const string& path, PureGaugeHMCParams& p)
{
  try {
    XMLReader top(xml, path);

    read(top, "IO_version/version", p.io_version.version);
    
    switch( p.io_version.version ) {
    case 1:
      {
	read(top, "MC_Iterations", p.MC_iters);

	read(top, "MC_StartUp", p.MC_startup);

	p.gauge_bc_handle = readGaugeBCParams(top, "GaugeBC");

	p.gauge_act_MC_handle = readGaugeActParams(top, "GaugeAction_MC");

	if( top.count("GaugeAction_MD") == 1 ) { 
	  p.gauge_act_MD_handle = readGaugeActParams(top, "GaugeAction_MD");
	}
	else { 
	  p.gauge_act_MD_handle = p.gauge_act_MC_handle;
	}

	p.MD_params_handle = readMDIntegratorParams(top, "MDIntegrator");
      }
      break;
    default:
      QDPIO::cerr << "Unknown IO-version: " << p.io_version.version << endl;
      QDP_abort(1);
    }
  }
  catch( const string& e ) {
    QDPIO::cerr << "Caught Exception while reading XML : " << e << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const PureGaugeHMCParams& p)
{
  try { 
    push(xml, path);
   
    push(xml,"IO_version");
    write(xml, "version", p.io_version.version);
    pop(xml);

    write(xml, "MC_Iterations", p.MC_iters);
    write(xml, "MC_StartUp", p.MC_startup);
    write(xml, "GaugeBC", *(p.gauge_bc_handle));
    write(xml, "GaugeAction_MC", *(p.gauge_act_MC_handle));
    write(xml, "GaugeAction_MD", *(p.gauge_act_MD_handle));
    write(xml, "MDIntegrator", *(p.MD_params_handle));
  }
  catch( const string& e ) {
    QDPIO::cerr << "Caught exception while writing XML : " << e << endl;
    QDP_abort(1);
  }
}
