#include "util/gauge/gauge_startup.h"

void gaugeStartup(XMLReader& gauge_file_xml,
		  XMLReader& gauge_xml,
		  multi1d<LatticeColorMatrix>& u,
		  Cfg_t& cfg)
{

  switch (cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, cfg.cfg_file);
    break;

  case CFG_TYPE_SZINQIO:
    readGauge(gauge_file_xml, gauge_xml, u, cfg.cfg_file, QDPIO_SERIAL);
    break;

  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, cfg.cfg_file);
    break;
  
  case CFG_TYPE_DISORDERED:
    QDPIO::cout << "Starting up disordered (random/hot) config" << endl;
    for(int dim=0; dim < Nd; dim++) { 
	random(u[dim]);
	reunit(u[dim]);
    }
    break;
  case CFG_TYPE_UNIT:
    QDPIO::cout << "Starting up unit gauge (free) config" << endl;
    for(int dim=0; dim < Nd; dim++) { 
	u[dim] = Real(1);
    }
    break; 
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

}
