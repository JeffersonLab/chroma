// $Id: gauge_startup.cc,v 1.3 2004-04-28 14:28:58 edwards Exp $
/*! \file
 *  \brief Initialize the gauge fields
 */

#include "chromabase.h"
#include "util/gauge/gauge_startup.h"

using namespace QDP;

//! Initialize the gauge fields
/*!
 * \ingroup gauge
 *
 * \param gauge_file_xml  File xml
 * \param gauge_xml       Record xml
 * \param u               Gauge fields
 * \param cfg             Configuration structure
 */
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
  
  case CFG_TYPE_MILC:
    readMILC(gauge_xml, u, cfg.cfg_file);
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
