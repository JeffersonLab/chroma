// $Id: gauge_startup.cc,v 1.11 2005-01-13 18:15:54 edwards Exp $
/*! \file
 *  \brief Initialize the gauge fields
 */

#include "chromabase.h"
#include "util/gauge/gauge_startup.h"

#include "qdp_iogauge.h"
#include "io/param_io.h"
#include "io/gauge_io.h"
#include "io/readszin.h"
#include "io/readmilc.h"
#include "io/kyugauge_io.h"
#include "io/readcppacs.h"

#include "util/gauge/hotst.h"


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
  u.resize(Nd);

  switch (cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, cfg.cfg_file);
    break;

  case CFG_TYPE_SZINQIO:
  case CFG_TYPE_SCIDAC:
    readGauge(gauge_file_xml, gauge_xml, u, cfg.cfg_file, QDPIO_SERIAL);
    break;

  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, cfg.cfg_file);
    break;
  
  case CFG_TYPE_MILC:
    readMILC(gauge_xml, u, cfg.cfg_file);
    break;

  case CFG_TYPE_CPPACS :
    readCPPACS(gauge_xml, u, cfg.cfg_file);
    break;

  case CFG_TYPE_KYU:
    readKYU(u, cfg.cfg_file);
    break;

  case CFG_TYPE_DISORDERED:
    QDPIO::cout << "Starting up disordered (random/hot) config" << endl;
    HotSt(u);
    break;

  case CFG_TYPE_UNIT:
    QDPIO::cout << "Starting up unit gauge (free) config" << endl;
    u = 1;
    break; 

  default:
    QDPIO::cerr << "gaugeStartup: Configuration type is unsupported." << endl;
    QDP_abort(1);
  }

}
