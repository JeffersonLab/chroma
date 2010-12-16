// $Id: gauge_startup.cc,v 3.5 2009-05-13 13:12:19 edwards Exp $
/*! \file
 *  \brief Initialize the gauge fields
 */

#include "chromabase.h"
#include "util/gauge/gauge_startup.h"
#include "util/gauge/reunit.h"

#include "qdp_iogauge.h"
#include "io/param_io.h"
#include "io/gauge_io.h"
#include "io/readszin.h"
#include "io/readwupp.h"
#include "io/readmilc.h"
#include "io/kyugauge_io.h"
#include "io/readcppacs.h"

#include "util/gauge/hotst.h"
#include "util/gauge/weak_field.h"

#include "actions/gauge/gaugebcs/schr_nonpert_gaugebc.h"

namespace Chroma 
{

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
    START_CODE();

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

    case CFG_TYPE_WUPP :
      readWupp(gauge_xml, u, cfg.cfg_file);
      break;


    case CFG_TYPE_KYU:
    {
      readKYU(u, cfg.cfg_file);

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      push(record_xml, "kentucky");
      pop(record_xml);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
    break;

    case CFG_TYPE_DISORDERED:
    {
      QDPIO::cout << "Starting up disordered (random/hot) config" << endl;
      HotSt(u);

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      push(record_xml, "disordered");
      pop(record_xml);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
    break;

    case CFG_TYPE_UNIT:
    {
      QDPIO::cout << "Starting up unit gauge (free) config" << endl;
      u = 1;

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      push(record_xml, "unit");
      pop(record_xml);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
    break; 

    case CFG_TYPE_WEAK_FIELD:
    {
      QDPIO::cout << "Starting up a weak field config" << endl;
      weakField(u);

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      push(record_xml, "weak_field");
      pop(record_xml);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
    break; 

    case CFG_TYPE_CLASSICAL_SF:
    {
      QDPIO::cout << "Starting up a classical Schroedinger functional config" << endl;
      SchrGaugeBCParams params;
      params.loop_extent = 1;
      params.SchrPhiMult = 1;
      params.decay_dir   = Nd-1;
      SchrNonPertGaugeBC gaugebc(params);
      
      u = gaugebc.SFBndFld();

      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      write(record_xml, "SF_classical", params);

      gauge_file_xml.open(file_xml);
      gauge_xml.open(record_xml);
    }
    break; 

    default:
      QDPIO::cerr << __func__ << ": Configuration type is unsupported." << endl;
      QDP_abort(1);
    }

    // Reunitarize gauge field (eg for Double prec?)
    // This should get rid of those pesky big Delta H's 
    // going from single to double and they can't possibly hurt
    // in either single or double
    for(int mu=0; mu < Nd; mu++){ 
      reunit(u[mu]);
    }


    END_CODE();
  }

}  // end namespace Chroma
