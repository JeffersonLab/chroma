#ifndef GAUGE_STARTUP_H
#define GAUGE_STARTUP_H

#include "chroma.h"



void gaugeStartup(XMLReader& gauge_file_xml,
		  XMLReader& gauge_xml,
		  multi1d<LatticeColorMatrix>& u,
		  Cfg_t& cfg);


#endif
