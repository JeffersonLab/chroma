#ifndef GAUGEACT_UTILS_H
#define GAUGEACT_UTILS_H

#include "chromabase.h"
#include "gaugebc.h"
#include "gaugeact.h"
#include "io/gaugeact_io.h"

GaugeAction* getGaugeActFromParams(Handle< GaugeBC > gbc, 
				   const GaugeActParamsBase& b);

#endif
