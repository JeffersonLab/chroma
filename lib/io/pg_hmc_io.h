#ifndef PG_HMC_IO
#define PG_HMC_IO

#include "chromabase.h"
#include "handle.h"
#include "io/param_io.h"
#include "io/mc_io.h"
#include "io/md_io.h"
#include "io/gaugebc_io.h"
#include "io/gaugeact_io.h"
#include <string>

using namespace QDP;
using namespace std;

struct PureGaugeHMCParams { 
  struct IO_version_t io_version;

  struct MCParams MC_iters;

  struct MCStartUpParams MC_startup;

  Handle< GaugeBCParamsBase > gauge_bc_handle;

  Handle< GaugeActParamsBase > gauge_act_MC_handle;
  Handle< GaugeActParamsBase > gauge_act_MD_handle;

  Handle< MDIntegratorParamsBase > MD_params_handle;
};

void read(XMLReader& xml, const string& path, PureGaugeHMCParams& p);
void write(XMLWriter& xml, const string& path, const PureGaugeHMCParams& p);

#endif
