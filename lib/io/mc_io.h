#ifndef MC_IO_BASE
#define MC_IO_BASE

#include "chromabase.h"
#include "io/param_io.h"

using namespace QDP;
using namespace std;

// A Base MC Param reader...
struct MCParams {
  multi1d<int> nrow;   // Lattice Size
  int n_warm_up;       // No of warm up updates
  int n_updates;       // No of Updates to do in this run
  int final_update_no; // The no of the final update
  int start_update_no; // The no of the start update in this run
  int save_interval;   // No of updates between saving configs
};

struct MCStartUpParams {
  bool  restartP;       // Are we restarting from a saved state

  Cfg_t start_cfg;      // Starting from HOT/COLD/Loaded CFG
  QDP::Seed seed;       // RNG Seed for start

  std::string restart_file; // LIME Filename for restart from state
};

void read(XMLReader& xml, const string& path, MCParams& p);
void read(XMLReader& xml, const string& path, MCStartUpParams& s);
void write(XMLWriter& xml, const string& path, const MCParams& p);
void write(XMLWriter& xml, const string& path, const MCStartUpParams& s);

#endif
