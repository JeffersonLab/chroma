#include "chromabase.h"
#include "io/mc_io.h"

using namespace QDP;
using namespace std;

void read(XMLReader& xml, const string& path, MCParams& p)
{
  XMLReader top(xml, path);

  try { 
    read(top, "nrow", p.nrow);
    read(top, "NWarmUp", p.n_warm_up);
    read(top, "NUpdates", p.n_updates);
    read(top, "FinalUpdateNo", p.final_update_no);
    read(top, "StartUpdateNo", p.start_update_no);
    read(top, "SaveInterval", p.save_interval);
  }
  catch(const string& e) {
    QDPIO::cerr << "Caught exception reading XML: " << e << endl;
    QDP_abort(1);
  }
}

void write(XMLWriter& xml, const string& path, const MCParams& p)
{
  try { 
    push(xml, path);
    write(xml, "nrow", p.nrow);
    write(xml, "NWarmUp", p.n_warm_up);
    write(xml, "NUpdates", p.n_updates);
    write(xml, "FinalUpdateNo", p.final_update_no);
    write(xml, "StartUpdateNo", p.start_update_no);
    write(xml, "SaveInterval", p.save_interval);
    pop(xml);
  }
  catch(const string& e) {
    QDPIO::cerr << "Caught exception writing XML" << e << endl;
    QDP_abort(1);
  }
}

void read(XMLReader& xml, const string& path, MCStartUpParams& s)
{
  XMLReader top(xml, path);
  try { 
    read(top, "RestartP", s.restartP);
  }
  catch(const string& e) { 
    QDPIO::cerr << "Caught exception reading XML" << e << endl;
    QDP_abort(1);
  }

  if ( s.restartP ) { 
    // We are restarting from a saved state -- read its LIME filename
    try { 
      read(top, "RestartFilename", s.restart_file);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught exception reading XML" << e << endl;
      QDP_abort(1);
    }
  }
  else { 
    try {
      read(top, "Cfg", s.start_cfg);
      read(top, "Seed", s.seed);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught exception reading XML" << e << endl;
      QDP_abort(1);
    }
  }
}

void write(XMLWriter& xml, const string& path, const MCStartUpParams& s)
{
  push(xml, path);
  write(xml, "RestartP", s.restartP);

  if( s.restartP ) { 
    write(xml, "RestartFilename", s.restart_file);
  }
  else { 
    write(xml, "Cfg", s.start_cfg);
    write(xml, "Seed", s.seed);
  }
  pop(xml);
}

