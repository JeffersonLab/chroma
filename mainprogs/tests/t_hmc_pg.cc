#include "chroma.h"

#include "io/pg_hmc_io.h"

#include <iostream>

using namespace QDP;
using namespace std;


int main(int argc, char *argv[])
{
  // Initialise QDP
  QDP_initialize(&argc, &argv);

  XMLReader xml_in("DATA");
  struct PureGaugeHMCParams input;
  read(xml_in, "PureGaugeHMC/params", input);

  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "PureGaugeHMC");
  write(xml_out, "params", input);
  pop(xml_out);

  QDP_finalize();
  exit(0);
}

