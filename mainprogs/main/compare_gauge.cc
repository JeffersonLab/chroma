/*! \file
 *  \brief Compare two gauge configurations (data only)
 */

#include "chroma.h"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace Chroma;

namespace
{
  void usage(const char* prog)
  {
    std::cerr << "Usage: " << prog
              << " --nrow n0 n1 n2 n3 [--tol tol] [--rtol rtol] <cfg1.lime> <cfg2.lime>\n";
  }
}

int main(int argc, char* argv[])
{
  Chroma::initialize(&argc, &argv);

  if (argc < 3) {
    usage(argv[0]);
    Chroma::finalize();
    return 1;
  }

  multi1d<int> nrow(Nd);
  bool have_nrow = false;
  double tol = 0.0;
  double rtol = 0.0;
  std::string cfg1;
  std::string cfg2;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--nrow") {
      if (i + Nd >= argc) {
        usage(argv[0]);
        Chroma::finalize();
        return 1;
      }
      for (int mu = 0; mu < Nd; ++mu) {
        nrow[mu] = std::atoi(argv[++i]);
      }
      have_nrow = true;
    } else if (arg == "--tol") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        Chroma::finalize();
        return 1;
      }
      tol = std::atof(argv[++i]);
    } else if (arg == "--rtol") {
      if (i + 1 >= argc) {
        usage(argv[0]);
        Chroma::finalize();
        return 1;
      }
      rtol = std::atof(argv[++i]);
    } else if (cfg1.empty()) {
      cfg1 = arg;
    } else if (cfg2.empty()) {
      cfg2 = arg;
    } else {
      usage(argv[0]);
      Chroma::finalize();
      return 1;
    }
  }

  if (!have_nrow || cfg1.empty() || cfg2.empty()) {
    usage(argv[0]);
    Chroma::finalize();
    return 1;
  }

  Layout::setLattSize(nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u1(Nd);
  multi1d<LatticeColorMatrix> u2(Nd);

  try {
    XMLReader file_xml_1;
    XMLReader gauge_xml_1;
    readGauge(file_xml_1, gauge_xml_1, u1, cfg1, QDPIO_SERIAL);

    XMLReader file_xml_2;
    XMLReader gauge_xml_2;
    readGauge(file_xml_2, gauge_xml_2, u2, cfg2, QDPIO_SERIAL);
  } catch (const std::string& e) {
    QDPIO::cerr << "compare_gauge: error reading gauge files: " << e << std::endl;
    Chroma::finalize();
    return 1;
  }

  double total_n2 = 0.0;
  double ref_n2 = 0.0;
  for (int mu = 0; mu < Nd; ++mu) {
    LatticeColorMatrix diff = u1[mu] - u2[mu];
    Double n2 = norm2(diff);
    total_n2 += toDouble(n2);
    Double ref = norm2(u1[mu]);
    ref_n2 += toDouble(ref);
  }

  double rms = std::sqrt(total_n2);
  double ref_rms = std::sqrt(ref_n2);
  double rel = (ref_rms > 0.0) ? (rms / ref_rms) : 0.0;
  QDPIO::cout << "Gauge diff: sqrt(sum_mu norm2) = " << rms << std::endl;
  QDPIO::cout << "Gauge diff (relative) = " << rel << std::endl;

  bool ok = (rms <= tol) && (rel <= rtol);
  QDPIO::cout << "Compare status: " << (ok ? "MATCH" : "DIFF")
              << " (tol=" << tol << ", rtol=" << rtol << ")" << std::endl;

  Chroma::finalize();
  return ok ? 0 : 2;
}
