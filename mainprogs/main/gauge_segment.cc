/*! \file
 *  \brief Extract a time segment from a QIO gauge configuration
 */

#include "chroma.h"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

using namespace Chroma;

namespace
{
  struct Params
  {
    multi1d<int> nrow;
    bool have_nrow;
    int t0;
    int t1;
    int t_dir;
    bool have_t0;
    bool have_t1;
    QDP_serialparallel_t read_serpar;
    QDP_serialparallel_t write_serpar;
    QDP_volfmt_t volfmt;
    std::string input;
    std::string output;

    Params()
      : nrow(Nd),
        have_nrow(false),
        t0(0),
        t1(-1),
        t_dir(Nd - 1),
        have_t0(false),
        have_t1(false),
        read_serpar(QDPIO_SERIAL),
        write_serpar(QDPIO_SERIAL),
        volfmt(QDPIO_SINGLEFILE)
    {
    }
  };

  void usage(const char* prog)
  {
    std::cerr
      << "Usage: " << prog << " --nrow n0 n1 n2 n3 --t0 t0 --t1 t1 [options] <input.lime> <output.lime>\n"
      << "\n"
      << "Extracts the inclusive time range [t0,t1] into a normal QIO gauge field\n"
      << "with output lattice size n0 n1 n2 (t1-t0+1). Gauge links are copied\n"
      << "from base sites in the selected time range.\n"
      << "\n"
      << "Options:\n"
      << "  --tdir dir          Time direction (default: Nd-1)\n"
      << "  --parallel-read     Use QDPIO_PARALLEL for input\n"
      << "  --parallel-write    Use QDPIO_PARALLEL for output\n"
      << "  --parallel-io       Use QDPIO_PARALLEL for input and output\n"
      << "  --multi-file        Write QDPIO_MULTIFILE output (default: single file)\n";
  }

  bool parse_int(const char* text, int& value)
  {
    char* end = 0;
    long parsed = std::strtol(text, &end, 10);
    if (end == text || *end != '\0')
      return false;
    if (parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max())
      return false;

    value = static_cast<int>(parsed);
    return true;
  }

  bool parse_args(int argc, char* argv[], Params& params)
  {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];

      if (arg == "--help" || arg == "-h") {
        return false;
      } else if (arg == "--nrow") {
        if (i + Nd >= argc)
          return false;
        for (int mu = 0; mu < Nd; ++mu) {
          if (!parse_int(argv[++i], params.nrow[mu]))
            return false;
        }
        params.have_nrow = true;
      } else if (arg == "--t0") {
        if (i + 1 >= argc || !parse_int(argv[++i], params.t0))
          return false;
        params.have_t0 = true;
      } else if (arg == "--t1") {
        if (i + 1 >= argc || !parse_int(argv[++i], params.t1))
          return false;
        params.have_t1 = true;
      } else if (arg == "--tdir") {
        if (i + 1 >= argc || !parse_int(argv[++i], params.t_dir))
          return false;
      } else if (arg == "--parallel-read") {
        params.read_serpar = QDPIO_PARALLEL;
      } else if (arg == "--parallel-write") {
        params.write_serpar = QDPIO_PARALLEL;
      } else if (arg == "--parallel-io") {
        params.read_serpar = QDPIO_PARALLEL;
        params.write_serpar = QDPIO_PARALLEL;
      } else if (arg == "--multi-file") {
        params.volfmt = QDPIO_MULTIFILE;
      } else if (params.input.empty()) {
        params.input = arg;
      } else if (params.output.empty()) {
        params.output = arg;
      } else {
        return false;
      }
    }

    return params.have_nrow && params.have_t0 && params.have_t1 &&
           !params.input.empty() && !params.output.empty();
  }

  bool validate_params(const Params& params)
  {
    if (params.t_dir < 0 || params.t_dir >= Nd) {
      QDPIO::cerr << "gauge_segment: tdir must be in [0," << (Nd - 1) << "]\n";
      return false;
    }

    for (int mu = 0; mu < Nd; ++mu) {
      if (params.nrow[mu] <= 0) {
        QDPIO::cerr << "gauge_segment: all lattice dimensions must be positive\n";
        return false;
      }
    }

    if (params.t0 < 0 || params.t1 < params.t0 ||
        params.t1 >= params.nrow[params.t_dir]) {
      QDPIO::cerr << "gauge_segment: require 0 <= t0 <= t1 < nrow[tdir]\n";
      return false;
    }

    return true;
  }

  int checked_volume(const multi1d<int>& nrow)
  {
    long long volume = 1;
    for (int mu = 0; mu < Nd; ++mu) {
      volume *= static_cast<long long>(nrow[mu]);
      if (volume > std::numeric_limits<int>::max()) {
        QDPIO::cerr << "gauge_segment: output volume is too large for this utility\n";
        QDP_abort(1);
      }
    }

    return static_cast<int>(volume);
  }

  multi1d<int> linear_to_coord(int linear, const multi1d<int>& nrow)
  {
    multi1d<int> coord(Nd);
    for (int mu = 0; mu < Nd; ++mu) {
      coord[mu] = linear % nrow[mu];
      linear /= nrow[mu];
    }

    return coord;
  }
}

int main(int argc, char* argv[])
{
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      return 0;
    }
  }

  Chroma::initialize(&argc, &argv);

  Params params;
  if (!parse_args(argc, argv, params)) {
    usage(argv[0]);
    Chroma::finalize();
    return 1;
  }

  if (!validate_params(params)) {
    Chroma::finalize();
    return 1;
  }

  const int segment_extent = params.t1 - params.t0 + 1;
  multi1d<int> out_nrow(params.nrow);
  out_nrow[params.t_dir] = segment_extent;
  const int out_volume = checked_volume(out_nrow);

  XMLReader input_file_xml;
  XMLReader input_record_xml;
  multi1d< multi1d<ColorMatrix> > segment(Nd);
  for (int mu = 0; mu < Nd; ++mu)
    segment[mu].resize(out_volume);

  Layout::setLattSize(params.nrow);
  Layout::create();

  {
    multi1d<LatticeColorMatrix> u(Nd);

    try {
      readGauge(input_file_xml, input_record_xml, u, params.input, params.read_serpar);
    } catch (const std::string& e) {
      QDPIO::cerr << "gauge_segment: error reading input gauge field: " << e << std::endl;
      Chroma::finalize();
      return 1;
    }

    QDPIO::cout << "Extracting time slices " << params.t0 << " through "
                << params.t1 << " from direction " << params.t_dir << std::endl;

    for (int site = 0; site < out_volume; ++site) {
      multi1d<int> out_coord = linear_to_coord(site, out_nrow);
      multi1d<int> in_coord(out_coord);
      in_coord[params.t_dir] += params.t0;

      for (int mu = 0; mu < Nd; ++mu)
        segment[mu][site] = peekSite(u[mu], in_coord);
    }
  }

  Layout::destroy();

  Layout::setLattSize(out_nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u_segment(Nd);
  for (int mu = 0; mu < Nd; ++mu)
    u_segment[mu] = zero;

  for (int site = 0; site < out_volume; ++site) {
    multi1d<int> out_coord = linear_to_coord(site, out_nrow);

    for (int mu = 0; mu < Nd; ++mu)
      pokeSite(u_segment[mu], segment[mu][site], out_coord);
  }

  XMLBufferWriter output_file_xml;
  push(output_file_xml, "gauge_segment");
  write(output_file_xml, "id", int(0));
  write(output_file_xml, "input_file", params.input);
  write(output_file_xml, "output_file", params.output);
  write(output_file_xml, "input_nrow", params.nrow);
  write(output_file_xml, "output_nrow", out_nrow);
  write(output_file_xml, "t_dir", params.t_dir);
  write(output_file_xml, "t0", params.t0);
  write(output_file_xml, "t1", params.t1);
  write(output_file_xml, "input_file_xml", input_file_xml);
  pop(output_file_xml);

  XMLBufferWriter output_record_xml;
  push(output_record_xml, "gauge_segment_record");
  write(output_record_xml, "input_record_xml", input_record_xml);
  write(output_record_xml, "segment_extent", segment_extent);
  write(output_record_xml, "copied_base_site_t0", params.t0);
  write(output_record_xml, "copied_base_site_t1", params.t1);
  pop(output_record_xml);

  writeGauge(output_file_xml, output_record_xml, u_segment, params.output,
             params.volfmt, params.write_serpar);

  QDPIO::cout << "Wrote segmented gauge field to " << params.output << std::endl;

  Chroma::finalize();
  return 0;
}
