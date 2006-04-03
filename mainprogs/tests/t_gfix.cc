/// $Id: t_gfix.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Read the lattice size
  multi1d<int> nrow(Nd);
  QDPIO::cout << "Enter lattice size" << endl;
  QDPIO::cin >> nrow;

  // Setup the layout
  Layout::setLattSize(nrow);
  Layout::create();

  int type;
  QDPIO::cout << "Enter Gauge field type\n"
	      << "  (1) Free field\n"
	      << "  (2) Random-transformed free field\n"
	      << "  (3) Hot start (call hotst)\n"
	      << "  (4) SZIN configuration\n"
	      << "  (5) NERSC configuration" << endl;
  QDPIO::cin >> type;

  int j_decay;
  QDPIO::cout << "Enter the direction of decay" << endl;
  QDPIO::cin >> j_decay;

  Real GFAccu;
  QDPIO::cout << "Enter the gauge fixing accuracy" << endl;
  QDPIO::cin >> GFAccu;

  int GFMax;
  QDPIO::cout << "Enter the maximum number of gauge fixing sweeps" << endl;
  QDPIO::cin >> GFMax;

  bool OrlxDo;
  QDPIO::cout << "Want over-relaxation? (yes=YES)" << endl;
  QDPIO::cin >> OrlxDo;

  Real OrPara = 1;
  if (OrlxDo)
  {
    QDPIO::cout << "Enter the over-relaxtion parameter" << endl;
    QDPIO::cin >> OrPara;
  }

  QDPIO::cout << "Now, I am running..." << endl;

  XMLFileWriter xml_out("t_gfix.xml");
  push(xml_out, "t_gfix");

  push(xml_out,"Lattice_dimensions");
  write(xml_out,"Nc", Nc);
  write(xml_out,"Nd", Nd);
  write(xml_out, "nrow", nrow);
  pop(xml_out);

  multi1d<LatticeColorMatrix> u(Nd);    // Gauge field
  LatticeColorMatrix g;

  switch (type)
  {
  case 1:
    push(xml_out,"Free_Field");
    write(xml_out, "type", type);
    pop(xml_out);
    QDPIO::cout << "Fill u with free field" << endl;
    u = 1;
    break;

  case 2:
    push(xml_out,"Free_Field_with_random_gauge_transformation");
    write(xml_out, "type", type);
    pop(xml_out);
    QDPIO::cout << "Fill u with random gauge transformed free field" << endl;
    u = 1;
    rgauge(u, g);
    break;

  case 3:
    push(xml_out,"Semi-Haar_measure");
    write(xml_out, "type", type);
    pop(xml_out);
    QDPIO::cout << "Fill u with semi-Haar" << endl;
    HotSt(u);
    break;

  case 4:
  {
    string cfg_file_in;
    QDPIO::cout << "Enter SZIN input file name" << endl;
    QDPIO::cin >> cfg_file_in;

    push(xml_out,"Configuration");
    write(xml_out, "type", type);
    pop(xml_out);
    QDPIO::cout << "Read SZIN config from " << cfg_file_in << endl;

    XMLReader gauge_xml;
    readSzin(gauge_xml, u, cfg_file_in);
  }
  break;

  case 5:
  {
    string cfg_file_in;
    QDPIO::cout << "Enter NERSC input file name" << endl;
    QDPIO::cin >> cfg_file_in;

    push(xml_out,"Configuration");
    write(xml_out, "type", type);
    pop(xml_out);
    QDPIO::cout << "Read NERSC config from " << cfg_file_in << endl;

    XMLReader gauge_xml;
    readArchiv(gauge_xml, u, cfg_file_in);
  }
  break;

  default:
    QDP_error_exit("unknown type", type);
  }

  string cfg_file_out;
  QDPIO::cout << "Enter NERSC output file name" << endl;
  QDPIO::cin >> cfg_file_out;


  push(xml_out,"Gfix_parameters");
  write(xml_out, "j_decay", j_decay);
  write(xml_out, "GFAccu", GFAccu);
  write(xml_out, "GFMax", GFMax);
  write(xml_out, "OrlxDo", OrlxDo);
  write(xml_out, "OrPara", OrPara);
  pop(xml_out);

  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << " Initial plaqettes and link: " << w_plaq
	      << " " << s_plaq << " " << t_plaq << " " << link << endl;

  push(xml_out,"Initial_gauge_invariant_observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  // Now gauge fix
  int nrl_gf;
  coulGauge(u, nrl_gf, j_decay, GFAccu, GFMax, OrlxDo, OrPara);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << " Final plaqettes and link: " << w_plaq
	      << " " << s_plaq << " " << t_plaq << " " << link << endl;

  push(xml_out,"Final_gauge_invariant_observables");
  write(xml_out, "w_plaq",w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  push(xml_out,"Relaxation_iterations_in_GFIX");
  write(xml_out, "nrl_gf", nrl_gf);
  pop(xml_out);

  // Now write the gauge field in NERSC format
  QDPIO::cout << "Trying to write NERSC Archive to file  " << cfg_file_out << endl;
  writeArchiv(u, cfg_file_out);
    
  pop(xml_out);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
