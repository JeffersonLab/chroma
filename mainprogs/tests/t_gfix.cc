/// $Id: t_gfix.cc,v 1.1 2004-01-02 21:19:39 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Read the lattice size
  multi1d<int> nrow(Nd);
  QDPIO::cout << "Enter lattice size" << endl;
  QDPIO::cin >> nrow;

  // Setup the layout
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_mesplq.xml");
  push(xml, "t_mesplq");

  push(xml,"lattis");
  Write(xml,Nd);
  Write(xml,Nc);
  Write(xml,nrow);
  pop(xml);

  int type;
  QDPIO::cout << "Enter Gauge field type\n"
	      << "  (1) Free field\n"
	      << "  (2) Random-transformed free field\n"
	      << "  (3) Hot start (call hotst)\n"
	      << "  (4) SZIN configuration" << endl;
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

  push(xml_out,"Lattice_dimensions");
  Write(xml_out, Nc);
  Write(xml_out, Nd);
  Write(xml_out, nrow);
  pop(xml_out);
  push(xml_out,"Boundary_conditions");
  Write(xml_out, boundary);
  pop(xml_out);

  switch (type)
  {
  case 1:
    push(xml_out,"Free_Field");
    Write(xml_out, type);
    pop(xml_out);
    QDPIO::cout << "Fill u with free field" << endl;
    u = 1;
    break;

  case 2:
    push(xml_out,"Free_Field_with_random_gauge_transformation");
    Write(xml_out, type);
    pop(xml_out);
    QDPIO::cout << "Fill u with random gauge transformed free field" << endl;
    u = 1;
    rgauge(u, g);
    break;

  case 3:
    push(xml_out,"Semi-Haar_measure");
    Write(xml_out, type);
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
    Write(xml_out, type);
    pop(xml_out);
    QDPIO::cout << "Read SZIN config from " << cfg_file_in << endl;

    XMLReader gauge_xml;
    readSzin(gauge_xml, u, cfg_file);
  }
  break;

  case 5:
  {
    string cfg_file_in;
    QDPIO::cout << "Enter NERSC input file name" << endl;
    QDPIO::cin >> cfg_file_in;

    push(xml_out,"Configuration");
    Write(xml_out, type);
    pop(xml_out);
    QDPIO::cout << "Read NERSC config from " << cfg_file_in << endl;

    XMLReader gauge_xml;
    readArchiv(gauge_xml, u, cfg_file);
  }
  break;

  default:
    QDP_error_exit("unknown type", type);
  }

  string cfg_file_out;
  QDPIO::cout << "Enter output file name" << endl;
  QDPIO::cin >> cfg_file_out;


  push(xml_out,"Gfix_parameters");
  Write(xml_out, j_decay);
  Write(xml_out, GFAccu);
  Write(xml_out, GFMax);
  Write(xml_out, OrlxDo);
  Write(xml_out, OrPara);
  pop(xml_out);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << " Initial plaqettes and link: " << w_plaq
	      << " " << s_plaq << " " << t_plaq << " " << link << endl;

  push(xml_out,"Initial_gauge_invariant_observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  // Now gauge fix
  gfix(u, j_decay, GFAccu, GFMax, nrl_gf, OrlxDo, OrPara);

  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  QDPIO::cout << " Final plaqettes and link: " << w_plaq
	      << " " << s_plaq << " " << t_plaq << " " << link << endl;

  push(xml_out,"Final_gauge_invariant_observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  push(xml_out,"Relaxation_iterations_in_GFIX");
  Write(xml_out, nrl_gf);
  pop(xml_out);

  // Now write the gauge field in NERSC format
  QDPIO::cout << "Trying to write NERSC Archive  t_nersc.cfg" << endl;
  writeArchiv(u, "t_nersc.cfg");
    
  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
