// $Id: t_invborici.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

#include "actions/ferm/invert/inv_borici_w.h"

using namespace Chroma;

struct App_input_t {
  Zolotarev4DFermActParams* zolo4D;
  Zolotarev5DFermActParams* zolo5D;
  multi1d<int> nrow;
  multi1d<int> boundary;
  InvertParam_t invParam;
  Cfg_t        cfg;
};

// Reader for input parameters
void read(XMLReader& xml, const string& path, App_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    read(inputtop, "nrow", input.nrow);
    read(inputtop, "boundary", input.boundary);
    input.zolo4D = dynamic_cast<Zolotarev4DFermActParams*>(read(inputtop, "Zolo4D"));
    input.zolo5D = dynamic_cast<Zolotarev5DFermActParams*>(read(inputtop, "Zolo5D"));
    read(inputtop, "InvertParam", input.invParam);
    read(inputtop, "Cfg", input.cfg);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  App_input_t input;
  XMLReader xml_in(Chroma::getXMLInputFileName());

  try {
    read(xml_in, "/BoriciTest", input);
  }
   catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }


  // Setup the lattice
  Layout::setLattSize(input.nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_borici");


  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Create a 4D FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.boundary));

  // Create a 5D FermBC
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc_a(new SimpleFermBC<multi1d<LatticeFermion> >(input.boundary));

 
  // Make the Actions
  Zolotarev5DFermActArray S5(fbc_a, fbc, *(input.zolo5D), xml_out);
  Zolotarev4DFermAct S4(fbc, *(input.zolo4D), xml_out);

  
  // Make the connect state(s)
  Handle<const ConnectState> state_4(S4.createState(u, input.zolo4D->StateInfo, xml_out, input.zolo4D->AuxFermActHandle->getMass()));

  Handle<const ConnectState> state_5(S5.createState(u, input.zolo5D->StateInfo, xml_out, input.zolo5D->AuxFermActHandle->getMass()));


  // Make the LinOps
  Handle<const LinearOperator<LatticeFermion> > D_4(S4.linOp(state_4));
  Handle<const LinearOperator< multi1d<LatticeFermion> > > D_5(S5.linOp(state_5));
  Handle<const LinearOperator< multi1d<LatticeFermion> > > D_dag_D_5(S5.lMdagM(state_5));

  LatticeFermion b;
  // Solve on chiral point source
  multi1d<int> coord(4);
  coord[0] = 0;
  coord[1] = 0;
  coord[2] = 0;
  coord[3] = 0;
  b=zero;
  srcfil(b, coord, 0, 0);

  LatticeFermion x = zero;

  int n_iters;
  QDP::StopWatch swatch;
  swatch.reset();
  swatch.start();

  InvBorici(*D_4, 
	    *D_5,
	    *D_dag_D_5,
	    b, 
	    x, 
	    input.invParam.RsdCG,
	    input.invParam.RsdCGPrec,
	    input.invParam.MaxCG, 
	    input.invParam.MaxCGPrec, 
	    input.zolo5D->AuxFermActHandle->getMass(),
	    n_iters);

  swatch.stop();

  QDPIO::cout << "InvBorici: " << n_iters << " iterations " << endl;

  LatticeFermion tmp;
  (*D_4)(tmp, x, PLUS);
  tmp -= b;
  Double tmpnorm = sqrt(norm2(tmp)/norm2(b));
  QDPIO::cout << "Final residue: " << tmpnorm << endl;
  QDPIO::cout << "Time: " << swatch.getTimeInSeconds() << " s" << endl;

  x = zero;

  swatch.reset();
  swatch.start();

  S4.qprop(x, state_4, b, REL_GMRESR_SUMR_INVERTER,   
	     input.invParam.RsdCG,
	     input.invParam.RsdCGPrec,
	     input.invParam.MaxCG, 
	     input.invParam.MaxCGPrec, 
	     n_iters);

  swatch.stop();
  QDPIO::cout << "Time: " << swatch.getTimeInSeconds() << " s" << endl;
  
  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
