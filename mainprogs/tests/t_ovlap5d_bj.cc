// $Id: t_ovlap5d_bj.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

#include "actions/ferm/invert/inv_gmresr_cg_array.h"
#include "actions/ferm/invert/inv_minres_array.h"
using namespace Chroma;

struct App_input_t {
  ChromaProp_t param;
  Cfg_t        cfg;
};

void read(XMLReader& xml, const string& path, App_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // The parameters holds the version number
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
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
    read(xml_in, "/ovlapTest", input);
  }
   catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }



  // Setup the lattice
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_ovlap5d_bj");

  
  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.

  // Create a FermBC
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));
 
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc5(new SimpleFermBC<multi1d<LatticeFermion> >(input.param.boundary));

  QDPIO::cout << "FERM_ACT_ZOLOTAREV_5D" << endl;
  const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(input.param.FermActHandle));
  
  //! N order Zolo approx, with wilson action.
  Zolotarev5DFermActArray   S(fbc5, fbc, zolo5d, xml_out);
			      
  Handle<const ConnectState> connect_state(S.createState(u, zolo5d.StateInfo, xml_out,zolo5d.AuxFermActHandle->getMass()));
  
  // Make me a linop (this callls the initialise function)
  Handle<const LinearOperator< multi1d< LatticeFermion > > > D_op(S.linOp(connect_state));

  Handle<const LinearOperator< multi1d< LatticeFermion > > > DD_op(S.lMdagM(connect_state));

  int N5 = S.size();
  int G5 = Ns*Ns-1;

  LatticeFermion chi4;
  LatticeFermion psi4;
  LatticeFermion tmp;
  int n_count;


  gaussian(chi4);
  chi4 /= sqrt(norm2(chi4));

 
  // We solve M_5d psi = gamma_5 chi
  multi1d<LatticeFermion> psi( N5 );
  multi1d<LatticeFermion> chi( N5 );


  QDP::StopWatch swatch;
  double t;


  for(int i = 0; i < N5; i++) { 
    psi[i] = zero;
    chi[i] = zero;
  }
 
  chi[N5-1] = Gamma(G5)*chi4;             // 4D source in last component
                                         // is gamma5 chi


  swatch.reset();
  swatch.start();
  //InvGMRESR_CG(*DD_op, *D_op, *D_op, chi, psi, 
  //	       input.param.invParam.RsdCG, 
  //	       input.param.invParam.RsdCGPrec,
  //	       input.param.invParam.MaxCG,
  //	       input.param.invParam.MaxCGPrec,
  //	       n_count);

  InvMINRES(*D_op, chi, psi, input.param.invParam.RsdCG, 
	    input.param.invParam.MaxCG, n_count);

  swatch.stop();
  t = swatch.getTimeInSeconds();

  multi1d<LatticeFermion> tmp5_1(N5);
  // Multiply back to check inverse
  (*D_op)(tmp5_1, psi, PLUS);

  
  multi1d<LatticeFermion> tmp5_2( N5 );
  Double dnorm = Double(0);
  Double g5chi_norm = Double(0);
  for(int i = 0; i < N5; i++) { 
    tmp5_2[i] = chi[i] - tmp5_1[i];
    dnorm += norm2(tmp5_2[i]);
    g5chi_norm += norm2(chi[i]);
  }
  
  
  Double r_norm5 = sqrt(dnorm/g5chi_norm);
  QDPIO::cout << "|| chi - D psi ||/ || chi || = " << r_norm5 << endl;
  QDPIO::cout << "time = " << t << " seconds " << endl;
  push(xml_out, "Inv5DCheck");
  write(xml_out, "r_norm", r_norm5);
  pop(xml_out);



  for(int i = 0; i < N5; i++) { 
    psi[i] = zero;
    chi[i] = zero;
    tmp5_1[i] = zero;
  }
 
  chi[N5-1] = Gamma(G5)*chi4;             // 4D source in last component
                                         // is gamma5 chi


  swatch.reset();
  swatch.start();

  // Part of solution process
  (*D_op)(tmp5_1, chi, MINUS);

  InvCG1(*DD_op, tmp5_1, psi, 
	input.param.invParam.RsdCG, 
	input.param.invParam.MaxCG,
	n_count);


  swatch.stop();
  t = swatch.getTimeInSeconds();

  // Multiply back to check inverse
  (*D_op)(tmp5_1, psi, PLUS);

  
  dnorm = Double(0);
  g5chi_norm = Double(0);
  for(int i = 0; i < N5; i++) { 
    tmp5_2[i] = chi[i] - tmp5_1[i];
    dnorm += norm2(tmp5_2[i]);
    g5chi_norm += norm2(chi[i]);
  }
  
  
  r_norm5 = sqrt(dnorm/g5chi_norm);
  QDPIO::cout << "|| chi - D psi ||/ || chi || = " << r_norm5 << endl;
  QDPIO::cout << "time = " << t << " seconds " << endl;
  push(xml_out, "Inv5DCheck");
  write(xml_out, "r_norm", r_norm5);
  pop(xml_out);

 
  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
 
