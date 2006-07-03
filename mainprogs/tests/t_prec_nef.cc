// $Id: t_prec_nef.cc,v 3.1 2006-07-03 15:26:11 edwards Exp $

#include "chroma.h"

using namespace Chroma;

struct App_input_t {
  GroupXML_t      invParam;   // Inverter parameters
  Cfg_t           cfg;
  multi1d<int>    nrow;
};

// Reader for input parameters
void read(XMLReader& xml, const string& path, App_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);
    read(inputtop, "nrow", input.nrow);
    input.invParam = readXMLGroup(paramtop, "InvertParam", "invType");
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
    read(xml_in, "/NEFTest", input);
  }
   catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }


  // Setup the lattice
  Layout::setLattSize(input.nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader gauge_file_xml, gauge_xml;
    gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);
  }

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"NEFTest");


  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Create actions
  string zpath = "/NEFTest/ZNEF";
  EvenOddPrecZoloNEFFermActArray  S_znef(WilsonTypeFermBCArrayEnv::reader(xml_in, zpath), 
					 EvenOddPrecZoloNEFFermActArrayParams(xml_in, zpath));

  string npath = "/NEFTest/NEF";
  EvenOddPrecNEFFermActArray  S_nef(WilsonTypeFermBCArrayEnv::reader(xml_in, npath), 
				    EvenOddPrecNEFFermActArrayParams(xml_in, npath));

  Handle<const ConnectState> state(S_znef.createState(u));

  int N5 = S_znef.size();
  QDPIO::cout << "Znef size = " << S_znef.size() << endl;
  QDPIO::cout << "Nef  size = " << S_nef.size() << endl;

  {
    // Make the znef linOp
    Handle< const LinearOperator< multi1d<LatticeFermion> > > M_z(S_znef.linOp(state));
  
    // Make the nef linOp
    Handle< const LinearOperator< multi1d<LatticeFermion> > > M_n(S_nef.linOp(state));

    multi1d<LatticeFermion> psi5a(N5);
    multi1d<LatticeFermion> chi5(N5);
    multi1d<LatticeFermion> psi5b(N5);
    multi1d<LatticeFermion> tmp5(N5);

    for(int n=0; n<N5; n++) 
    {
      gaussian(psi5a[n]);
      gaussian(psi5b[n]);
      gaussian(chi5[n]);
    }

    (*M_z)(psi5a, chi5, PLUS);
    (*M_n)(psi5b, chi5, PLUS);
    for(int n=0; n < N5; ++n)
    {
      tmp5[n][M_z->subset()] = psi5a[n];
      tmp5[n][M_z->subset()] -= psi5b[n];
    }
    QDPIO::cout << "PLUS: norm2(psi5a-psi5b)=" << norm2(tmp5,M_z->subset()) << endl;

    (*M_z)(psi5a, chi5, MINUS);
    (*M_n)(psi5b, chi5, MINUS);
    for(int n=0; n < N5; ++n)
    {
      tmp5[n][M_z->subset()] = psi5a[n];
      tmp5[n][M_z->subset()] -= psi5b[n];
    }
    QDPIO::cout << "MINUS: norm2(psi5a-psi5b)=" << norm2(tmp5,M_z->subset()) << endl;


    // Make the znef unprec linOp
    Handle< const LinearOperator< multi1d<LatticeFermion> > > 
      U_z(S_znef.unprecLinOp(S_znef.createState(u),Real(1)));
  
    // Make the nef unprec linOp
    Handle< const LinearOperator< multi1d<LatticeFermion> > > 
      U_n(S_nef.unprecLinOp(S_nef.createState(u),Real(1)));

    (*U_z)(psi5a, chi5, PLUS);
    (*U_n)(psi5b, chi5, PLUS);
    for(int n=0; n < N5; ++n)
    {
      tmp5[n][U_z->subset()] = psi5a[n];
      tmp5[n][U_z->subset()] -= psi5b[n];
    }
    QDPIO::cout << "PLUS: norm2(psi5a-psi5b)=" << norm2(tmp5,U_z->subset()) << endl;

    (*U_z)(psi5a, chi5, MINUS);
    (*U_n)(psi5b, chi5, MINUS);
    for(int n=0; n < N5; ++n)
    {
      tmp5[n][U_z->subset()] = psi5a[n];
      tmp5[n][U_z->subset()] -= psi5b[n];
    }
    QDPIO::cout << "MINUS: norm2(psi5a-psi5b)=" << norm2(tmp5,U_z->subset()) << endl;
  }

  {
    Handle< const SystemSolver< multi1d<LatticeFermion> > > ZQ(S_znef.qpropT(S_znef.createState(u), input.invParam));
    Handle< const SystemSolver< multi1d<LatticeFermion> > > NQ(S_nef.qpropT(S_nef.createState(u), input.invParam));

    multi1d<LatticeFermion> psi5a(N5);
    multi1d<LatticeFermion> chi5(N5);
    multi1d<LatticeFermion> psi5b(N5);
    multi1d<LatticeFermion> tmp5(N5);

    for(int n=0; n<N5; n++) 
    {
      gaussian(psi5a[n]);
      gaussian(psi5b[n]);
      gaussian(chi5[n]);
    }
   
    int n_counta = (*ZQ)(psi5a, chi5);
    int n_countb = (*NQ)(psi5b, chi5);

    for(int n=0; n<N5; n++)
      tmp5[n] = psi5b[n] - psi5a[n];

    QDPIO::cout << "norm(qpropT_diff)=" << Real(norm2(tmp5)) 
		<< "  norm2(psi5a)=" << Real(norm2(psi5a)) 
		<< "  norm2(psi5b)=" << Real(norm2(psi5b)) 
		<< endl;
    for(int n=0; n < N5; ++n)
    {
      QDPIO::cout << "QpropT:" 
		  << " norm2(tmp5[" << n << "])= " << Real(norm2(tmp5[n])) 
		  << " norm2(psi5a[])= " << Real(norm2(psi5a[n])) 
		  << " norm2(psi5b[])= " << Real(norm2(psi5b[n])) 
		  << endl;
    }

#if 1
    push(xml_out,"QpropT");
    write(xml_out,"psi5a",psi5a);
    write(xml_out,"psi5b",psi5b);
    write(xml_out,"tmp5",tmp5);
    pop(xml_out);
#endif
  }

  {
    Handle< const SystemSolver<LatticeFermion> > ZQ(S_znef.qprop(S_znef.createState(u), input.invParam));
    Handle< const SystemSolver<LatticeFermion> > NQ(S_nef.qprop(S_nef.createState(u), input.invParam));

    LatticeFermion psia, psib, chi;
    gaussian(chi);
    gaussian(psia);
    gaussian(psib);

    int n_counta = (*ZQ)(psia, chi);
    int n_countb = (*NQ)(psib, chi);

    QDPIO::cout << "norm(qprop_diff)=" << Real(norm2(psib - psia))
		<< "  norm2(psia)=" << Real(norm2(psia)) 
		<< "  norm2(psib)=" << Real(norm2(psib)) 
		<< endl;
  }

  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
