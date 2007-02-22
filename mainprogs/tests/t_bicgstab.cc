// $Id: t_bicgstab.cc,v 3.1 2007-02-22 21:11:50 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "actions/ferm/invert/invsumr.h"

using namespace Chroma;
using namespace std;

struct App_input_t {
  ChromaProp_t param;
  Cfg_t        cfg;
};

// Reader for input parameters
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
    read(xml_in, "/BiCGStabTest", input);
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
  push(xml_out,"BiCGStabTest");


  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();
  
  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

  // Create a FermBC
  Handle<FermBC<T,P,Q> >  fbc(new SimpleFermBC<T,P,Q>(input.param.boundary));

  // Initialize fermion action
  //
  FermionAction<LatticeFermion>* S_f_ptr = 0;
  FermionAction< multi1d<LatticeFermion> >* S_f_a_ptr = 0;

  switch (input.param.FermActHandle->getFermActType() ) {
  case FERM_ACT_WILSON:
    {
      const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_WILSON" << endl;
      S_f_ptr = new EvenOddPrecWilsonFermAct(fbc, wils.Mass,
					     wils.anisoParam);
    }
    break;
    
  case FERM_ACT_UNPRECONDITIONED_WILSON:
    {
      const WilsonFermActParams& wils = dynamic_cast<const WilsonFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_UNPRECONDITIONED_WILSON" << endl;
      S_f_ptr = new UnprecWilsonFermAct(fbc, wils.Mass);
    }
    break;
    
  case FERM_ACT_ZOLOTAREV_4D:
    {
      QDPIO::cout << "FERM_ACT_ZOLOTAREV_4D" << endl;
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      
      // Construct Fermact -- now uses constructor from the zolo4d params
      // struct
      S_f_ptr = new Zolotarev4DFermAct(fbc, zolo4d, xml_out);
    }
    break;
  
  case FERM_ACT_ZOLOTAREV_5D:
    {
      QDPIO::cout << "FERM_ACT_ZOLOTAREV_5D" << endl;
      const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(input.param.FermActHandle));
    
      // Construct Fermact -- now uses constructor from the zolo4d params
      // struct
      S_f_a_ptr = new Zolotarev5DFermActArray(fbc, fbc, zolo5d, xml_out);
    }
    break;

  case FERM_ACT_DWF:
    {
      const DWFFermActParams& dwf = dynamic_cast<const DWFFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_DWF" << endl;
      S_f_a_ptr = new EvenOddPrecDWFermActArray(fbc,
						dwf.chiralParam.OverMass, 
						dwf.Mass, 
						dwf.chiralParam.N5);
    }
  break;

  case FERM_ACT_UNPRECONDITIONED_DWF:
    {
      const DWFFermActParams& dwf = dynamic_cast<const DWFFermActParams&>(*(input.param.FermActHandle));
      
      QDPIO::cout << "FERM_ACT_UNPRECONDITONED_DWF" << endl;
      S_f_a_ptr = new UnprecDWFermActArray(fbc,
					   dwf.chiralParam.OverMass, 
					   dwf.Mass, 
					   dwf.chiralParam.N5);
    }
    break;
  default:
    QDPIO::cerr << "Unsupported fermion action" << endl;
    QDP_abort(1);
  }

  QDPIO::cout << "Ferm Act Created " << endl;

  // Create a useable handle on the action
  // The handle now owns the pointer
  Handle< FermionAction<T,P,Q> > S_f(S_f_ptr);
  Handle< FermionAction<T,P,Q> > S_f_a(S_f_a_ptr);
  
  // FIrst we have to set up the state -- this is fermact dependent
  const FermState< *state_ptr;

  switch(input.param.FermActHandle->getFermActType()) {
  case FERM_ACT_WILSON:
    state_ptr = S_f->createState(u);
    break;
  case FERM_ACT_UNPRECONDITIONED_WILSON:
    state_ptr = S_f->createState(u);
    break;
  case FERM_ACT_DWF:
    state_ptr = S_f_a->createState(u);
    break;
  case FERM_ACT_UNPRECONDITIONED_DWF:
    state_ptr = S_f_a->createState(u);
    break;
    
  case FERM_ACT_ZOLOTAREV_4D:
    {    
      const Zolotarev4DFermActParams& zolo4d = dynamic_cast<const Zolotarev4DFermActParams& > (*(input.param.FermActHandle));
      const Zolotarev4DFermAct& S_zolo4 = dynamic_cast<const Zolotarev4DFermAct&>(*S_f);
      
      state_ptr = S_zolo4.createState(u, zolo4d.StateInfo, xml_out,zolo4d.AuxFermActHandle->getMass());
      
    }
    break;
  case FERM_ACT_ZOLOTAREV_5D:
    {
      const Zolotarev5DFermActParams& zolo5d = dynamic_cast<const Zolotarev5DFermActParams& > (*(input.param.FermActHandle));
      const Zolotarev5DFermActArray& S_zolo5 = dynamic_cast<const Zolotarev5DFermActArray&>(*S_f_a);
      
      
      state_ptr = S_zolo5.createState(u, zolo5d.StateInfo, xml_out, zolo5d.AuxFermActHandle->getMass());
    }
    break;
    
  default:
    QDPIO::cerr << "Unsupported fermion action (state creation)" << endl;
    QDP_abort(1);
  }
  
  // Now do the quarkprop here again depending on which of the
  // action pointers is not null
  Handle<const ConnectState> state(state_ptr);  // inserts any BC

  QDPIO::cout << "ConnectState Created" << endl;

  switch(input.param.FermActHandle->getFermActType()) {
  case FERM_ACT_WILSON:
  case FERM_ACT_UNPRECONDITIONED_WILSON:
  case FERM_ACT_ZOLOTAREV_4D:
    {
      // 4D BiCGStab test
      Handle<const LinearOperator<LatticeFermion> > D(S_f->linOp(state));
      
      const Subset& s = D->subset();

      LatticeFermion psi, chi;
      psi[s] = zero;
      chi[s] = zero;

      gaussian(chi[s]);
      Double chi_norm = sqrt(norm2(chi,s));

      chi[s] /= chi_norm;
      psi[s] = zero;
      
      int n_count=0;
      InvBiCGStab( *D, chi, psi, input.param.invParam.RsdCG,
		   input.param.invParam.MaxCG,n_count);

      // Check back solution
      LatticeFermion r;
      (*D)(r, psi, PLUS);
      r[s] -= chi;
      
      QDPIO::cout << " || chi - D psi || / || chi || = " << sqrt(norm2(r,s))/chi_norm << endl;
    }
    break;
  case FERM_ACT_DWF:
  case FERM_ACT_UNPRECONDITIONED_DWF:
  case FERM_ACT_ZOLOTAREV_5D:
  {
      // 4D BiCGStab test
      Handle<const LinearOperator<multi1d< LatticeFermion> > > D(S_f_a->linOp(state));
      
      const Subset& s = (*D).subset();
      int N = (*D).size();


      multi1d<LatticeFermion> psi(N);
      multi1d<LatticeFermion> chi(N);

      for(int n=0; n < N; n++) {
	psi[n][s] = zero;
	chi[n][s] = zero;
      }

      gaussian(chi[N-1][s]);
      Double chi_norm = sqrt(norm2(chi[N-1],s));

      chi[N-1][s] /= chi_norm;
      
      
      int n_count=0;
      InvBiCGStab( *D, chi, psi, input.param.invParam.RsdCG,
		   input.param.invParam.MaxCG,n_count);

      // Check back solution
      multi1d<LatticeFermion> r(N);
      (*D)(r, psi, PLUS);

      for(int n=0; n < N; n++) { 
	r[n][s] -= chi[n];
      }

      QDPIO::cout << " || chi - D psi || / || chi || = " << sqrt(norm2(r,s))/chi_norm << endl;

    }
    break;
  default:
    QDPIO::cerr << "Unsupported fermion action (state creation)" << endl;
    QDP_abort(1);
  }



  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
