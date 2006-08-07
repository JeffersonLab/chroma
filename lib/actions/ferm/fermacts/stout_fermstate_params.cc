// -*- C++ -*-
// $Id: stout_fermstate_params.cc,v 1.4 2006-08-07 18:13:28 edwards Exp $

#include "actions/ferm/fermacts/stout_fermstate_params.h"


namespace Chroma 
{

  StoutFermStateParams::StoutFermStateParams() 
  {
    rho.resize(Nd, Nd);
    smear_in_this_dirP.resize(Nd);

    n_smear = 0;
    rho = sm_fact = zero;
    orthog_dir = Nd;
    smear_in_this_dirP = true;
  }


  StoutFermStateParams::StoutFermStateParams(XMLReader& in, const std::string& path) 
  {
    rho.resize(Nd, Nd);
    smear_in_this_dirP.resize(Nd);

    sm_fact = zero;
    orthog_dir = Nd;

    try 
    { 
      XMLReader paramtop(in, path);

      read(paramtop, "rho", sm_fact);
      read(paramtop, "n_smear", n_smear);
      if( paramtop.count("orthog_dir") == 1 ) { 
	read(paramtop, "orthog_dir", orthog_dir);
      }
      else { 
	// default value for orthog dir is 3 -- smear in space only
	QDPIO::cout << "Using Default value: orthog_dir = 3, spatial only smearing" << endl;
      }
    }
    catch(const std::string& e) 
    { 
      QDPIO::cout << "Failed to read stout action XML:" << e << endl;
    }
    

    // For each (mu,nu) set sm_fact_array(mu,nu)=sm_fact
    // (Isotropy). Since mu != nu ever, we set those
    // to zero for safety
    for(int mu=0; mu < Nd; mu++) 
    { 
      for(int nu=0; nu < Nd; nu++) 
      { 
	if( mu != nu ) {
	  rho[mu][nu] = sm_fact;
	}
	else {
	  // Set the rho to 0 if mu==nu
	  rho[mu][nu] = 0;
	}
      }
      
      // Mask out the orthog dir
      if( mu == orthog_dir ) {  // Direction is same as orthog dir
	smear_in_this_dirP[mu]=false;
      }
      else 
      {                    // Direction orthogonal to orthog dir
	smear_in_this_dirP[mu]=true;
      }
    }

  }

  void read(XMLReader& xml, const std::string& path, StoutFermStateParams& p)
  {
    StoutFermStateParams tmp_p(xml, path);
    p = tmp_p;
  }

  void write(XMLWriter& xml, const std::string& path, const StoutFermStateParams& p) 
  {
    push(xml, path);
    write(xml, "rho", p.sm_fact);
    write(xml, "n_smear", p.n_smear);
    write(xml, "orthog_dir", p.orthog_dir);
    pop(xml);
  }


};    

