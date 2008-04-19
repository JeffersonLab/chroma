// $Id: stout_fermstate_params.cc,v 1.4 2008-04-19 03:12:49 edwards Exp $
/*! \file
 *  \brief Stout fermstate params
 */

#include "actions/ferm/fermstates/stout_fermstate_params.h"


namespace Chroma 
{

  StoutFermStateParams::StoutFermStateParams() 
  {
    rho.resize(Nd, Nd);
    smear_in_this_dirP.resize(Nd);

    n_smear = 0;
    rho = sm_fact = zero;
    smear_in_this_dirP = true;
  }


  StoutFermStateParams::StoutFermStateParams(XMLReader& in, const std::string& path) 
  {
    try
    {
      XMLReader paramtop(in, path);

      rho.resize(Nd, Nd);

      int version = 1;
      if( paramtop.count("version") == 1 )
	read(paramtop, "version", version);

      read(paramtop, "rho", sm_fact);
      read(paramtop, "n_smear", n_smear);

      switch (version) 
      {
      case 1:
      {
	int orthog_dir = Nd;

	if( paramtop.count("orthog_dir") == 1 ) { 
	  read(paramtop, "orthog_dir", orthog_dir);
	}
	else { 
	  // default value for orthog dir is 3 -- smear in space only
	  QDPIO::cout << "Using Default value: orthog_dir = 3, spatial only smearing" << endl;
	}

	if (paramtop.count("smear_in_this_dirP") > 0)
	{
	  QDPIO::cerr << __func__ << ": found a smear_in_this_dirP in version 1. You need version 2 or higher" << endl;
	  QDP_abort(1);
	}

	smear_in_this_dirP.resize(Nd);
	smear_in_this_dirP = true;

	if (orthog_dir >= 0 && orthog_dir < Nd)
	{
	  smear_in_this_dirP[orthog_dir] = false;
	}
      }
      break;

      case 2:
      {
	read(paramtop, "smear_in_this_dirP", smear_in_this_dirP);
      }
      break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }
    }
    catch(const std::string& e) 
    { 
      QDPIO::cout << "Failed to read stout action XML:" << e << endl;
    }
    

    // Sanity check
    if (smear_in_this_dirP.size() != Nd)
    {
      QDPIO::cerr << __func__ << ": invalid size of smear_in_this_dirP, expecting size=Nd" << endl;
      QDP_abort(1);
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
    }

    // Zero out any directions that are not smeared
    for(int mu=0; mu < Nd; mu++) 
    { 
      if( ! smear_in_this_dirP[mu] )
      {
	for(int nu=0; nu < Nd; nu++) { 
	  rho[mu][nu] = 0;
	  rho[nu][mu] = 0;
	}
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
    int version = 2;
    write(xml, "version", version);
    write(xml, "rho", p.sm_fact);
    write(xml, "n_smear", p.n_smear);
    write(xml, "smear_in_this_dirP", p.smear_in_this_dirP);
    pop(xml);
  }


};    

