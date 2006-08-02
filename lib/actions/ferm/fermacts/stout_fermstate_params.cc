// -*- C++ -*-
// $Id: stout_fermstate_params.cc,v 1.1 2006-08-02 04:10:22 edwards Exp $

#include "actions/ferm/fermacts/stout_fermstate_params.h"


namespace Chroma 
{

  StoutFermStateParams::StoutFermStateParams(XMLReader& in, const std::string& path) 
  {
    try 
    { 
      XMLReader paramtop(in, path);

      read(paramtop, "rho", rho);
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
  }

  void read(XMLReader& xml, const std::string& path, StoutFermStateParams& p)
  {
    StoutFermStateParams tmp_p(xml, path);
    p = tmp_p;
  }

  void write(XMLWriter& xml, const std::string& path, const StoutFermStateParams& p) 
  {
    push(xml, path);
    write(xml, "rho", p.rho);
    write(xml, "n_smear", p.n_smear);
    write(xml, "orthog_dir", p.orthog_dir);
    pop(xml);
  }


};    

