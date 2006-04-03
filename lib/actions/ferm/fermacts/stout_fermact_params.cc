// -*- C++ -*-
// $Id: stout_fermact_params.cc,v 3.0 2006-04-03 04:58:46 edwards Exp $

#include "actions/ferm/fermacts/stout_fermact_params.h"
#include <sstream>


namespace Chroma { 


  StoutFermActParams::StoutFermActParams(XMLReader& in, const std::string& path) {

    try { 
      XMLReader paramtop(in, path);

      XMLReader internal_fermact_reader(paramtop, "InternalFermionAction");
      std::ostringstream os; 
      internal_fermact_reader.print(os);
      internal_fermact = os.str();
      
      QDPIO::cout << internal_fermact << endl;
      
      read(paramtop, "./rho", rho);
      read(paramtop, "./n_smear", n_smear);
      if( paramtop.count("./orthog_dir") == 1 ) { 
	read(paramtop, "./orthog_dir", orthog_dir);
      }
      else { 
	// default value for orthog dir is 3 -- smear in space only
	QDPIO::cout << "Using Default value: orthog_dir = 3, spatial only smearing" << endl;
      }
    }
    catch(const std::string& e) { 
      QDPIO::cout << "Failed to read stout action XML:" << e << endl;
    }
  }

  void read(XMLReader& xml, const std::string& path, StoutFermActParams& p)
  {
    StoutFermActParams tmp_p(xml, path);
    p = tmp_p;
  }

  void write(XMLWriter& xml, const std::string& path, const StoutFermActParams& p) 
  {
    push(xml, path);
    std::istringstream is(p.internal_fermact);
    XMLReader r(is);
    xml << r;
    write(xml, "rho", p.rho);
    write(xml, "n_smear", p.n_smear);
    write(xml, "orthog_dir", p.orthog_dir);
    pop(xml);
  }


};    

