// -*- C++ -*-
// $Id: stout_fermact_params.cc,v 2.1 2005-10-02 03:08:49 bjoo Exp $

#include "actions/ferm/fermacts/stout_fermact_params.h"
#include <sstream>


namespace Chroma { 


  StoutFermActParams::StoutFermActParams(XMLReader& in, const std::string& path) {
    XMLReader paramtop(in, path);

    XMLReader internal_fermact_reader(paramtop, "./InternalFermionAction");
    std::ostringstream os; 
    internal_fermact_reader.print(os);
    internal_fermact = os.str();

    QDPIO::cout << internal_fermact << endl;

    read(paramtop, "./rho", rho);
    read(paramtop, "./n_smear", n_smear);
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
    pop(xml);
  }


};    

