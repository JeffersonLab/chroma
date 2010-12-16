// $Id: stout_fermstate_params.cc,v 1.4 2008-04-19 03:12:49 edwards Exp $
/*! \file
 *  \brief Hex fermstate params
 */

#include "actions/ferm/fermstates/hex_fermstate_params.h"


namespace Chroma 
{

  HexFermStateParams::HexFermStateParams() 
  {
    n_smear = 0;

  }


  HexFermStateParams::HexFermStateParams(XMLReader& in, const std::string& path) 
  {
    try
    {
      XMLReader paramtop(in, path);

      int version = 1;
      if( paramtop.count("version") == 1 )
	read(paramtop, "version", version);

      read(paramtop, "n_smear", n_smear);


    }
    catch(const std::string& e) 
    { 
      QDPIO::cout << "Failed to read hex action XML:" << e << endl;
    }
    

  }

  void read(XMLReader& xml, const std::string& path, HexFermStateParams& p)
  {
    HexFermStateParams tmp_p(xml, path);
    p = tmp_p;
  }

  void write(XMLWriter& xml, const std::string& path, const HexFermStateParams& p) 
  {
    push(xml, path);
    int version = 1;
    write(xml, "version", version);
    write(xml, "n_smear", p.n_smear);
    pop(xml);
  }


};    

