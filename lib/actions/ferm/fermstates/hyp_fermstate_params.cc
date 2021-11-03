/*! \file
 *  \brief Hyp fermstate params
 */

#include "actions/ferm/fermstates/hyp_fermstate_params.h"


namespace Chroma 
{

  HypFermStateParams::HypFermStateParams() 
  {
    n_smear = 0;

  }


  HypFermStateParams::HypFermStateParams(XMLReader& in, const std::string& path) 
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
      QDPIO::cout << "Failed to read hyp action XML:" << e << std::endl;
    }
    

  }

  void read(XMLReader& xml, const std::string& path, HypFermStateParams& p)
  {
    HypFermStateParams tmp_p(xml, path);
    p = tmp_p;
  }

  void write(XMLWriter& xml, const std::string& path, const HypFermStateParams& p) 
  {
    push(xml, path);
    int version = 1;
    write(xml, "version", version);
    write(xml, "n_smear", p.n_smear);
    pop(xml);
  }


}    

