// $Id: wilson_fermact_params_w.cc,v 3.0 2006-04-03 04:58:47 edwards Exp $
/*! \file
 *  \brief Wilson fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"

#include "io/param_io.h"

namespace Chroma
{
  //! Default constructor
  WilsonFermActParams::WilsonFermActParams()
  {
    Mass = 0.0;
  }


  //! Read parameters
  WilsonFermActParams::WilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }

    //  Read optional anisoParam.
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, WilsonFermActParams& param)
  {
    WilsonFermActParams tmp(xml, path);
    param = tmp;
  }

  //! Writer parameters
  void write(XMLWriter& xml, const string& path, const WilsonFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    if (param.anisoParam.anisoP)
      write(xml, "AnisoParam", param.anisoParam);

    pop(xml);
  }
}
