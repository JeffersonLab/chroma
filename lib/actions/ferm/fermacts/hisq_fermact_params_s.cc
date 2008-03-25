// $Id: hisq_fermact_params_s.cc,v 1.2 2008-03-25 10:53:36 mcneile Exp $
/*! \file
 *  \brief Hisq fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/hisq_fermact_params_s.h"

#include "io/param_io.h"

namespace Chroma
{
  //! Default constructor
  HisqFermActParams::HisqFermActParams()
  {
    Mass = 0.0;
    u0  = 1.0;
    epsilon = 0.0 ;
  }


  //! Read parameters
  HisqFermActParams::HisqFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "Mass", Mass);
    read(paramtop, "u0", u0);

    if( paramtop.count("epsilon") > 0 ) {
      read(paramtop, "epsilon", epsilon);
    }    
    else
      {
	epsilon = 0.0 ;
      }


    //  Read optional anisoParam.
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, HisqFermActParams& param)
  {
    HisqFermActParams tmp(xml, path);
    param = tmp;
  }

  //! Writer parameters
  void write(XMLWriter& xml, const string& path, const HisqFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    write(xml, "u0", param.u0);
    write(xml, "epsilon", param.epsilon);
    pop(xml);
  }
}
