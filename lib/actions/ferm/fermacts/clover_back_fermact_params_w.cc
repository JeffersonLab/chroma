/*! \file
 *  \brief Clover fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_back_fermact_params_w.h"
#include "io/param_io.h"

namespace Chroma
{

  //! Default constructor
  CloverBackFermActParams::CloverBackFermActParams():CloverFermActParams(),gamma(0)
  {
  }

  //! Read parameters
  CloverBackFermActParams::CloverBackFermActParams(XMLReader& xml, const std::string& path):CloverFermActParams(xml, path)
  {
    XMLReader paramtop(xml, path);
    
    if( paramtop.count("gamma") != 0) 
      read(paramtop,"gamma",gamma);
    else
      gamma=0 ;
      
  }

  //! Read parameters
  void read(XMLReader& xml, const std::string& path, CloverBackFermActParams& param)
  {
    CloverBackFermActParams tmp(xml, path);
    param = tmp;
  }


  //! Write parameters
  void write(XMLWriter& xml, const std::string& path, const CloverBackFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    write(xml, "u0", param.u0);

    if (param.anisoParam.anisoP)
    {
      write(xml, "clovCoeffR", param.clovCoeffR);
      write(xml, "clovCoeffT", param.clovCoeffT);
    }
    else
    {
      write(xml, "clovCoeff", param.clovCoeffR);
    }

    if (param.max_norm_usedP){
      write(xml, "MaxNorm", param.max_norm);
    }

    if (param.sub_zero_usedP == true ) {
      write(xml, "ZeroEnergy", param.sub_zero);
    }

    if (param.twisted_m_usedP == true ) { 
      write(xml, "TwistedM", param.twisted_m);
    }

    write(xml,"gamma",param.gamma) ;

    pop(xml);
  }


}

