// $Id: syssolver_eigcg_params.cc,v 1.3 2007-09-28 02:09:42 kostas Exp $
/*! \file
 *  \brief Params of EigCG inverter
 */

#include "actions/ferm/invert/syssolver_eigcg_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverEigCGParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);
    read(paramtop, "eigen_id", param.eigen_id);
    
    if( paramtop.count("RsdCGRestart") > 0 ) { 
      read(paramtop, "RsdCGRestart", param.RsdCGRestart);
    }
    else {
      param.RsdCGRestart = 33.0*param.RsdCG;
    }
    read(paramtop, "Nmax", param.Nmax);
    read(paramtop, "Neig", param.Neig);
    param.Neig_max = 0 ;
    if(paramtop.count("Neig_max")!=0){
      read(paramtop, "Neig_max", param.Neig_max);
    }

    param.vPrecCGvecs = 0 ;
    param.vPrecCGvecStart = 0 ;
    if(paramtop.count("vPrecCGvecs")!=0){
      read(paramtop, "vPrecCGvecs", param.vPrecCGvecs);
      read(paramtop, "vPrecCGvecStart", param.vPrecCGvecStart);
    }

    
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverEigCGParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "EIG_CG_INVERTER");
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "Nmax", param.Nmax);
    write(xml, "Neig", param.Neig);
    write(xml, "Neig_max", param.Neig_max);
    write(xml, "RsdCGRestart", param.RsdCGRestart);
    write(xml, "vPrecCGvecs", param.vPrecCGvecs);
    write(xml, "vPrecCGvecs", param.vPrecCGvecStart);
    write(xml, "eigen_id", param.eigen_id);

    pop(xml);
  }

  //! Default constructor
  SysSolverEigCGParams::SysSolverEigCGParams()
  {
    RsdCG = zero;
    MaxCG = 0;
    RsdCGRestart = RsdCG;
    Neig =0 ;
    Neig_max =0 ;
    vPrecCGvecs=0 ;
    eigen_id="NULL";
  }

  //! Read parameters
  SysSolverEigCGParams::SysSolverEigCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
