// $Id: syssolver_OPTeigcg_params.cc,v 1.3 2008-04-03 15:58:43 kostas Exp $
/*! \file
 *  \brief Params of EigCG inverter
 */

#include "actions/ferm/invert/syssolver_OPTeigcg_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverOptEigCGParams& param)
  {
    XMLReader paramtop(xml, path);

    param.defaults();

    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);

    if( paramtop.count("PrintLevel") > 0 ) { 
      read(paramtop, "PrintLevel", param.PrintLevel);
    }

    read(paramtop, "Nmax", param.Nmax);
    read(paramtop, "Neig", param.Neig);
    if(paramtop.count("Neig_max")!=0){
      read(paramtop, "Neig_max", param.Neig_max);
    }

    if( paramtop.count("esize") > 0 ) { 
      read(paramtop, "esize", param.esize);
    }

    read(paramtop, "eigen_id", param.eigen_id);
    
    if( paramtop.count("restartTol") > 0 ) { 
      read(paramtop, "restartTol", param.restartTol);
    }

    if( paramtop.count("NormAest") > 0 ) { 
      read(paramtop, "NormAest", param.NormAest);
    }

    if( paramtop.count("updateRestartTol") > 0 ) { 
      read(paramtop, "updateRestartTol", param.updateRestartTol);
    }

    read(paramtop, "cleanUpEvecs", param.cleanUpEvecs);
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverOptEigCGParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "OPT_EIG_CG_INVERTER");
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "Nmax", param.Nmax);
    write(xml, "Neig", param.Neig);
    write(xml, "Neig_max", param.Neig_max);
    write(xml, "restartTol", param.restartTol);
    write(xml, "updateRestartTol", param.updateRestartTol);
    write(xml, "NormAest", param.NormAest);
    write(xml, "cleanUpEvecs", param.cleanUpEvecs);
    write(xml, "eigen_id", param.eigen_id);

    pop(xml);
  }

  //! Default constructor
  SysSolverOptEigCGParams::SysSolverOptEigCGParams()
  {
    defaults();
  }

  //! Read parameters
  SysSolverOptEigCGParams::SysSolverOptEigCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
