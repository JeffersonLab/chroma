// $Id: syssolver_OPTeigbicg_params.cc,v 1.2 2009-10-22 20:57:26 kostas Exp $
/*! \file
 *  \brief Params of EigCG inverter
 */

#include "actions/ferm/invert/syssolver_OPTeigbicg_params.h"
#include "io/qprop_io.h"

namespace Chroma
{

  //! File output                                                                                        
  void write(XMLWriter& xml, const string& path, const SysSolverOptEigBiCGParams::File_t& input){
    push(xml, path);

    write(xml, "read", input.read);
    write(xml, "write", input.write);
    write(xml, "file_name", input.file_name);
    write(xml, "file_volfmt", input.file_volfmt);

    pop(xml);
  }


  //! File output                                                                                             
  void read(XMLReader& xml, const string& path, SysSolverOptEigBiCGParams::File_t& input){
    XMLReader inputtop(xml, path);

    read(inputtop, "read", input.read);
    read(inputtop, "write", input.write);
    read(inputtop, "file_name", input.file_name);
    read(inputtop, "file_volfmt", input.file_volfmt);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverOptEigBiCGParams& param)
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

    read(paramtop, "cleanUpEvecs", param.cleanUpEvecs);

    if(paramtop.count("FileIO")!=0){
      read(paramtop, "FileIO", param.file);
    }

    if( paramtop.count("sort_option") > 0 ) { 
      read(paramtop, "sort_option", param.sort_option);
    }

    if( paramtop.count("epsi") > 0 )
      read(paramtop, "epsi", param.epsi);

    if( paramtop.count("ConvTestOpt")>0 ) 
      read(paramtop, "ConvTestOpt", param.ConvTestOpt);
    
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverOptEigBiCGParams& param)
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

    write(xml, "NormAest", param.NormAest);
    write(xml, "cleanUpEvecs", param.cleanUpEvecs);
    write(xml, "eigen_id", param.eigen_id);
    
    write(xml, "sort_option", param.sort_option);
    write(xml, "epsi", param.epsi);
    write(xml, "ConvTestOpt", param.ConvTestOpt);
    
    write(xml, "FileIO",param.file);

    pop(xml);
  }

  //! Default constructor
  SysSolverOptEigBiCGParams::SysSolverOptEigBiCGParams()
  {
    defaults();
  }

  //! Read parameters
  SysSolverOptEigBiCGParams::SysSolverOptEigBiCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
