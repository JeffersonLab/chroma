// $Id: syssolver_eigcg_params.cc,v 1.9 2008-12-15 05:02:06 kostas Exp $
/*! \file
 *  \brief Params of EigCG inverter
 */

#include "actions/ferm/invert/syssolver_eigcg_params.h"
#include "io/qprop_io.h"

namespace Chroma
{

  //! File output                                                                       
  void write(XMLWriter& xml, const string& path, const SysSolverEigCGParams::File_t& input){
    push(xml, path);

    write(xml, "read", input.read);
    write(xml, "write", input.write);
    write(xml, "file_name", input.file_name);
    write(xml, "file_volfmt", input.file_volfmt);

    pop(xml);
  }


  //! File output
  void read(XMLReader& xml, const string& path, SysSolverEigCGParams::File_t& input){
    XMLReader inputtop(xml, path);

    read(inputtop, "read", input.read);
    read(inputtop, "write", input.write);
    read(inputtop, "file_name", input.file_name);
    read(inputtop, "file_volfmt", input.file_volfmt);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverEigCGParams& param)
  {
    XMLReader paramtop(xml, path);
    
    param.defaults();

    read(paramtop, "invType", param.invType);
    
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

     if( paramtop.count("restartTol") > 0 ) { 
      read(paramtop, "restartTol", param.restartTol);
    }

    if( paramtop.count("NormAest") > 0 ) { 
      read(paramtop, "NormAest", param.NormAest);
    }

    if( paramtop.count("updateRestartTol") > 0 ) { 
      read(paramtop, "updateRestartTol", param.updateRestartTol);
    }

    if(paramtop.count("vPrecCGvecs")!=0){
      read(paramtop, "vPrecCGvecs", param.vPrecCGvecs);
      read(paramtop, "vPrecCGvecStart", param.vPrecCGvecStart);
    }

    read(paramtop, "cleanUpEvecs", param.cleanUpEvecs);
    read(paramtop, "eigen_id", param.eigen_id);

    if(paramtop.count("FileIO")!=0){
      read(paramtop, "FileIO", param.file);
    }

  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverEigCGParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
  
    write(xml, "invType", param.invType);
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "Nmax", param.Nmax);
    write(xml, "Neig", param.Neig);
    write(xml, "Neig_max", param.Neig_max);
    write(xml, "restartTol", param.restartTol);
    write(xml, "updateRestartTol", param.updateRestartTol);
    write(xml, "NormAest", param.NormAest);
    write(xml, "vPrecCGvecs", param.vPrecCGvecs);
    write(xml, "vPrecCGvecs", param.vPrecCGvecStart);
    write(xml, "cleanUpEvecs", param.cleanUpEvecs);
    write(xml, "eigen_id", param.eigen_id);

    write(xml, "FileIO",param.file);

    pop(xml);
  }

  //! Default constructor
  SysSolverEigCGParams::SysSolverEigCGParams()
  {
    defaults();
  }

  //! Read parameters
  SysSolverEigCGParams::SysSolverEigCGParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
