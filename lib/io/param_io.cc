// $Id: param_io.cc,v 1.37 2005-02-22 02:49:59 edwards Exp $
/*! \file
 *  \brief Various parameter readers/writers for main programs
 */

#include "chromabase.h"
#include "io/param_io.h"


namespace Chroma { 

//! Convert a Kappa to a mass
Real kappaToMass(const Real& Kappa)
{
  return 1.0/(2*Kappa) - Nd;
}


//! Convert a Kappa to a mass
multi1d<Real> kappaToMass(const multi1d<Real>& Kappa)
{
  multi1d<Real> Mass(Kappa.size());

  for(int i=0; i < Kappa.size(); ++i)
    Mass[i] = 1.0/(2*Kappa[i]) - Nd;

  return Mass;
}


//! Convert a Kappa to a mass
Real massToKappa(const Real& Mass)
{
  return 0.5/(Nd + Mass);
}


//! Convert a mass to a Kappa
multi1d<Real> massToKappa(const multi1d<Real>& Mass)
{
  multi1d<Real> Kappa(Mass.size());

  for(int i=0; i < Kappa.size(); ++i)
    Kappa[i] = 0.5/(Nd + Mass[i]);

  return Kappa;
}


//! Read the input version
void read(XMLReader& xml, const string& path, IO_version_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "version", param.version);
}


//! Configuration input
void read(XMLReader& xml, const string& path, Cfg_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "cfg_type", input.cfg_type);
  read(inputtop, "cfg_file", input.cfg_file);
}


#if 1
//! Initialize a smearing param struct
SmearingParam_t::SmearingParam_t()
{
  wvf_kind    = WVF_KIND_GAUSSIAN;
  wvf_param   = 0;
  wvfIntPar   = 0;
}

//! Initialize a anisotropy param struct
AnisoParam_t::AnisoParam_t()
{
  anisoP = false;
  t_dir  = Nd-1;   // doesn't matter - should not be used
  xi_0   = 1;
  nu     = 1;
}

//! Initialize a chiral param struct
ChiralParam_t::ChiralParam_t()
{
  OverMass = 0;
  N5       = 0;
  a5       = 1;
  NWilsVec = 0;
}
#endif


//! Read a anisotropy param struct
void read(XMLReader& xml, const string& path, AnisoParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "anisoP", param.anisoP);
  read(paramtop, "t_dir", param.t_dir);
  read(paramtop, "xi_0", param.xi_0);
  read(paramtop, "nu", param.nu);
}


//! Read a smearing param struct
void read(XMLReader& xml, const string& path, SmearingParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "wvf_kind", param.wvf_kind);
  read(paramtop, "wvf_param", param.wvf_param);
  read(paramtop, "wvfIntPar", param.wvfIntPar);
}


//! Read chiral action like parameters
void read(XMLReader& xml, const string& path, ChiralParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "OverMass", param.OverMass);
  read(paramtop, "N5", param.N5);

  string xpath;
  xpath = "a5";
  if (paramtop.count(xpath) != 0)
    read(paramtop, xpath, param.a5);
  else
    param.a5 = 1;

  xpath = "NWilsVec";
  if (paramtop.count(xpath) != 0)
    read(paramtop, xpath, param.NWilsVec);
  else
    param.NWilsVec = 0;
}


//! Read inverter parameters
void read(XMLReader& xml, const string& path, InvertParam_t& param)
{
  XMLReader paramtop(xml, path);

  try {
    read(paramtop, "invType", param.invType);
    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);
    param.MROver = 1;
    
    if( paramtop.count("RsdCGPrec") == 1 ) {
      read(paramtop, "RsdCGPrec", param.RsdCGPrec);
    }
    else {
      param.RsdCGPrec = param.RsdCG;
    }

    if( paramtop.count("MaxCGPrec") == 1 ) {
      read(paramtop, "MaxCGPrec", param.MaxCGPrec);
    }
    else {
      param.MaxCGPrec = param.MaxCG;
    }

  }
  catch( const string& e ) { 
    QDPIO::cerr << "Caught exception : " << e << endl;
    QDP_abort(1);
  }

}

//! Read inverter parameters
void read(XMLReader& xml, const string& path, MultiInvertParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "invType", param.invType);
  read(paramtop, "RsdCG", param.RsdCG);
  read(paramtop, "MaxCG", param.MaxCG);

  param.MROver = 1;
}


//---------------------------- Writers -----------------------------
//! Write a config struct
void write(XMLWriter& xml, const string& path, const Cfg_t& cfg)
{
  push(xml, "Cfg");
  write(xml, "cfg_type", cfg.cfg_type);
  write(xml, "cfg_file", cfg.cfg_file);
  pop(xml);
}


//! Write a anisotropy param struct
void write(XMLWriter& xml, const string& path, const AnisoParam_t& param)
{
  push(xml, path);

  write(xml, "anisoP", param.anisoP);
  write(xml, "t_dir", param.t_dir);
  write(xml, "xi_0", param.xi_0);
  write(xml, "nu", param.nu);

  pop(xml);
}


//! Write a smearing param struct
void write(XMLWriter& xml, const string& path, const SmearingParam_t& param)
{
  push(xml, path);

  write(xml, "wvf_kind", param.wvf_kind);
  write(xml, "wvf_param", param.wvf_param);
  write(xml, "wvfIntPar", param.wvfIntPar);

  pop(xml);
}


//! Write chiral action like parameters
void write(XMLWriter& xml, const string& path, const ChiralParam_t& param)
{
  push(xml, path);

  write(xml, "OverMass", param.OverMass);
  write(xml, "N5", param.N5);
  write(xml, "a5", param.a5);
  write(xml, "NWilsVec", param.NWilsVec);

  pop(xml);
}


//! Write inverter parameters
void write(XMLWriter& xml, const string& path, const InvertParam_t& param)
{
  push(xml, path);

  write(xml, "invType", param.invType);
  write(xml, "RsdCG", param.RsdCG);
  write(xml, "MaxCG", param.MaxCG);
  write(xml, "MROver", param.MROver);
  write(xml, "RsdCGPrec", param.RsdCGPrec);
  write(xml, "MaxCGPrec", param.MaxCGPrec);
  pop(xml);
}

//! Write inverter parameters
void write(XMLWriter& xml, const string& path, const MultiInvertParam_t& param)
{
  push(xml, path);

  write(xml, "invType", param.invType);
  write(xml, "RsdCG", param.RsdCG);
  write(xml, "MaxCG", param.MaxCG);
  write(xml, "MROver", param.MROver);

  pop(xml);
}

}  // end namespace Chroma
