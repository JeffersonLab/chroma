// $Id: param_io.cc,v 1.4 2004-01-06 04:51:06 edwards Exp $
/*! \file
 *  \brief Various parameter readers/writers for main programs
 */

#include "chromabase.h"
#include "io/param_io.h"

using namespace QDP;

//! Read a fermion type enum
void read(XMLReader& xml, const string& path, FermType& param)
{
  string ferm_type_str;
  read(xml, path, ferm_type_str);
  if (ferm_type_str == "WILSON")
    param = FERM_TYPE_WILSON;
  else if (ferm_type_str == "STAGGERED")
    param = FERM_TYPE_STAGGERED;
  else 
  {
    QDPIO::cerr << "Unsupported fermion type" << endl;
    QDP_abort(1);
  }
}

//! Read a fermion action type enum
void read(XMLReader& xml, const string& path, FermActType& param)
{
  string ferm_type_str;
  read(xml, path, ferm_type_str);
  if (ferm_type_str == "WILSON")
    param = FERM_ACT_WILSON;
  else if (ferm_type_str == "DWF")
    param = FERM_ACT_DWF;
  else 
  {
    QDPIO::cerr << "Unsupported fermion action type" << endl;
    QDP_abort(1);
  }
}

//! Read a configuration type enum
void read(XMLReader& xml, const string& path, CfgType& param)
{
  string cfg_type_str;
  read(xml, path, cfg_type_str);
  if (cfg_type_str == "SZIN")
    param = CFG_TYPE_SZIN;
  else if (cfg_type_str == "NERSC")
    param = CFG_TYPE_NERSC;
  else 
  {
    QDPIO::cerr << "Unsupported configuration type" << endl;
    QDP_abort(1);
  }
}


//! Read a propagator type enum
void read(XMLReader& xml, const string& path, PropType& param)
{
  string prop_type_str;
  read(xml, path, prop_type_str);
  if (prop_type_str == "SZIN")
    param = PROP_TYPE_SZIN;
  else if (prop_type_str == "SCIDAC")
    param = PROP_TYPE_SCIDAC;
  else 
  {
    QDPIO::cerr << "Unsupported propagator type" << endl;
    QDP_abort(1);
  }
}


//! Read a wave-function type enum
void read(XMLReader& xml, const string& path, WvfKind& param)
{
  string wvf_kind_str;
  read(xml, path, wvf_kind_str);
  if (wvf_kind_str == "GAUGE_INV_GAUSSIAN")
    param = WVF_KIND_GAUGE_INV_GAUSSIAN;
  else 
  {
    QDPIO::cerr << "Unsupported gauge-invariant wvf_kind" << endl;
    QDP_abort(1);
  }
}


//! Read a inverter type enum
void read(XMLReader& xml, const string& path, InvType& param)
{
  string inv_type_str;
  read(xml, path, inv_type_str);
  if (inv_type_str == "CG_INVERTER")
    param = CG_INVERTER;
  else if (inv_type_str == "MR_INVERTER")
    param = MR_INVERTER;
  else if (inv_type_str == "BICG_INVERTER")
    param = BICG_INVERTER;
  else 
  {
    QDPIO::cerr << "Unsupported inverter type" << endl;
    QDP_abort(1);
  }
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

  read(inputtop, "cfg_file", input.cfg_file);
}


//! Read a smearing param struct
void read(XMLReader& xml, const string& path, SmearingParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "Wvf_kind", param.Wvf_kind);
  read(paramtop, "wvf_param", param.wvf_param);
  read(paramtop, "WvfIntPar", param.WvfIntPar);
}


//! Read chiral action like parameters
void read(XMLReader& xml, const string& path, ChiralParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "overMass", param.overMass);
  read(paramtop, "N5", param.N5);

  string xpath;
  xpath = "a5";
  if (paramtop.count(xpath) != 0)
    read(paramtop, xpath, param.a5);
  else
    param.a5 = 1;

  xpath = "nWilsVec";
  if (paramtop.count(xpath) != 0)
    read(paramtop, xpath, param.nWilsVec);
  else
    param.nWilsVec = 0;
}


//! Read inverter parameters
void read(XMLReader& xml, const string& path, InvertParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "invType", param.invType);
  read(paramtop, "RsdCG", param.RsdCG);
  read(paramtop, "MaxCG", param.MaxCG);

  param.MROver = 1;
}


