// $Id: param_io.cc,v 1.11 2004-01-13 03:57:32 edwards Exp $
/*! \file
 *  \brief Various parameter readers/writers for main programs
 */

#include "chromabase.h"
#include "io/param_io.h"

using namespace QDP;

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
  else if (ferm_type_str == "UNPRECONDITIONED_WILSON")
    param = FERM_ACT_UNPRECONDITIONED_WILSON;
  else if (ferm_type_str == "PARITY_BREAKING_WILSON")
    param = FERM_ACT_PARITY_BREAKING_WILSON;
  else if (ferm_type_str == "CLOVER")
    param = FERM_ACT_CLOVER;
  else if (ferm_type_str == "UNPRECONDITONED_CLOVER")
    param = FERM_ACT_UNPRECONDITIONED_CLOVER;
  else if (ferm_type_str == "DWF")
    param = FERM_ACT_DWF;
  else if (ferm_type_str == "UNPRECONDITIONED_DWF")
    param = FERM_ACT_UNPRECONDITIONED_DWF;
  else if (ferm_type_str == "PROJECTED_DWF")
    param = FERM_ACT_PROJECTED_DWF;
  else if (ferm_type_str == "ZOLOTAREV_4D")
    param = FERM_ACT_ZOLOTAREV_4D;
  else if (ferm_type_str == "OVERLAP_DWF")
    param = FERM_ACT_OVERLAP_DWF;
  else if (ferm_type_str == "EXTENDED_OVERLAP")
    param = FERM_ACT_EXTENDED_OVERLAP;
  else if (ferm_type_str == "UNPRECONDITIONED_EXTENDED_OVERLAP")
    param = FERM_ACT_UNPRECONDITIONED_EXTENDED_OVERLAP;
  else if (ferm_type_str == "SMEARED_LAPLACIAN_WILSON")
    param = FERM_ACT_SMEARED_LAPLACIAN_WILSON;
  else if (ferm_type_str == "PLANAR_WILSON")
    param = FERM_ACT_PLANAR_WILSON;
  else if (ferm_type_str == "HAMBER_WU")
    param = FERM_ACT_HAMBER_WU;
  else if (ferm_type_str == "STAGGERED")
    param = FERM_ACT_STAGGERED;
  else if (ferm_type_str == "NAIK")
    param = FERM_ACT_NAIK;
  else if (ferm_type_str == "ASQTAD")
    param = FERM_ACT_ASQTAD;
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
  if (wvf_kind_str == "GAUSSIAN")
    param = WVF_KIND_GAUSSIAN;
  else if (wvf_kind_str == "EXPONENTIAL")
    param = WVF_KIND_EXPONENTIAL;
  else if (wvf_kind_str == "GAUGE_INV_GAUSSIAN")
    param = WVF_KIND_GAUGE_INV_GAUSSIAN;
  else if (wvf_kind_str == "WUPPERTAL")
    param = WVF_KIND_WUPPERTAL;
  else if (wvf_kind_str == "JACOBI")
    param = WVF_KIND_JACOBI;
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


//! Read a source type enum
void read(XMLReader& xml, const string& path, SourceType& param)
{
  string src_type_str;
  read(xml, path, src_type_str);
  if (src_type_str == "POINT_SOURCE")
    param = SRC_TYPE_POINT_SOURCE;
  else if (src_type_str == "WALL_SOURCE")
    param = SRC_TYPE_WALL_SOURCE;
  else if (src_type_str == "SHELL_SOURCE")
    param = SRC_TYPE_SHELL_SOURCE;
  else if (src_type_str == "BNDST_SOURCE")
    param = SRC_TYPE_BNDST_SOURCE;
  else 
  {
    QDPIO::cerr << "Unsupported SourceType" << endl;
    QDP_abort(1);
  }
}


//! Read a sink type enum
void read(XMLReader& xml, const string& path, SinkType& param)
{
  string src_type_str;
  read(xml, path, src_type_str);
  if (src_type_str == "POINT_SINK")
    param = SNK_TYPE_POINT_SINK;
  else if (src_type_str == "WALL_SINK")
    param = SNK_TYPE_WALL_SINK;
  else if (src_type_str == "SHELL_SINK")
    param = SNK_TYPE_SHELL_SINK;
  else if (src_type_str == "BNDST_SINK")
    param = SNK_TYPE_BNDST_SINK;
  else 
  {
    QDPIO::cerr << "Unsupported SinkType" << endl;
    QDP_abort(1);
  }
}


//! Read a wave type enum
void read(XMLReader& xml, const string& path, WaveStateType& param)
{
  string wave_type_str;
  read(xml, path, wave_type_str);
  if (wave_type_str == "S_WAVE")
    param = WAVE_TYPE_S_WAVE;
  else if (wave_type_str == "P_WAVE")
    param = WAVE_TYPE_P_WAVE;
  else if (wave_type_str == "D_WAVE")
    param = WAVE_TYPE_D_WAVE;
  else 
  {
    QDPIO::cerr << "Unsupported particle wave-state type" << endl;
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

  read(paramtop, "invType", param.invType);
  read(paramtop, "RsdCG", param.RsdCG);
  read(paramtop, "MaxCG", param.MaxCG);

  param.MROver = 1;
}


