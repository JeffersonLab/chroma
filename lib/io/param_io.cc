// $Id: param_io.cc,v 1.18 2004-04-14 12:53:21 bjoo Exp $
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


//! Read a QDP volume format type
void read(XMLReader& xml, const string& path, QDP_volfmt_t& param)
{
  string volfmt_str;
  read(xml, path, volfmt_str);
  if (volfmt_str == "SINGLEFILE")
    param = QDPIO_SINGLEFILE;
  else if (volfmt_str == "MULTIFILE")
    param = QDPIO_MULTIFILE;
  else 
  {
    QDPIO::cerr << "Unsupported QDP volume format type" << endl;
    QDP_abort(1);
  }
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
  if (cfg_type_str == "MILC")
    param = CFG_TYPE_MILC;
  else if (cfg_type_str == "NERSC")
    param = CFG_TYPE_NERSC;
  else if (cfg_type_str == "SCIDAC")
    param = CFG_TYPE_SCIDAC;
  else if (cfg_type_str == "SZIN")
    param = CFG_TYPE_SZIN;
  else if (cfg_type_str == "SZINQIO")
    param = CFG_TYPE_SZINQIO;
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

  read(inputtop, "cfg_type", input.cfg_type);
  read(inputtop, "cfg_file", input.cfg_file);
}


//! Initialize a smearing param struct
void initHeader(SmearingParam_t& param)
{
  param.wvf_kind    = WVF_KIND_GAUSSIAN;
  param.wvf_param   = 0;
  param.wvfIntPar   = 0;
}

//! Initialize a anisotropy param struct
void initHeader(AnisoParam_t& param)
{
  param.anisoP = false;
  param.t_dir  = Nd-1;   // doesn't matter - should not be used
  param.xi_0   = 1;
  param.nu     = 1;
}

//! Initialize a chiral param struct
void initHeader(ChiralParam_t& param)
{
  param.OverMass = 0;
  param.N5       = 0;
  param.a5       = 1;
  param.NWilsVec = 0;
}


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

  read(paramtop, "invType", param.invType);
  read(paramtop, "RsdCG", param.RsdCG);
  read(paramtop, "MaxCG", param.MaxCG);

  param.MROver = 1;
}


//---------------------------- Writers -----------------------------
//! Write a QDP volume format type
void write(XMLWriter& xml, const string& path, QDP_volfmt_t param)
{
  string volfmt_str;
  if (param == QDPIO_SINGLEFILE)
    volfmt_str = "SINGLEFILE";
  else if (param == QDPIO_MULTIFILE)
    volfmt_str = "MULTIFILE";
  else 
  {
    QDPIO::cerr << "Unsupported QDP volume format type" << endl;
    QDP_abort(1);
  }
  write(xml, path, volfmt_str);
}

//! Write a fermion type enum
void write(XMLWriter& xml, const string& path, FermType param)
{
  string ferm_type_str;
  if (param == FERM_TYPE_WILSON)
    ferm_type_str = "WILSON";
  else if (param == FERM_TYPE_STAGGERED)
    ferm_type_str = "STAGGERED";
  else 
  {
    QDPIO::cerr << "Unsupported fermion type" << endl;
    QDP_abort(1);
  }
  write(xml, path, ferm_type_str);
}

//! Write a fermion action type enum
void write(XMLWriter& xml, const string& path, FermActType param)
{
  string ferm_type_str;
  if (param == FERM_ACT_WILSON)
    ferm_type_str = "WILSON";
  else if (param == FERM_ACT_UNPRECONDITIONED_WILSON)
    ferm_type_str = "UNPRECONDITIONED_WILSON";
  else if (param == FERM_ACT_PARITY_BREAKING_WILSON)
    ferm_type_str = "PARITY_BREAKING_WILSON";
  else if (param == FERM_ACT_CLOVER)
    ferm_type_str = "CLOVER";
  else if (param == FERM_ACT_UNPRECONDITIONED_CLOVER)
    ferm_type_str = "UNPRECONDITONED_CLOVER";
  else if (param == FERM_ACT_DWF)
    ferm_type_str = "DWF";
  else if (param == FERM_ACT_UNPRECONDITIONED_DWF)
    ferm_type_str = "UNPRECONDITIONED_DWF";
  else if (param == FERM_ACT_PROJECTED_DWF)
    ferm_type_str = "PROJECTED_DWF";
  else if (param == FERM_ACT_ZOLOTAREV_4D)
    ferm_type_str = "ZOLOTAREV_4D";
  else if (param == FERM_ACT_OVERLAP_DWF)
    ferm_type_str = "OVERLAP_DWF";
  else if (param == FERM_ACT_EXTENDED_OVERLAP)
    ferm_type_str = "EXTENDED_OVERLAP";
  else if (param == FERM_ACT_UNPRECONDITIONED_EXTENDED_OVERLAP)
    ferm_type_str = "UNPRECONDITIONED_EXTENDED_OVERLAP";
  else if (param == FERM_ACT_SMEARED_LAPLACIAN_WILSON)
    ferm_type_str = "SMEARED_LAPLACIAN_WILSON";
  else if (param == FERM_ACT_PLANAR_WILSON)
    ferm_type_str = "PLANAR_WILSON";
  else if (param == FERM_ACT_HAMBER_WU)
    ferm_type_str = "HAMBER_WU";
  else if (param == FERM_ACT_STAGGERED)
    ferm_type_str = "STAGGERED";
  else if (param == FERM_ACT_NAIK)
    ferm_type_str = "NAIK";
  else if (param == FERM_ACT_ASQTAD)
    ferm_type_str = "ASQTAD";
  else 
  {
    QDPIO::cerr << "Unsupported fermion action type" << endl;
    QDP_abort(1);
  }
  write(xml, path, ferm_type_str);
}

//! Write a configuration type enum
void write(XMLWriter& xml, const string& path, CfgType param)
{
  string cfg_type_str;
  if (param == CFG_TYPE_MILC)
    cfg_type_str = "MILC";
  else if (param == CFG_TYPE_NERSC)
    cfg_type_str = "NERSC";
  else if (param == CFG_TYPE_SCIDAC)
    cfg_type_str = "SCIDAC";
  else if (param == CFG_TYPE_SZIN)
    cfg_type_str = "SZIN";
  else if (param == CFG_TYPE_SZINQIO)
    cfg_type_str = "SZINQIO";
  else 
  {
    QDPIO::cerr << "Unsupported configuration type" << endl;
    QDP_abort(1);
  }
  write(xml, path, cfg_type_str);
}

//! Write a config struct
void write(XMLWriter& xml, const string& path, const Cfg_t& cfg)
{
  push(xml, "Cfg");
  write(xml, "cfg_type", cfg.cfg_type);
  write(xml, "cfg_file", cfg.cfg_file);
  pop(xml);
}


//! Write a propagator type enum
void write(XMLWriter& xml, const string& path, PropType param)
{
  string prop_type_str;
  if (param == PROP_TYPE_SZIN)
    prop_type_str = "SZIN";
  else if (param == PROP_TYPE_SCIDAC)
    prop_type_str = "SCIDAC";
  else 
  {
    QDPIO::cerr << "Unsupported propagator type" << endl;
    QDP_abort(1);
  }
  write(xml, path, prop_type_str);
}


//! Write a wave-function type enum
void write(XMLWriter& xml, const string& path, WvfKind param)
{
  string wvf_kind_str;
  if (param == WVF_KIND_GAUSSIAN)
    wvf_kind_str = "GAUSSIAN";
  else if (param == WVF_KIND_EXPONENTIAL)
    wvf_kind_str = "EXPONENTIAL";
  else if (param == WVF_KIND_GAUGE_INV_GAUSSIAN)
    wvf_kind_str = "GAUGE_INV_GAUSSIAN";
  else if (param == WVF_KIND_WUPPERTAL)
    wvf_kind_str = "WUPPERTAL";
  else if (param == WVF_KIND_JACOBI)
    wvf_kind_str = "JACOBI";
  else 
  {
    QDPIO::cerr << "Unsupported gauge-invariant wvf_kind" << endl;
    QDP_abort(1);
  }
  write(xml, path, wvf_kind_str);
}


//! Write a inverter type enum
void write(XMLWriter& xml, const string& path, InvType param)
{
  string inv_type_str;
  if (param == CG_INVERTER)
    inv_type_str = "CG_INVERTER";
  else if (param == MR_INVERTER)
    inv_type_str = "MR_INVERTER";
  else if (param == BICG_INVERTER)
    inv_type_str = "BICG_INVERTER";
  else 
  {
    QDPIO::cerr << "Unsupported inverter type" << endl;
    QDP_abort(1);
  }
  write(xml, path, inv_type_str);
}


//! Write a source type enum
void write(XMLWriter& xml, const string& path, SourceType param)
{
  string src_type_str;
  if (param == SRC_TYPE_POINT_SOURCE)
    src_type_str = "POINT_SOURCE";
  else if (param == SRC_TYPE_WALL_SOURCE)
    src_type_str = "WALL_SOURCE";
  else if (param == SRC_TYPE_SHELL_SOURCE)
    src_type_str = "SHELL_SOURCE";
  else if (param == SRC_TYPE_BNDST_SOURCE)
    src_type_str = "BNDST_SOURCE";
  else 
  {
    QDPIO::cerr << "Unsupported SourceType" << endl;
    QDP_abort(1);
  }
  write(xml, path, src_type_str);
}


//! Write a sink type enum
void write(XMLWriter& xml, const string& path, SinkType param)
{
  string src_type_str;
  if (param == SNK_TYPE_POINT_SINK)
    src_type_str = "POINT_SINK";
  else if (param == SNK_TYPE_WALL_SINK)
    src_type_str = "WALL_SINK";
  else if (param == SNK_TYPE_SHELL_SINK)
    src_type_str = "SHELL_SINK";
  else if (param == SNK_TYPE_BNDST_SINK)
    src_type_str = "BNDST_SINK";
  else 
  {
    QDPIO::cerr << "Unsupported SinkType" << endl;
    QDP_abort(1);
  }
  write(xml, path, src_type_str);
}


//! Write a wave type enum
void write(XMLWriter& xml, const string& path, WaveStateType param)
{
  string wave_type_str;
  if (param == WAVE_TYPE_S_WAVE)
    wave_type_str = "S_WAVE";
  else if (param == WAVE_TYPE_P_WAVE)
    wave_type_str = "P_WAVE";
  else if (param == WAVE_TYPE_D_WAVE)
    wave_type_str = "D_WAVE";
  else 
  {
    QDPIO::cerr << "Unsupported particle wave-state type" << endl;
    QDP_abort(1);
  }
  write(xml, path, wave_type_str);
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

  pop(xml);
}


