// $Id: param_io.cc,v 1.1 2004-01-06 01:30:48 edwards Exp $
/*! \file
 *  \brief Various parameter readers/writers for main programs
 */

#include "chromabase.h"

using namespace QDP;

// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType         FermTypeP;
  int              FermAct;
  multi1d<Real>    mass;       // Quark mass and **NOT** kappa
 
  ChiralParam_t    chiralParam;
  SrceSinkParam_t  srceSinkParam;
  InvertParam_t    invParam;

  CfgType          cfg_type;   // storage order for stored gauge configuration
  PropType         prop_type;  // storage order for stored propagator

  int              j_decay;    // decay direction

  multi1d<int>     seq_src;    // integer array holding sequential source numbers

  multi1d<int>     nrow;
  multi1d<int>     boundary;
};

struct Cfg_t
{
  string       cfg_file;
};

struct Prop_t
{
  string       source_file;
  string       prop_file;
};


//! Read a fermion type enum
void read(XMLReader& xml, const string& path, FermType& param)
{
  string ferm_type_str;
  read(xml, path, ferm_type_str);
  if (ferm_type_str == "WILSON") {
    param = FERM_TYPE_WILSON;
  } else if (ferm_type_str == "STAGGERED") {
    param = FERM_TYPE_STAGGERED;
  else 
  {
    QDPIO::cerr << "Unsupported fermion type" << endl;
    QDP_abort(1);
  }
}

//! Read a configuration type enum
void read(XMLReader& xml, const string& path, CfgType& param)
{
  string cfg_type_str;
  read(xml, path, cfg_type_str);
  if (cfg_type_str == "SZIN") {
    param = CFG_TYPE_SZIN;
  } else if (cfg_type_str == "NERSC") {
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
  string cfg_type_str;
  read(xml, path, cfg_type_str);
  if (cfg_type_str == "SZIN") {
    param = CFG_TYPE_SZIN;
  } else if (cfg_type_str == "NERSC") {
    param = CFG_TYPE_NERSC;
  else 
  {
    QDPIO::cerr << "Unsupported propagator type" << endl;
    QDP_abort(1);
  }
}


//! Read a wave-function type enum
void read(XMLReader& xml, const string& path, WvfType& param)
{
  string wvf_type_str;
  read(xml, path, wvf_type_str);
  if (wvf_type_str == "GAUGE_INV_GAUSSIAN") {
    param = WVF_TYPE_GAUGE_INV_GAUSSIAN;
  } 
  else 
  {
    QDPIO::cerr << "Unsupported gauge-invariant wvf_type" << endl;
    QDP_abort(1);
  }
}


//! Read a inverter type enum
void read(XMLReader& xml, const string& path, InvType& param)
{
  string inv_type_str;
  read(xml, path, inv_type_str);
  if (inv_type_str == "CG_INVERTER") {
    param = CG_INVERTER;
  } else if (inv_type_str == "MR_INVERTER") {
    param = MR_INVERTER;
  } else if (inv_type_str == "BICG_INVERTER") {
    param = BICG_INVERTER;
  else 
  {
    QDPIO::cerr << "Unsupported inverter type" << endl;
    QDP_abort(1);
  }
}


//
void read(XMLReader& xml, const string& path, Cfg_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "cfg_file", input.cfg_file);
}


//
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
}


//! Read a smearing param struct
void read(XMLReader& xml, const string& path, SmearingParam_t& param)
{
  XMLReader paramtop(xml, path);

  {
    string wvf_type_str;
    read(paramtop, "wvf_type", wvf_type_str);
    if (wvf_type_str == "GAUGE_INV_GAUSSIAN") {
      param.wvf_type = WVF_TYPE_GAUGE_INV_GAUSSIAN;
    } 
    else 
    {
      QDPIO::cerr << "Unsupported gauge-invariant wvf_type." << endl;
      QDPIO::cerr << "  wvf_type = " << wvf_type_str << endl;
      QDP_abort(1);
    }
  }

  read(paramtop, "wvf_param", param.wvf_param);
  read(paramtop, "wvfIntPar", param.wvfIntPar);
}



//
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



//
void read(XMLReader& xml, const string& path, SrceSinkParam_t& param)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "Pt_src", param.Pt_src);
  read(paramtop, "Sl_src", param.Sl_src);
  read(paramtop, "Pt_snk", param.Pt_snk);
  read(paramtop, "Sl_snk", param.Sl_snk);

  if (param.Sl_src || param.Sl_snk)
    read(paramtop, "smearParam", param.smearParam);

  read(paramtop, "t_srce", param.t_srce);
  read(paramtop, "t_sink", param.t_sink);
  read(paramtop, "sink_mom", param.sink_mom);
}

//
void read(XMLReader& xml, const string& path, InvertParam_t& param)
{
  XMLReader paramtop(xml, path);

//  read(paramtop, "invType", param.param.invType);
  param.invType = CG_INVERTER;   //need to fix this
  read(paramtop, "RsdCG", param.RsdCG);
  read(paramtop, "MaxCG", param.MaxCG);

  param.MROver = 1;
}


