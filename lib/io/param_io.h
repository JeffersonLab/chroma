// -*- C++ -*-
// $Id: param_io.h,v 1.1 2004-01-06 01:30:48 edwards Exp $
/*! \file
 *  \brief Reunitarize (to a SU(N)) inplace the matrix A under some option
 */

#ifndef __param_io_h__
#define __param_io_h__


/*
 *  Here we have various temporary definitions
 */
enum CfgType {
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_UNKNOWN
} ;

enum PropType {
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
  PROP_TYPE_UNKNOWN
} ;

enum FermType {
  FERM_TYPE_WILSON,
  FERM_TYPE_STAGGERED,
  FERM_TYPE_UNKNOWN
};


/*
 * Input 
 */
struct IO_version_t
{
  int version;
};


//! Parameters for chiral fermion actions
struct ChiralParam_t
{
  Real       overMass;
  Real       N5;
  Real       a5;
  int        nWilsVec;
};


//! Parameters for anisotropy
/*! NOT USED YET */
struct AnisoParam_t
{
  bool       anisoP;
  int        t_dir;
  Real       xi_0;
  Real       xiF_0;
  Real       Wilsr_s;
};


//! Parameters for sources and sinks
struct SmearingParam_t
{
  WvfType       wvf_type;
  multi1d<Real> wvf_param;
  multi1d<int>  wvfIntPar;
};


//! Parameters for sources and sinks
struct SrceSinkParam_t
{
  bool             Pt_src;   // point source
  bool             Sl_src;   // shell source
  bool             Pt_snk;   // point sink
  bool             Sl_snk;   // shell sink

  SmearingParam_t  smearParam;

  multi1d<int>     t_srce;
  multi1d<int>     sink_mom;
  int              t_sink;
};


//! Parameters for inverter
struct InvertParam_t
{
  InvType       invType;   // Inverter type
  Real          MROver;
  Real          RsdCG;
  int           MaxCG;	   // Iteration parameters
};

#endif
