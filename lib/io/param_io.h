// -*- C++ -*-
// $Id: param_io.h,v 1.2 2004-01-06 02:09:01 edwards Exp $
/*! \file
 *  \brief Reunitarize (to a SU(N)) inplace the matrix A under some option
 */

#ifndef __param_io_h__
#define __param_io_h__

#include "meas/smear/sink_smear2_w.h"
#include "invtype.h"

/*
 *  Here we have various temporary definitions
 */
enum CfgType {
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
};

enum PropType {
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
};

enum FermType {
  FERM_TYPE_WILSON,
  FERM_TYPE_STAGGERED,
};


/*
 * Input 
 */
struct IO_version_t
{
  int version;
};

struct Cfg_t
{
  string       cfg_file;
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


//! Parameters for inverter
struct InvertParam_t
{
  InvType       invType;   // Inverter type
  Real          MROver;
  Real          RsdCG;
  int           MaxCG;	   // Iteration parameters
};

#endif
