// -*- C++ -*-
// $Id: param_io.h,v 1.27 2005-01-14 20:13:07 edwards Exp $
/*! \file
 *  \brief Reunitarize (to a SU(N)) inplace the matrix A under some option
 */

#ifndef __param_io_h__
#define __param_io_h__
#include "invtype.h"

//! reading enums
#include "io/enum_io/enum_io.h"

// ! now needed

namespace Chroma { 

//! Convert a Kappa to a mass
Real kappaToMass(const Real& Kappa);

//! Convert a Kappa to a mass
multi1d<Real> kappaToMass(const multi1d<Real>& Kappa);

//! Convert a Kappa to a mass
Real massToKappa(const Real& Mass);

//! Convert a mass to a Kappa
multi1d<Real> massToKappa(const multi1d<Real>& Mass);


/*!
 * Types and structures
 *
 * \ingroup io
 *
 * @{
 */


/*
 * Input 
 */
struct IO_version_t
{
  int version;
};


struct Cfg_t
{
  CfgType      cfg_type;   // storage order for stored gauge configuration
  string       cfg_file;
};


//! Parameters for chiral fermion actions
struct ChiralParam_t
{
  Real       OverMass;
  int        N5;
  Real       a5;
  int        NWilsVec;
};


//! Parameters for anisotropy
struct AnisoParam_t
{
  bool       anisoP;
  int        t_dir;
  Real       xi_0;
  Real       nu;
};


//! Parameters for sources and sinks
struct SmearingParam_t
{
  WvfKind       wvf_kind;
  Real          wvf_param;
  int           wvfIntPar;
};


//---------------------------- Initializers -----------------------------
//! Initialize a smearing param struct
void initHeader(SmearingParam_t& param);

//! Initialize a anisotropy param struct
void initHeader(AnisoParam_t& param);

//! Initialize a chiral param struct
void initHeader(ChiralParam_t& param);


//---------------------------- Readers -----------------------------
//! Configuration input
void read(XMLReader& xml, const string& path, Cfg_t& input);

//! Read a anisotropy param struct
void read(XMLReader& xml, const string& path, AnisoParam_t& param);

//! Read a smearing param struct
void read(XMLReader& xml, const string& path, SmearingParam_t& param);

//! Read chiral action like parameters
void read(XMLReader& xml, const string& path, ChiralParam_t& param);

//! Read inverter parameters
void read(XMLReader& xml, const string& path, InvertParam_t& param);

//! Read inverter parameters
void read(XMLReader& xml, const string& path, MultiInvertParam_t& param);


//---------------------------- Writers -----------------------------
//! Configuration input
void write(XMLWriter& xml, const string& path, const Cfg_t& input);

//! Write a anisotropy param struct
void write(XMLWriter& xml, const string& path, const AnisoParam_t& param);

//! Write a smearing param struct
void write(XMLWriter& xml, const string& path, const SmearingParam_t& param);

//! Write chiral action like parameters
void write(XMLWriter& xml, const string& path, const ChiralParam_t& param);

//! Write inverter parameters
void write(XMLWriter& xml, const string& path, const InvertParam_t& param);

//! Write inverter parameters
void write(XMLWriter& xml, const string& path, const MultiInvertParam_t& param);
/*! @} */  // end of group io

}; //end namespace chroma
#endif
