// -*- C++ -*-
// $Id: param_io.h,v 1.20 2004-04-26 11:19:13 bjoo Exp $
/*! \file
 *  \brief Reunitarize (to a SU(N)) inplace the matrix A under some option
 */

#ifndef __param_io_h__
#define __param_io_h__

#include "invtype.h"
#include "meas/smear/wvfkind.h"
#include "meas/sources/srcsnktype.h"
#include "meas/sources/wavetype.h"

//! Convert a Kappa to a mass
Real kappaToMass(const Real& Kappa);

//! Convert a Kappa to a mass
multi1d<Real> kappaToMass(const multi1d<Real>& Kappa);

//! Convert a Kappa to a mass
Real massToKappa(const Real& Mass);

//! Convert a mass to a Kappa
multi1d<Real> massToKappa(const multi1d<Real>& Mass);


/*
 * Types and structures
 *
 * \ingroup io
 *
 * @{
 */

//! Configuration type
enum CfgType 
{
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_SZINQIO,
  CFG_TYPE_DISORDERED,
  CFG_TYPE_UNIT, 
};

//! Propagator type
enum PropType 
{
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
  PROP_TYPE_KYU,
};

//! Eigenvector type
enum EigenVecType {
  EVEC_TYPE_SCIDAC = 2,
  EVEC_TYPE_SZIN
};

//! Types of fermion
enum FermType 
{
  FERM_TYPE_WILSON,
  FERM_TYPE_STAGGERED,
};


//! Types of fermion actions
enum FermActType 
{
  FERM_ACT_WILSON,
  FERM_ACT_UNPRECONDITIONED_WILSON,
  FERM_ACT_PARITY_BREAKING_WILSON,
  FERM_ACT_CLOVER,
  FERM_ACT_UNPRECONDITIONED_CLOVER,
  FERM_ACT_DWF,                           // Precond. Shamir DWF
  FERM_ACT_UNPRECONDITIONED_DWF,          // Unprec. Shamir DWF
  FERM_ACT_PROJECTED_DWF,                 // Shamir precond. DWF with E&H projection
  FERM_ACT_ZOLOTAREV_4D,                  // Overlap pole with Zolotarev coeffs
  FERM_ACT_ZOLOTAREV_5D,                  // 5D overlap op. ZOlotarev coeffs
  FERM_ACT_OVERLAP_DWF,                   // Borici
  FERM_ACT_EXTENDED_OVERLAP,              // Unprecond. N&N 5D overlap
  FERM_ACT_UNPRECONDITIONED_EXTENDED_OVERLAP,  // Precond. N&N 5D overlap
  FERM_ACT_SMEARED_LAPLACIAN_WILSON,
  FERM_ACT_PLANAR_WILSON,
  FERM_ACT_HAMBER_WU,
  FERM_ACT_STAGGERED,
  FERM_ACT_NAIK,
  FERM_ACT_ASQTAD
};


//! Sequential source type
enum SeqSourceType 
{
  SEQ_SOURCE_TYPE_NUCL_U_UNPOL  = 0,
  SEQ_SOURCE_TYPE_NUCL_D_UNPOL  = 1,
  SEQ_SOURCE_TYPE_NUCL_U_POL    = 2,   // same as polarized in 2pt baryon case
  SEQ_SOURCE_TYPE_NUCL_D_POL    = 3,   // same as polarized in 2pt baryon case
  SEQ_SOURCE_TYPE_DELTA_U_UNPOL = 4,
  SEQ_SOURCE_TYPE_DELTA_D_UNPOL = 5,
  SEQ_SOURCE_TYPE_NUCL_U_UNPOL_NONREL = 6,
  SEQ_SOURCE_TYPE_NUCL_D_UNPOL_NONREL = 7,
  SEQ_SOURCE_TYPE_NUCL_U_POL_NONREL   = 8,
  SEQ_SOURCE_TYPE_NUCL_D_POL_NONREL   = 9,
  SEQ_SOURCE_TYPE_NUCL_U_MIXED_NONREL = 21,
  SEQ_SOURCE_TYPE_NUCL_D_MIXED_NONREL = 22,

  SEQ_SOURCE_TYPE_PION = 10,
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


//! Parameters for inverter
struct InvertParam_t
{
  InvType       invType;   // Inverter type
  Real          MROver;
  Real          RsdCG;
  int           MaxCG;	   // Iteration parameters
};

//! Parameters for MultiMass Inverter
struct MultiInvertParam_t
{
  InvType      invType;
  Real         MROver;
  int          MaxCG;
  multi1d<Real> RsdCG;
};

//---------------------------- Initializers -----------------------------
//! Initialize a smearing param struct
void initHeader(SmearingParam_t& param);

//! Initialize a anisotropy param struct
void initHeader(AnisoParam_t& param);

//! Initialize a chiral param struct
void initHeader(ChiralParam_t& param);


//---------------------------- Readers -----------------------------
//! Read a QDP volume format type
void read(XMLReader& xml, const string& path, QDP_volfmt_t& param);

//! Read a fermion type enum
void read(XMLReader& xml, const string& path, FermType& param);

//! Read a fermion action type enum
void read(XMLReader& xml, const string& path, FermActType& param);

//! Read a configuration type enum
void read(XMLReader& xml, const string& path, CfgType& param);

//! Read a propagator type enum
void read(XMLReader& xml, const string& path, PropType& param);

//! Read an Eigenvector type enum
void read(XMLReader& xml, const string& path, EigenVecType& param);


//! Read a wave-function type enum
void read(XMLReader& xml, const string& path, WvfKind& param);

//! Read a inverter type enum
void read(XMLReader& xml, const string& path, InvType& param);

//! Read a multi inverter type enum
//! Read a source type enum
void read(XMLReader& xml, const string& path, SourceType& param);

//! Read a sink type enum
void read(XMLReader& xml, const string& path, SinkType& param);

//! Read a wave type enum
void read(XMLReader& xml, const string& path, WaveStateType& param);

//! Read a sequential source type enum
void read(XMLReader& xml, const string& path, SeqSourceType& param);

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
//! Write a QDP volume format type
void write(XMLWriter& xml, const string& path, QDP_volfmt_t param);

//! Write a fermion type enum
void write(XMLWriter& xml, const string& path, FermType param);

//! Write a fermion action type enum
void write(XMLWriter& xml, const string& path, FermActType param);

//! Write a configuration type enum
void write(XMLWriter& xml, const string& path, CfgType param);

//! Write a propagator type enum
void write(XMLWriter& xml, const string& path, PropType param);

//! Write an eigenvector type enum
void write(XMLWriter& xml, const string& path, EigenVecType param);

//! Write a wave-function type enum
void write(XMLWriter& xml, const string& path, WvfKind param);

//! Write a inverter type enum
void write(XMLWriter& xml, const string& path, InvType param);

//! Write a source type enum
void write(XMLWriter& xml, const string& path, SourceType param);

//! Write a sink type enum
void write(XMLWriter& xml, const string& path, SinkType param);

//! Write a wave type enum
void write(XMLWriter& xml, const string& path, WaveStateType param);

//! Write a seqsource type enum
void write(XMLWriter& xml, const string& path, SeqSourceType param);

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

#endif
