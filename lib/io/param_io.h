// -*- C++ -*-
// $Id: param_io.h,v 2.0 2005-09-25 21:04:32 edwards Exp $
/*! \file
 *  \brief Various parameter structs and reader/writers
 */

#ifndef __param_io_h__
#define __param_io_h__
#include "invtype.h"

//! reading enums
#include "io/enum_io/enum_io.h"

// ! now needed

namespace Chroma 
{

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */

  //! Convert a Kappa to a mass
  Real kappaToMass(const Real& Kappa);

  //! Convert a Kappa to a mass
  multi1d<Real> kappaToMass(const multi1d<Real>& Kappa);

  //! Convert a Kappa to a mass
  Real massToKappa(const Real& Mass);

  //! Convert a mass to a Kappa
  multi1d<Real> massToKappa(const multi1d<Real>& Mass);



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
    ChiralParam_t();  // default constructor
    ~ChiralParam_t() {}

    Real       OverMass;
    int        N5;
    Real       a5;
    int        NWilsVec;
  };


  //---------------------------- Readers -----------------------------
  //! Configuration input
  void read(XMLReader& xml, const string& path, Cfg_t& input);

  //! Read chiral action like parameters
  void read(XMLReader& xml, const string& path, ChiralParam_t& param);

  //! Read inverter parameters
  void read(XMLReader& xml, const string& path, InvertParam_t& param);

  //! Read inverter parameters
  void read(XMLReader& xml, const string& path, MultiInvertParam_t& param);


  //---------------------------- Writers -----------------------------
  //! Configuration input
  void write(XMLWriter& xml, const string& path, const Cfg_t& input);

  //! Write chiral action like parameters
  void write(XMLWriter& xml, const string& path, const ChiralParam_t& param);

  //! Write inverter parameters
  void write(XMLWriter& xml, const string& path, const InvertParam_t& param);

  //! Write inverter parameters
  void write(XMLWriter& xml, const string& path, const MultiInvertParam_t& param);

  /*! @} */  // end of group io

} //end namespace chroma

#endif
