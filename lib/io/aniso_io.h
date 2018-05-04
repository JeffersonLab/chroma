// -*- C++ -*-
/*! \file
 *  \brief Anisotropy parameters
 */

#ifndef __aniso_io_h__
#define __aniso_io_h__

#include "chromabase.h"

namespace Chroma 
{ 

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */

  //! Parameters for anisotropy
  struct AnisoParam_t
  {
    AnisoParam_t();  // default constructor
    ~AnisoParam_t() {}

    bool       anisoP;
    int        t_dir;
    Real       xi_0;
    Real       nu;
  };

  //! Read a anisotropy param struct
  void read(XMLReader& xml, const std::string& path, AnisoParam_t& param);

  //! Write a anisotropy param struct
  void write(XMLWriter& xml, const std::string& path, const AnisoParam_t& param);

  //! Make fermion coefficients
  multi1d<Real> makeFermCoeffs(const AnisoParam_t& aniso);

  /*! @} */  // end of group io

} //end namespace chroma
#endif
