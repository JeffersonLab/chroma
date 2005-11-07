// -*- C++ -*-
// $Id: gaus_quark_smearing.h,v 2.1 2005-11-07 06:40:55 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color vector and propagator
 */

#ifndef __gaus_quark_smearing_h__
#define __gaus_quark_smearing_h__

#include "meas/smear/quark_smearing.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup smear */
  namespace GausPropSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  
  //! Name and registration
  /*! @ingroup smear */
  namespace GausFermSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  
  //! Name and registration
  /*! @ingroup smear */
  namespace GausColorVecSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Params for Gauss quark smearing
  /*! @ingroup smear */
  struct GausQuarkSmearingParams
  {
    GausQuarkSmearingParams() {}
    GausQuarkSmearingParams(XMLReader& in, const std::string& path);
    
    struct Param_t
    {
      Real wvf_param;                   /*!< Smearing width */
      int  wvfIntPar;                   /*!< Number of smearing hits */
      int  no_smear_dir;		/*!< No smearing in this direction */
    } param;
  };


  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, GausQuarkSmearingParams& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const GausQuarkSmearingParams& param);


  //! Gaussian quark smearing
  /*! @ingroup smear
   *
   * Gaussian quark smearing object
   */
  template<typename T>
  class GausQuarkSmearing : public QuarkSmearing<T>
  {
  public:
    //! Full constructor
    GausQuarkSmearing(const GausQuarkSmearingParams& p) : params(p) {}

    //! Smear the quark
    void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    GausQuarkSmearing() {}

  private:
    GausQuarkSmearingParams  params;   /*!< smearing params */
  };

}  // end namespace Chroma

#endif
