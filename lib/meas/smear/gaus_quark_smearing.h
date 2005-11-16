// -*- C++ -*-
// $Id: gaus_quark_smearing.h,v 2.3 2005-11-16 02:34:58 edwards Exp $
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
  namespace GausQuarkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  

    //! Params for Gauss quark smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Real wvf_param;                   /*!< Smearing width */
      int  wvfIntPar;                   /*!< Number of smearing hits */
      int  no_smear_dir;		/*!< No smearing in this direction */
    };


    //! Gaussian quark smearing
    /*! @ingroup smear
     *
     * Gaussian quark smearing object
     */
    template<typename T>
    class QuarkSmear : public QuarkSmearing<T>
    {
    public:
      //! Full constructor
      QuarkSmear(const Params& p) : params(p) {}
      
      //! Smear the quark
      void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      QuarkSmear() {}

    private:
      Params  params;   /*!< smearing params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, GausQuarkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const GausQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif
