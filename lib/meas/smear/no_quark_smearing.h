// -*- C++ -*-
// $Id: no_quark_smearing.h,v 3.3 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief No quark smearing
 */

#ifndef __no_quark_smearing_h__
#define __no_quark_smearing_h__

#include "meas/smear/quark_smearing.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup smear */
  namespace NoQuarkSmearingEnv
  {
    bool registerAll();
  
    //! Return the name
    std::string getName();

    //! Params for No quark smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    };


    //! No quark smearing
    /*! @ingroup smear
     *
     * No quark smearing object
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
  void read(XMLReader& xml, const string& path, NoQuarkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const NoQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif
