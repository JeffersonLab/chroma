// -*- C++ -*-
// $Id: no_link_smearing.h,v 3.1 2006-05-19 21:34:06 edwards Exp $
/*! \file
 *  \brief No link smearing
 */

#ifndef __no_link_smearing_h__
#define __no_link_smearing_h__

#include "meas/smear/link_smearing.h"

namespace Chroma
{

  //! Name and registration
  namespace NoLinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;

  
    //! Params for No link smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    };


    //! No link smearing
    /*! @ingroup smear
     *
     * No link smearing object
     */
    class LinkSmear : public LinkSmearing
    {
    public:
      //! Full constructor
      LinkSmear(const Params& p) : params(p) {}

      //! Smear the links
      void operator()(multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      LinkSmear() {}

    private:
      Params  params;   /*!< smearing params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, NoLinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const NoLinkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
