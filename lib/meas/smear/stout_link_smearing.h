// -*- C++ -*-
// $Id: stout_link_smearing.h,v 3.4 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Stout link smearing
 */

#ifndef __stout_link_smearing_h__
#define __stout_link_smearing_h__

#include "meas/smear/link_smearing.h"

namespace Chroma
{

  //! Name and registration
  namespace StoutLinkSmearingEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();

    //! Params for Stout link smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Real           link_smear_fact;    /*!< Smearing parameters */
      int            link_smear_num;     /*!< Number of smearing hits */
      multi1d<bool>  smear_dirs;         /*!< Only allow smearing and staples in these directions */

      multi2d<Real>  rho;                /*!< Parameters actually used by stout::smear_links */
    };



    //! Stout link smearing
    /*! @ingroup smear
     *
     * Stout link smearing object
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
  void read(XMLReader& xml, const string& path, StoutLinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const StoutLinkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
