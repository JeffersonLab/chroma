// -*- C++ -*-
// $Id: stout_link_smearing_old.h,v 3.1 2008-01-20 03:07:51 edwards Exp $
/*! \file
 *  \brief Stout link smearing using the old (non-gauge covariant) method
 */

#ifndef __stout_link_smearing_old_h__
#define __stout_link_smearing_old_h__

#include "meas/smear/link_smearing.h"

namespace Chroma
{

  //! Name and registration
  namespace StoutLinkSmearingOldEnv
  {
    extern const std::string name;
    bool registerAll();
  

    //! Params for Stout link smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Real link_smear_fact;		/*!< Smearing parameters */
      int  link_smear_num;              /*!< Number of smearing hits */
      multi1d<bool> smear_dirs;         /*!< Only allow smearing and staples in these directions */
    };



    //! Stout link smearing
    /*! @ingroup smear
     *
     * Stout link smearing using the old (non-gauge covariant) method 
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
  void read(XMLReader& xml, const string& path, StoutLinkSmearingOldEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const StoutLinkSmearingOldEnv::Params& param);

}  // end namespace Chroma


#endif
