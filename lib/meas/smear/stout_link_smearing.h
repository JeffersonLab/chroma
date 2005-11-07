// -*- C++ -*-
// $Id: stout_link_smearing.h,v 1.1 2005-11-07 18:05:42 edwards Exp $
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
    extern const std::string name;
    extern const bool registered;
    //! Name to be used
  }
  
  //! Params for Stout link smearing
  /*! @ingroup smear */
  struct StoutLinkSmearingParams
  {
    StoutLinkSmearingParams() {}
    StoutLinkSmearingParams(XMLReader& in, const std::string& path);
    
    struct Param_t
    {
      Real link_smear_fact;		/*!< Smearing parameters */
      int  link_smear_num;              /*!< Number of smearing hits */
      int  no_smear_dir;		/*!< Decay direction */
    } param;

  };


  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, StoutLinkSmearingParams& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const StoutLinkSmearingParams& param);


  //! Stout link smearing
  /*! @ingroup smear
   *
   * Stout link smearing object
   */
  class StoutLinkSmearing : public LinkSmearing
  {
  public:
    //! Full constructor
    StoutLinkSmearing(const StoutLinkSmearingParams& p) : params(p) {}

    //! Smear the links
    void operator()(multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    StoutLinkSmearing() {}

  private:
    StoutLinkSmearingParams  params;   /*!< smearing params */
  };

}  // end namespace Chroma


#endif
