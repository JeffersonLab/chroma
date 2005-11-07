// -*- C++ -*-
// $Id: ape_link_smearing.h,v 1.2 2005-11-07 06:40:55 edwards Exp $
/*! \file
 *  \brief APE link smearing
 */

#ifndef __ape_link_smearing_h__
#define __ape_link_smearing_h__

#include "meas/smear/link_smearing.h"

namespace Chroma
{

  //! Name and registration
  namespace APELinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
    //! Name to be used
  }
  
  //! Params for APE link smearing
  /*! @ingroup smear */
  struct APELinkSmearingParams
  {
    APELinkSmearingParams() {}
    APELinkSmearingParams(XMLReader& in, const std::string& path);
    
    struct Param_t
    {
      int link_smear_num;               /*!< Smearing hits */
      Real link_smear_fact;		/*!< Smearing parameters */
      int no_smear_dir;			/*!< Direction to not smear */
    } param;
  };


  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, APELinkSmearingParams& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const APELinkSmearingParams& param);


  //! APE link smearing
  /*! @ingroup smear
   *
   * APE link smearing object
   */
  class APELinkSmearing : public LinkSmearing
  {
  public:
    //! Full constructor
    APELinkSmearing(const APELinkSmearingParams& p) : params(p) {}

    //! Smear the links
    void operator()(multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    APELinkSmearing() {}

  private:
    APELinkSmearingParams  params;   /*!< smearing params */
  };

}  // end namespace Chroma


#endif
