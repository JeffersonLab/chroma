// -*- C++ -*-
// $Id: hyp_link_smearing.h,v 1.1 2005-11-07 18:05:42 edwards Exp $
/*! \file
 *  \brief HYP link smearing
 */

#ifndef __hyp_link_smearing_h__
#define __hyp_link_smearing_h__

#include "meas/smear/link_smearing.h"

namespace Chroma
{

  //! Name and registration
  namespace HypLinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  
  //! Params for Hyp link smearing
  /*! @ingroup smear */
  struct HypLinkSmearingParams
  {
    HypLinkSmearingParams() {}
    HypLinkSmearingParams(XMLReader& in, const std::string& path);
    
    struct Param_t
    {
      Real alpha1;	                /*!< staple coefficient "1" */
      Real alpha2;	                /*!< staple coefficient "2" */
      Real alpha3;	                /*!< staple coefficient "3" */
      int num_smear;			/*!< Number of smearing hits */
      int no_smear_dir;			/*!< Direction to not smear */
    } param;
  };


  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, HypLinkSmearingParams& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const HypLinkSmearingParams& param);


  //! Hyp link smearing
  /*! @ingroup smear
   *
   * Hyp link smearing object
   */
  class HypLinkSmearing : public LinkSmearing
  {
  public:
    //! Full constructor
    HypLinkSmearing(const HypLinkSmearingParams& p) : params(p) {}

    //! Smear the links
    void operator()(multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    HypLinkSmearing() {}

  private:
    HypLinkSmearingParams  params;   /*!< smearing params */
  };

}  // end namespace Chroma


#endif
