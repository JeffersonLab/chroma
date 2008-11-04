// -*- C++ -*-
// $Id: hyp_link_smearing.h,v 3.3 2008-11-04 18:43:57 edwards Exp $
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
    bool registerAll();

    //! Return the name
    std::string getName();
  
    //! Params for Hyp link smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Real alpha1;	                /*!< staple coefficient "1" */
      Real alpha2;	                /*!< staple coefficient "2" */
      Real alpha3;	                /*!< staple coefficient "3" */
      int num_smear;			/*!< Number of smearing hits */
      int no_smear_dir;			/*!< Direction to not smear */
      int BlkMax;                       /*!< Max number of iterations */
      Real BlkAccu;                     /*!< Relative error to maximize trace */
    };


    //! Hyp link smearing
    /*! @ingroup smear
     *
     * Hyp link smearing object
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
  void read(XMLReader& xml, const string& path, HypLinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const HypLinkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
