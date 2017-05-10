// -*- C++ -*-
/*! \file
 *  \brief Stout link smearing
 */

#ifndef __phase_stout_link_smearing_h__
#define __phase_stout_link_smearing_h__

#include "meas/smear/link_smearing.h"
#include "meas/smear/stout_link_smearing.h"

namespace Chroma
{

  //! Name and registration
  namespace PhaseStoutLinkSmearingEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();

    //! Params for Stout link smearing
    /*! @ingroup smear */
    class Params : public StoutLinkSmearingEnv::Params
    {
    public:
      Params(): StoutLinkSmearingEnv::Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      multi1d<int> k;
      Real zeta ;
    };



    //! Stout link smearing
    /*! @ingroup smear
     *
     * Stout link smearing object
     */
    class LinkSmear : public  LinkSmearing
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
  void read(XMLReader& xml, const std::string& path, PhaseStoutLinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const std::string& path, const PhaseStoutLinkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
