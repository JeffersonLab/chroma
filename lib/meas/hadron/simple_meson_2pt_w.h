// -*- C++ -*-
// $Id: simple_meson_2pt_w.h,v 1.3 2007-06-10 14:49:06 edwards Exp $
/*! \file
 *  \brief Construct meson 2pt correlators.
 */

#ifndef __simple_meson_2pt_w_h__
#define __simple_meson_2pt_w_h__

#include "meas/hadron/hadron_2pt.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleMeson2PtEnv
  {
    bool registerAll();

  
    //! Simple meson 2pt parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int           mom2_max;           /*!< (mom - mom_origin)^2 <= mom2_max */
      multi1d<int>  mom_origin;         /*!< Origin for the momentum */
      bool          avg_equiv_mom;      /*!< average over equivalent momenta */

      std::string   first_id;        /*!< First/light quark id */
      std::string   second_id;       /*!< Second/heavy quark id */
    };


    //! Simple meson 2pt construction - all simple mesons
    /*! @ingroup hadron
     *
     * Create all the 16 (Ns=4) simple diagonal meson 2pt correlators
     */
    class DiagGammaMesonCorrs : public Hadron2PtCorr
    {
    public:
      //! Full constructor
      DiagGammaMesonCorrs(const Params& p) : params(p) {}

      //! Default destructor
      ~DiagGammaMesonCorrs() {}
      
      //! Construct the correlators
      std::list< Handle<HadronContractResult_t> > operator()(const multi1d<LatticeColorMatrix>& u,
							     const std::string& xml_group,
							     const std::string& id_tag);

    private:
      //! Hide partial constructor
      DiagGammaMesonCorrs() {}

    private:
      Params   params;     /*!< The common params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleMeson2PtEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleMeson2PtEnv::Params& param);


}  // end namespace Chroma

#endif
