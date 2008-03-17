// -*- C++ -*-
// $Id: meson_spec_2pt_w.h,v 1.1 2008-03-17 15:23:58 edwards Exp $
/*! \file
 *  \brief Construct meson 2pt correlators leaving all spin indices open
 */

#ifndef __meson_spec_2pt_w_h__
#define __meson_spec_2pt_w_h__

#error "STILL WORKING ON THIS"

#include "meas/hadron/hadron_2pt.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace MesonSpec2PtEnv
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


    //! 2pt-mesons but with all 4 spin indices open
    /*! @ingroup hadron
     *
     * Contract mesons but leave source and sink spin indices open
     */
    class MesonSpecCorrs : public Hadron2PtCorr
    {
    public:
      //! Full constructor
      MesonSpecCorrs(const Params& p) : params(p) {}

      //! Default destructor
      ~MesonSpecCorrs() {}
      
      //! Construct the correlators
      std::list< Handle<HadronContractResult_t> > operator()(const multi1d<LatticeColorMatrix>& u,
							     const std::string& xml_group,
							     const std::string& id_tag);

    protected:
    private:
      //! Hide partial constructor
      MesonSpecCorrs() {}

    private:
      Params   params;     /*!< The common params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, MesonSpec2PtEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const MesonSpec2PtEnv::Params& param);


}  // end namespace Chroma

#endif
