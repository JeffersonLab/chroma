// -*- C++ -*-
// $Id: delta_2pt_w.h,v 3.3 2008-05-22 18:50:00 caubin Exp $
/*! \file
 *  \brief Construct delta 2pt correlators.
 */

#ifndef __delta_2pt_w_h__
#define __delta_2pt_w_h__

#include "meas/hadron/hadron_2pt.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace Delta2PtEnv
  {
    bool registerAll();

  
    //! Simple meson 2pt parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      int           mom2_max;        /*!< (mom - mom_origin)^2 <= mom2_max */
      multi1d<int>  mom_origin;      /*!< Origin for the momentum */
      bool          avg_equiv_mom;   /*!< average over equivalent momenta */
	  
	  bool          min_contractions;/*!< Only do minimal number of contractions*/
	  std::string   parity;          /*!< Which parity/parities to do?*/
      std::string   first_id;        /*!< First/light quark id */
      std::string   second_id;       /*!< Second/heavy quark id */
    };


    //! Decuplet baryon 2pt construction 
    /*! @ingroup hadron
     *
     * Create all possible Decuplet polarization states
     */
    class DeltaCorrs : public Hadron2PtCorr
    {
    public:
      //! Full constructor
      DeltaCorrs(const Params& p) : params(p) {}

      //! Default destructor
      ~DeltaCorrs() {}
      
      //! Construct the correlators
      std::list< Handle<HadronContractResult_t> > operator()(const multi1d<LatticeColorMatrix>& u,
							     const std::string& xml_group,
							     const std::string& id_tag);

    private:
      //! Hide partial constructor
      DeltaCorrs() {}

    private:
      Params   params;     /*!< The common params */
    };
  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, Delta2PtEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const Delta2PtEnv::Params& param);


}  // end namespace Chroma

#endif
