// -*- C++ -*-
// $Id: stoch_cond_cont_w.h,v 3.1 2007-08-11 22:43:21 edwards Exp $
/*! \file
 * \brief Stoch quark condensates
 *
 * Stoch quark condensates via the hadron contraction interface.
 * This code is using the Hadron2PtCorr interface because the
 * condensates are projected onto a fixed 3-space momentum. Thus,
 * they are time-slice dependent and similar then to the 2-pt
 * correlators. However, the origin is not an issue like that in 
 * the 2pt corrs., so that part can be avoided.
 */

#ifndef __stoch_cond_cont_w_h__
#define __stoch_cond_cont_w_h__

#include "meas/hadron/hadron_2pt.h"

namespace Chroma 
{ 
  //! Stochastic quark condensates
  /*! \ingroup hadron */
  namespace StochCondContEnv 
  {
    bool registerAll();


    //! Parameter structure
    /*! \ingroup hadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;

      int           mom2_max;           /*!< (mom)^2 <= mom2_max */
      multi1d<int>  mom_origin;         /*!< Origin for the momentum */
      bool          avg_equiv_mom;      /*!< average over equivalent momenta */

      multi1d<std::string> soln_files;  /*!< Stochastic quark prop solutions */
    };


    //! Stochastic quark condensates
    /*! \ingroup hadron */
    class StochCondCont : public Hadron2PtCorr
    {
    public:
      //! Full constructor
      StochCondCont(const Params& p) : params(p) {}

      //! Default destructor
      ~StochCondCont() {}
      
      //! Construct the correlators
      std::list< Handle<HadronContractResult_t> > operator()(const multi1d<LatticeColorMatrix>& u,
							     const std::string& xml_group,
							     const std::string& id_tag);

    private:
      Params params;
    };

  } // namespace StochCondContEnv


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, StochCondContEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const StochCondContEnv::Params& param);


} // namespace Chroma

#endif
