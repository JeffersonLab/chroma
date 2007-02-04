// -*- C++ -*-
// $Id: sf_gauge_init.h,v 3.1 2007-02-04 22:06:42 edwards Exp $
/*! \file
 *  \brief Initialize a Schroedinger BC config
 */

#ifndef __sf_gauge_init_h__
#define __sf_gauge_init_h__

#include "util/gauge/gauge_init.h"
#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma
{

  //! Name and registration
  namespace SFGaugeInitEnv
  {
    extern const std::string name;
    bool registerAll();
  

    //! Params for initializing config
    /*! @ingroup gauge */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      SchrGaugeBCParams sf;  	/*!< File name */
    };



    //! Gauge initialization
    /*! @ingroup gauge
     *
     * SF reader
     */
    class GaugeIniter : public GaugeInit
    {
    public:
      //! Full constructor
      GaugeIniter(const Params& p) : params(p) {}

      //! Initialize the gauge field
      void operator()(XMLReader& gauge_file_xml,
		      XMLReader& gauge_xml,
		      multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      GaugeIniter() {}

    private:
      Params  params;
    };

  }  // end namespace


  //! Reader
  /*! @ingroup gauge */
  void read(XMLReader& xml, const string& path, SFGaugeInitEnv::Params& param);

  //! Writer
  /*! @ingroup gauge */
  void write(XMLWriter& xml, const string& path, const SFGaugeInitEnv::Params& param);

}  // end namespace Chroma


#endif
