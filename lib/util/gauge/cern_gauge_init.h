// -*- C++ -*-
/*! \file
 *  \brief Read a CERN config
 */

#ifndef __cern_gauge_init_h__
#define __cern_gauge_init_h__

#include "util/gauge/gauge_init.h"

namespace Chroma
{

  //! Name and registration
  namespace CERNGaugeInitEnv
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
    
      std::string cfg_file;		/*!< File name */
    };



    //! Gauge initialization
    /*! @ingroup gauge
     *
     * CERN reader
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
  void read(XMLReader& xml, const std::string& path, CERNGaugeInitEnv::Params& param);

  //! Writer
  /*! @ingroup gauge */
  void write(XMLWriter& xml, const std::string& path, const CERNGaugeInitEnv::Params& param);

}  // end namespace Chroma


#endif
