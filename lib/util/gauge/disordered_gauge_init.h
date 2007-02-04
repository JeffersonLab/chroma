// -*- C++ -*-
// $Id: disordered_gauge_init.h,v 3.1 2007-02-04 22:06:42 edwards Exp $
/*! \file
 *  \brief Create a disordered config
 */

#ifndef __disordered_gauge_init_h__
#define __disordered_gauge_init_h__

#include "util/gauge/gauge_init.h"

namespace Chroma
{

  //! Name and registration
  namespace DisorderedGaugeInitEnv
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
    };



    //! Gauge initialization
    /*! @ingroup gauge
     *
     * Disordered reader
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

}  // end namespace Chroma


#endif
