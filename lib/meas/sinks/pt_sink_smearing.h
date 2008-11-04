// -*- C++ -*-
// $Id: pt_sink_smearing.h,v 3.4 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Point sink smearing
 */

#ifndef __pt_sink_smearing_h__
#define __pt_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sinks */
  namespace PointQuarkSinkSmearingEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();


    //! Point sink parameters
    /*! @ingroup sinks */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      GroupXML_t       quark_displacement;   /*!< displacement xml */
      GroupXML_t       link_smearing;        /*!< link smearing xml */
    };


    //! Point sink smearing
    /*! @ingroup sinks
     *
     * Create a point propagator sink
     */
    template<typename T>
    class SinkSmear : public QuarkSourceSink<T>
    {
    public:
      //! Full constructor
      SinkSmear(const Params& p, const multi1d<LatticeColorMatrix>& u) :
	params(p), u_smr(u) 
	{
	  this->create(u_smr, params.link_smearing);
	}

      //! Smear the sink
      void operator()(T& obj) const;

    private:
      //! Hide partial constructor
      SinkSmear() {}

    private:
      Params  params;                         /*!< sink params */
      multi1d<LatticeColorMatrix>   u_smr;    /*!< hold a smeared copy for efficiency */
    };

  } // end namespace


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, PointQuarkSinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const PointQuarkSinkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
