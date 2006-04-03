// -*- C++ -*-
// $Id: pt_sink_smearing.h,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Point sink smearing
 */

#ifndef __pt_sink_smearing_h__
#define __pt_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sinks */
  namespace PointQuarkSinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;


    //! Point sink parameters
    /*! @ingroup sinks */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      std::string      link_smearing;        /*!< link smearing xml */
      std::string      link_smearing_type;   /*!< link smearing type name */

      std::string      quark_displacement;      /*!< displacement xml */
      std::string      quark_displacement_type; /*!< displacement type name */
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
	  this->create(u_smr, params.link_smearing, params.link_smearing_type);
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
