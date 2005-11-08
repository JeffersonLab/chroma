// -*- C++ -*-
// $Id: pt_sink_smearing.h,v 1.3 2005-11-08 05:29:37 edwards Exp $
/*! \file
 *  \brief Point sink smearing
 */

#ifndef __pt_sink_smearing_h__
#define __pt_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Point sink parameters
  /*! @ingroup sinks */
  struct PointQuarkSinkSmearingParams
  {
    PointQuarkSinkSmearingParams();
    PointQuarkSinkSmearingParams(XMLReader& in, const std::string& path);

    std::string      link_smearing;        /*!< link smearing xml */
    std::string      link_smearing_type;   /*!< link smearing type name */

    int              disp_length;          /*!< displacement length */
    int              disp_dir;             /*!< x(0), y(1), z(2) */   
  };


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, PointQuarkSinkSmearingParams& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const PointQuarkSinkSmearingParams& param);



  //! Name and registration
  /*! @ingroup sinks */
  namespace PointQuarkSinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Point sink smearing
  /*! @ingroup sinks
   *
   * Create a point propagator sink
   */
  template<typename T>
  class PointQuarkSinkSmearing : public QuarkSourceSink<T>
  {
  public:
    //! Full constructor
    PointQuarkSinkSmearing(const PointQuarkSinkSmearingParams& p, const multi1d<LatticeColorMatrix>& u) :
      params(p), u_smr(u) {create();}

    //! Smear the sink
    void operator()(T& obj) const;

  private:
    //! Hide partial constructor
    PointQuarkSinkSmearing() {}

    //! Creator
    void create();

  private:
    PointQuarkSinkSmearingParams  params;   /*!< sink params */
    multi1d<LatticeColorMatrix>   u_smr;    /*!< hold a smeared copy for efficiency */
  };

}  // end namespace Chroma


#endif
