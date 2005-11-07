// -*- C++ -*-
// $Id: sh_sink_smearing.h,v 1.1 2005-11-07 06:24:09 edwards Exp $
/*! \file
 *  \brief Shell sink smearing
 */

#ifndef __sh_sink_smearing_h__
#define __sh_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Point sink parameters
  /*! @ingroup sinks */
  struct ShellQuarkSinkSmearingParams
  {
    ShellQuarkSinkSmearingParams() {}
    ShellQuarkSinkSmearingParams(XMLReader& in, const std::string& path);
    
    std::string      sink_type;            /*!< sink smearing type */

    std::string      quark_smearing;       /*!< xml string holding smearing params */
    std::string      quark_smearing_type;  /*!< quark smearing type name */

    int              disp_length;          /*!< displacement length */
    int              disp_dir;             /*!< x(0), y(1), z(2) */

    std::string      link_smearing;        /*!< link smearing xml */
    std::string      link_smearing_type;   /*!< link smearing type name */
  };


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, ShellQuarkSinkSmearingParams& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSinkSmearingParams& param);



  //! Name and registration
  /*! @ingroup sinks */
  namespace ShellPropSinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Name and registration
  /*! @ingroup sinks */
  namespace ShellFermSinkSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Shell sink smearing
  /*! @ingroup sinks
   *
   * Sink smeared quark
   */
  template<typename T>
  class ShellQuarkSinkSmearing : public QuarkSourceSink<T>
  {
  public:
    //! Full constructor
    ShellQuarkSinkSmearing(const ShellQuarkSinkSmearingParams& p, const multi1d<LatticeColorMatrix>& u) :
      params(p), u_smr(u) {create();}

    //! Smear the sink
    void operator()(T& obj) const;

  private:
    //! Hide partial constructor
    ShellQuarkSinkSmearing() {}

    //! Creator
    void create();

  private:
    ShellQuarkSinkSmearingParams  params;   /*!< sink params */
    multi1d<LatticeColorMatrix>   u_smr;    /*!< hold a smeared copy for efficiency */
  };

}  // end namespace Chroma


#endif
