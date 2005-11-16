// -*- C++ -*-
// $Id: sh_sink_smearing.h,v 1.6 2005-11-16 02:34:58 edwards Exp $
/*! \file
 *  \brief Shell sink smearing
 */

#ifndef __sh_sink_smearing_h__
#define __sh_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sinks */
  namespace ShellQuarkSinkSmearingEnv
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
    
      std::string      sink_type;            /*!< sink smearing type */

      std::string      quark_smearing;       /*!< xml string holding smearing params */
      std::string      quark_smearing_type;  /*!< quark smearing type name */

      int              disp_length;          /*!< displacement length */
      int              disp_dir;             /*!< x(0), y(1), z(2) */

      std::string      link_smearing;        /*!< link smearing xml */
      std::string      link_smearing_type;   /*!< link smearing type name */
    };



    //! Shell sink smearing
    /*! @ingroup sinks
     *
     * Sink smeared quark
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


  }  // end namespace


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, ShellQuarkSinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSinkSmearingEnv::Params& param);


}  // end namespace Chroma


#endif
