// -*- C++ -*-
// $Id: sh_sink_smearing.h,v 3.6 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Shell sink smearing
 */

#ifndef __sh_sink_smearing_h__
#define __sh_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sinks */
  namespace ShellQuarkSinkSmearingEnv
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
    
      bool             quark_smear_firstP;   /*!< Flag controlling order of smearing */

      GroupXML_t       quark_smearing;       /*!< xml string holding smearing params */
      GroupXML_t       quark_displacement;   /*!< displacement xml */
      GroupXML_t       link_smearing;        /*!< link smearing xml */
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


  }  // end namespace


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, ShellQuarkSinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSinkSmearingEnv::Params& param);


}  // end namespace Chroma


#endif
