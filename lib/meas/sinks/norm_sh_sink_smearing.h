// -*- C++ -*-
// $Id: norm_sh_sink_smearing.h,v 3.1 2009-04-12 03:45:00 kostas Exp $
/*! \file
 *  \brief NormShell sink smearing
 */

#ifndef __norm_sh_sink_smearing_h__
#define __norm_sh_sink_smearing_h__

#include "meas/smear/quark_source_sink.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sinks */
  namespace NormShellQuarkSinkSmearingEnv
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
      int Nhits ;
      bool site_normalized ; // not implemented  
    };



    //! NormShell sink smearing
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
  void read(XMLReader& xml, const string& path, NormShellQuarkSinkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const NormShellQuarkSinkSmearingEnv::Params& param);


}  // end namespace Chroma


#endif
