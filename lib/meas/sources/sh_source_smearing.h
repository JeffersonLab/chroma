// -*- C++ -*-
// $Id: sh_source_smearing.h,v 3.7 2008-11-04 18:43:59 edwards Exp $
/*! \file
 *  \brief Shell source smearing
 */

#ifndef __sh_source_smearing_h__
#define __sh_source_smearing_h__

#include "meas/smear/quark_source_sink.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace ShellQuarkSourceSmearingEnv
  {
    bool registerAll();
  
    //! Return the name
    std::string getName();

    //! Point source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      bool             quark_smear_lastP;    /*!< Flag controlling order of smearing */

      GroupXML_t       quark_smearing;       /*!< xml string holding smearing params */
      GroupXML_t       quark_displacement;   /*!< displacement xml */
      GroupXML_t       link_smearing;        /*!< link smearing xml */

      int              j_decay;              /*!< Decay direction */
    };


    //! Shell source smearing
    /*! @ingroup sources
     *
     * A shell source smearing of a quark.
     */
    template<typename T>
    class SourceSmearing : public QuarkSourceSink<T>
    {
    public:
      //! Full constructor
      SourceSmearing(const Params& p, const multi1d<LatticeColorMatrix>& u) :
	params(p), u_smr(u) 
	{
	  this->create(u_smr, params.link_smearing);
	}

      //! Construct the source
      void operator()(T& obj) const;

    private:
      //! Hide partial constructor
      SourceSmearing() {}

    private:
      Params  params;                           /*!< source params */
      multi1d<LatticeColorMatrix>     u_smr;    /*!< hold a smeared copy for efficiency */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, ShellQuarkSourceSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
