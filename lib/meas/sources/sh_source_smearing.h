// -*- C++ -*-
// $Id: sh_source_smearing.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Shell source smearing
 */

#ifndef __sh_source_smearing_h__
#define __sh_source_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace ShellQuarkSourceSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  

    //! Point source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      std::string      source_type;          /*!< source smearing type */

      std::string      quark_smearing;       /*!< xml string holding smearing params */
      std::string      quark_smearing_type;  /*!< quark smearing type name */

      std::string      quark_displacement;      /*!< displacement xml */
      std::string      quark_displacement_type; /*!< displacement type name */

      std::string      link_smearing;        /*!< link smearing xml */
      std::string      link_smearing_type;   /*!< link smearing type name */
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
	  this->create(u_smr, params.link_smearing, params.link_smearing_type);
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
