// -*- C++ -*-
// $Id: pt_source_smearing.h,v 2.6 2005-11-16 02:34:58 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#ifndef __pt_source_smearing_h__
#define __pt_source_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace PointQuarkSourceSmearingEnv
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

      int              disp_length;          /*!< displacement length */
      int              disp_dir;             /*!< x(0), y(1), z(2) */   
    };


    //! Point source smearing
    /*! @ingroup sources
     *
     * Point source smearing. Really not a smearing; however, there can be
     * displacements which use smeared links
     */
    template<typename T>
    class SourceSmear : public QuarkSourceSink<T>
    {
    public:
      //! Full constructor
      SourceSmear(const Params& p, const multi1d<LatticeColorMatrix>& u) :
	params(p), u_smr(u) 
	{
	  this->create(u_smr, params.link_smearing, params.link_smearing_type);
	}

      //! Construct the source
      void operator()(T& obj) const;

    private:
      //! Hide partial constructor
      SourceSmear() {}

    private:
      Params  params;                           /*!< source params */
      multi1d<LatticeColorMatrix>     u_smr;    /*!< hold a smeared copy for efficiency */
    };


  }  // end namespace Chroma


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, PointQuarkSourceSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
