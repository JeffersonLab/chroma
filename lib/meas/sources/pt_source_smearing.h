// -*- C++ -*-
// $Id: pt_source_smearing.h,v 3.5 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#ifndef __pt_source_smearing_h__
#define __pt_source_smearing_h__

#include "meas/smear/quark_source_sink.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace PointQuarkSourceSmearingEnv
  {
    bool registerAll();
  
    //! Return the name
    std::string getName();

    //! Point sink parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      GroupXML_t       quark_displacement;   /*!< displacement xml */
      GroupXML_t       link_smearing;        /*!< link smearing xml */

      int              j_decay;              /*!< decay direction */
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
	  this->create(u_smr, params.link_smearing);
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
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, PointQuarkSourceSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
