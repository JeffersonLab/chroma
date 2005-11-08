// -*- C++ -*-
// $Id: pt_source_smearing.h,v 2.4 2005-11-08 18:41:53 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#ifndef __pt_source_smearing_h__
#define __pt_source_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Point sink parameters
  /*! @ingroup sinks */
  struct PointQuarkSourceSmearingParams
  {
    PointQuarkSourceSmearingParams();
    PointQuarkSourceSmearingParams(XMLReader& in, const std::string& path);

    std::string      link_smearing;        /*!< link smearing xml */
    std::string      link_smearing_type;   /*!< link smearing type name */

    int              disp_length;          /*!< displacement length */
    int              disp_dir;             /*!< x(0), y(1), z(2) */   
  };


  //! Reader
  /*! @ingroup sinks */
  void read(XMLReader& xml, const string& path, PointQuarkSourceSmearingParams& param);

  //! Writer
  /*! @ingroup sinks */
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceSmearingParams& param);



  //! Name and registration
  /*! @ingroup sources */
  namespace PointQuarkSourceSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Point source smearing
  /*! @ingroup sources
   *
   * Point source smearing. Really not a smearing; however, there can be
   * displacements which use smeared links
   */
  template<typename T>
  class PointQuarkSourceSmearing : public QuarkSourceSink<T>
  {
  public:
    //! Full constructor
    PointQuarkSourceSmearing(const PointQuarkSourceSmearingParams& p, 
			     const multi1d<LatticeColorMatrix>& u) :
      params(p), u_smr(u) 
      {
	create(u_smr, params.link_smearing, params.link_smearing_type);
      }

    //! Construct the source
    void operator()(T& obj) const;

  private:
    //! Hide partial constructor
    PointQuarkSourceSmearing() {}

    //! Potentially smear the gauge field
    /*!
     * \param u                   Gauge field to smear ( Modify )
     * \param link_smearing       XML of link smearing ( Read )
     * \param link_smearing_type  link smearing type ( Read )
     */
    void create(multi1d<LatticeColorMatrix>& u,
		std::string link_smearing,
		std::string link_smearing_type);

  private:
    PointQuarkSourceSmearingParams  params;   /*!< source params */
    multi1d<LatticeColorMatrix>     u_smr;    /*!< hold a smeared copy for efficiency */
  };

}  // end namespace Chroma


#endif
