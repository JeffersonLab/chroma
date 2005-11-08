// -*- C++ -*-
// $Id: sh_source_smearing.h,v 2.3 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Shell source smearing
 */

#ifndef __sh_source_smearing_h__
#define __sh_source_smearing_h__

#include "meas/smear/quark_source_sink.h"

namespace Chroma
{

  //! Point source parameters
  /*! @ingroup sources */
  struct ShellQuarkSourceSmearingParams
  {
    ShellQuarkSourceSmearingParams();
    ShellQuarkSourceSmearingParams(XMLReader& in, const std::string& path);
    
    std::string      source_type;          /*!< source smearing type */

    std::string      quark_smearing;       /*!< xml string holding smearing params */
    std::string      quark_smearing_type;  /*!< quark smearing type name */

    int              disp_length;          /*!< displacement length */
    int              disp_dir;             /*!< x(0), y(1), z(2) */

    std::string      link_smearing;        /*!< link smearing xml */
    std::string      link_smearing_type;   /*!< link smearing type name */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, ShellQuarkSourceSmearingParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceSmearingParams& param);



  //! Name and registration
  /*! @ingroup sources */
  namespace ShellQuarkSourceSmearingEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Shell source smearing
  /*! @ingroup sources
   *
   * A shell source smearing of a quark.
   */
  template<typename T>
  class ShellQuarkSourceSmearing : public QuarkSourceSink<T>
  {
  public:
    //! Full constructor
    ShellQuarkSourceSmearing(const ShellQuarkSourceSmearingParams& p, const multi1d<LatticeColorMatrix>& u) :
      params(p), u_smr(u) {create();}

    //! Construct the source
    void operator()(T& obj) const;

  private:
    //! Hide partial constructor
    ShellQuarkSourceSmearing() {}

    //! Creator
    void create();

  private:
    ShellQuarkSourceSmearingParams  params;   /*!< source params */
    multi1d<LatticeColorMatrix>     u_smr;    /*!< hold a smeared copy for efficiency */
  };

}  // end namespace Chroma


#endif
