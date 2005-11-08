// -*- C++ -*-
// $Id: sh_source_const.h,v 2.3 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#ifndef __sh_source_const_h__
#define __sh_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Point source parameters
  /*! @ingroup sources */
  struct ShellQuarkSourceConstParams
  {
    ShellQuarkSourceConstParams();
    ShellQuarkSourceConstParams(XMLReader& in, const std::string& path);
    
    std::string      source_type;          /*!< source smearing type */

    std::string      quark_smearing;       /*!< xml string holding smearing params */
    std::string      quark_smearing_type;  /*!< quark smearing type name */

    int              j_decay;              /*!< Decay direction */
    int              disp_length;          /*!< displacement length */
    int              disp_dir;             /*!< x(0), y(1), z(2) */
    multi1d<int>     t_srce;               /*!< source location */

    std::string      link_smearing;        /*!< link smearing xml */
    std::string      link_smearing_type;   /*!< link smearing type name */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, ShellQuarkSourceConstParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceConstParams& param);



  //! Name and registration
  /*! @ingroup sources */
  namespace ShellQuarkSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Shell source construction
  /*! @ingroup sources
   *
   * Create a shell quark source
   */
  template<typename T>
  class ShellQuarkSourceConst : public QuarkSourceConstruction<T>
  {
  public:
    //! Full constructor
    ShellQuarkSourceConst(const ShellQuarkSourceConstParams& p) : params(p) {}

    //! Construct the source
    T operator()(const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    ShellQuarkSourceConst() {}

  private:
    ShellQuarkSourceConstParams  params;   /*!< source params */
  };

}  // end namespace Chroma


#endif
