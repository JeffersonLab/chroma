// -*- C++ -*-
// $Id: pt_source_const.h,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#ifndef __pt_source_const_h__
#define __pt_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Point source parameters
  /*! @ingroup sources */
  struct PointQuarkSourceConstParams
  {
    PointQuarkSourceConstParams() {}
    PointQuarkSourceConstParams(XMLReader& in, const std::string& path);
    
    std::string      link_smearing;        /*!< link smearing xml */
    std::string      link_smearing_type;   /*!< link smearing type name */

    int              disp_length;          /*!< displacement length */
    int              disp_dir;             /*!< x(0), y(1), z(2) */   

    int              j_decay;              /*!< decay direction */
    multi1d<int>     t_srce;               /*!< source location */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, PointQuarkSourceConstParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const PointQuarkSourceConstParams& param);


  //! Name and registration
  /*! @ingroup sources */
  namespace PointPropSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Name and registration
  /*! @ingroup sources */
  namespace PointFermSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

  
  //! Point source construction
  /*! @ingroup sources
   *
   * Create a point propagator source
   */
  template<typename T>
  class PointQuarkSourceConst : public QuarkSourceConstruction<T>
  {
  public:
    //! Full constructor
    PointQuarkSourceConst(const PointQuarkSourceConstParams& p) : params(p) {}

    //! Construct the source
    T operator()(const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    PointQuarkSourceConst() {}

  private:
    PointQuarkSourceConstParams  params;   /*!< source params */
  };

}  // end namespace Chroma


#endif
