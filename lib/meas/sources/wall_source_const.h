// -*- C++ -*-
// $Id: wall_source_const.h,v 2.1 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#ifndef __wall_source_const_h__
#define __wall_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Wall source parameters
  /*! @ingroup sources */
  struct WallQuarkSourceConstParams
  {
    WallQuarkSourceConstParams();
    WallQuarkSourceConstParams(XMLReader& in, const std::string& path);
    
    int              j_decay;              /*!< decay direction */
    int              t_source;             /*!< source time slice location */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, WallQuarkSourceConstParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const WallQuarkSourceConstParams& param);


  //! Name and registration
  /*! @ingroup sources */
  namespace WallQuarkSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

  
  //! Wall source construction
  /*! @ingroup sources
   *
   * Create a wall propagator source
   */
  template<typename T>
  class WallQuarkSourceConst : public QuarkSourceConstruction<T>
  {
  public:
    //! Full constructor
    WallQuarkSourceConst(const WallQuarkSourceConstParams& p) : params(p) {}

    //! Construct the source
    T operator()(const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    WallQuarkSourceConst() {}

  private:
    WallQuarkSourceConstParams  params;   /*!< source params */
  };

}  // end namespace Chroma


#endif
