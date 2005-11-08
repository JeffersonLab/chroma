// -*- C++ -*-
// $Id: rndz2wall_source_const.h,v 2.1 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief Random Z2 wall source construction
 */

#ifndef __rndz2wall_source_const_h__
#define __rndz2wall_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Random Z2 wall source parameters
  /*! @ingroup sources */
  struct RandZ2WallQuarkSourceConstParams
  {
    RandZ2WallQuarkSourceConstParams();
    RandZ2WallQuarkSourceConstParams(XMLReader& in, const std::string& path);
    
    Seed             ran_seed;             /*!< Set the seed to this value */

    int              j_decay;              /*!< decay direction */
    int              t_source;             /*!< source time slice location */
  };


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, RandZ2WallQuarkSourceConstParams& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const RandZ2WallQuarkSourceConstParams& param);


  //! Name and registration
  /*! @ingroup sources */
  namespace RandZ2WallQuarkSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
  }

  
  //! Random Z2 wall source construction
  /*! @ingroup sources
   *
   * Create a random z2 wall source
   */
  template<typename T>
  class RandZ2WallQuarkSourceConst : public QuarkSourceConstruction<T>
  {
  public:
    //! Full constructor
    RandZ2WallQuarkSourceConst(const RandZ2WallQuarkSourceConstParams& p) : params(p) {}

    //! Construct the source
    T operator()(const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    RandZ2WallQuarkSourceConst() {}

  private:
    RandZ2WallQuarkSourceConstParams  params;   /*!< source params */
  };

}  // end namespace Chroma


#endif
