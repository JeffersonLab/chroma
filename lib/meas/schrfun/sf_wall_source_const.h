// -*- C++ -*-
// $Id: sf_wall_source_const.h,v 1.1 2007-08-25 04:07:41 edwards Exp $
/*! \file
 *  \brief Wall source construction for Schroedinger Functional
 */

#ifndef __sf_wall_source_const_h__
#define __sf_wall_source_const_h__

#include "meas/schrfun/sf_source_construction.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup schrfun */
  namespace WallSFSourceConstEnv
  {
    extern const std::string name;
    bool registerAll();

  
    //! Wall source parameters
    /*! @ingroup schrfun */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    };


    //! Wall source construction
    /*! @ingroup schrfun
     *
     * Create a wall propagator source
     */
    template<typename T>
    class SourceConst : public SFSourceConstruction<T>
    {
    public:
      //! Full constructor
      SourceConst(const Params& p) : params(p) {}

      //! Construct the source
      T operator()(const multi1d<LatticeColorMatrix>& u, int t0, int decay_dir) const;

    private:
      //! Hide partial constructor
      SourceConst() {}

    private:
      Params  params;   /*!< source params */
    };

  }  // end namespace

}  // end namespace Chroma


#endif
