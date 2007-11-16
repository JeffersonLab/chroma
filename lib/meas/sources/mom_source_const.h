// -*- C++ -*-
// $Id: mom_source_const.h,v 3.2 2007-11-16 22:27:33 kostas Exp $
/*! \file
 *  \brief Fixed momentum (wall) source construction
 */

#ifndef __mom_source_const_h__
#define __mom_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace MomWallQuarkSourceConstEnv
  {
    extern const std::string name;
    bool registerAll();

  
    //! MomWall source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      multi1d<int>     mom;                  /*<! prototype momentum */
      int              t_dir;                /*<! time direction */
      multi1d<int>     t_srce;               /*<! the origin for the FT */
    };


    //! MomWall source construction
    /*! @ingroup sources
     *
     * Create a momentum wall propagator source
     */
    template<typename T>
    class SourceConst : public QuarkSourceConstruction<T>
    {
    public:
      //! Full constructor
      SourceConst(const Params& p) : params(p) {}

      //! Construct the source
      T operator()(const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      SourceConst() {}

    private:
      Params  params;   /*!< source params */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, MomWallQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const MomWallQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
