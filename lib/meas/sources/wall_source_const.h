// -*- C++ -*-
// $Id: wall_source_const.h,v 3.2 2008-11-04 18:43:59 edwards Exp $
/*! \file
 *  \brief Wall source construction
 */

#ifndef __wall_source_const_h__
#define __wall_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace WallQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();
  
    //! Wall source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      int              j_decay;              /*!< decay direction */
      int              t_source;             /*!< source time slice location */
    };


    //! Wall source construction
    /*! @ingroup sources
     *
     * Create a wall propagator source
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
  void read(XMLReader& xml, const string& path, WallQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const WallQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
