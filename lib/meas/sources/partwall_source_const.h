// -*- C++ -*-
// $Id: partwall_source_const.h,v 3.2 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Partial wall source construction
 */

#ifndef __partial_wall_source_const_h__
#define __partial_wall_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Structure holding directions
  /*! @ingroup sources */
  struct FixedDir_t
  {
    int            dir;                  /*!< Direction to hold fixed */
    int            coord;                /*!< Coordinate along dir held fixed */
  };


  //! Name and registration
  /*! @ingroup sources */
  namespace PartialWallQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();
  
    //! PartialWall source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
//      int              j_decay;              /*!< decay direction */
//      int              t_source;             /*!< source time slice location */

      multi1d<FixedDir_t>  fixed_dirs;
    };


    //! PartialWall source construction
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
  void read(XMLReader& xml, const string& path, PartialWallQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const PartialWallQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
