// -*- C++ -*-
// $Id: rndz2wall_source_const.h,v 3.2 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Random Z2 wall source construction
 */

#ifndef __rndz2wall_source_const_h__
#define __rndz2wall_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace RandZ2WallQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();

    //! Random Z2 wall source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Seed             ran_seed;             /*!< Set the seed to this value */

      int              j_decay;              /*!< decay direction */
      int              t_source;             /*!< source time slice location */
    };

  
    //! Random Z2 wall source construction
    /*! @ingroup sources
     *
     * Create a random z2 wall source
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
  void read(XMLReader& xml, const string& path, RandZ2WallQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const RandZ2WallQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
