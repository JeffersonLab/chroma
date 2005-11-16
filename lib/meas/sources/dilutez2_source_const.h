// -*- C++ -*-
// $Id: dilutez2_source_const.h,v 2.2 2005-11-16 02:34:58 edwards Exp $
/*! \file
 *  \brief Random complex Z(2) source construction using dilution
 *
 * Uses the Dublin "dilution" stochastic source strategy of hep-lat/0505023
 */

#ifndef __dilutez2_source_const_h__
#define __dilutez2_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{

  //! Dilute Z2 quark source namespace, parameters, and classes
  /*! @ingroup sources */
  namespace DiluteZ2QuarkSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;

  
    //! Random complex Z(2) sources using dilution
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Seed                     ran_seed;             /*!< Set the seed to this value */

      multi1d<int>             spatial_mask_size;    /*!< Spatial size of periodic mask */
      multi1d< multi1d<int> >  spatial_mask;         /*!< Sites included in site mask */
      multi1d<int>             color_mask;           /*!< Color size of periodic mask */
      multi1d<int>             spin_mask;            /*!< Spin size of periodic mask */

      int                      j_decay;              /*!< decay direction */
      int                      t_source;             /*!< source time slice location */
    };


    //! Random complex Z(2) sources using dilution
    /*! @ingroup sources
     *
     * Create a random Z(2) using dilution
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

  }  // end namespace DiluteZ2QuarkSourceConstEnv


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DiluteZ2QuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DiluteZ2QuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
