// -*- C++ -*-
// $Id: dilutezN_source_const.h,v 3.3 2008-11-10 00:04:58 kostas Exp $
/*! \file
 *  \brief Random Z(N) source construction using dilution
 *
 * Uses the Dublin "dilution" stochastic source strategy of hep-lat/0505023
 */

#ifndef __dilutezN_source_const_h__
#define __dilutezN_source_const_h__

#include "meas/sources/source_construction.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Dilute Z(N) quark source namespace, parameters, and classes
  /*! @ingroup sources */
  namespace DiluteZNQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();
  
    //! Random complex Z(N) sources using dilution
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      GroupXML_t  smr; /*!< xml holding smearing params */
      GroupXML_t  displace; /*!< xml holding displacement params */
      GroupXML_t  link_smear;  /*!< link smearing xml */

      bool smear ; // a flag that tells me to smear or not to smear

      Seed                     ran_seed;             /*!< Set the seed to this value */
      int                      N;                    /*!< Z(N) */
      
      multi1d<int>             spatial_mask_size;    /*!< Spatial size of periodic mask */
      multi1d< multi1d<int> >  spatial_mask;         /*!< Sites included in site mask */
      multi1d<int>             color_mask;           /*!< Color size of periodic mask */
      multi1d<int>             spin_mask;            /*!< Spin size of periodic mask */

      int                      j_decay;              /*!< decay direction */
      int                      t_source;             /*!< source time slice location */
    };


    //! Random complex Z(N) sources using dilution
    /*! @ingroup sources
     *
     * Create a random Z(N) using dilution
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

  }  // end namespace DiluteZNQuarkSourceConstEnv


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DiluteZNQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DiluteZNQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
