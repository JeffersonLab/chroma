// -*- C++ -*-
// $Id: diluteGrid_source_const.h,v 3.2 2008-11-27 03:22:17 kostas Exp $
/*! \file
 *  \brief Random Z(M) source construction using dilution
 *
 * 
 */

#ifndef __diluteGrid_source_const_h__
#define __diluteGrid_source_const_h__

#include "meas/sources/source_construction.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Dilute Z(N) quark source namespace, parameters, and classes
  /*! @ingroup sources */
  namespace DiluteGridQuarkSourceConstEnv
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

      multi1d<int>             spatial_mask_size;    /*!< Spatial size of periodic mask */
      multi1d< multi1d<int> >  spatial_mask;         /*!< Sites included in site mask */
      //multi1d<int>             color_mask;           /*!< Color size of periodic mask */
      //multi1d<int>             spin_mask;            /*!< Spin size of periodic mask */
      int color ;                                    /*!< the color */
      int spin ;                                     /*!< the spin index */
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

  }  // end namespace DiluteGridQuarkSourceConstEnv


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, DiluteGridQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const DiluteGridQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
