// -*- C++ -*-
// $Id: rndzNwall_source_const.h,v 1.2 2009-05-01 22:41:13 kostas Exp $
/*! \file
 *  \brief Random Z2 wall source construction
 */

#ifndef __rndzNwall_source_const_h__
#define __rndzNwall_source_const_h__

#include "meas/sources/source_construction.h"
#include "io/xml_group_reader.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace RandZNWallQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();

    //! Random ZN wall source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Seed     ran_seed;             /*!< Set the seed to this value */

      int      j_decay;              /*!< decay direction */
      int      t_source;             /*!< source time slice location */
      int      N ;                   /*!< allow arbitrary Z_N */

      bool     quark_smear_lastP;    /*!< Flag controlling order of smearing */
      
      GroupXML_t quark_smearing;     /*!< xml string holding smearing params */
      GroupXML_t quark_displacement; /*!< displacement xml */
      GroupXML_t link_smearing;      /*!< link smearing xml */
    };

  
    //! Random ZN wall source construction
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
  void read(XMLReader& xml, const string& path, RandZNWallQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const RandZNWallQuarkSourceConstEnv::Params& param);

}  // end namespace Chroma


#endif
