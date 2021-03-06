// -*- C++ -*-
/*! \file
 *  \brief Wave source construction fpr Schroedinger Functional
 */

#ifndef __sf_wave_source_const_h__
#define __sf_wave_source_const_h__

#include "meas/sources/source_construction.h"
#include "io/xml_group_reader.h"
#include "io/enum_io/enum_plusminus_io.h"

namespace Chroma
{

  //! Name and registration
  /*! @ingroup sources */
  namespace SFWaveQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();
  
    //! Wave source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      GroupXML_t       fermbc;               /*!< xml std::string holding smearing params */
      PlusMinus        direction;            /*!< direction for decay */

      int              j_decay;              /*!< decay direction */
      int              t_source;             /*!< source time slice location */
      multi1d<int>     mom ;                 /*!< plane wave momentum */
    };


    //! Wave source construction for Schroedinger Functional
    /*! @ingroup sources
     *
     * Create a wave propagator source
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

}  // end namespace Chroma


#endif
