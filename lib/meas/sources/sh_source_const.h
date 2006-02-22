// -*- C++ -*-
// $Id: sh_source_const.h,v 2.5 2006-02-22 04:34:05 edwards Exp $
/*! \file
 *  \brief Shell source construction
 */

#ifndef __sh_source_const_h__
#define __sh_source_const_h__

#include "meas/sources/source_construction.h"

namespace Chroma
{


  //! Name and registration
  /*! @ingroup sources */
  namespace ShellQuarkSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;



    //! Point source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& out, const std::string& path) const;
    
      std::string      source_type;          /*!< source smearing type */

      std::string      quark_smearing;       /*!< xml string holding smearing params */
      std::string      quark_smearing_type;  /*!< quark smearing type name */

      std::string      quark_displacement;      /*!< displacement xml */
      std::string      quark_displacement_type; /*!< displacement type name */

      int              j_decay;              /*!< Decay direction */
      multi1d<int>     t_srce;               /*!< source location */

      std::string      link_smearing;        /*!< link smearing xml */
      std::string      link_smearing_type;   /*!< link smearing type name */
    };


    //! Shell source construction
    /*! @ingroup sources
     *
     * Create a shell quark source
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

  }  // end namespace ShellQuarkSourceConstEnv


  //! Reader
  /*! @ingroup sources */
  void read(XMLReader& xml, const string& path, ShellQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const string& path, const ShellQuarkSourceConstEnv::Params& param);


}  // end namespace Chroma


#endif
