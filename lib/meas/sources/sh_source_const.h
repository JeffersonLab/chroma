// -*- C++ -*-
/*! \file
 *  \brief Shell source construction
 */

#ifndef __sh_source_const_h__
#define __sh_source_const_h__

#include "meas/sources/source_construction.h"
#include "io/xml_group_reader.h"

namespace Chroma
{


  //! Name and registration
  /*! @ingroup sources */
  namespace ShellQuarkSourceConstEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();


    //! Point source parameters
    /*! @ingroup sources */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& out, const std::string& path) const;
    
      bool             quark_smear_lastP;    /*!< Flag controlling order of smearing */

      GroupXML_t       quark_smearing;       /*!< xml std::string holding smearing params */
      GroupXML_t       quark_displacement;   /*!< displacement xml */
      GroupXML_t       link_smearing;        /*!< link smearing xml */

      int              j_decay;              /*!< Decay direction */
      multi1d<int>     t_srce;               /*!< source location */
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
  void read(XMLReader& xml, const std::string& path, ShellQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup sources */
  void write(XMLWriter& xml, const std::string& path, const ShellQuarkSourceConstEnv::Params& param);


}  // end namespace Chroma


#endif
