// -*- C++ -*-
// $Id: sf_sh_source_const.h,v 1.1 2007-08-25 04:07:40 edwards Exp $
/*! \file
 *  \brief Shell source construction with Schroedinger Functional
 */

#ifndef __sf_sh_source_const_h__
#define __sf_sh_source_const_h__

#include "meas/schrfun/sf_source_construction.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Name and registration
  /*! @ingroup schrfun */
  namespace ShellSFSourceConstEnv
  {
    extern const std::string name;
    bool registerAll();


    //! Point source parameters
    /*! @ingroup schrfun */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& out, const std::string& path) const;
    
      GroupXML_t       quark_smearing;       /*!< xml string holding smearing params */
      multi1d<int>     t_srce;               /*!< 3d source location */
    };


    //! Shell source construction
    /*! @ingroup schrfun
     *
     * Create a shell quark source
     */
    template<typename T>
    class SourceConst : public SFSourceConstruction<T>
    {
    public:
      //! Full constructor
      SourceConst(const Params& p) : params(p) {}

      //! Construct the source
      T operator()(const multi1d<LatticeColorMatrix>& u, int t0, int decay_dir) const;

    private:
      //! Hide partial constructor
      SourceConst() {}

    private:
      Params  params;   /*!< source params */
    };

  }  // end namespace ShellSFSourceConstEnv

}  // end namespace Chroma


#endif
