// -*- C++ -*-
// $Id: jacobi_quark_smearing.h,v 3.1 2009-02-20 15:10:24 edwards Exp $
/*! \file
 *  \brief Jacobi smearing of color vector and propagator
 */

#ifndef __jacobi_quark_smearing_h__
#define __jacobi_quark_smearing_h__

#include "meas/smear/quark_smearing.h"

namespace Chroma 
{
  //! Name and registration
  /*! @ingroup smear */
  namespace JacobiQuarkSmearingEnv
  {
    bool registerAll();

    //! Return the name
    std::string getName();

    //! Params for Jacobi quark smearing
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      Real kappa;			/*!< Hopping parameter */
      int  iter;			/*!< Number of smearing hits */
      int  no_smear_dir;		/*!< No smearing in this direction */
    };


    //! Jacobi quark smearing
    /*! @ingroup smear
     *
     * Jacobi quark smearing object
     */
    template<typename T>
    class QuarkSmear : public QuarkSmearing<T>
    {
    public:
      //! Full constructor
      QuarkSmear(const Params& p) : params(p) {}
      
      //! Smear the quark
      void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      QuarkSmear() {}

    private:
      Params  params;   /*!< smearing params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, JacobiQuarkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const JacobiQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif
