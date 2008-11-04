// -*- C++ -*-
// $Id: vector_quark_smearing.h,v 3.3 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Vector Smearing: Use an outerproduct of vectors as the 
 *  smearing scheme. 
 */

#ifndef __vector_smearing_h__
#define __vector_smearing_h__

#include "util/ferm/subset_vectors.h"
#include "meas/smear/quark_smearing.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{

  //! Name and registration
  namespace VectorQuarkSmearingEnv
  {
    bool registerAll();
  
    //! Return the name
    std::string getName();

    //! Params for Vector Smearing 
    /*! @ingroup smear */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
     
      std::string vecs_id;
      Real sigma;                /*!< exponential smearing wieght */
      int no_smear_dir;         /*!< don't allow smearing in this direction */
    };


    //! Vector smearing
    /*! @ingroup smear
     *
     * vector quark smearing object
     */
    template<typename T> 
    class VectorQuarkSmear : public QuarkSmearing<T>
    {
    public:
      //! Full constructor
      VectorQuarkSmear(const Params& p) :
	params(p), 
	vecs(TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.vecs_id))
	{}

      //! Smear the quark
      void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      VectorQuarkSmear() {}

    private:
      Params  params;   /*!< smearing params */
      SubsetVectors<LatticeColorVector> vecs; /*!< vectors to be used */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup smear */
  void read(XMLReader& xml, const string& path, VectorQuarkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write(XMLWriter& xml, const string& path, const VectorQuarkSmearingEnv::Params& param);

}  // end namespace Chroma


#endif
