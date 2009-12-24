// -*- C++ -*-
// $Id: vector_quark_smearing.h,v 3.4 2008-11-10 22:05:54 jbulava Exp $
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
    class QuarkSmear : public QuarkSmearing<T>
    {
    public:
      //! Full constructor
      QuarkSmear(const Params& p) : params(p) 
	{
	  try{
	    vecs = TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.vecs_id);
	  }
	  catch(std::string & e)
	  {
	    QDPIO::cerr << " Caught Exception reading vecs: " << e << endl;
	    QDP_abort(1);
	  }
	}

      //! Smear the quark
      void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      QuarkSmear() {}

    private:
      Params  params;   /*!< smearing params */
      Handle< MapObject<int,EVPair<LatticeColorVector> > >  vecs; /*!< vectors to be used */
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
