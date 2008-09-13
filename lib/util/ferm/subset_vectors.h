// -*- C++ -*-
// $Id: subset_vectors.h,v 1.1 2008-09-13 19:57:05 edwards Exp $
/*! \file
 *  \brief Holds of vectors and weights
 */

#ifndef __subset_vectors_h__
#define __subset_vectors_h__

#include "chromabase.h"

namespace Chroma
{
  //! Weights for subset of vectors
  /*! \ingroup ferm */
  struct SubsetVectorWeight_t
  {
    multi1d<Real> weights;
  };

  //! Reader
  /*! \ingroup ferm */
  void read(XMLReader& xml, const std::string& path, SubsetVectorWeight_t& param);
 
  //! Writer
  /*! \ingroup ferm */
  void write(XMLWriter& xml, const std::string& path, const SubsetVectorWeight_t& param);
 

  //! Holds of vectors and weights
  /*!
   * \ingroup ferm
   */
  template<typename T>
  class SubsetVectors
  {
  public:
    //! Partial constructor
    SubsetVectors() {}

    //! Full constructor
    SubsetVectors(const multi1d<SubsetVectorWeight_t>& eval, int decay_dir_, const multi1d<T>& evec) 
      : evalues(eval), decay_dir(decay_dir_), evectors(evec) {
      if (eval.size() != evec.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }
    }

    //! Destructor
    ~SubsetVectors() {}
      
    //! Number of vectors
    int getNumVectors() const {return evectors.size();}

    //! Extent of decay direction
    int getDecayExtent() const {return evectors.size();}

    //! Get decay direction
    int getDecayDir() const {return decay_dir;}

    //! Set decay direction
    int& getDecayDir() {return decay_dir;}

    //! Getter
    const multi1d<SubsetVectorWeight_t>& getEvalues() const {return evalues;}
    //! Setter
    multi1d<SubsetVectorWeight_t>& getEvalues() {return evalues;}

    //! Getter
    const multi1d<T>& getEvectors() const {return evectors;}
    //! Setter
    multi1d<T>& getEvectors() {return evectors;}

  private:
    int decay_dir;
    multi1d<SubsetVectorWeight_t> evalues;
    multi1d<T> evectors;
  };
}

















#endif
