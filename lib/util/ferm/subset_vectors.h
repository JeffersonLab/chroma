// -*- C++ -*-
/*! \file
 *  \brief Holds of vectors and weights
 */

#ifndef __subset_vectors_h__
#define __subset_vectors_h__

#include "chromabase.h"
#include "util/ferm/map_obj/map_obj_memory.h"

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
    typedef T Type_t;

    //! Partial constructor
    SubsetVectors() {}

    //! Full constructor
    SubsetVectors(const multi1d<SubsetVectorWeight_t>& eval, int decay_dir_, const multi1d<T>& evec) 
      : evalues(eval), decay_dir(decay_dir_) {
      if (eval.size() != evec.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }

      evectors.openWrite();
      for(int i=0; i < evec.size(); i++) { 
	evectors.insert(i,evec[i]);
      }
      evectors.openRead();

    }

    //! Destructor
    ~SubsetVectors() {}
      
    //! Number of vectors
    int getNumVectors() const {return evectors.size();}

    //! Extent of decay direction
    int getDecayExtent() const {return QDP::Layout::lattSize()[decay_dir];}

    //! Get decay direction
    int getDecayDir() const {return decay_dir;}

    //! Set decay direction
    int& getDecayDir() {return decay_dir;}


    //! Get all ev-s as an array ('hide' internal representation from client).
    void getEvalues(multi1d<SubsetVectorWeight_t>& evs) const 
    {
      evs.resize(evaluesSize());
      for(int i=0; i < evaluesSize(); i++) {
	evs[i] = getEvalue(i);
      }
    }

    //! Indexed e-Value Getter
    const SubsetVectorWeight_t& getEvalue(int i) const { return evalues[i];}

    //! Indexed e-Value Setter
    SubsetVectorWeight_t& getEvalue(int i) { return evalues[i];}

    //! Indexed e-Vector Getter
#if 0
    const T& getEvector(int i) const { return evectors[i]; }
    T& getEvector(int i)  { return evectors[i]; }
#endif

    void lookup(int i, T& v) const { evectors.lookup(i,v); }
    void insert(int i, const T& v) { evectors.insert(i,v); }

    //! Lookup vector 

    //! Resize Evalues array 
    void resizeEvalues(int size) { evalues.resize(size); }
    
    //! Resize Evectors array (to be removed post refactoring)
    void resizeEvectors(int size) { 
      // evectors.resize(size); 
    }


    //! set write mode
    void openWrite() { evectors.openWrite(); }

    //! set read mode
    void openRead() { evectors.openRead(); }

    //! get size of Evalues array 
    int evaluesSize() const { return evalues.size(); }

    //! get size of Evectors array (to be removed post refactoring) 
    int evectorsSize() const { return evectors.size(); }

  private:
    int decay_dir;
    multi1d<SubsetVectorWeight_t> evalues;
    MapObjectMemory<int,T> evectors;
  };
}

















#endif
