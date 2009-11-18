// -*- C++ -*-
/*! \file
 *  \brief Holds of vectors and weights
 */

#ifndef __subset_vectors_h__
#define __subset_vectors_h__

#include "chromabase.h"
#include "handle.h"
#include "util/ferm/map_obj/map_obj_memory.h"
#include "util/ferm/map_obj/map_obj_disk.h"

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

    //! Constructor -- No Filename: Use memory based map
    SubsetVectors(): evectors(new MapObjectMemory<int,T>()) {}

    //! Constructor -- Filename: Use disk pased map
    SubsetVectors(const std::string& filename) : evectors(new MapObjectDisk<int,T>(filename)) {} 

    //! Full constructor
    SubsetVectors(const multi1d<SubsetVectorWeight_t>& eval, int decay_dir_, const multi1d<T>& evec) 
      : evalues(eval), evectors(new MapObjectMemory<int,T>()), decay_dir(decay_dir_) {
      if (eval.size() != evec.size() ) {
	QDPIO::cout << "Eval array size not equal to evec array size" << endl;
	QDP_abort(1);
      }

      evectors->openWrite();
      for(int i=0; i < evec.size(); i++) { 
	evectors->insert(i,evec[i]);
      }
      evectors->openRead();

    }

    //! Destructor
    ~SubsetVectors() {}
      
    //! Number of vectors
    int getNumVectors() const {return evectors->size();}

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

    //! Lookup vector -- must be in Read Mode
    void lookup(int i, T& v) const { evectors->lookup(i,v); }

    //! Insert vector -- must be in Write Mode
    void insert(int i, const T& v) { evectors->insert(i,v); }

    //! Update vector -- must be in Update mode
    void update(int i, const T& v) { evectors->update(i,v); }


    //! Resize Evalues array 
    void resizeEvalues(int size) { evalues.resize(size); }
    

    //! set write mode
    void openWrite() { evectors->openWrite(); }

    //! set read mode
    void openRead() { evectors->openRead(); }

    //! seet update mode 
    void openUpdate() { evectors->openUpdate();}

    //! get size of Evalues array 
    int evaluesSize() const { return evalues.size(); }

    //! get size of Evectors array (to be removed post refactoring) 
    int evectorsSize() const { return evectors->size(); }

  private:
    int decay_dir;
    multi1d<SubsetVectorWeight_t> evalues;
    Handle< MapObject<int,T> > evectors;
  };
}

















#endif
