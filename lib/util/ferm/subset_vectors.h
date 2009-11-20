// -*- C++ -*-
/*! \file
 *  \brief Holds of vectors and weights
 */

#ifndef __subset_vectors_h__
#define __subset_vectors_h__

#include "chromabase.h"
#include "handle.h"
#include "util/ferm/subset_ev_pair.h"
#include "util/ferm/map_obj/map_obj_memory.h"
#include "util/ferm/map_obj/map_obj_disk.h"


namespace Chroma
{
 

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
    SubsetVectors(): ev_pairs(new MapObjectMemory<int,EVPair<T> >()) {}

    //! Constructor -- Filename: Use disk pased map
    SubsetVectors(const std::string& filename) : ev_pairs(new MapObjectDisk<int,EVPair<T> >(filename)) {} 

    //! Full constructor
    SubsetVectors(const multi1d<SubsetVectorWeight_t>& eval, int decay_dir_, const multi1d<T>& evec) 
      : ev_pairs(new MapObjectMemory<int,EVPair<T> >()), decay_dir(decay_dir_) {

      ev_pairs->openWrite();
      for(int i=0; i < evec.size(); i++) { 
	EVPair<T> pair;
	pair.eigenValue= eval[i];
	pair.eigenVector = evec[i];

	ev_pairs->insert(i,pair);
      }
      ev_pairs->openRead();

    }

    //! Destructor
    ~SubsetVectors() {}
      
    //! Number of vectors
    int getNumVectors() const {return ev_pairs->size();}

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
	lookup(i, evs[i]);
      }
    }

#if 0
    //! Indexed e-Value Getter
    const SubsetVectorWeight_t& getEvalue(int i) const { return evalues[i];}

    //! Indexed e-Value Setter
    SubsetVectorWeight_t& getEvalue(int i) { return evalues[i];}
#endif


    //! Lookup ev_pair -- must be in Read Mode
    void lookup(int i, EVPair<T>& val) const { ev_pairs->lookup(i,val); }

    //! Insert ev_pair -- must be in Write Mode
    void insert(int i, const EVPair<T>& val) { ev_pairs->insert(i,val); }

    //! Update ev_pair -- must be in Update mode
    void update(int i, const EVPair<T>& val) { ev_pairs->update(i,val); }



    //! Convenience Lookup of just EVal
    void lookup(int i, SubsetVectorWeight_t& w) const {
      EVPair<T> tmp;
      ev_pairs->lookup(i,tmp);
      w = tmp.eigenValue;
    }

    //! Convenience Lookup of just Evector
    void lookup(int i, T& vec) const {
      EVPair<T> tmp;
      ev_pairs->lookup(i,tmp);
      vec = tmp.eigenVector;
    }

    //! Convenience Update of just EVal
    void update(int i, const SubsetVectorWeight_t& w) {
      EVPair<T> tmp;

      // Read 
      ev_pairs->lookup(i,tmp);

      // Modify 
      tmp.eigenValue = w;

      // Write 
      ev_pairs->update(i,tmp);
    }

    //! Convenience Update of just Evector
    void update(int i, const T& vec) {
      EVPair<T> tmp;

      // Read 
      ev_pairs->lookup(i,tmp);

      // Modify 
      tmp.eigenVector = vec;

      // Write 
      ev_pairs->update(i,tmp);
    }

    //! Resize Evalues array 
    // void resizeEvalues(int size) { }
    

    //! set write mode
    void openWrite() { ev_pairs->openWrite(); }

    //! set read mode
    void openRead() { ev_pairs->openRead(); }

    //! seet update mode 
    void openUpdate() { ev_pairs->openUpdate();}

    //! get size of Evalues array 
    int evaluesSize() const { return ev_pairs->size(); }

    //! get size of Evectors array (to be removed post refactoring) 
    int evectorsSize() const { return ev_pairs->size(); }

    //! get size of map 
    int size() const { return ev_pairs->size(); } 

  private:
    int decay_dir;
    // multi1d<SubsetVectorWeight_t> evalues;
    Handle< MapObject<int, EVPair<T> > > ev_pairs;
  };
}

















#endif
