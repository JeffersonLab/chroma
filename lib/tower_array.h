// -*- C++ -*-

#ifndef __ARRAY_TOWER_H__
#define __ARRAY_TOWER_H__

#include "chromabase.h"

#ifndef __TOWER_H__
#include "tower.h"
#endif


using namespace QDP;

namespace Chroma { 

template<typename T>
class TowerArray { 
public:

  //! Constructor
  /*!
    Initializes an array of towers. In this constructor
    only the height is given, and the width of the tower array
    will be Nd.
    @param[in] height The height of the tower
  */
  TowerArray(int height) : towerHeight(height)  {
    t_array.resize(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      t_array[mu].resize(height);
    }
  }
  
  //! Constructor
  /*!
    Initializes an array of towers. In this constructor
    only the height is given, and the width of the tower array
    will be Nd.
    @param[in] width  The width of the tower
    @param[in] height The height of the tower
  */
  TowerArray(int width, int height) : towerHeight(height) {
    t_array.resize(width);
    for(int mu=0; mu < t_array.size(); mu++) { 
      t_array[mu].resize(height);
    }
  }
  
  
  //! Copy Constructor
  /*!
    Copy a tower array from another tower array.
    Uses operator=() of constituent towers

    @param[in] t   The tower to copy from
  */
  TowerArray(const TowerArray<T>& t) : towerHeight(t.getHeight())
  {
    t_array.resize(t.size());
    for(int mu=0; mu < t_array.size(); mu++) { 
      // Tower assignment will resize
      t_array[mu] = t[mu];
    }
  }

  //! Destructor
  /*!
    Default destructor
  */
  ~TowerArray() { }

  
  //! Assignment
  /*!
    Assigns another TowerArray to this one
    @param[in] t  The input tower to assign
    @return  A const ref to the self of the assigned tower
  */

  TowerArray<T>& operator=(const TowerArray<T>& t) {
    // Resize my t_array
    t_array.resize(t.size());
    towerHeight = t.height();

    for(int mu=0; mu < t_array.size(); mu++) { 
      t_array[mu] = t[mu];
    }

    return *this;
  }

  //! Zero Assignment
  /*!
    Assigns another TowerArray to this one
    @param[in] t  The zero object...
    @return  A const ref to the self of the assigned tower
  */

  TowerArray<T>& operator=(const QDP::Zero& rhs) {
    // Resize my t_array
    for(int mu=0; mu < t_array.size(); mu++) { 
      t_array[mu] = zero;
    }
    return *this;
  }

  //! Size function
  /*!
    Return the size (width) of the tower
    @return The width of the tower
  */
  int size() const {
    return t_array.size();
  }

  //! Height of tower
  /*!
    Return the height of component towers
    @return an int representing the height of the tower
  */

  int getHeight() const { 
    return towerHeight;
  }


  //! Indexing: Getter
  /*!
    Returns one of the component towers.
    @param[in]  mu The index of the tower to return
    @return a const Tower<T>& corresponding to the tower for element m
  */
  const Tower<T>& operator[](int mu) const { 
    return t_array[mu];
  }

  //! Indexing Setter
  /*!
    Returns a 'live' reference to one of the component towers.
    @param[in] mu  The index of the tower to return
    @return a Tower<T>& corresponding to the tower for element mu
  */
  Tower<T>& operator[](int mu) {
    return t_array[mu];
  }

  //! Resize Width function
  /*!
    Resizes width, keeping current height
    @param[in] width
  */
  void resize(int width) 
  {
    t_array.resize(width);
    for(int mu=0; mu < width; mu++) {
      t_array[mu].resize(towerHeight);
    }
  }


  //! Resize Height function
  /*!
    Resize height, keeping current width
    @param[in] h  The tower height
  */
  void resizeHeight(int h)
  {
    for(int mu=0; mu < t_array.size(); mu++){ 
      t_array[mu].resize(h);
    }
    towerHeight=h;
  }

  //! Resize Width & Height function
  /*!
    Resize height, keeping current width
    @param[in] w  The tower width
    @param[in] h  The tower height
  */
  void resizeWidthAndHeight(int w,int h)
  {
    t_array.resize(w);
    for(int mu=0; mu < t_array.size(); mu++){ 
      t_array[mu].resize(h);
    }
    towerHeight=h;
  }

private:
  int towerHeight;
  multi1d< Tower<T> > t_array;

}; // Class


}; // Namespace 




#endif
