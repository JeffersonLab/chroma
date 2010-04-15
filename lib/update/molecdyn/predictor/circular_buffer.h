// -*- C++ -*-
// $Id: circular_buffer.h,v 3.1 2006-08-26 02:08:42 edwards Exp $
/*! \file
 * \brief Circular buffers
 *
 * Circular buffers for predictors
 */

#ifndef CIRCULAR_BUFFER_H
#define CIRCULAR_BUFFER_H

#include "chromabase.h"

using namespace QDP;

namespace Chroma 
{ 

  //! Circular Buffer
  /*! 
   * @ingroup predictor
   * Definitions for a templated circular buffer class to be used in Chronological Predictor.
   *
   * Usage: construct with CircularBuffer(int n_max)  
   *        with n_max the maximum number of elements -- this allocates storage
  
   * add elements with push(elem);
  
   * reset with reset()
  
   * Access elements using operator[unsigned int i]
   *
   * i=0 most recent
   * i=1 next most rcent 
   * ...
   * i=size()-1 // least recent
   */
  template<typename T>
  class CircularBuffer 
  {
  public:
    
    // Exception struct 
    struct OutOfBoundsException { 
      OutOfBoundsException(const std::string& s, unsigned int i_, unsigned int size_) : error_string(s), i(i_), size(size_) {}
      const std::string error_string;
      const unsigned int i;
      const unsigned int size;
    };
    
    
    
    //! Constructor
    CircularBuffer(int n_max) : size_max(n_max), size_internal(0), start(n_max), end(n_max-1) {
      q.resize(n_max);
    }
    
    //! Destructor is automatic
    ~CircularBuffer() {};
    
    //! Copy constructor
    CircularBuffer(const CircularBuffer<T>& c) : size_max(c.size_max),
						 size_internal(c.size), q(c.q), start(c.start), end(c.end) {}
    
    
    // Operations:
    
    //! Nuke it:
    void reset(void) { // Discard all circular buffer elements
      size_internal = 0;
      start= size_max - 1; 
      end = size_max - 1;
    }
    
    
    // Accessors   This should maybe be replaced with iterators
    // at some stage
    
    //! Get the maximum number of data items one can store
    int sizeMax(void) const { // Maximum number of elements
      return size_max;
    }
    
    //! Get the current number of items stored
    int size(void) const {  // THe current number of elements
      return size_internal;
    }
    
    //! get the ith most recent item
    void get(const unsigned int i, T& x) const 
    {
      START_CODE();

      if( i >= size_internal ) { 
	throw OutOfBoundsException(std::string("Index Out of bounds"), i, size_internal);
      }
      else { 
	// Get index of ith element from start with wraparound
	unsigned int index = (start + i) % size_max;
	x= q[index];
      }
    
      END_CODE();
    }
    
    //! Is empty check
    bool isEmpty(void) const { 
      return (size_internal == 0);
    }
    
    
    //! push in an item as most recent.
    void push(const T& e) 
    {
      START_CODE();
      
      if( size_internal == 0 ) { 
	// First element -- to the end of the list
	
	start = size_max -1; // end of list for first element
	end = size_max -1;   // first element is also last element
	
	q[start]=e;          // store element
	
	size_internal = 1;   // size is now 1
      }
      else {
	
	// Not empty: 
	
	// if adding this element overruns buffer size
	// then we drop the least recent element,
	// and we don't increment the size of the array
	
	// Otherwise if adding the element doesn't exhaust buffer
	// we decrease start(with wraparound) and increase the size.
	
	
	if( size_internal+1 > size_max ) {
	  // We would exceed max size
	  
	  // decrease end pointer
	  // with wraparound -- this drops last element
	  if (end > 0) { 
	    end--;
	  }
	  else { 
	    end = size_max - 1;
	  }
	  
	  // decrease start pointer
	  // with wraparound
	  // the result of this that the start pointer can point
	  // to the slot just vacated by decreasing the end pointer
	  if( start > 0 ) { 
	    start--;
	  }
	  else {
	    start = size_max -1;
	  } 
	  
	  // copy element to start
	  q[start] = e;
	  
	  // don't increase size
	}
	else { 
	  // We still have room before we hit max size
	  // decrease start (with wraparound of course).
	  // If there is still space we should not need to 
	  // worry about overtaking the end pointer.
	  if ( start > 0 ) { 
	    start--;
	  }
	  else {
	    start = size_max -1;
	  }
	  
	  // Store element
	  q[start] = e;
	  
	  // Increase size count
	  size_internal++;
	}
      }
    
      END_CODE();
    }
    
    void replaceHead(const T& x)
    {
      q[start]=x;
    }

    
  private:
    unsigned int size_max;
    unsigned int size_internal;
    multi1d<T> q;
    unsigned int start; // Start index -- these move as things are added
    unsigned int end; // End index -- these move as things are added
  };


  //! Circular buffer of arrays
  /*! @ingroup predictor */
  template<typename T>
  class CircularBufferArray 
  {
  public:
    
    // Exception struct 
    struct OutOfBoundsException { 
      OutOfBoundsException(const std::string& s, unsigned int i_, unsigned int size_, unsigned int n5_) : error_string(s), i(i_), size(size_), n5(n5_) {}
      const std::string error_string;
      const unsigned int i;
      const unsigned int size;
      const unsigned int n5;
    };
    

    struct SizeMismatchException {
      SizeMismatchException(const std::string& s, unsigned int N5_buf_, unsigned int N5_vec_) : error_string(s), N5_buf(N5_buf_), N5_vec(N5_vec_) {}
      const std::string error_string;
      const unsigned int N5_buf;
      const unsigned int N5_vec;
    };
    
    //! Constructor
    CircularBufferArray(unsigned int n_max, unsigned int n5) : size_max(n_max), size_internal(0), start(n_max), end(n_max-1), N5(n5) {
      q.resize(n_max, n5);
    }
    
    //! Destructor is automatic
    ~CircularBufferArray() {};
    
    //! Copy constructor
    CircularBufferArray(const CircularBufferArray<T>& c) : size_max(c.size_max),
							   size_internal(c.size), q(c.q), start(c.start), end(c.end), N5(c.N5) {}
    
    
    // Operations:
    
    //! Nuke it:
    void reset(void) { // Discard all circular buffer elements
      size_internal = 0;
      start= size_max - 1; 
      end = size_max - 1;
    }
    
    
    // Accessors   This should maybe be replaced with iterators
    // at some stage
    
    //! Get the maximum number of data items one can store
    int sizeMax(void) const { // Maximum number of elements
      return size_max;
    }
    
    //! Get the current number of items stored
    int size(void) const {  // THe current number of elements
      return size_internal;
    }
    
    //! get the ith most recent item
    void  get(const unsigned int i, multi1d<T>& x) const 
    {
      START_CODE();

      if( i >= size_internal ) { 
	throw OutOfBoundsException(std::string("Index Out of bounds"), i, size_internal, N5);
      }
      else { 
	// Get index of ith element from start with wraparound
	unsigned int index = (start + i) % size_max;

	x.resize(N5);
	for(int s=0; s < N5; s++) { 
	  x[s] = q[index][s];
	}
      }
    
      END_CODE();
    }
    
    //! Is empty check
    bool isEmpty(void) const { 
      return (size_internal == 0);
    }
    
    
    //! push in an item as most recent.
    void push(const multi1d<T>& e) 
    {
      START_CODE();

      if( e.size() != N5) { 
	throw SizeMismatchException("Attempting to push vector of wrong size into circular buffer", N5, e.size());
      }

      if( size_internal == 0 ) { 
	// First element -- to the end of the list
	
	start = size_max -1; // end of list for first element
	end = size_max -1;   // first element is also last element

	for(int s = 0; s < N5; s++) {
	  q[start][s]=e[s];          // store element
	}

	size_internal = 1;   // size is now 1
      }
      else {
	
	// Not empty: 
	
	// if adding this element overruns buffer size
	// then we drop the least recent element,
	// and we don't increment the size of the array
	
	// Otherwise if adding the element doesn't exhaust buffer
	// we decrease start(with wraparound) and increase the size.
	
	
	if( size_internal+1 > size_max ) {
	  // We would exceed max size
	  
	  // decrease end pointer
	  // with wraparound -- this drops last element
	  if (end > 0) { 
	    end--;
	  }
	  else { 
	    end = size_max - 1;
	  }
	  
	  // decrease start pointer
	  // with wraparound
	  // the result of this that the start pointer can point
	  // to the slot just vacated by decreasing the end pointer
	  if( start > 0 ) { 
	    start--;
	  }
	  else {
	    start = size_max -1;
	  } 
	  
	  // copy element to start
	  for(int s = 0; s < N5; s++) { 
	    q[start][s] = e[s];
	  }
	  // don't increase size
	}
	else { 
	  // We still have room before we hit max size
	  // decrease start (with wraparound of course).
	  // If there is still space we should not need to 
	  // worry about overtaking the end pointer.
	  if ( start > 0 ) { 
	    start--;
	  }
	  else {
	    start = size_max -1;
	  }
	  
	  // Store element
	  for(int s = 0; s < N5; s++) { 
	    q[start][s] = e[s];
	  }
	  // Increase size count
	  size_internal++;
	}
      }
    
      END_CODE();
    }
    
    
    
  private:
    unsigned int size_max;
    unsigned int size_internal;
    multi2d<T> q;
    unsigned int start; // Start index -- these move as things are added
    unsigned int end; // End index -- these move as things are added
    unsigned int N5;
  };

} // End namespace Chroma 
#endif
