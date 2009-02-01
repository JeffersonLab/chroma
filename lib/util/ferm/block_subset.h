// -*- C++ -*-
// $Id: block_subset.h,v 1.3 2009-02-01 05:34:54 kostas Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#ifndef __FUNC_block_h__
#define __FUNC_block_h__

#include "chromabase.h"

namespace Chroma
{

  //! Function object used for constructing a block  set
  class BlockFunc : public SetFunc
  {
  public:
    BlockFunc(const multi1d<int>& blk ): time_dir(-1){
      init(blk);
    }
    BlockFunc(int dir, const multi1d<int>& blk ): time_dir(dir){
      init(blk);
    }

    int operator() (const multi1d<int>& coordinate) const{
      int b = 0 ;
      
      for(int d(Nd-1); d>-1 ;d--)
	b = blockNum[d]*b + coordinate[d]/block[d] ;

      return b ;      
    }

    int numSubsets() const{
      return Nblocks ;
    }

  private:
    BlockFunc() {}  // hide default constructor
    
    void init(const multi1d<int>& blk)  {
      blockNum.resize(Nd) ;
      block.resize(Nd) ;
      Nblocks = 1 ;
      int k(0);
      for(int d(0); d < Nd; d++){
	if (d == time_dir)
	  block[d] = Layout::lattSize()[time_dir] ;
	else{
	  if(k<blk.size()){
	    block[d] = blk[k] ;
	    k++ ;
	  }
	  else
	    block[d] = Layout::lattSize()[d] ;
	}
	blockNum[d] = Layout::lattSize()[d] / block[d];
	Nblocks *= blockNum[d] ;
      }// d

      /** DEBUG
      QDPIO::cout<<"BlockNums: ";
      for(int d(0); d < Nd; d++){
	QDPIO::cout<<blockNum[d]<<" ";
      }
      QDPIO::cout<<endl;

      QDPIO::cout<<"BlockSize: ";
      for(int d(0); d < Nd; d++){
	QDPIO::cout<<block[d]<<" ";
      }
      QDPIO::cout<<endl;
      END DEBUG **/
    }


    int time_dir;
    multi1d<int> block ;
    multi1d<int> blockNum;
    int Nblocks ;
  };

} // namespace Chroma

#endif
