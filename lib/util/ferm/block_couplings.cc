// -*- C++ -*-
// $Id: block_couplings.cc,v 1.2 2009-01-30 21:07:21 kostas Exp $
/*! \file
 *  \brief Caclulates the couplings between neighboring blocks given a displacement path
 */

#include "chromabase.h"
#include "util/ferm/block_couplings.h"
namespace Chroma
{

  LatticeInteger displace(const LatticeInteger& in,
			  const int len, const int dir){
    LatticeInteger out = in ;
    LatticeInteger tmp  ;
    if (len > 0)
      for(int n = 0; n < len; ++n){
	tmp = shift(out, FORWARD, dir);
	out = tmp ;
      }
    else // If length = or < 0.  If length == 0, does nothing.
      for(int n = 0; n > len; --n){
	tmp = shift(out, BACKWARD, dir);
	out = tmp;
      }
    return out;    
  }

  vector<int> block_couplings(const int b,
			      const Set& S,const multi1d<int>& disp,
			      const int len){
    vector<int> b_list ;
    LatticeInteger blk = zero ;
    LatticeInteger One = 1 ;
    blk[S[b]] = One ;
    for(int b2(0) ; b2<S.numSubsets();b2++){
      LatticeInteger blk2 = zero ;
      blk[S[b2]] = One ;
      for(int i=0; i < disp.size(); ++i)
	if (disp[i] > 0){
	  int disp_dir = disp[i] - 1;
	  int disp_len = len;
	  blk2 = displace(blk2, disp_len, disp_dir);
	}
	else if (disp[i] < 0){
	  int disp_dir = -disp[i] - 1;
	  int disp_len = -len;
	  blk2 = displace(blk2, disp_len, disp_dir);
	}

      int flag = toInt(sum(blk2*blk));
      if(flag>0)
	b_list.push_back(b2);
    }

    return b_list ;
      
  }

} // namespace Chroma
