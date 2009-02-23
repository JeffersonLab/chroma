// -*- C++ -*-
// $Id: block_couplings.cc,v 1.5 2009-02-23 19:50:15 edwards Exp $
/*! \file
 *  \brief Caclulates the couplings between neighboring blocks given a displacement path
 */

#include "chromabase.h"
#include "util/ferm/block_couplings.h"

namespace Chroma
{
  LatticeInteger displace(const LatticeInteger& in,
			  const int len, 
			  const int dir)
  {
    LatticeInteger out = in;
    LatticeInteger tmp ;
    if (len > 0)
      for(int n = 0; n < len; ++n){
	tmp = shift(out, FORWARD, dir);
	out = tmp;
      }
    else // If length = or < 0.  If length == 0, does nothing.
      for(int n = 0; n > len; --n){
	tmp = shift(out, BACKWARD, dir);
	out = tmp;
      }
    return out;    
  }

  void displace(LatticeInteger& blk2,
		const int len, 
		const multi1d<int>& disp)
  {
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
  }


  vector<int> block_couplings(const int b,
			      const Set& S,
			      const multi1d<int>& disp,
			      const int len)
  {
    vector<int> b_list;
    LatticeInteger blk = zero;
    LatticeInteger One = 1;
    blk[S[b]] = One;

    for(int b2(0); b2 < S.numSubsets(); b2++)
    {
      LatticeInteger blk2 = zero;
      blk2[S[b2]] = One;
      //QDPIO::cout<<"Points in block "<<b2<<" "<<sum(blk2)<<endl;
      displace(blk2,len,disp);

      int flag = toInt(sum(blk2*blk));
      //QDPIO::cout<<"b2 ="<<b2<<" flag "<<flag << endl;
      if(flag > 0)
	b_list.push_back(b2);
    }

    return b_list;
  }

  bool blocks_couple(const multi1d<DisplacedBlock>& b, 
		     const Set& S,
		     const int len,
		     int blocks_to_check)
  {
    if (blocks_to_check <= 0)
      return false;

    LatticeInteger One  = 1;
    LatticeInteger flag = 1;

    if (blocks_to_check > b.size())
      blocks_to_check = b.size();

    for(int k(0); k < blocks_to_check; k++)
    {
      LatticeInteger blk = zero;
      blk[S[b[k].blk]] = One;
      displace(blk,len,b[k].disp);
      flag *= blk;
    }

    //QDPIO::cout<<" Coupling flag "<<sum(flag)<<endl;
    if (toBool(sum(flag) > 0))
      return true;
    else
      return false;
  }

  bool blocks_couple(const multi1d<DisplacedBlock>& b, 
		     const Set& S,
		     const int len)
  {
    return blocks_couple(b,S,len,b.size());
  }

} // namespace Chroma
