// $Id: shift2.cc,v 1.1 2006-08-24 02:34:42 edwards Exp $
/*! \file
 *  \brief Shift by a power of 2
 */

#include "chromabase.h"
#include "util/gauge/shift2.h"

namespace Chroma 
{ 

  //! A simple not-fancy power of 2 shift
  /*! \ingroup gauge */
  LatticeColorMatrix shift2(const LatticeColorMatrix& s1, int isign, int dir, int level)
  {
    LatticeColorMatrix d;

    if (level == 0)
    {
      d = shift(s1,isign,dir);
    }
    else
    {
      LatticeColorMatrix tmp = shift(s1,isign,dir);
      d = shift(tmp,isign,dir);

      int cnt = 1 << (level -1);
      for (; cnt-- > 0; )
      {
	tmp = shift(d,isign,dir);
	d   = shift(tmp,isign,dir);
      }
    }

    return d;
  }

}
