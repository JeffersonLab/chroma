#ifndef __coarsenable_op_w_h__
#define __coarsenable_op_w_h__

#include "linearop.h"


namespace Chroma 
{ 


template<class T>
class MGCoarsenableOperator
{
public:
  virtual ~MGCoarsenableOperator() {}

  virtual void applyLocal(
      T& chi,
      const T& psi,
      PlusMinus isign) const
  {
    QDPIO::cout << "MGCoarsenableOperator: applyLocal not implemented\n";
  };

  virtual void applyDirection(
      T& chi,
      const T& psi,
      PlusMinus isign,
      int dir ) const
  {
    QDPIO::cout << "MGCoarsenableOperator: applyDirection not implemented\n";
  }

  //virtual int numDirections() const = 0;
};


} // namespace

#endif
