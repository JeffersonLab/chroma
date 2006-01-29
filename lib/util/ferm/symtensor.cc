// $Id: symtensor.cc,v 1.1 2006-01-29 05:21:55 edwards Exp $
/*! \file
 *  \brief Compute symmetric tensors
 */

#include "chromabase.h"
#include "util/ferm/symtensor.h"

namespace Chroma 
{

  //! Anonymous namespace
  /*! \ingroup ferm */
  namespace
  {
    multi3d<int> sym_tensor3d;
    bool initP = false;


    //! Init gamma matrix multiplication table factors
    /*!
     * \ingroup ferm
     *
     * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
     */
    void symTensor3dInit()
    {
      START_CODE();

      sym_tensor3d.resize(3,3,3);
      sym_tensor3d = 0;

      // d = |\epsilon^{i,j,k}|
      // Permutations: +(0,1,2)+(1,2,0)+(2,0,1)+(1,0,2)+(0,2,1)+(2,1,0)

      sym_tensor3d(0,1,2) = 1.0;
      sym_tensor3d(1,2,0) = 1.0;
      sym_tensor3d(2,0,1) = 1.0;

      sym_tensor3d(1,0,2) = 1.0;
      sym_tensor3d(0,2,1) = 1.0;
      sym_tensor3d(2,1,0) = 1.0;
  
      initP = true;

      END_CODE();
    }
  } // end namespace


  //! Return 3d symmetric tensor
  /*!
   * \ingroup ferm
   *
   * \return  \f$s_{ijk} = |\epsilon_{ijk}|\f$
   */
  int symTensor3d(int i, int j, int k)
  {
    if (! initP)
      symTensor3dInit();

    return sym_tensor3d(i,j,k);
  }

}  // end namespace Chroma
