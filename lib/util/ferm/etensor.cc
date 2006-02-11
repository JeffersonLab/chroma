// $Id: etensor.cc,v 1.1 2006-02-11 17:22:38 edwards Exp $
/*! \file
 *  \brief Tensor used for E representations
 */

#include "chromabase.h"
#include "util/ferm/etensor.h"

namespace Chroma 
{

  //! Anonymous namespace
  /*! \ingroup ferm */
  namespace
  {
    multi3d<int> E_tensor3d;
    bool initP = false;


    //! Init E table factors
    /*!
     * \ingroup ferm
     */
    void ETensor3dInit()
    {
      START_CODE();

      E_tensor3d.resize(2,3,3);
      E_tensor3d = 0;

      // \f$S_{\alpha jk} = 0\quad j\ne k, S_{111}=S_{222}=+1, S_{122}=S_{233}=-1\f$
      E_tensor3d(0,0,0) = 1.0;
      E_tensor3d(1,1,1) = 1.0;

      E_tensor3d(0,1,1) = -1.0;
      E_tensor3d(1,2,2) = -1.0;
  
      initP = true;

      END_CODE();
    }
  } // end namespace


  //! Return E antisymmetric tensor
  /*!
   * \ingroup ferm
   *
   * \return  \f$S_{\alpha jk} = 0\quad j\ne k, S_{111}=S_{222}=+1, S_{122}=S_{233}=-1\f$
   */
  int ETensor3d(int alpha, int j, int k)
  {
    if (! initP)
      ETensor3dInit();

    return E_tensor3d(alpha,j,k);
  }

}  // end namespace Chroma
