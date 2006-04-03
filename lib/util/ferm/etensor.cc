// $Id: etensor.cc,v 3.0 2006-04-03 04:59:11 edwards Exp $
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
      E_tensor3d(0,0,0) = 1;
      E_tensor3d(1,1,1) = 1;

      E_tensor3d(0,1,1) = -1;
      E_tensor3d(1,2,2) = -1;
  
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

    if (alpha < 0 || alpha > 1 || j < 0 || j > 2 || k < 0 || k > 2)
    {
      QDPIO::cerr << __func__ << ": indices out of bounds: alpha,j,k=" 
		  << alpha << " " << j << " " << k << endl;
      QDP_abort(1);
    }

    return E_tensor3d(alpha,j,k);
  }

}  // end namespace Chroma
