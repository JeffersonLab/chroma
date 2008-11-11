// $Id: etensor.cc,v 3.1 2008-11-11 21:27:42 edwards Exp $
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
    multi3d<Real> E_tensor3d;
    bool initP = false;


    //! Init E table factors
    /*!
     * \ingroup ferm
     */
    void ETensor3dInit()
    {
      START_CODE();

      E_tensor3d.resize(2,3,3);
      E_tensor3d = zero;
      
      // \f$Q_{\alpha jk} = 0\quad j\ne k, Q_{111}=1/sqrt(2),Q_{122}=-1/sqrt(2), Q_{211}=-1/sqrt(6), Q_{222}=-1/sqrt(6), Q_{233}=2/sqrt(6)\f$
      E_tensor3d(0,0,0) =  1.0/sqrt(Real(2));
      E_tensor3d(0,1,1) = -1.0/sqrt(Real(2));

      E_tensor3d(1,0,0) = -1.0/sqrt(Real(6));
      E_tensor3d(1,1,1) = -1.0/sqrt(Real(6));
      E_tensor3d(1,2,2) =  2.0/sqrt(Real(6));
  
      initP = true;

      END_CODE();
    }
  } // end namespace


  //! Return E antisymmetric tensor
  /*!
   * \ingroup ferm
   *
   * \return  \f$Q_{\alpha jk} = 0\quad j\ne k, Q_{111}=1/sqrt(2),Q_{122}=-1/sqrt(2), Q_{211}=-1/sqrt(6), Q_{222}=-1/sqrt(6), Q_{233}=2/sqrt(6)\f$
   */
  Real ETensor3d(int alpha, int j, int k)
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
