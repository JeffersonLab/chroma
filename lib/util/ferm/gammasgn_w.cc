// $Id: gammasgn_w.cc,v 3.1 2007-10-29 05:04:08 edwards Exp $
/*! \file
 *  \brief Compute gamma matrix multiplication table factors
 */

#include "chromabase.h"
#include "util/ferm/gammasgn_w.h"

namespace Chroma 
{

  //! Anonymous namespace
  /*! \ingroup ferm */
  namespace
  {
    multi2d<int> meson_eta2;
    bool initP = false;


    //! Init gamma matrix multiplication table factors
    /*!
     * \ingroup ferm
     *
     * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
     */
    void gammaSgnInit()
    {
      START_CODE();

      if (Ns != 4)
      {
	QDPIO::cerr << __func__ << ": only supports Ns=4 currently" << endl;
	QDP_abort(1);
      }

      meson_eta2.resize(Ns*Ns,Ns*Ns);

      for(int n = 0; n < Ns*Ns; ++n)
      {
	int nn = n;
	int n1 = nn & 1;
	nn >>= 1;
	int n2 = nn & 1;
	nn >>= 1;
	int n3 = nn & 1;
	int n4 = nn >> 1;

//    int ssum = (n1 + n2 + n3 + n4) & 1;

	for(int m = 0; m < Ns*Ns; ++m)
	{
	  int mm = m;
	  int m1 = mm & 1;
	  mm >>= 1;
	  int m2 = mm & 1;
	  mm >>= 1;
	  int m3 = mm & 1;
	  int m4 = mm >> 1;

	  int ssum = ((n2+n3+n4)*m1 + (n3+n4)*m2 + n4*m3) & 1;
	  meson_eta2(n,m) = (ssum == 0) ? 1 : -1;
	}
      }
  
      initP = true;

      END_CODE();
    }
  } // end namespace

  //! Return gamma matrix multiplication table factors
  /*!
   * \ingroup ferm
   *
   * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
   */
  int gammaSgn(int n, int m)
  {
    if (! initP)
      gammaSgnInit();

    return meson_eta2(n,m);
  }

}  // end namespace Chroma
