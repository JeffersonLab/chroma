//  $Id: diquark_w.cc,v 1.3 2007-02-28 03:28:16 edwards Exp $
/*! \file
 *  \brief Construct a diquark object
 */

#include "meas/hadron/diquark_w.h"

namespace Chroma 
{
 
  //! Unpack a quark
  /*!
   * \ingroup hadron
   *
   * We need this fast, so at the expense of a lot of memory we will
   * expose all the color/spin indices of each propagator into a temporary
   */
  multi2d< multi2d<LatticeComplex> > unpackQuark(const LatticePropagator& quark_propagator)
  {
    multi2d< multi2d<LatticeComplex> > qc(Ns,Ns);

    for(int scol=0; scol < Ns; ++scol)             // spin col
    {
      for(int srow=0; srow < Ns; ++srow)           // spin row
      {
	LatticeColorMatrix color_mat = peekSpin(quark_propagator,
						srow,scol);  // (srow,scol)

	qc(srow,scol).resize(Nc,Nc);

	for(int ccol=0; ccol < Nc; ++ccol)         // color col
	  for(int crow=0; crow < Nc; ++crow)       // color row
	  {
	    qc(srow,scol)(crow,ccol) = peekColor(color_mat,
						 crow,ccol);  // (crow,ccol)
	  }
      }
    }

    return qc;
  }
  

  //! Construct a QQ diquark object
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all baryons the colour components are contracted with the totally
   * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   *
   * \param diquark                  diquark object (in real space) ( Write )
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   */

  void QQDiquark(QQDiquarkContract_t& diquark, 
		 const LatticePropagator& quark_propagator_1, 
		 const LatticePropagator& quark_propagator_2)
  {
    START_CODE();

    // We need this fast, so at the expense of a lot of memory we will
    // expose all the color/spin indices of each propagator into a temporary
    multi2d< multi2d<LatticeComplex> > qc_1(unpackQuark(quark_propagator_1));
    multi2d< multi2d<LatticeComplex> > qc_2(unpackQuark(quark_propagator_2));

    // Final diquark contractions
    multi1d<int> ranks(6);
    ranks = Ns;
    ranks[4] = Nc;  // cf
    ranks[5] = Nc;  // ci

    diquark.comp.resize(ranks);

    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])             // sf_2
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])           // sf_1
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])         // si_2
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])       // si_1
	  {
	    const multi2d<LatticeComplex>& q_1 = qc_1(ranks[3],ranks[1]);  // Note: oddball transpose
	    const multi2d<LatticeComplex>& q_2 = qc_2(ranks[2],ranks[0]);  // Note: oddball transpose

	    // Contract over color indices with antisym tensors
	    // Permutations: +(0,1,2)+(1,2,0)+(2,0,1)-(1,0,2)-(0,2,1)-(2,1,0)
	    // d = (\epsilon^{i1,j1,k1}\epsilon^{i2,j2,k2} Q1^{i1,i2} Q2^{j1,j2})^{k1,k2}
	    ranks[4] = 2; ranks[5] = 2;
	    diquark.comp[ranks] = q_1(0,0)*q_2(1,1) - q_1(1,0)*q_2(0,1) - q_1(0,1)*q_2(1,0) + q_1(1,1)*q_2(0,0);

	    ranks[4] = 0; ranks[5] = 0;
	    diquark.comp[ranks] = q_1(1,1)*q_2(2,2) - q_1(2,1)*q_2(1,2) - q_1(1,2)*q_2(2,1) + q_1(2,2)*q_2(1,1);

	    ranks[4] = 1; ranks[5] = 1;
	    diquark.comp[ranks] = q_1(2,2)*q_2(0,0) - q_1(0,2)*q_2(2,0) - q_1(2,0)*q_2(0,2) + q_1(0,0)*q_2(2,2);

	    ranks[4] = 0; ranks[5] = 2;
	    diquark.comp[ranks] = q_1(1,0)*q_2(2,1) - q_1(2,0)*q_2(1,1) - q_1(1,1)*q_2(2,0) + q_1(2,1)*q_2(1,0);

	    ranks[4] = 0; ranks[5] = 1;
	    diquark.comp[ranks] = q_1(1,2)*q_2(2,0) - q_1(2,2)*q_2(1,0) - q_1(1,0)*q_2(2,2) + q_1(2,0)*q_2(1,2);
		  
	    ranks[4] = 1; ranks[5] = 2;
	    diquark.comp[ranks] = q_1(2,0)*q_2(0,1) - q_1(0,0)*q_2(2,1) - q_1(2,1)*q_2(0,0) + q_1(0,1)*q_2(2,0);

	    ranks[4] = 1; ranks[5] = 0;
	    diquark.comp[ranks] = q_1(2,1)*q_2(0,2) - q_1(0,1)*q_2(2,2) - q_1(2,2)*q_2(0,1) + q_1(0,2)*q_2(2,1);

	    ranks[4] = 2; ranks[5] = 0;
	    diquark.comp[ranks] = q_1(0,1)*q_2(1,2) - q_1(1,1)*q_2(0,2) - q_1(0,2)*q_2(1,1) + q_1(1,2)*q_2(0,1);

	    ranks[4] = 2; ranks[5] = 1;
	    diquark.comp[ranks] = q_1(0,2)*q_2(1,0) - q_1(1,2)*q_2(0,0) - q_1(0,0)*q_2(1,2) + q_1(1,0)*q_2(0,2);
	  }

    END_CODE();
  }


}  // end namespace Chroma
