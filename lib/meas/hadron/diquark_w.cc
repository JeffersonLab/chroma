//  $Id: diquark_w.cc,v 1.1 2007-02-22 06:58:55 edwards Exp $
/*! \file
 *  \brief Construct a diquark object
 */

#include "meas/hadron/diquark_w.h"

namespace Chroma 
{
 
  //! Serialize generalized object
  multi1d<LatticeComplex> QQDiquarkContract_t::serialize()
  {
    multi1d<int> ranks(6);

    multi1d<LatticeComplex> barprop_1d(Ns*Ns*Ns*Ns*Nc*Nc);

    int cnt = 0;
    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])             // sf_2
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])           // sf_1
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])         // si_2
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])       // si_1
	    for(ranks[4] = 0; ranks[4] < Nc; ++ranks[4])   // cf
	      for(ranks[5] = 0; ranks[5] < Nc; ++ranks[5]) // ci
	      {
		diquark_1d[cnt++] = comp[ranks];
	      }

    return diquark_1d;
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
    multi2d< multi2d<LatticeComplex> > qc_1(Ns,Ns);
    multi2d< multi2d<LatticeComplex> > qc_2(Ns,Ns);

    for(int sf=0; sf < Ns; ++sf)             // sf
    {
      for(int si=0; si < Ns; ++si)           // si
      {
	LatticeColorMatrix color_mat1 = peekSpin(quark_propagator_1,
						 si,sf);  // (si,sf)
	LatticeColorMatrix color_mat2 = peekSpin(quark_propagator_2,
						 si,sf);  // (si,sf)

	qc_1(sf,si).resize(Nc,Nc);
	qc_2(sf,si).resize(Nc,Nc);

	for(int cf=0; cf < Nc; ++cf)         // cf
	  for(int ci=0; ci < Nc; ++ci)       // ci
	  {
	    qc_1(sf,si)(cf,ci) = peekColor(color_mat1,
					   ci,cf);  // (ci,cf)
	    
	    qc_2(sf,si)(cf,ci) = peekColor(color_mat2,
					   ci,cf);  // (ci,cf)
	  }
      }
    }


    // Final diquark contractions
    multi1d<int> ranks(6);
    ranks = Ns;
    ranks[4] = Nc;
    ranks[5] = Nc;

    diquark.comp.resize(ranks);

    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])             // sf_2
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])           // sf_1
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])         // si_2
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])       // si_1
	  {
	    const multi2d<LatticeComplex>& q_1 = qc_1(ranks[1],ranks[3]);
	    const multi2d<LatticeComplex>& q_2 = qc_2(ranks[0],ranks[2]);

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
      }

    END_CODE();
  }


}  // end namespace Chroma
