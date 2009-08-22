// $Id: qactden.cc,v 3.2 2009-08-22 20:05:49 edwards Exp $
/*! \file
 *  \brief Measure the lattice density of the lattice energy and the naive topological charge.
 */

#include "chromabase.h"
#include "meas/glue/qnaive.h"

namespace Chroma 
{

  //! Measure the lattice density of the lattice energy and the naive topological charge.
  /*!
   * \ingroup glue
   *
   * \param lrqtop  topological charge density (Write)
   * \param lract   action to continuum instanton action density (Write) 
   * \param u       gauge field (Read)
   */

  void qactden(LatticeReal& lract, LatticeReal& lrqtop, const multi1d<LatticeColorMatrix>& u)
  {
    START_CODE();
  
    LatticeColorMatrix u_clov1;
    LatticeColorMatrix u_clov2;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeColorMatrix tmp_3;
    LatticeReal qtop_tmp;
    LatticeReal plaq_tmp;
  
    /* Lattice version of S_ratio */
    lract = Real(2*Nd*(Nd-1)*Nc);
  
    /* Lattice version of Qtop */
    lrqtop = zero;
  
    /* Loop over checkerboards and triplet of perpendicular planes */
    int mu1 = 0;
    for(int nu1=1; nu1 < Nd; ++nu1)
    {
      /* First "plus-plus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1 = shift(u[nu1], FORWARD, mu1);

      /* tmp_2(x) = u(x+nu1,mu1) */
      tmp_2 = shift(u[mu1], FORWARD, nu1);

      /* tmp_3(x) = tmp_1 * tmp_2^dag = u(x+mu1,nu1)*u_dag(x+nu1,mu1) */
      tmp_3 = tmp_1 * adj(tmp_2);

      /* tmp_1(x) = tmp_3 * u_dag(x,nu1)= u(x+mu1,nu1)*u_dag(x+nu1,mu1)*u_dag(x,nu1) */
      tmp_1 = tmp_3 * adj(u[nu1]);

      /* u_clov1(x) = u(x,mu1) * tmp_1 = u(x,mu1)*u(x+mu1,nu1)* */
      /*                                 u_dag(x+nu1,mu1)*u_dag(x,nu1) */
      u_clov1 = u[mu1] * tmp_1;
      
      plaq_tmp = real(trace(u_clov1));
      lract -= plaq_tmp;
      
      /* First "plus-minus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1 = shift(u[nu1], FORWARD, mu1);

      /* tmp_2(x) = u(x,mu1) * tmp_1 = u(x,mu1)*u(x+mu1,nu1) */
      tmp_2 = u[mu1] * tmp_1;

      /* tmp_1(x) = tmp_2_dag * u(x,nu1) = u_dag(x+mu1,nu1)*u_dag(x,mu1)*u(x,nu1) */
      tmp_1 = adj(tmp_2) * u[nu1];

      /* tmp_2(x) = tmp_1(x-nu1) */
      tmp_2 = shift(tmp_1, BACKWARD, nu1);

      /* tmp_1(x) = u(x,mu1) * tmp_2 = u(x,mu1)*u_dag(x-nu1+mu1,nu1)* */
      /*                               u_dag(x-nu1,mu1)*u(x-nu1,nu1) */
      tmp_1 = u[mu1] * tmp_2;

      u_clov1 -= tmp_1;

      plaq_tmp = real(trace(tmp_1));
      lract -= plaq_tmp;
      
      /* First "minus-minus" plaquette */
      /* tmp_1(x) = u(x+mu1,nu1) */
      tmp_1 = shift(u[nu1], FORWARD, mu1);

      /* tmp_2(x) = u(x,mu1) * tmp_1 = u(x,mu1)*u(x+mu1,nu1) */
      tmp_2 = u[mu1] * tmp_1;

      /* tmp_1(x) = u_dag(x,nu1) * tmp_2 = u_dag(x,nu1)*u(x,mu1)*u(x+mu1,nu1) */
      tmp_1 = adj(u[nu1]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-nu1) */
      tmp_2 = shift(tmp_1, BACKWARD, nu1);

      /* tmp_1(x) = u_dag(x,mu1) * tmp_2 = u_dag(x,mu1)*u_dag(x-nu1,nu1)* */
      /*                                   u(x-nu1,mu1)*u(x-nu1+mu1,nu1) */
      tmp_1 = adj(u[mu1]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-mu1) */
      tmp_2 = shift(tmp_1, BACKWARD, mu1);

      u_clov1 += tmp_2;

      plaq_tmp = real(trace(tmp_2));
      lract -= plaq_tmp;
      
      /* First "minus-plus" plaquette */
      /* tmp_1(x) = u_dag(x,mu1) * u(x,nu1) */
      tmp_1 = adj(u[mu1]) * u[nu1];

      /* tmp_2(x) = u(x+nu1,mu1) */
      tmp_2 = shift(u[mu1], FORWARD, nu1);

      /* tmp_3(x) = tmp_1 * tmp_2 = u_dag(x,mu1)*u(x,nu1)*u(x+nu1,mu1) */
      tmp_3 = tmp_1 * tmp_2;

      /* tmp_1(x) = tmp_3(x-mu1) */
      tmp_1 = shift(tmp_3, BACKWARD, mu1);

      /* tmp_2(x) = tmp_1 * u_dag(x,nu1) = u_dag(x-mu1,mu1)*u(x-mu1,nu1)* */
      /*                                   u(x-mu1+nu1,mu1)*u_dag(x,nu1) */
      tmp_2 = tmp_1 * adj(u[nu1]);

      u_clov1 -= tmp_2;

      plaq_tmp = real(trace(tmp_2));
      lract -= plaq_tmp;
      
      int mu2 = (nu1 % 3) + 1;
      int nu2 = (mu2 % 3) + 1;

      /* Second "plus-plus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1 = shift(u[nu2], FORWARD, mu2);

      /* tmp_2(x) = u(x+nu2,mu2) */
      tmp_2 = shift(u[mu2], FORWARD, nu2);

      /* tmp_3(x) = tmp_1 * tmp_2^dag = u(x+mu2,nu2)*u_dag(x+nu2,mu2) */
      tmp_3 = tmp_1 * adj(tmp_2);

      /* tmp_1(x) = tmp_3 * u_dag(x,nu2)= u(x+mu2,nu2)*u_dag(x+nu2,mu2)*u_dag(x,nu2) */
      tmp_1 = tmp_3 * adj(u[nu2]);

      /* u_clov2(x) = u(x,mu2) * tmp_1 = u(x,mu2)*u(x+mu2,nu2)* */
      /*                                 u_dag(x+nu2,mu2)*u_dag(x,nu2) */
      u_clov2 = u[mu2] * tmp_1;

      plaq_tmp = real(trace(u_clov2));
      lract -= plaq_tmp;
      
      /* Second "plus-minus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1 = shift(u[nu2], FORWARD, mu2);

      /* tmp_2(x) = u(x,mu2) * tmp_1 = u(x,mu2)*u(x+mu2,nu2) */
      tmp_2 = u[mu2] * tmp_1;

      /* tmp_1(x) = tmp_2_dag * u(x,nu2) = u_dag(x+mu2,nu2)*u_dag(x,mu2)*u(x,nu2) */
      tmp_1 = adj(tmp_2) * u[nu2];

      /* tmp_2(x) = tmp_1(x-nu2) */
      tmp_2 = shift(tmp_1, BACKWARD, nu2);

      /* tmp_1(x) = u(x,mu2) * tmp_2 = u(x,mu2)*u_dag(x-nu2+mu2,nu2)* */
      /*                               u_dag(x-nu2,mu2)*u(x-nu2,nu2) */
      tmp_1 = u[mu2] * tmp_2;

      u_clov2 -= tmp_1;

      plaq_tmp = real(trace(tmp_1));
      lract -= plaq_tmp;
      
      /* Second "minus-minus" plaquette */
      /* tmp_1(x) = u(x+mu2,nu2) */
      tmp_1 = shift(u[nu2], FORWARD, mu2);

      /* tmp_2(x) = u(x,mu2) * tmp_1 = u(x,mu2)*u(x+mu2,nu2) */
      tmp_2 = u[mu2] * tmp_1;

      /* tmp_1(x) = u_dag(x,nu2) * tmp_2 = u_dag(x,nu2)*u(x,mu2)*u(x+mu2,nu2) */
      tmp_1 = adj(u[nu2]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-nu2) */
      tmp_2 = shift(tmp_1, BACKWARD, nu2);

      /* tmp_1(x) = u_dag(x,mu2) * tmp_2 = u_dag(x,mu2)*u_dag(x-nu2,nu2)* */
      /*                                   u(x-nu2,mu2)*u(x-nu2+mu2,nu2) */
      tmp_1 = adj(u[mu2]) * tmp_2;

      /* tmp_2(x) = tmp_1(x-mu2) */
      tmp_2 = shift(tmp_1, BACKWARD, mu2);

      u_clov2 += tmp_2;

      plaq_tmp = real(trace(tmp_2));
      lract -= plaq_tmp;
      
      /* Second "minus-plus" plaquette */
      /* tmp_1(x) = u_dag(x,mu2) * u(x,nu2) */
      tmp_1 = adj(u[mu2]) * u[nu2];

      /* tmp_2(x) = u(x+nu2,mu2) */
      tmp_2 = shift(u[mu2], FORWARD, nu2);

      /* tmp_3(x) = tmp_1 * tmp_2 = u_dag(x,mu2)*u(x,nu2)*u(x+nu2,mu2) */
      tmp_3 = tmp_1 * tmp_2;

      /* tmp_1(x) = tmp_3(x-mu2) */
      tmp_1 = shift(tmp_3, BACKWARD, mu2);

      /* tmp_2(x) = tmp_1 * u_dag(x,nu2) = u_dag(x-mu2,mu2)*u(x-mu2,nu2)* */
      /*                                   u(x-mu2+nu2,mu2)*u_dag(x,nu2) */
      tmp_2 = tmp_1 * adj(u[nu2]);

      u_clov2 -= tmp_2;

      plaq_tmp = real(trace(tmp_2));
      lract -= plaq_tmp;
      
      /* Now comes the contribution to the topological charge */
      tmp_1 = 1;
      tmp_2 = adj(u_clov1) * tmp_1;
      u_clov1 -= tmp_2;
      tmp_2 = adj(u_clov2) * tmp_1;
      u_clov2 -= tmp_2;
      tmp_2 = u_clov1 * u_clov2;

      qtop_tmp = real(trace(tmp_2));
      lrqtop -= qtop_tmp;
    }  
            
    /* Lattice version of S_ratio */
    lract /= (4*Chroma::twopi*Chroma::twopi);
  
    /* Lattice version of qtop */
    lrqtop /= (64*Chroma::twopi*Chroma::twopi);
  
    END_CODE();
  }

} // namespace Chroma
