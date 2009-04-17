// -*- C++ -*-
// $Id: mesfield.h,v 3.1 2009-04-17 02:05:37 bjoo Exp $
/*! \file
 *  \brief Calculates the antihermitian field strength tensor  iF(mu,nu)
 */

#ifndef __mesfield_h__
#define __mesfield_h__

namespace Chroma 
{

  //! Calculates the antihermitian field strength tensor  iF(mu,nu)
  /* 
   * \ingroup glue
   *
   *    F(mu,nu) =  (1/4) sum_p (1/2) [ U_p(x) - U^dag_p(x) ]
   *
   *  where
   *    U_1 = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
   *    U_2 = u(x,nu)*u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu)
   *    U_3 = u_dag(x-mu,mu)*u_dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
   *    U_4 = u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u_dag(x,mu)

   *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

   * Arguments:

   *  \param f   field strength tensor f(cb,mu,nu) (Write)
   *  \param u   gauge field (Read)
   */


  void mesField(multi1d<LatticeColorMatrixF>& f,
		const multi1d<LatticeColorMatrixF>& u);

  void mesField(multi1d<LatticeColorMatrixD>& f,
		const multi1d<LatticeColorMatrixD>& u);

}  // end namespace Chroma

#endif
