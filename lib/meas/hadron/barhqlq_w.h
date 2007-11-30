// -*- C++ -*-
// $Id: barhqlq_w.h,v 3.2 2007-11-30 06:38:20 kostas Exp $
/*! \file
 *  \brief Heavy-light baryon 2-pt functions
 */

#ifndef __barhqlq_w_h__
#define __barhqlq_w_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Baryon 2pt contractions
  /*! \ingroup hadron */
  namespace Baryon2PtContractions
  {
    //! Sigma 2-pt
    /*! \ingroup hadron */
    LatticeComplex sigma2pt(const LatticePropagator& quark_propagator_1,
			    const LatticePropagator& quark_propagator_2,
			    const SpinMatrix& T, const SpinMatrix& sp);

    //! Sigma 2-pt
    /*! \ingroup hadron */
    LatticeComplex xi2pt(const LatticePropagator& quark_propagator_1,
			 const LatticePropagator& quark_propagator_2,
			 const SpinMatrix& T, const SpinMatrix& sp);

    //! Lambda 2-pt
    /*! \ingroup hadron */
    LatticeComplex lambda2pt(const LatticePropagator& quark_propagator_1,
			     const LatticePropagator& quark_propagator_2,
			     const SpinMatrix& T, const SpinMatrix& sp);

    //! Lambda 2-pt
    /*! \ingroup hadron */
    LatticeComplex lambdaNaive2pt(const LatticePropagator& quark_propagator_1,
				  const LatticePropagator& quark_propagator_2,
				  const SpinMatrix& T, const SpinMatrix& sp);

    //! Delta 2-pt
    /*! \ingroup hadron */
    LatticeComplex sigmast2pt(const LatticePropagator& quark_propagator_1,
			      const LatticePropagator& quark_propagator_2,
			      const SpinMatrix& T, const SpinMatrix& sp);

    //! Delta 2-pt
    /*! \ingroup hadron */
    LatticeComplex sigmast2pt(const LatticePropagator& quark_propagator_1,
			      const LatticePropagator& quark_propagator_2,
			      const SpinMatrix& T, 
			      const SpinMatrix& spSRC, 
			      const SpinMatrix& spSNK );

  }  // namespace  Baryon2PtContractions


  //! Heavy-light baryon 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * Construct baryon propagators for the Proton and the Delta^+ with
   * degenerate "u" and "d" quarks, as well as the Lambda for, in
   * addition, a degenerate "s" quark. For these degenerate quarks, the
   * Lambda is degenerate with the Proton, but we keep it for compatibility
   * with the sister routine that treats non-degenerate quarks.

   * The routine optionally computes time-charge reversed baryons and adds them
   * in for increased statistics.

   * \param propagator_1   "s" quark propagator ( Read )
   * \param propagator_2   "u" quark propagator ( Read )
   * \param t0             cartesian coordinates of the source ( Read )
   * \param bc_spec        boundary condition for spectroscopy ( Read )
   * \param time_rev       add in time reversed contribution if true ( Read )
   * \param phases         object holds list of momenta and Fourier phases ( Read )
   * \param xml            xml file object ( Read )
   * \param xml_group      group name for xml data ( Read )
   *
   */

  void barhqlq(const LatticePropagator& propagator_1, 
	       const LatticePropagator& propagator_2, 
	       const SftMom& phases,
	       int t0, int bc_spec, bool time_rev,
	       XMLWriter& xml,
	       const string& xml_group);



  //! Heavy-light baryon 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   *###########################################################################
   * WARNING: No symmetrization over the spatial part of the wave functions   #
   *          is performed. Therefore, if this routine is called with         #
   *          "shell-sink" quark propagators of different widths the          #
   *          resulting octet baryons may have admixters of excited           #
   *          decouplet baryons with mixed symmetric spatial wave functions,  #
   *          and vice-versa!!!                                               #
   *                                                                          #
   * WARNING 2: The time reversal will be wrong if both quark propagators     #
   *            are identical!                                                #
   *###########################################################################

   * Construct heavy-light baryon propagators with two "u" quarks and
   * one separate "s" quark for the Sigma^+, the Lambda and the Sigma^{*+}.
   * In the Lambda we take the "u" and "d" quark as degenerate!

   * The routine also computes time-charge reversed baryons and adds them
   * in for increased statistics.

   * \param quark_propagator_1   "s" quark propagator ( Read )
   * \param quark_propagator_2   "u" quark propagator ( Read )
   * \param barprop              baryon propagator ( Modify )
   * \param phases               object holds list of momenta and Fourier phases ( Read )
   *
   *        ____
   *        \
   * b(t) =  >  < b(t_source, 0) b(t + t_source, x) >
   *        /                    
   *        ----
   *          x

   * For the Sigma^+ we take

   * |S_1, s_z=1/2> = (s C gamma_5 u) "u_up"

   * for the Lambda

   * |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
   *                  + (u C gamma_5 s) "d_up"

   * and for the Sigma^{*+}

   * |S*_1, s_z=3/2> = 2*(s C gamma_- u) "u_up" + (u C gamma_- u) "s_up".

   * We have put "q_up" in quotes, since this is meant in the Dirac basis,
   * not in the 'DeGrand-Rossi' chiral basis used in the program!
   * In gamma_- we ignore a factor sqrt(2).

   * For all baryons we compute a 'B_2' that differs from the 'B_1' above
   * by insertion of a gamma_4 between C and the gamma_{5,-}.
   * And finally, we also compute the non-relativistic baryons, 'B_3',
   * which up to a factor 1/2 are just the difference B_1 - B_2, as can
   * be seen by projecting to the "upper" components in the Dirac basis,
   * achieved by (1 + gamma_4)/2 q, for quark q.

   * The Sigma^+_k is baryon 3*(k-1), the Lambda_k is baryon 3*(k-1)+1
   * and the Sigma^{*+}_k is baryon 3*(k-1)+2.

   * We are using a chiral basis for the Dirac matrices (gamma_5 diagonal).
   * Therefore a spin-up quark in the Dirac basis corresponds to
   * 1/sqrt(2) * ( - q_1 - q_3 ) in this chiral basis. We shall neglect
   * the sign and the 1/sqrt(2) here.
   * The projection on "spin_up" is done with S_proj. 
   */

  void barhqlq(const LatticePropagator& propagator_1,
	       const LatticePropagator& propagator_2,
	       const SftMom& phases,
	       multi3d<DComplex>& barprop);

}  // end namespace Chroma


#endif
