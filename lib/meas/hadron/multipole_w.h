// -*- C++ -*-
// $Id: multipole_w.h,v 1.3 2005-03-27 18:10:17 edwards Exp $
/*! \file
 *  \brief Multipole moments
 *
 *  Compute multipole moment within two hadron states using spherical
 *  harmonic expansion
 */

#ifndef __multipole_h__
#define __multipole_h__

namespace Chroma 
{

  //! Storage structure to hold electric and magnetic multipole moments
  /*!
   * \ingroup hadron
   *
   * Here I mean  array_L<array_M<array_T>> 
   */
  struct Multipole_t
  {
    struct ElecMag_t
    {
      multi1d<Complex>  electric;
      multi1d<Complex>  magnetic;
    };

    int  j_decay;
    multi1d< multi1d< ElecMag_t > >  corr;
  };


  //! Write a Multipole_t
  /*!
   * \ingroup hadron
   */
  void write(BinaryWriter& bin, const Multipole_t& pole);

  
  //! Compute contractions for multipole moments
  /*!
   * \ingroup hadron
   *
   * \param pole               structures holding formfactors ( Write )
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param GammaInsertion     extra gamma matrix insertion ( Read )
   * \param max_power          max value of L ( Read )
   * \param j_decay            direction of decay ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   */

  void multipole(Multipole_t& pole,
		 const LatticePropagator& quark_propagator,
		 const LatticePropagator& seq_quark_prop, 
		 int   GammaInsertion,
		 int   max_power,
		 int   j_decay,
		 int   t0);

}  // end namespace Chroma

#endif
