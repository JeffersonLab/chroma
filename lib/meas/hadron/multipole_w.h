// -*- C++ -*-
// $Id: multipole_w.h,v 3.0 2006-04-03 04:59:00 edwards Exp $
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
      int  L;
      int  M;
      multi1d<Complex>  electric;
      multi1d<Complex>  magnetic;
    };

    int  j_decay;
    multi1d< multi1d< ElecMag_t > >  corr;
  };


  //! Read a Multipole_t
  /*!
   * \ingroup hadron
   */
  void read(XMLReader& xml, const string& path, Multipole_t& pole);

  //! Write a Multipole_t
  /*!
   * \ingroup hadron
   */
  void write(XMLWriter& xml, const string& path, const Multipole_t& pole);

  
  //! Compute contractions for multipole moments
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param GammaInsertion     extra gamma matrix insertion ( Read )
   * \param max_power          max value of L ( Read )
   * \param j_decay            direction of decay ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   * \param xml                xml file object ( Write )
   * \param xml_group          string used for writing xml data ( Read )
   */

  void multipole(const LatticePropagator& quark_propagator,
		 const LatticePropagator& seq_quark_prop, 
		 int   GammaInsertion,
		 int   max_power,
		 int   j_decay,
		 int   t0,
		 XMLWriter& xml,
		 const string& xml_group);

}  // end namespace Chroma

#endif
