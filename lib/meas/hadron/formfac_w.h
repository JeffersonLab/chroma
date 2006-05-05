// -*- C++ -*-
// $Id: formfac_w.h,v 3.1 2006-05-05 03:07:20 edwards Exp $
/*! \file
 *  \brief Form-factors 
 *
 *  Form factors constructed from a quark and a sequential quark propagator
 */

#ifndef __formfac_h__
#define __formfac_h__

#include "util/ft/sftmom.h"

namespace Chroma 
{

 /*
  * Structures for hadron parts
  *
  * \ingroup hadron
  *
  * @{
  */
  struct FormFac_momenta_t
  {
    int              magic;     // magic number for sanity checks
    multi1d<int>     inser_mom;
    multi1d<Complex> local_current;
    multi1d<Complex> nonlocal_current;
  };

  struct FormFac_insertion_t
  {
    int              gamma_value;
    multi1d<FormFac_momenta_t> momenta;
  };

  struct FormFac_insertions_t
  {
    int  output_version;   // Unique id for each output version of the structures
    multi1d<FormFac_insertion_t>  formFac;
  };


  // Readers and writers
  void read(BinaryReader& bin, FormFac_momenta_t& mom);
  void read(BinaryReader& bin, FormFac_insertion_t& mes);
  void read(BinaryReader& bin, FormFac_insertions_t& form);
  void write(BinaryWriter& bin, const FormFac_momenta_t& mom);
  void write(BinaryWriter& bin, const FormFac_insertion_t& mes);
  void write(BinaryWriter& bin, const FormFac_insertions_t& form);


  /*! @} */  // end of group hadron



  //! Compute contractions for current insertion 3-point functions.
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * \param form               structures holding formfactors ( Write )
   * \param u                  gauge fields (used for non-local currents) ( Read )
   * \param quark_propagator   quark propagator ( Read )
   * \param seq_quark_prop     sequential quark propagator ( Read )
   * \param gamma_insertion    extra gamma insertion at source ( Read )
   * \param phases             fourier transform phase factors ( Read )
   * \param t0                 cartesian coordinates of the source ( Read )
   */

  void FormFac(FormFac_insertions_t& form,
	       const multi1d<LatticeColorMatrix>& u,
	       const LatticePropagator& quark_propagator,
	       const LatticePropagator& seq_quark_prop, 
	       int gamma_insertion,
	       const SftMom& phases,
	       int t0);

}  // end namespace Chroma

#endif
