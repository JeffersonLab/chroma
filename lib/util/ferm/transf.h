// -*- C++ -*-
// $Id: transf.h,v 3.0 2006-04-03 04:59:12 edwards Exp $

#ifndef __transf_h__
#define __transf_h__

namespace Chroma
{
 //! Convert (insert) a LatticeColorVector into a LatticeFermion
 /*!
  * \ingroup ferm
  *
  * \param a      Source color vector
  * \param b      Destination fermion
  * \param spin_index   Spin index
  */
  void CvToFerm(const LatticeColorVectorF& a, LatticeFermionF& b, 
		int spin_index);

  //! Convert (insert) a LatticeColorVector into a LatticeFermion
  /*!
   * \ingroup ferm
   *
   * \param a      Source color vector
   * \param b      Destination fermion
   * \param spin_index   Spin index
   */
  void CvToFerm(const LatticeColorVectorD& a, LatticeFermionD& b, 
		int spin_index);



 //! Convert (insert) a LatticeColorVector into a LatticeStaggeredFermion
 /*!
  * \ingroup ferm
  *
  * \param a      Source color vector
  * \param b      Destination fermion
  * \param spin_index   Spin index
  */
  void CvToFerm(const LatticeColorVectorF& a, LatticeStaggeredFermionF& b);

  //! Convert (insert) a LatticeColorVector into a LatticeStaggeredFermion
  /*!
   * \ingroup ferm
   *
   * \param a      Source color vector
   * \param b      Destination fermion
   * \param spin_index   Spin index
   */
  void CvToFerm(const LatticeColorVectorD& a, LatticeStaggeredFermionD& b);



  //! Insert a LatticeFermion into a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source fermion
   * \param b      Destination propagator
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void FermToProp(const LatticeFermionF& a, LatticePropagatorF& b, 
		  int color_index, int spin_index);

  //! Insert a LatticeFermion into a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source fermion
   * \param b      Destination propagator
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void FermToProp(const LatticeFermionD& a, LatticePropagatorD& b, 
		  int color_index, int spin_index);



  //! Insert a LatticeFermion into a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source fermion
   * \param b      Destination propagator
   * \param color_index  Color index
   */
  void FermToProp(const LatticeStaggeredFermionF& a, LatticeStaggeredPropagatorF& b, 
		  int color_index);

  //! Insert a LatticeFermion into a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source fermion
   * \param b      Destination propagator
   * \param color_index  Color index
   */
  void FermToProp(const LatticeStaggeredFermionD& a, LatticeStaggeredPropagatorD& b, 
		  int color_index);



  //! Extract a LatticeFermion from a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source propagator
   * \param b      Destination fermion
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void PropToFerm(const LatticePropagatorF& b, LatticeFermionF& a, 
		  int color_index, int spin_index);

  //! Extract a LatticeFermion from a LatticePropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source propagator
   * \param b      Destination fermion
   * \param color_index  Color index
   * \param spin_index   Spin index
   */
  void PropToFerm(const LatticePropagatorD& b, LatticeFermionD& a, 
		  int color_index, int spin_index);



  //! Extract a LatticeStaggeredFermion from a LatticeStaggeredPropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source propagator
   * \param b      Destination fermion
   * \param color_index  Color index
   */
  void PropToFerm(const LatticeStaggeredPropagatorF& b, LatticeStaggeredFermionF& a, 
		  int color_index);

  //! Extract a LatticeStaggeredFermion from a LatticeStaggeredPropagator
  /*!
   * \ingroup ferm
   *
   * \param a      Source propagator
   * \param b      Destination fermion
   * \param color_index  Color index
   */
  void PropToFerm(const LatticeStaggeredPropagatorD& b, LatticeStaggeredFermionD& a, 
		  int color_index);
}


#endif


