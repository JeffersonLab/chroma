// -*- C++ -*-
/*! \file
 *  \brief Jacobi smearing of color std::vector
 */

#ifndef __jacobi_smear_h__
#define __jacobi_smear_h__

namespace Chroma 
{

  //! Do a covariant Jacobi smearing of a lattice color std::vector field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u             gauge field ( Read )
   *  \param chi           color std::vector field ( Modify )
   *  \param kappa         hopping parameter ( Read )
   *  \param iter          number of iterations ( Read )
   *  \param no_smear_dir  no smearing in this direction ( Read )
   */
  void jacobiSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeColorVector& chi, 
		 const Real& kappa, int iter, int no_smear_dir);


  //! Do a covariant Jacobi smearing of a lattice fermion field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u             gauge field ( Read )
   *  \param chi           fermion field ( Modify )
   *  \param kappa         hopping parameter ( Read )
   *  \param iter          number of iterations ( Read )
   *  \param no_smear_dir  no smearing in this direction ( Read )
   */
  void jacobiSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeFermion& chi, 
		 const Real& kappa, int iter, int no_smear_dir);


  //! Do a covariant Jacobi smearing of a lattice propagator field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u             gauge field ( Read )
   *  \param chi           propagator field ( Modify )
   *  \param kappa         hopping parameter ( Read )
   *  \param iter          number of iterations ( Read )
   *  \param no_smear_dir  no smearing in this direction ( Read )
   */
  void jacobiSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeStaggeredPropagator& chi, 
		 const Real& kappa, int iter, int no_smear_dir);


  //! Do a covariant Jacobi smearing of a lattice propagator field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u             gauge field ( Read )
   *  \param chi           propagator field ( Modify )
   *  \param kappa         hopping parameter ( Read )
   *  \param iter          number of iterations ( Read )
   *  \param no_smear_dir  no smearing in this direction ( Read )
   */
  void jacobiSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticePropagator& chi, 
		 const Real& kappa, int iter, int no_smear_dir);

}  // end namespace Chroma

#endif
