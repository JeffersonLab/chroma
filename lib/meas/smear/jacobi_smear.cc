// $Id: jacobi_smear.cc,v 3.1 2009-02-20 15:10:24 edwards Exp $
/*! \file
 *  \brief Jacobi smearing of color vector
 */

#include "chromabase.h"
#include "meas/smear/jacobi_smear.h"

namespace Chroma 
{

    //! Do a covariant Jacobi smearing of a lattice field
    /*!
     * Arguments:
     *
     *  \param u             gauge field ( Read )
     *  \param chi           propagator field ( Modify )
     *  \param kappa         hopping parameter ( Read )
     *  \param iter          number of iterations ( Read )
     *  \param no_smear_dir  no smearing in this direction ( Read )
     */

    template<typename T>
    void jacobiSmear(const multi1d<LatticeColorMatrix>& u, 
		     T& chi, 
		     const Real& kappa, int iter, int no_smear_dir)
    {
	T psi;
	Real norm;

	T s_0,h_smear;
	s_0 = chi;

	for(int n = 0; n < iter; ++n)
	    {
		psi = chi;
		bool first = true;

		for(int mu = 0; mu < Nd; ++mu )
		    if( mu != no_smear_dir )
			{
			    if (first)
				h_smear =  u[mu]*shift(psi, FORWARD, mu) + shift(adj(u[mu])*psi, BACKWARD, mu);
			    else
				h_smear += u[mu]*shift(psi, FORWARD, mu) + shift(adj(u[mu])*psi, BACKWARD, mu);
			    first = false;
			}
		chi = s_0 + kappa * h_smear;
	    }
    }


    //! Do a covariant Jacobi smearing of a lattice color vector field
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
		     LatticeColorVector& chi, 
		     const Real& kappa, int iter, int no_smear_dir)
    {
	jacobiSmear<LatticeColorVector>(u, chi, kappa, iter, no_smear_dir);
    }


    //! Do a covariant Jacobi smearing of a lattice fermion field
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
		     LatticeFermion& chi, 
		     const Real& kappa, int iter, int no_smear_dir)
    {
	jacobiSmear<LatticeFermion>(u, chi, kappa, iter, no_smear_dir);
    }


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
		     const Real& kappa, int iter, int no_smear_dir)
    {
	jacobiSmear<LatticeStaggeredPropagator>(u, chi, kappa, iter, no_smear_dir);
    }


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
		     const Real& kappa, int iter, int no_smear_dir)
    {
	jacobiSmear<LatticePropagator>(u, chi, kappa, iter, no_smear_dir);
    }


}  // end namespace Chroma
