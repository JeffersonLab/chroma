/*! \file
 *  \brief Gaussian smearing of color std::vector
 */

#include "chromabase.h"
#include "meas/smear/gaus_smear.h"
#include "actions/boson/operator/klein_gord.h"

namespace Chroma 
{

  //! Do a covariant Gaussian smearing of a lattice field
  /*!
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color std::vector field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  template<typename T>
  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 T& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    T psi;

    Real ftmp = - (width*width) / Real(4*ItrGaus);
    /* The Klein-Gordon operator is (Lapl + mass_sq), where Lapl = -d^2/dx^2.. */
    /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */
    Real ftmpi = Real(1) / ftmp;
  
    for(int n = 0; n < ItrGaus; ++n)
    {
      psi = chi * ftmp;
      klein_gord(u, psi, chi, ftmpi, j_decay);
    }
  }


  template<>
  void gausSmear(const multi1d<LatticeColorMatrix>& u,
                 LatticePropagator& chi,
                 const Real& width, int ItrGaus, int j_decay)
  {
    Real ftmp = - (width*width) / Real(4*ItrGaus);
    Real ftmpi = Real(1) / ftmp;

    for (int sa = 0 ; sa < Ns ; sa++ )
      for (int sb = 0 ; sb < Ns ; sb++ )
	{
          LatticeColorMatrix psi_a;
          LatticeColorMatrix chi_a = peekSpin( chi , sa , sb );

          for(int n = 0; n < ItrGaus; ++n)
            {
              psi_a = chi_a * ftmp;
              klein_gord(u, psi_a, chi_a, ftmpi, j_decay);
            }

          pokeSpin( chi , chi_a , sa , sb );
        }
  }
  

  //! Do a covariant Gaussian smearing of a lattice color std::vector field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color std::vector field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeColorVector& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    gausSmear<LatticeColorVector>(u, chi, width, ItrGaus, j_decay);
  }


  //! Do a covariant Gaussian smearing of a lattice fermion field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      fermion field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeFermion& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    gausSmear<LatticeFermion>(u, chi, width, ItrGaus, j_decay);
  }


  //! Do a covariant Gaussian smearing of a lattice propagator field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      propagator field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticeStaggeredPropagator& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    gausSmear<LatticeStaggeredPropagator>(u, chi, width, ItrGaus, j_decay);
  }


  //! Do a covariant Gaussian smearing of a lattice propagator field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      propagator field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void gausSmear(const multi1d<LatticeColorMatrix>& u, 
		 LatticePropagator& chi, 
		 const Real& width, int ItrGaus, int j_decay)
  {
    gausSmear<LatticePropagator>(u, chi, width, ItrGaus, j_decay);
  }


}  // end namespace Chroma
