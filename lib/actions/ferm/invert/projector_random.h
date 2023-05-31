// -*- C++ -*-
/*! \file
 *  \brief Oblique projector to a random vector (useful for testing)
 */

#ifndef __projector_random_h__
#define __projector_random_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"

namespace Chroma
{

  //! Solve a M*psi=chi linear system by BICGSTAB
  /*! \ingroup invert
   */
  template<typename T>
  class ProjectorRandom : public Projector<T>
  {
  public:
    using Ts = const std::vector<std::shared_ptr<T>>;
    using const_Ts = const std::vector<std::shared_ptr<const T>>;

    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     */
    ProjectorRandom(Handle<LinearOperator<T>> A_) : A(A_)
    {
      T u0 = zero, v0 = zero;
      random(v0, A->subset());
      v = v0 / sqrt(norm2(v0, A->subset()));
      (*A)(u0, v, PLUS);
      l = sqrt(norm2(u0, A->subset()));
      u = u0 / l;
    }

    //! Destructor is automatic
    ~ProjectorRandom() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const
    {
      return A->subset();
    }

    //! Apply the oblique projector A*V*inv(U^H*A*V)*U^H
    /*! 
     *! Returns A*V*inv(U^H*A*V)*U^H*chi = psi
     */
    void AVUObliqueProjector(Ts& psi, const_Ts& chi) const override
    {
      for (int i = 0; i < psi.size(); ++i)
      {
	// Return u * (u^* * chi)
	*psi[i] = u * innerProduct(u, *chi[i], A->subset());
      }
    }

    //! Apply the oblique projector V*inv(U^H*A*V)*U^H*A
    /*! 
     * Returns V*inv(U^H*A*V)*U^H*A*chi = psi
     */
    void VUAObliqueProjector(Ts& psi, const_Ts& chi) const override
    {
      for (int i = 0; i < psi.size(); ++i)
      {
	T A_chi = zero;
	(*A)(A_chi, *chi[i], PLUS);
	// Return v * (u^* A * chi) / l
	*psi[i] = v * (innerProduct(u, A_chi, A->subset()) / l);
      }
    }

    //! Rank of the projector, which is the rank of U and V also
    unsigned int rank() const override
    {
      return 1;
    }

    //! Return U[i]
    void U(unsigned int i, T& psi) const override
    {
      psi = u;
    }

    //! Return V[i]
    void V(unsigned int i, T& psi) const override
    {
      psi = v;
    }

    //! Return U[i]^H*A*V[i]
    void lambda(unsigned int i, DComplex& lambda) const override
    {
      lambda = l;
    }

  private:
    Handle<LinearOperator<T>> A;
    T u, v;
    DComplex l; ///< = u^* * A * v
  };

  //! Random projector namespace
  namespace ProjectorRandomEnv
  {
    Projector<LatticeFermion>* createProjector(
      XMLReader&, const std::string&,
      Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>,
      Handle<LinearOperator<LatticeFermion>> A)
    {
      return new ProjectorRandom<LatticeFermion>(A);
    }

    //! Register the projector
    inline bool registerAll() {
      static bool registered = false;
      if (registered) return true;
      registered = true;
      return Chroma::TheLinOpFermProjectorFactory::Instance().registerObject("RANDOM_PROJECTOR",
									     createProjector);
    }
  }
} // End namespace

#endif 
