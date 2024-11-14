// -*- C++ -*-
/*! \file
 *  \brief Oblique projector to a random vector (useful for testing)
 */

#ifndef __projector_null_h__
#define __projector_null_h__

#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"

namespace Chroma
{

  //! Return a null projector
  /*! \ingroup invert
   */
  template<typename T>
  class ProjectorNull : public Projector<T>
  {
  public:
    using Ts = const std::vector<std::shared_ptr<T>>;
    using const_Ts = const std::vector<std::shared_ptr<const T>>;

    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     */
    ProjectorNull(Handle<LinearOperator<T>> A_) : A(A_)
    {
    }

    //! Destructor is automatic
    ~ProjectorNull() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const
    {
      return A->subset();
    }

    //! Apply the oblique projector A*V*inv(U^H*A*V)*U^H
    /*! 
     *! Returns A*V*inv(U^H*A*V)*U^H*chi = psi
     */
    void AVUObliqueProjector(Ts& psi, const_Ts&) const override
    {
      for (int i = 0; i < psi.size(); ++i)
      {
	*psi[i] = zero;
      }
    }

    //! Apply the oblique projector V*inv(U^H*A*V)*U^H*A
    /*! 
     * Returns V*inv(U^H*A*V)*U^H*A*chi = psi
     */
    void VUAObliqueProjector(Ts& psi, const_Ts&) const override
    {
      for (int i = 0; i < psi.size(); ++i)
      {
	*psi[i] = zero;
      }
    }

    //! Rank of the projector, which is the rank of U and V also
    unsigned int rank() const override
    {
      return 0;
    }

    //! Return U[i]
    void U(unsigned int, T&) const override
    {
      throw std::runtime_error("ProjectorNull: rank of the projector is null");
    }

    //! Return V[i]
    void V(unsigned int, T&) const override
    {
      throw std::runtime_error("ProjectorNull: rank of the projector is null");
    }

    //! Return U[i]^H*A*V[i]
    void lambda(unsigned int, DComplex&) const override
    {
      throw std::runtime_error("ProjectorNull: rank of the projector is null");
    }

  private:
    Handle<LinearOperator<T>> A;
  };

  //! Null projector namespace
  namespace ProjectorNullEnv
  {
    Projector<LatticeFermion>* createProjector(
      XMLReader&, const std::string&,
      Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>>,
      Handle<LinearOperator<LatticeFermion>> A)
    {
      return new ProjectorNull<LatticeFermion>(A);
    }

    //! Register the projector
    inline bool registerAll() {
      static bool registered = false;
      if (registered) return true;
      registered = true;
      return Chroma::TheLinOpFermProjectorFactory::Instance().registerObject("NULL_PROJECTOR",
									     createProjector);
    }
  }
} // End namespace

#endif 
