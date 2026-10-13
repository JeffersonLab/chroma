#ifndef __coarsenable_op_w_h__
#define __coarsenable_op_w_h__

#include "linearop.h"


namespace Chroma 
{


  enum class MGOperatorForm
    {
      UNPRECONDITIONED,         // M
      EVEN_ODD_PRECONDITIONED   // S_o = M_oo - M_oe M_ee^{-1} M_eo
    };



  template<class T>
  class MGCoarsenableOperator
  {
  public:
    virtual ~MGCoarsenableOperator() {}

    virtual MGOperatorForm mgForm() const = 0;

    virtual const Subset& paritySubset(int p) const = 0;

    virtual bool hasLocal(int p, int q) const = 0;
    virtual bool hasDirection(int p, int q) const = 0;

    virtual void applyLocal(T& chi, const T& psi, enum PlusMinus isign,
                            int p, int q) const = 0;

    virtual void applyDirection(T& chi, const T& psi, enum PlusMinus isign,
                                int dir, int p, int q) const = 0;

    virtual void liftToEven(T& v_e, const T& v_o) const
    {
      QDPIO::cerr << "MGCoarsenableOperator::liftToEven: "
                  << "not available for this operator" << std::endl;
      QDP_abort(1);
    }


    void applyLocalFull(T& chi, const T& psi, enum PlusMinus isign) const
    {
      assertCheckerboarded();
      chi = zero;
      T tmp;
      for (int p = 0; p < 2; ++p)
        for (int q = 0; q < 2; ++q)
          if (hasLocal(p, q))
            {
              tmp = zero;
              applyLocal(tmp, psi, isign, p, q);
              chi[paritySubset(p)] += tmp;
            }
    }

    void applyDirectionFull(T& chi, const T& psi, enum PlusMinus isign,
                            int dir) const
    {
      assertCheckerboarded();
      chi = zero;
      T tmp;
      for (int p = 0; p < 2; ++p)
        for (int q = 0; q < 2; ++q)
          if (hasDirection(p, q))
            {
              tmp = zero;
              applyDirection(tmp, psi, isign, dir, p, q);
              chi[paritySubset(p)] += tmp;
            }
    }

  private:
    void assertCheckerboarded() const
    {
      if (paritySubset(0).numSiteTable() + paritySubset(1).numSiteTable()
          != all.numSiteTable())
        {
          QDPIO::cerr << "MGCoarsenableOperator: parity spaces are separate "
                      << "fields; use the block interface" << std::endl;
          QDP_abort(1);
        }
    }
  };



  template<class T>
  const MGCoarsenableOperator<T>&
  asCoarsenable(const LinearOperator<T>& op, int level)
  {
    const MGCoarsenableOperator<T>* c = dynamic_cast<const MGCoarsenableOperator<T>*>(&op);

    if (c == 0)
      {
        QDPIO::cerr << "MG level " << level << ": operator does not implement "
                    << "MGCoarsenableOperator" << std::endl;
        QDP_abort(1);
      }

    const Subset& expect =
      (c->mgForm() == MGOperatorForm::EVEN_ODD_PRECONDITIONED)
      ? c->paritySubset(1) : all;

    if (op.subset().numSiteTable() != expect.numSiteTable())
      {
        QDPIO::cerr << "MG level " << level << ": operator subset does not "
                    << "match its declared MGOperatorForm" << std::endl;
        QDP_abort(1);
      }

    return *c;
  }



  
} // namespace

#endif
