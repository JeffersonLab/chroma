// -*- C++ -*-

#ifndef UNIT_CHECK_INCLUDE
#define UNIT_CHECK_INCLUDE

namespace Chroma 
{
  //! Check the unitarity of color matrix in SU(N)
  /*!
   * \ingroup gauge
   *
   * \param  u  The multi1d<LatticeColorMatrix> to be tested
   */
  template <typename Q>
    void unitarityCheck(const multi1d<Q>& u);

}

#endif
